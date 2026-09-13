#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
runS1interferogram - single command to produce one GrIMP unwrapped-
interferogram product from two Sentinel-1 SLC SAFE files.

    runS1interferogram REFERENCE.zip SECONDARY.zip [--project uwproject.yaml]

Run it from your base environment. The Python workflow runs in base; the
isce-dependent commands (topsApp.py, convert_isce) are dispatched into the
`isce2grimp` conda env with `conda run` behind the scenes (see isceenv.py), so
you never have to activate it.

It consolidates the old prep_pair -> run_isce -> azPhaseCorrect ->
convert_isce/SETide/setupisceuw chain into one step:

  1. read the project config (scratch dir per machine, DEM, orbit dir, region
     data, ISCE parameters)
  2. build a topsApp.xml for the pair (no downloading; the SAFE zips are
     assumed to be on disk in dataDir; .zip and .zip.1 are both accepted)
  3. run topsApp through the burst interferogram step (--end=burstifg) in a
     machine-local scratch dir
  4. apply the azimuth burst correction, remove simulated phase, unwrap once,
     and remap to the GrIMP product (azPhaseCorrect.runAzPhaseCorrect)
  5. copy the low-volume final product to outputRoot and, unless --debug,
     delete the large ISCE working files from scratch.

The heavy ISCE files (~2 GB/pair) live only in the machine-local scratch dir to
avoid NFS disk contention; only the final ~150 MB product is written to the
shared output volume.

Selection of which two SLCs to process is handled separately - here a valid
pair is assumed to be given on the command line.
"""
import argparse
import datetime
import glob
import os
import re
import shutil
import signal
import sys
import threading
import time
import traceback
import zipfile


_LOGFILE = None

# Burst-time -> frame constant, matching findS1InsarPairs and
# isce2grimp.convert_isce.get_frame_number
BTIME = 2.759


def setLogFile(path):
    ''' Route log() lines (and isce command output) to this file. '''
    global _LOGFILE
    _LOGFILE = path
    isceenv.setLog(path)


def log(msg):
    ''' Print a timestamped line to the console and the run log. '''
    line = f'{datetime.datetime.now():%H:%M:%S} {msg}'
    print(line, flush=True)
    if _LOGFILE is not None:
        with open(_LOGFILE, 'a') as fp:
            fp.write(line + '\n')


def timed(label, func, *a, **kw):
    ''' Run func, log its wall-clock duration, and return its result. '''
    t0 = time.time()
    result = func(*a, **kw)
    log(f'[TIME] {label}: {time.time() - t0:.1f} s')
    return result

import utilities as u
import isce2grimp.util.dinosar as dinosar

import uwproject as uwp
import azPhaseCorrect as apc
import isceenv


def cmdLineParse():
    ''' Command line parser. '''
    parser = argparse.ArgumentParser(
        description='Produce one GrIMP unwrapped-interferogram product from '
                    'two Sentinel-1 SAFE zips (runs isce2grimp behind the '
                    'scenes).')
    parser.add_argument('reference', type=str,
                        help='reference SLC SAFE zip (.zip or .zip.1); path or '
                        'basename in dataDir')
    parser.add_argument('secondary', type=str,
                        help='secondary SLC SAFE zip (.zip or .zip.1)')
    parser.add_argument('--project', type=str, default='uwproject.yaml',
                        help='project yaml [uwproject.yaml]')
    parser.add_argument('--track', type=int, default=None,
                        help='relative-orbit/track number (else read from SAFE '
                        'manifest)')
    parser.add_argument('--cpus', type=int, default=None,
                        help='OMP thread count [threads[host] or isce.cpus]')
    parser.add_argument('--stallMinutes', type=int, default=None,
                        help='abort if the run makes no progress for this '
                        'many minutes [project stallMinutes, or 30]')
    parser.add_argument('--azThreads', type=int, default=None,
                        help='bursts corrected concurrently; lower spreads '
                        'the allocation spike that has drawn oomd kills '
                        '[project azThreads, or 3]')
    parser.add_argument('--pressureLimit', type=int, default=None,
                        help='release the scratch page cache when this run\'s '
                        'cgroup memory pressure (PSI avg10) exceeds this '
                        'percent; 0 disables [project pressureLimit, or 20]')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='keep the ISCE and intermediate files in scratch')
    parser.add_argument('--keepScratch', action='store_true', default=False,
                        help='do not delete the scratch working dir')
    parser.add_argument('--remapOnly', action='store_true', default=False,
                        help='reuse the existing scratch ISCE dir and redo '
                        'only the remap to GrIMP (recovery after a failure in '
                        'convert_isce / SETide / setupUW)')
    parser.add_argument('--noAzCorrect', action='store_true', default=False,
                        help='skip the azimuth burst correction')
    parser.add_argument('--noPhaseRemove', action='store_true', default=False,
                        help='skip simulated-phase removal before unwrapping')
    return parser.parse_args()


# ---- SAFE-file handling


def resolveSafe(name, dataDir, archiveDir=None):
    ''' Absolute path to a SAFE zip. Looks for an existing path, then in
    dataDir, then across the archive's YYYY-MM subdirs (so cross-month pairs
    resolve).

    The .zip / .zip.1 suffix is not part of the granule identity, so both are
    tried at every location: .1 is the marker fileS1 writes once a granule has
    been unpacked into the assembly tree, so the same granule answers to either
    name depending on whether the nightly filing run has reached it yet. A
    queue line naming one must still resolve if filing renamed it in the
    meantime. '''
    stem = os.path.basename(name).split('.zip')[0]
    directory = os.path.dirname(name)
    for candidateName in (f'{stem}.zip', f'{stem}.zip.1'):
        asGiven = os.path.join(directory, candidateName)
        if os.path.exists(asGiven):
            return os.path.abspath(asGiven)
        candidate = os.path.join(dataDir, candidateName)
        if os.path.exists(candidate):
            return candidate
        if archiveDir:
            hits = glob.glob(os.path.join(
                archiveDir, '[0-9][0-9][0-9][0-9]-[0-9][0-9]', candidateName))
            if hits:
                return hits[0]
    u.myerror(f'resolveSafe: SAFE file not found: {stem}[.zip|.zip.1] (also '
              f'tried {dataDir} and {archiveDir}/*/)')


def safeStem(safePath):
    ''' Strip directory and any extension from .zip onward (handles .zip and
    the .zip.1 partial-download suffix). '''
    return os.path.basename(safePath).split('.zip')[0]


def parseSafe(safePath):
    ''' Parse mission / absolute orbit / acquisition date from the SAFE name.
    e.g. S1C_IW_SLC__1SDH_20251108T095750_..._004921_009BB9_B629 '''
    stem = safeStem(safePath)
    pieces = stem.split('_')
    return {'mission': pieces[0],
            'absOrbit': int(pieces[7]),
            'date': datetime.datetime.strptime(pieces[5].split('T')[0],
                                               '%Y%m%d'),
            'start': pieces[5],   # e.g. 20251108T095750 (unique per frame)
            'stop': pieces[6],
            'stem': stem}


def readManifest(safePath):
    ''' Text of manifest.safe from inside the zip, or None if it cannot be
    read. Works regardless of the .zip / .zip.1 extension. '''
    try:
        with zipfile.ZipFile(safePath) as z:
            manifest = [n for n in z.namelist() if n.endswith('manifest.safe')]
            if not manifest:
                return None
            return z.read(manifest[0]).decode('utf-8', errors='ignore')
    except zipfile.BadZipFile:
        return None


def trackFromManifest(manifest):
    ''' Relative-orbit (track) number from the manifest text. '''
    if manifest is None:
        return None
    match = re.search(r'relativeOrbitNumber[^>]*>(\d+)<', manifest)
    return int(match.group(1)) if match else None


def frameFromManifest(manifest, start):
    ''' Burst-derived GrIMP frame number for an acquisition starting at
    `start`, from the manifest's ascending-node time. Same calculation as
    findS1InsarPairs and isce2grimp.convert_isce.get_frame_number, so the log
    names line up with the <orbit1>_<frame1> product dirs. '''
    if manifest is None:
        return None
    match = re.search(r'ascendingNodeTime>([^<]+)<', manifest)
    if match is None:
        return None
    ascNode = datetime.datetime.strptime(match.group(1)[:26],
                                         '%Y-%m-%dT%H:%M:%S.%f')
    startTime = datetime.datetime.strptime(start, '%Y%m%dT%H%M%S')
    return int((startTime - ascNode).seconds / BTIME + 0.5)


def linkSafe(safePath, workDir):
    ''' Expose the SAFE zip to topsApp under a .zip name in workDir (ISCE's
    reader dispatches on the .zip extension, so a .zip.1 must be linked). '''
    zipName = f'{safeStem(safePath)}.zip'
    linkPath = os.path.join(workDir, zipName)
    if not os.path.exists(linkPath):
        os.symlink(safePath, linkPath)
    return zipName


def frameBbox(safePath, margin=0.3):
    ''' [S, N, W, E] bounding box of the SAFE footprint (+margin degrees),
    from the manifest.safe gml:coordinates (lat,lon pairs). '''
    with zipfile.ZipFile(safePath) as z:
        m = [n for n in z.namelist() if n.endswith('manifest.safe')][0]
        text = z.read(m).decode('utf-8', errors='ignore')
    coords = re.search(r'<gml:coordinates>(.*?)</gml:coordinates>', text)
    lats, lons = [], []
    for pair in coords.group(1).split():
        la, lo = pair.split(',')
        lats.append(float(la))
        lons.append(float(lo))
    return (min(lats) - margin, max(lats) + margin,
            min(lons) - margin, max(lons) + margin)


# ---- Precise orbits


def orbitFile(orbitDir, safe):
    ''' The POEORB .EOF in orbitDir whose validity window spans this
    acquisition, or None. Name form:
    S1A_OPER_AUX_POEORB_OPOD_<produced>_V<validStart>_<validStop>.EOF '''
    fmt = '%Y%m%dT%H%M%S'
    start = datetime.datetime.strptime(safe['start'], fmt)
    stop = datetime.datetime.strptime(safe['stop'], fmt)
    for eof in glob.glob(os.path.join(
            orbitDir, f'{safe["mission"]}_OPER_AUX_POEORB_*.EOF')):
        pieces = os.path.basename(eof).split('_')
        try:
            validStart = datetime.datetime.strptime(pieces[6][1:], fmt)
            validStop = datetime.datetime.strptime(pieces[7].split('.')[0],
                                                   fmt)
        except (IndexError, ValueError):
            continue
        if validStart <= start and stop <= validStop:
            return eof
    return None


def checkOrbits(project, safes):
    ''' Verify a precise orbit covers each acquisition before any processing
    starts - topsApp would otherwise fail well into the run. Raises rather than
    calling u.myerror so the failure still goes through the fail-file path. '''
    orbitDir = project['orbitDir']
    if not os.path.exists(orbitDir):
        raise RuntimeError(f'checkOrbits: orbitDir does not exist: {orbitDir}')
    for safe in safes:
        eof = orbitFile(orbitDir, safe)
        if eof is None:
            raise RuntimeError(
                f'checkOrbits: no POEORB orbit in {orbitDir} covering '
                f'{safe["stem"]} ({safe["start"]}-{safe["stop"]})')
        log(f'orbit {safe["mission"]} {safe["start"]}: '
            f'{os.path.basename(eof)}')


# ---- ISCE DEM


def prepareIsceDem(project, isceDir, bbox):
    ''' topsApp needs a lat/lon (EPSG:4326) DEM. If the project sets an explicit
    `dem`, use it; otherwise derive one for this frame by cropping+reprojecting
    the polar-stereo `demTiff` to the footprint bbox. '''
    dem = project.get('dem')
    if dem and dem != 'TBD' and os.path.exists(dem):
        return dem
    demTiff = project.get('demTiff')
    if not demTiff or not os.path.exists(demTiff):
        u.myerror('prepareIsceDem: need `demTiff` (or an explicit `dem`) to '
                  'build the topsApp DEM')
    S, N, W, E = bbox
    tr = project.get('isce', {}).get('demPosting', 0.000277)
    demOut = os.path.join(isceDir, 'dem.wgs84')
    print(f'building topsApp DEM {demOut} for bbox S{S:.3f} N{N:.3f} '
          f'W{W:.3f} E{E:.3f}')
    # Run each command as its own isceShell call (a single `conda run` only
    # dispatches the first command of a ;-chain into the env).
    isceenv.isceShell(f"gdalwarp -q -overwrite -t_srs EPSG:4326 "
                      f"-te {W} {S} {E} {N} -tr {tr} {tr} -r bilinear "
                      f"-of ISCE '{demTiff}' '{demOut}'")
    isceenv.isceShell(f"gdal2isce_xml.py -i '{demOut}'")
    return demOut


# ---- topsApp.xml


def writeTopsAppXML(project, demPath, refZip, secZip, workDir):
    ''' Build topsApp.xml in workDir from the isce2grimp template, overriding
    the DEM, orbit dir, SAFE names, looks, polarization and ionosphere flag
    with the project config. '''
    templatePath = os.path.abspath(os.path.join(
        os.path.dirname(dinosar.__file__), '..', 'data', 'template.yml'))
    inputDict = dinosar.read_yaml_template(templatePath)
    isce = project.get('isce', {})
    tops = inputDict['topsinsar']
    tops['demfilename'] = demPath
    tops['azimuthlooks'] = isce.get('azimuthlooks', tops.get('azimuthlooks'))
    tops['rangelooks'] = isce.get('rangelooks', tops.get('rangelooks'))
    tops['unwrappername'] = isce.get('unwrappername',
                                     tops.get('unwrappername'))
    tops['doionospherecorrection'] = isce.get(
        'doionospherecorrection', tops.get('doionospherecorrection'))
    pol = isce.get('polarization', 'hh')
    for role, safeZip in [('reference', refZip), ('secondary', secZip)]:
        tops[role]['safe'] = safeZip
        tops[role]['orbit directory'] = project['orbitDir']
        tops[role]['polarization'] = pol
    tops['reference']['output directory'] = 'referencedir'
    tops['secondary']['output directory'] = 'secondarydir'
    xml = dinosar.dict2xml(inputDict)
    dinosar.write_xml(xml, outname=os.path.join(workDir, 'topsApp.xml'))


# ---- Product handling


def copyProduct(productDir, outputRoot, track):
    ''' Copy the final GrIMP product to outputRoot/track-N/, dropping only the
    links that point outside it (the back-link to the intermediate dir in
    scratch) so the product is self-contained. Relative intra-product links are
    kept as links - tieScript locates the phase through the fixed-name
    phase.uw.vrt. Copies to a temp dir then swaps in, tolerating a pre-existing
    dest (incl. NFS .nfs* leftovers from an open handle). '''
    destTrackDir = os.path.join(outputRoot, f'track-{track}')
    os.makedirs(destTrackDir, exist_ok=True)
    dest = os.path.join(destTrackDir, os.path.basename(productDir))
    productRoot = os.path.realpath(productDir) + os.sep

    def ignoreLinks(directory, names):
        return [n for n in names
                if os.path.islink(os.path.join(directory, n))
                and not os.path.realpath(
                    os.path.join(directory, n)).startswith(productRoot)]

    tmp = dest + '.new'
    shutil.rmtree(tmp, ignore_errors=True)
    shutil.copytree(productDir, tmp, ignore=ignoreLinks, symlinks=True)
    # Move any existing product aside by rename (works even when a plain
    # rmtree would fail on non-empty NFS dirs), then swap the new one in.
    if os.path.exists(dest):
        old = dest + '.old'
        shutil.rmtree(old, ignore_errors=True)
        os.rename(dest, old)
        shutil.rmtree(old, ignore_errors=True)
    os.rename(tmp, dest)
    return dest


class KilledBySignal(Exception):
    ''' A fatal signal, turned into an exception so it reaches the fail-record
    path instead of ending the run silently. returncode follows the shell
    convention of 128 + signal number. '''

    def __init__(self, signalNumber):
        self.signalNumber = signalNumber
        self.returncode = 128 + signalNumber
        super().__init__(f'killed by '
                         f'{signal.Signals(signalNumber).name} '
                         f'(signal {signalNumber})')


def installSignalHandlers():
    ''' Convert the catchable fatal signals into KilledBySignal, so a job that
    is terminated still logs why and writes a Fail record.

    SIGKILL cannot be caught. That is the point: with these installed, a run
    that dies leaving no record at all was SIGKILLed - in practice the OOM
    killer - and the absence of a record becomes the diagnosis rather than an
    ambiguity. SIGINT is left alone so Ctrl-C keeps its usual behaviour; it is
    caught as KeyboardInterrupt in main instead. '''
    def handler(signalNumber, frame):
        raise KilledBySignal(signalNumber)

    for name in ('SIGHUP', 'SIGTERM'):
        signal.signal(getattr(signal, name), handler)


class StalledRun(RuntimeError):
    ''' The run stopped making progress. Carries EX_TEMPFAIL (75) so a stall
    is distinguishable from a real processing failure in exitStatus.log, and
    so a queue driver can treat it as retryable. '''

    def __init__(self, minutes, where):
        self.returncode = 75
        super().__init__(f'stalled: no CPU time used and nothing written to '
                         f'the run log for {minutes:.1f} min; {where}')


def stackDump():
    ''' Innermost frame of every live thread, to name what the run was doing
    when it stopped. This is what was missing from the hung Aug-24 jobs: the
    logs ended mid-stage with no indication of which read never returned. '''
    lines = []
    for threadId, frame in sorted(sys._current_frames().items()):
        if threadId == threading.get_ident():   # the watchdog itself
            continue
        stack = traceback.extract_stack(frame)
        inner = stack[-1]
        lines.append(f'  thread {threadId}: '
                     f'{os.path.basename(inner.filename)}:{inner.lineno} '
                     f'in {inner.name}(): {(inner.line or "").strip()}')
    return lines


def ownCgroupDir():
    ''' This process's cgroup v2 directory, or None if unavailable. '''
    try:
        with open('/proc/self/cgroup') as fp:
            for line in fp:
                field = line.strip().split(':')
                if field[0] != '0':        # 0:: is the unified hierarchy
                    continue
                path = os.path.join('/sys/fs/cgroup', field[2].lstrip('/'))
                if os.path.exists(os.path.join(path, 'memory.current')):
                    return path
    except OSError:
        pass
    return None


def cgroupBytes(cgroupDir):
    ''' Current charge against this cgroup: anon plus its share of the page
    cache. This is the number systemd-oomd weighs, which is why the page cache
    a run leaves behind can get it killed. '''
    try:
        with open(os.path.join(cgroupDir, 'memory.current')) as fp:
            return int(fp.read())
    except (OSError, ValueError):
        return 0


def cgroupPressure(cgroupDir):
    ''' PSI memory-pressure avg10 for this cgroup, as a percentage.

    This is the signal systemd-oomd actually acts on -- not a byte count. It
    measures the share of time this cgroup's tasks spent stalled waiting on
    memory reclaim, so a run streaming tens of GB through the page cache
    registers pressure even on a machine with hundreds of GB free. Shedding
    well below oomd's limit keeps the run clear of it.

    Prefers the "full" line (every task stalled), which is what oomd compares
    against; falls back to "some" where "full" is absent. '''
    avg10 = {}
    try:
        with open(os.path.join(cgroupDir, 'memory.pressure')) as fp:
            for line in fp:
                field = line.split()
                for item in field[1:]:
                    key, _, value = item.partition('=')
                    if key == 'avg10':
                        avg10[field[0]] = float(value)
    except (OSError, ValueError, IndexError):
        return 0.
    return avg10.get('full', avg10.get('some', 0.))


def dropFileCache(root, minBytes=64 * 1024 ** 2):
    ''' Release the page cache holding the large files under root, returning
    the bytes released.

    A pair writes ~86 GB of ISCE scratch, and those file pages stay charged to
    the cgroup that touched them long after the stage that wrote them is done.
    They are clean and reclaimable, so nothing is lost by handing them back --
    the kernel re-reads on demand -- but while they sit there they count
    toward the memory-pressure limit that killed the Aug-24 runs.

    fsync first: only clean pages can be dropped. Symlinks are skipped so this
    stays inside the scratch tree and does not touch the source SAFE archive.
    '''
    released = 0
    for dirPath, _, fileNames in os.walk(root):
        for name in fileNames:
            path = os.path.join(dirPath, name)
            try:
                if os.path.islink(path) or os.path.getsize(path) < minBytes:
                    continue
                fd = os.open(path, os.O_RDONLY)
            except OSError:
                continue
            try:
                os.fsync(fd)
                os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
                released += os.path.getsize(path)
            except OSError:
                pass
            finally:
                os.close(fd)
    return released


def startWatchdog(logPath, failPath, scratchPath, limitMinutes,
                  pressureLimit=None, interval=60):
    ''' Abort the run if it stops making progress, instead of hanging until
    someone notices days later.

    Progress is either CPU time used by this process and its children, or a
    write to the run log. Between them they cover every stage: topsApp runs
    with our own CPU idle but tees continuously to the log, and the in-process
    stages burn CPU without touching the log. A stalled run does neither.

    The watchdog thread reports and exits the process itself rather than
    signalling the main thread, because the main thread may be blocked in a C
    extension where a signal would not be delivered until the call returns.
    Note the limit of this, and of any in-process approach: a thread stuck in
    an uninterruptible kernel wait (a hard NFS mount that never answers) dies
    only when that wait ends -- os._exit cannot preempt it either. '''
    def cpuSeconds():
        t = os.times()
        return t.user + t.system + t.children_user + t.children_system

    def logMtime():
        try:
            return os.path.getmtime(logPath)
        except OSError:
            return 0.

    cgroupDir = ownCgroupDir() if pressureLimit else None

    def governor(selfCpu, governorHold):
        ''' Release the scratch page cache when this cgroup's memory pressure
        rises, before it reaches the limit oomd acts on. Returns the updated
        (selfCpu, governorHold). '''
        pressure = cgroupPressure(cgroupDir)
        if pressure <= pressureLimit:
            return selfCpu, governorHold
        used = cgroupBytes(cgroupDir)
        t0 = time.thread_time()
        dropFileCache(scratchPath)
        selfCpu += time.thread_time() - t0
        freed = used - cgroupBytes(cgroupDir)
        # Back off when little comes back: the pressure is then coming from
        # something other than our scratch (another job sharing the scope),
        # and re-walking the tree every tick would cost more than it saves.
        quiet = freed < 1024 ** 3
        if quiet:
            governorHold = time.time() + 10 * interval
        note = ' (little left here, backing off)' if quiet else ''
        log(f'governor: memory pressure {pressure:.0f}%, cgroup at '
            f'{used / 1024 ** 3:.0f} GB, released '
            f'{max(freed, 0) / 1024 ** 3:.1f} GB{note}')
        return selfCpu, governorHold

    def monitor():
        lastCpu, lastLog, lastProgress = cpuSeconds(), logMtime(), time.time()
        selfCpu, governorHold, lastCheck = 0., 0., time.time()
        # avg10 is a 10 s average and oomd acts after 20 s over its limit, so
        # the governor is polled on a short tick; the stall check, which is
        # comparing minutes, stays on the full interval.
        tick = min(interval, 10)
        while True:
            time.sleep(tick)
            if cgroupDir is not None and time.time() >= governorHold:
                selfCpu, governorHold = governor(selfCpu, governorHold)
            if time.time() - lastCheck < interval:
                continue
            lastCheck = time.time()
            # Discount the governor's own CPU, or shedding cache would
            # register as progress and mask a genuine stall.
            cpu, mtime = cpuSeconds() - selfCpu, logMtime()
            # 1 s of CPU over the interval: enough to clear timer noise, far
            # below what any real stage uses.
            if cpu > lastCpu + 1.0 or mtime > lastLog:
                lastCpu, lastLog, lastProgress = cpu, mtime, time.time()
                continue
            stalled = (time.time() - lastProgress) / 60.
            if stalled < limitMinutes:
                continue
            threads = stackDump()
            exc = StalledRun(stalled, f'{len(threads)} threads, innermost '
                                      f'frames follow')
            log(f'WATCHDOG: {exc}')
            for line in threads:
                log(line)
            log(f'scratch kept for debugging: {scratchPath}')
            log(f'fail record: '
                f'{writeFailFile(failPath, exc.returncode, exc, logPath)}')
            sys.stdout.flush()
            os._exit(exc.returncode)

    threading.Thread(target=monitor, daemon=True,
                     name='watchdog').start()


def exitCode(exc):
    ''' Exit code to report for a failure: the ISCE command's own code or a
    signal's 128+n when we have it, the code from a u.myerror abort
    (SystemExit) otherwise, else 1. '''
    code = getattr(exc, 'returncode', None)
    if code is None and isinstance(exc, SystemExit):
        code = exc.code
    if code is None and isinstance(exc, KeyboardInterrupt):
        code = 128 + signal.SIGINT
    return code if isinstance(code, int) and code != 0 else 1


def failMessage(exc):
    ''' One-line description of a failure. A u.myerror abort raises a bare
    SystemExit carrying no message (myerror prints it to the console), so fall
    back to the source line that aborted - skipping myerror itself. '''
    if isinstance(exc, KeyboardInterrupt):
        return 'interrupted from the terminal (Ctrl-C / SIGINT)'
    message = str(exc)
    if message:
        return message
    for frame in reversed(traceback.extract_tb(exc.__traceback__)):
        if os.path.basename(frame.filename) != 'myerror.py':
            return (f'{type(exc).__name__} at '
                    f'{os.path.basename(frame.filename)}:{frame.lineno}: '
                    f'{frame.line}')
    return type(exc).__name__


def writeFailFile(failPath, code, exc, logPath):
    ''' Record the failure and its exit code beside the run log, so a queue
    driver can spot a failed pair without parsing the log. '''
    with open(failPath, 'w') as fp:
        fp.write(f'exitCode: {code}\n')
        fp.write(f'error: {failMessage(exc)}\n')
        fp.write(f'log: {logPath}\n')
    return failPath


def main():
    args = cmdLineParse()
    # Setup is guarded separately from the run below: u.myerror calls a bare
    # sys.exit(), which is SystemExit(None) and so exits 0 - a missing SAFE or
    # an unresolvable track would otherwise be recorded as a success by the
    # caller. No Fail record is possible this early (its name is built from the
    # very values being computed here), so this reports and exits non-zero; the
    # queue's exitStatus.log is what captures it.
    try:
        project = uwp.loadProject(args.project)
        isceenv.setEnv(project.get('isceEnv'))
        regionDefs = uwp.regionDefsFromProject(project)
        scratch = uwp.scratchForHost(project)
        cpus = args.cpus if args.cpus is not None \
            else uwp.threadsForHost(project)

        # Resolve and parse the two SAFE files
        refPath = resolveSafe(args.reference, project['dataDir'],
                              project.get('archiveDir'))
        secPath = resolveSafe(args.secondary, project['dataDir'],
                              project.get('archiveDir'))
        ref = parseSafe(refPath)
        sec = parseSafe(secPath)
        manifest = readManifest(refPath)
        track = args.track if args.track is not None \
            else trackFromManifest(manifest)
        if track is None:
            u.myerror('main: could not determine track; pass --track')
        # Frame is for naming only, so an unreadable manifest is not fatal here
        frame = frameFromManifest(manifest, ref['start'])
        if frame is None:
            print('main: could not determine frame from manifest; '
                  'naming logs 0')
            frame = 0
    except (Exception, SystemExit, KeyboardInterrupt) as exc:
        code = exitCode(exc)
        print(f'FAILED during setup (exit {code}): {failMessage(exc)}')
        sys.exit(code)
    o1, o2 = ref['absOrbit'], sec['absOrbit']
    # Frames of one datatake share absolute orbits, so tag scratch/log with the
    # reference acquisition start time to keep frames from colliding.
    frameTag = ref['start']

    # Per-pair scratch working area (holds isce/, gimp/, intermediate/)
    pairScratch = os.path.join(scratch, f'{track}-{o1}-{o2}-{frameTag}')
    isceDir = os.path.join(pairScratch, 'isce')
    gimpDir = os.path.join(pairScratch, 'gimp', f'track-{track}')
    intermediatePath = os.path.join(pairScratch, 'intermediate', f'{o1}-{o2}')

    # Persistent run log (survives scratch cleanup), in a per-run-date dir so
    # one directory holds everything processed that day. The run timestamp in
    # the name keeps repeat runs of a pair from colliding.
    runTime = datetime.datetime.now()
    runStamp = f'{runTime:%Y%m%dT%H%M%S}'
    logDir = os.path.join(project['outputRoot'], 'logs',
                          f'{runTime:%Y-%m-%d}')
    os.makedirs(logDir, exist_ok=True)
    logPath = os.path.join(logDir, f'{track}.{o1}.{frame}.{runStamp}.log')
    open(logPath, 'w').close()   # fresh log per run
    setLogFile(logPath)
    failPath = os.path.join(logDir, f'Fail.{track}.{o1}.{frame}.{runStamp}')
    log(f'==== runS1interferogram {ref["stem"]} / {sec["stem"]} ====')
    log(f'track {track} frame {frame}: reference {o1} ({ref["date"].date()}), '
        f'secondary {o2} ({sec["date"].date()})')
    log(f'host scratch {scratch}; cpus {cpus}; log {logPath}')

    # From here on a fatal signal is recorded rather than ending the run
    # silently, as happened to an overnight job that left no trace of why.
    installSignalHandlers()
    stallMinutes = args.stallMinutes if args.stallMinutes is not None \
        else project.get('stallMinutes', 30)
    pressureLimit = args.pressureLimit if args.pressureLimit is not None \
        else project.get('pressureLimit', 20)
    startWatchdog(logPath, failPath, pairScratch, stallMinutes,
                  pressureLimit=pressureLimit)
    governed = (f'; release scratch page cache above {pressureLimit}% '
                f'memory pressure' if pressureLimit else '')
    log(f'watchdog: abort after {stallMinutes} min without progress{governed}')

    tStart = time.time()
    try:
        if args.remapOnly:
            # Recover a run that got as far as the unwrap: reuse the scratch
            # ISCE dir and redo only convert_isce -> SETide -> setupUW.
            if not os.path.exists(f'{isceDir}/merged'):
                raise RuntimeError(f'--remapOnly: no unwrapped ISCE product '
                                   f'in {isceDir}; run the pair normally')
            log(f'--remapOnly: reusing {isceDir}, skipping topsApp and the '
                f'burst correction')
        else:
            # Fail now, not an hour in, if a precise orbit is missing
            checkOrbits(project, [ref, sec])
            os.makedirs(isceDir, exist_ok=True)
            # topsApp inputs: link SAFE zips, build the lat/lon DEM for the
            # frame footprint, and write topsApp.xml
            refZip = linkSafe(refPath, isceDir)
            secZip = linkSafe(secPath, isceDir)
            demPath = timed('build DEM', prepareIsceDem, project, isceDir,
                            frameBbox(refPath))
            writeTopsAppXML(project, demPath, refZip, secZip, isceDir)

            # Run topsApp through burst interferograms only; the single unwrap
            # happens in azPhaseCorrect after the burst correction.
            tStart = time.time()
            startDir = os.getcwd()
            os.chdir(isceDir)
            log(f'running topsApp in {isceDir} (--end=burstifg, {cpus} cpus)')
            timed('topsApp --end=burstifg', isceenv.isceShell,
                  'topsApp.py --end=burstifg', ompThreads=cpus,
                  ompPlaces='sockets(1)')
            os.chdir(startDir)

        # Burst correction, unwrap, and remap to the GrIMP product
        azThreads = args.azThreads if args.azThreads is not None \
            else project.get('azThreads', 3)
        productDir = timed('azPhaseCorrect (total)', apc.runAzPhaseCorrect,
                           isceDir, regionDefs, gimpDir, intermediatePath,
                           noAzCorrect=args.noAzCorrect,
                           noPhaseRemove=args.noPhaseRemove, cpus=cpus,
                           azThreads=azThreads,
                           gimpConvertOnly=args.remapOnly)

        # Copy the low-volume product out
        dest = timed('copy product', copyProduct, productDir,
                     project['outputRoot'], track)
        log(f'final product: {dest}')
        log(f'[TIME] TOTAL pipeline (excl. DEM): {time.time() - tStart:.1f} s')
    # SystemExit too: azPhaseCorrect / setupisceuw abort via u.myerror, which
    # calls sys.exit and would otherwise skip the fail record. KeyboardInterrupt
    # and KilledBySignal cover Ctrl-C and SIGHUP/SIGTERM.
    except (Exception, SystemExit, KeyboardInterrupt) as exc:
        code = exitCode(exc)
        log(f'FAILED (exit {code}): {failMessage(exc)}')
        log(traceback.format_exc())
        log(f'scratch kept for debugging: {pairScratch}')
        log(f'fail record: {writeFailFile(failPath, code, exc, logPath)}')
        sys.exit(code)

    # Clean up the large ISCE working files unless debugging
    if args.debug or args.keepScratch:
        log(f'--debug/--keepScratch: leaving scratch at {pairScratch}')
    else:
        shutil.rmtree(pairScratch)
        log(f'removed scratch {pairScratch}')


if __name__ == "__main__":
    main()
