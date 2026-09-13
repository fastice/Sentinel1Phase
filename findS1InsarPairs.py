#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
findS1InsarPairs - scan the Sentinel-1 archive and list interferometric pairs
that match the project temporal baseline.

    findS1InsarPairs TRACK [TRACK ...] [--firstDate YYYY-MM-DD] [--lastDate ...]
    findS1InsarPairs all
    findS1InsarPairs 112 --firstDate 2025-01-01 --lastDate 2025-12-31

For each pair it prints, space separated:

    track  orbit1  orbit2  frame  firstImage.zip  secondImage.zip

where frame is the burst-derived GrIMP frame (same numbering the products use),
orbit1/firstImage are the earlier acquisition. Pairs whose separation equals
`temporalBaseline` days are listed; acquisitions (either image) in
`excludeMonths` are dropped (unless --noExcludeMonths), as are those outside
[--firstDate, --lastDate].
Acquisitions younger than --minAge days (default 25) are also dropped, since
their precise orbit will not have been downloaded yet. --maxPairs N trims the
sorted list to the first N pairs, for a quick look at a large scan.

With --queue the pairs are instead printed as ready-to-run runS1interferogram
commands, skipping any pair whose product already exists under
outputRoot/track-N/<orbit1>_<frame1>.

Config keys (uwproject.yaml): archiveDir, temporalBaseline, excludeMonths,
maximumDT (per-frame baseline limits), outputRoot (--queue only),
searchArea / searchRegion (--onlinePairs search extent).
"""
import argparse
import datetime
import glob
import os
import re
import signal
import sys
import zipfile
from subprocess import call

import yaml
import utilities as u
import uwproject as uwp
from asfsearchdownload.fileS1 import computeTrack

# Burst-time -> frame constant, matching isce2grimp.convert_isce.get_frame_number
BTIME = 2.759


def cmdLineParse():
    parser = argparse.ArgumentParser(
        description='List Sentinel-1 interferometric pairs at the project '
                    'temporal baseline.')
    parser.add_argument('tracks', nargs='+',
                        help="track number(s), or 'all' for every track")
    parser.add_argument('--project', type=str, default='uwproject.yaml',
                        help='project yaml [uwproject.yaml]')
    parser.add_argument('--firstDate', type=str, default=None,
                        help='earliest acquisition date YYYY-MM-DD')
    parser.add_argument('--lastDate', type=str, default=None,
                        help='latest acquisition date YYYY-MM-DD')
    parser.add_argument('--excludedPairs', action='store_true', default=False,
                        help='inverse: list only pairs whose frame is OUTSIDE '
                        'the allowed phasePairRanges')
    parser.add_argument('--summary', action='store_true', default=False,
                        help='print only the pair count per frame '
                        '(track burstFrame asfFrame count)')
    parser.add_argument('--queue', action='store_true', default=False,
                        help='print runS1interferogram commands for the pairs '
                        'that have no product under outputRoot yet')
    parser.add_argument('--onlinePairs', action='store_true', default=False,
                        help='list granules ASF holds that would complete a '
                        'pair with one already downloaded. Runs searchASF over '
                        'the region unless --asfMeta gives an earlier result')
    parser.add_argument('--asfMeta', type=str, default=None, metavar='FILE',
                        help='reuse granule metadata from an earlier `searchASF --metadata` '
                        '(the <output>.meta file), used by --onlinePairs')
    parser.add_argument('--searchArea', type=str, default=None, metavar='FILE',
                        help='polygon file (lon,lat or geojson/shapefile) '
                        'forwarded to searchASF --searchArea; overrides '
                        '--greenland/--antarctica')
    parser.add_argument('--greenland', action='store_true', default=False,
                        help='search the bundled Greenland polygon '
                        '(forwarded to searchASF)')
    parser.add_argument('--antarctica', action='store_true', default=False,
                        help='search Antarctica (forwarded to searchASF)')
    parser.add_argument('--noExcludeMonths', action='store_true',
                        default=False,
                        help='ignore the project excludeMonths list')
    parser.add_argument('--minAge', type=int, default=25,
                        help='skip acquisitions younger than this many days, '
                        'so their precise orbit has been downloaded [25]')
    parser.add_argument('--maxPairs', type=int, default=None, metavar='N',
                        help='print only the first N pairs (after sorting by '
                        'track/frame/date), for a quick look at a large scan. '
                        'With --queue, N counts the pairs that still need '
                        'running -- already-done pairs do not use up the '
                        'budget')
    parser.add_argument('--parRun', type=str, default=None, metavar='DIR',
                        help='par_run directory whose queue/running/done tasks are also '
                             'treated as already handled, so --queue does not re-emit them. '
                             'Default: <project dir>/par_run when it exists')
    parser.add_argument('--noParRun', action='store_true', default=False,
                        help='ignore par_run entirely; --queue then skips only pairs whose '
                             'product already exists')
    parser.add_argument('--burstTolerance', type=int, default=1,
                        help='treat frames whose first burst differs by up to '
                        'this many bursts as the same frame [1]')
    return parser.parse_args()


def loadPhasePairRanges(project):
    ''' {track: [[lo, hi], ...]} from the phasePairRanges file, or None if not
    configured / missing (then no frame-range filtering is applied). '''
    f = project.get('phasePairRanges')
    if not f or not os.path.exists(f):
        return None
    with open(f) as fp:
        return yaml.safe_load(fp)


def loadMaximumDT(project):
    ''' {track: [[firstFrame, lastFrame, maxDays], ...]} from the maximumDT
    file, or None if not configured / missing (then no baseline limit is
    applied). '''
    f = project.get('maximumDT')
    if not f or not os.path.exists(f):
        return None
    with open(f) as fp:
        return yaml.safe_load(fp)


def maximumDTfor(limits, track, frame):
    ''' Shortest maximum baseline (days) covering this frame, or None if the
    frame is unlimited. Ranges are inclusive and carry their own burst buffer,
    so the raw and canonical frame both fall inside. '''
    if limits is None:
        return None
    covering = [maxDays for lo, hi, maxDays in limits.get(track, [])
                if lo <= frame <= hi]
    return min(covering) if covering else None


def loadAsfFrameLookup(project):
    ''' {track: {burstFrame: ASFframe}} from the asfFrameLookup file, or None.
    Maps our burst frame to ASF's frame number (their internal scheme). '''
    f = project.get('asfFrameLookup')
    if not f or not os.path.exists(f):
        return None
    with open(f) as fp:
        return yaml.safe_load(fp)


def invertAsfFrames(asfLookup):
    ''' {track: {asfFrame: [burstFrame, ...]}} from the burstFrame -> asfFrame
    lookup. One ASF frame covers the two or three burst frames the first burst
    drifts between across orbits, so the burst frame comes back only to within
    that spread - which is the same slop --burstTolerance already allows. '''
    inverse = {}
    for track, frames in (asfLookup or {}).items():
        perTrack = inverse.setdefault(track, {})
        for burstFrame, asfFrame in frames.items():
            perTrack.setdefault(asfFrame, []).append(burstFrame)
    for perTrack in inverse.values():
        for burstFrames in perTrack.values():
            burstFrames.sort()
    return inverse


def burstFrameFor(inverse, track, asfFrame):
    ''' Representative burst frame for an ASF frame: the smallest of the ones
    it covers, matching the canonical-frame convention used elsewhere here.
    None when this ASF frame has never been downloaded, so has no entry. '''
    burstFrames = inverse.get(track, {}).get(asfFrame)
    return burstFrames[0] if burstFrames else None


def regionArgs(args, project):
    ''' searchASF region flags, in precedence order:

      1. the command line (--searchArea / --antarctica / --greenland)
      2. the project's searchArea (a polygon file) then searchRegion (a name)
      3. the project region name, when it says antarctic
      4. Greenland

    Returned as a list so it drops straight into the command. '''
    if args.searchArea:
        return ['--searchArea', args.searchArea]
    if args.antarctica:
        return ['--antarctica']
    if args.greenland:
        return ['--greenland']
    if project.get('searchArea'):
        return ['--searchArea', project['searchArea']]
    named = str(project.get('searchRegion', '')).strip().lower()
    if named:
        if named not in ('greenland', 'antarctica'):
            u.myerror(f'regionArgs: searchRegion in the project must be '
                      f'greenland or antarctica, not "{named}"; use '
                      f'searchArea for a polygon file')
        return [f'--{named}']
    return ['--antarctica'] if 'antarctic' in \
        str(project.get('region', '')).lower() else ['--greenland']


def archiveExcludedTracks(archiveDir):
    ''' Tracks listed in the archive's own autoupdate.yaml tracksToExclude,
    as a sorted list. The config that governs an archive lives in the archive,
    so this stays right without duplicating the list into uwproject.yaml.
    Never fatal - an archive with no autoupdate.yaml excludes nothing. '''
    try:
        with open(os.path.join(archiveDir, 'autoupdate.yaml')) as fp:
            raw = (yaml.safe_load(fp) or {}).get('tracksToExclude', '')
    except Exception:
        return []
    if isinstance(raw, str):
        raw = raw.split()
    tracks = []
    for value in raw or []:
        try:
            tracks.append(int(value))
        except (TypeError, ValueError):
            pass
    return sorted(tracks)


def runSearchAsf(args, project, firstDate, lastDate):
    ''' Run searchASF over the project region for this window and return the
    path to the metadata file it writes.

    Shelling out to searchASF rather than querying ASF here keeps one place
    doing the search, so the region outline and the exclusions stay consistent
    with what actually gets downloaded. Output lands beside the archive's other
    search results, so a repeat run can be fed back with --asfMeta instead of
    searching again. '''
    archiveDir = project['archiveDir']
    resultDir = os.path.join(archiveDir, 'searchResults')
    if not os.path.isdir(resultDir):
        resultDir = os.getcwd()
    output = os.path.join(resultDir, f'onlinePairs.{firstDate}_{lastDate}')
    command = ['searchASF', str(firstDate), str(lastDate), output,
               '--sensor', 'SENTINEL1', '--products', 'SLC'] \
        + regionArgs(args, project) \
        + ['--metadata', '--archiveDir', f'{archiveDir}/*']
    # Tracks the archive deliberately does not hold (their zips were deleted)
    # would otherwise be quietly re-downloaded by the queue. Excluding at
    # search time keeps them out of the metadata entirely.
    exclude = archiveExcludedTracks(archiveDir)
    if exclude:
        command += ['--excludeTracks', ' '.join(str(t) for t in exclude)]
    print(f'# {" ".join(command)}', file=sys.stderr)
    # searchASF's own progress (Authenticated/Archive/Found/Volume/Metadata)
    # goes to stderr with everything else informational, so a redirected stdout
    # holds nothing but the runnable command lines.
    status = call(command, stdout=sys.stderr)
    if status != 0:
        u.myerror(f'runSearchAsf: searchASF failed (exit {status})')
    return output + '.meta'


def readAsfMetadata(metaFile):
    ''' Granule records from a `searchASF --metadata` file: one line of
    `granule track frame status sizeBytes url`, covering everything the search
    matched whether or not it is archived. Deduped by granule.

    Using searchASF rather than querying ASF here keeps one place doing the
    search, so the region outline, track and frame exclusions all stay
    consistent with what actually gets downloaded. '''
    if not os.path.exists(metaFile):
        u.myerror(f'readAsfMetadata: no such file: {metaFile}\n'
                  f'\tgenerate it with: searchASF ... --metadata')
    records = {}
    with open(metaFile) as fp:
        for line in fp:
            if line.startswith('#') or not line.strip():
                continue
            pieces = line.split()
            if len(pieces) < 6:
                continue
            granule, track, frame = pieces[0], pieces[1], pieces[2]
            if track == 'None' or frame == 'None':
                continue
            try:
                size = int(pieces[4])
            except ValueError:
                size = 0
            records[granule] = {'granule': granule, 'track': int(track),
                                'frame': int(frame), 'status': pieces[3],
                                'sizeBytes': size, 'url': pieces[5]}
    return list(records.values())


def archiveStems(archiveDir):
    ''' {stem: path} for every SAFE already on disk, from the YYYY-MM subdirs.
    Stems drop the .zip / .zip.1 suffix so they compare directly with ASF scene
    names, while the path preserves which suffix it actually has - a queue line
    naming an on-disk partner should name the file that is really there. '''
    return {os.path.basename(p).split('.zip')[0]: p
            for p in glob.glob(f'{archiveDir}/[0-9][0-9][0-9][0-9]-[0-9][0-9]/'
                               f'S1*_IW_SLC__*.zip*')}


def predictedPath(archiveDir, stem):
    ''' Where pullASF will put a granule that is not on disk yet:
    archiveDir/<YYYY-MM>/<stem>.zip, using the same first-date-token rule as
    autoupdateS1.monthDirFor. resolveSafe tolerates either suffix, so it still
    resolves if a later filing run renames it to .zip.1. '''
    match = re.search(r'(\d{8})T\d{6}', stem)
    if match is None:
        return None
    ymd = match.group(1)
    return os.path.join(archiveDir, f'{ymd[:4]}-{ymd[4:6]}', f'{stem}.zip')


def parseName(path):
    ''' sat, absolute orbit, acquisition-start datetime, and stem from a SAFE
    zip name (.zip or .zip.1). '''
    stem = os.path.basename(path).split('.zip')[0]
    f = stem.split('_')
    sat = f[0]
    absOrbit = int(f[7])
    start = datetime.datetime.strptime(f[5], '%Y%m%dT%H%M%S')
    return sat, absOrbit, start, stem


def polCode(stem):
    ''' The 1SDH/1SDV/1SSH/1SSV mode-and-polarisation token from a SAFE name,
    or None. Matched by pattern, not position: the name has a double
    underscore after SLC, so the token is field 4 and a positional split on
    field 3 silently yields an empty string. '''
    for field in stem.split('_'):
        if re.fullmatch(r'1S[SD][HV]', field):
            return field
    return None


def hasPolarization(stem, polarization):
    ''' Does this SAFE carry the polarization the project processes with?

    The 4th name field is the mode/polarisation code -- 1SDH (HH+HV), 1SDV
    (VV+VH), 1SSH (HH only), 1SSV (VV only). Its last character is the
    transmit polarisation and D/S says whether the cross-pol channel is there.

    Worth filtering here because the failure downstream is unrecognisable: a
    granule without the configured channel has no matching annotation xml, so
    topsApp extracts no swaths, says only "Could not extract swath N", carries
    on, and dies two steps later in runTopo with "IndexError: too many indices
    for array" from an empty bounding-box list. Six track-17 frame-396 pairs
    burned a DEM build and a topsApp startup each that way. '''
    pol = (polarization or 'hh').lower()
    code = polCode(stem)
    if code is None or len(pol) != 2:
        return True                       # unrecognised: leave it alone
    if code[3].upper() != pol[0].upper():
        return False                      # wrong transmit polarisation
    return pol[0] == pol[1] or code[2].upper() == 'D'   # cross-pol needs dual


def projectPolarization(project):
    ''' Polarisation topsApp will be run with (isce.polarization, else hh). '''
    return project.get('isce', {}).get('polarization', 'hh')


def ascNodeTime(path):
    ''' Ascending-node time from the SAFE manifest (one per absolute orbit). '''
    with zipfile.ZipFile(path) as z:
        man = [n for n in z.namelist() if n.endswith('manifest.safe')][0]
        text = z.read(man).decode('utf-8', errors='ignore')
    m = re.search(r'ascendingNodeTime>([^<]+)<', text)
    return datetime.datetime.strptime(m.group(1)[:26], '%Y-%m-%dT%H:%M:%S.%f')


def frameNumber(start, asc):
    ''' Burst-derived GrIMP frame number (matches convert_isce). '''
    return int((start - asc).seconds / BTIME + 0.5)


def parseDate(s):
    return datetime.datetime.strptime(s, '%Y-%m-%d').date() if s else None


def existingFrames(outputRoot, track):
    ''' {orbit1: [frame1, ...]} for the products already under
    outputRoot/track-N/, read from their <orbit1>_<frame1> dir names
    (setupisceuw.setupOutputDir). '''
    frames = {}
    for d in glob.glob(f'{outputRoot}/track-{track}/[0-9]*_[0-9]*'):
        if not os.path.isdir(d):
            continue
        try:
            orbit, frame = [int(v) for v in os.path.basename(d).split('_')]
        except ValueError:
            continue
        frames.setdefault(orbit, []).append(frame)
    return frames


def resolveParRun(args):
    """par_run directory to honour, or None.  Defaults to the one beside the project file."""
    if args.noParRun:
        return None
    if args.parRun:
        return args.parRun
    d = os.path.join(os.path.dirname(os.path.abspath(args.project)), 'par_run')
    return d if os.path.isdir(d) else None


def parRunTasks(parRunDir, states=('queue', 'running', 'done')):
    """(orbit1, orbit2, frame) for every task already sitting in par_run.

    A product existing under outputRoot is not the only reason to skip a pair: it may be
    queued, running, or finished as a par_run task without its product being visible yet
    (still running, or failed after the task was filed).  Re-emitting those duplicates work
    and, for a queued-but-not-started pair, produces two workers racing on the same frame.

    Each task file is a csh script whose echo trailer carries the same identifiers the queue
    itself emits: "track=N frame=F orbits=O1-O2".  Parsing that rather than the SAFE paths
    keeps this tied to what the emitter writes.
    """
    tasks = []
    if not parRunDir or not os.path.isdir(parRunDir):
        return tasks
    for state in states:
        for f in glob.glob(os.path.join(parRunDir, state, 'task*')):
            try:
                txt = open(f).read()
            except OSError:
                continue
            orb = re.search(r'orbits=(\d+)-(\d+)', txt)
            if not orb:
                continue
            fr = re.search(r'frame=(\d+)', txt)
            tasks.append((int(orb.group(1)), int(orb.group(2)),
                          int(fr.group(1)) if fr else None))
    return tasks


def alreadyQueued(tasks, orbit1, orbit2, frame, tolerance):
    """True if this pair is already a par_run task.  Frame is compared with the same burst
    slop as alreadyDone; a task or candidate with no frame matches on the orbit pair alone,
    which errs toward skipping rather than duplicating."""
    for a, b, f in tasks:
        if a != orbit1 or b != orbit2:
            continue
        if f is None or frame is None or abs(f - frame) <= tolerance:
            return True
    return False


def alreadyDone(frames, orbit, frame, tolerance):
    ''' True if a product exists for this reference orbit at this frame. The
    product frame comes from the merged sensingStart while ours comes from the
    SAFE-name start, so allow the same burst slop used to canonicalize frames.
    '''
    return any(abs(f - frame) <= tolerance
               for f in frames.get(orbit, []))


def onlinePairs(args, project, firstDate, lastDate, excludeMonths, trackSet):
    ''' List granules ASF holds that would complete a pair with something
    already downloaded.

    Pairing is done on the ASF frame rather than our burst frame: ASF's frame
    is stable across repeats (burst frames 390 and 391 of one ground frame are
    both ASF 215), so it needs no ascending-node time - which is only available
    inside the zip, and so cannot be had for a granule not yet downloaded. The
    burst frame is recovered by inverting asfFrameLookup, to within the burst
    drift.

    The search is geographic - the whole Greenland outline, every track - so
    frames and tracks never downloaded before are found too. The burst frame is
    reported where asfFrameLookup knows it and as '-' otherwise; the filters
    keyed on burst frames (phasePairRanges, maximumDT) are only applied when it
    is known, so an unseen frame is listed rather than silently dropped. '''
    archiveDir = project['archiveDir']
    baseline = int(project['temporalBaseline'])
    inverse = invertAsfFrames(loadAsfFrameLookup(project))
    if firstDate is None:
        firstDate = datetime.date(2014, 1, 1)
    # --asfMeta reuses an earlier search; without it, run one now
    metaFile = args.asfMeta or runSearchAsf(args, project, firstDate,
                                            lastDate)
    found = readAsfMetadata(metaFile)
    print(f'# {len(found)} granules in {metaFile}', file=sys.stderr)

    onDisk = archiveStems(archiveDir)
    limits = loadMaximumDT(project)
    ranges = loadPhasePairRanges(project)
    # Group by (track, ASF frame): the ASF frame is repeat stable, so it does
    # the job our burst-frame clustering does, without an ascending node time.
    byFrame, nWrongPol = {}, 0
    polarization = projectPolarization(project)
    for record in found:
        stem = record['granule']
        if trackSet and record['track'] not in trackSet:
            continue
        try:
            # everything else needed is in the granule name
            sat, orbit, start, stem = parseName(stem)
        except (ValueError, IndexError):
            continue
        if start.month in excludeMonths:
            continue
        if not (firstDate <= start.date() <= lastDate):
            continue
        if not hasPolarization(stem, polarization):
            nWrongPol += 1
            continue
        byFrame.setdefault((record['track'], record['frame']), []).append(
            {'stem': stem, 'start': start, 'orbit': orbit,
             'url': record['url'], 'sizeBytes': record.get('sizeBytes', 0),
             # mission and polarisation token both have to match: one relative
             # orbit carries 1SDH here and 1SDV elsewhere
             'kind': (sat, polCode(stem)),
             'have': stem in onDisk, 'path': onDisk.get(stem)})

    nBothMissing, nUnknownFrame, rows = 0, 0, []
    for (track, asfFrame), members in sorted(byFrame.items()):
        burstFrame = burstFrameFor(inverse, track, asfFrame)
        if burstFrame is None:
            nUnknownFrame += 1
        else:
            maxDT = maximumDTfor(limits, track, burstFrame)
            if maxDT is not None and baseline > maxDT:
                continue
            if ranges is not None and ranges.get(track) and not any(
                    lo <= burstFrame <= hi for lo, hi in ranges[track]):
                continue
        members.sort(key=lambda m: m['start'])
        for i, first in enumerate(members):
            for second in members[i + 1:]:
                days = (second['start'].date() - first['start'].date()).days
                if days > baseline:
                    break
                if days != baseline or first['kind'] != second['kind']:
                    continue
                if first['have'] and second['have']:
                    continue
                if not first['have'] and not second['have']:
                    nBothMissing += 1
                # first/second are already in acquisition order, which is the
                # order runS1interferogram wants (reference then secondary).
                rows.append((track, asfFrame, burstFrame, first, second))

    rows.sort(key=lambda r: (r[0], r[1], r[3]['start']))
    toPull = {m['stem']: m for _, _, _, a, b in rows for m in (a, b)
              if not m['have']}
    volume = sum(m['sizeBytes'] for m in toPull.values())
    if args.queue:
        emitOnlineQueue(args, project, rows)
    else:
        print('# track asfFrame burstFrame refDate refGranule secDate '
              'secGranule missing url(s)', file=sys.stderr)
        for track, asfFrame, burstFrame, first, second in rows:
            urls = ' '.join(m['url'] for m in (first, second)
                            if not m['have'])
            missingCount = sum(1 for m in (first, second) if not m['have'])
            print(f"{track} {asfFrame} {burstFrame if burstFrame else '-'} "
                  f"{first['start'].date()} {first['stem']}.zip "
                  f"{second['start'].date()} {second['stem']}.zip "
                  f"{missingCount} {urls}")
    wrongPol = (f'; skipped {nWrongPol} granules without {polarization}'
                if nWrongPol else '')
    print(f'# {len(rows)} pairs need a download ({len(rows) - nBothMissing} '
          f'complete one already on disk, {nBothMissing} need both); '
          f'{len(toPull)} distinct granules, {volume / 1e12:.2f} TB to fetch; '
          f'{nUnknownFrame} frames had no asfFrameLookup entry{wrongPol}',
          file=sys.stderr)


def emitOnlineQueue(args, project, rows):
    ''' csh lines that fetch what each pair needs and then run it.

    One pullASF call per pair, not one per granule: with a single download slot
    chaining two calls interleaves at the granule level, so a job can win the
    slot for its first image early and its second last and nothing starts until
    nearly every download is done. Taking the slot once per pair means each job
    leaves with something runnable. '''
    archiveDir = project['archiveDir']
    outputRoot = project.get('outputRoot')
    if not outputRoot:
        u.myerror('--onlinePairs --queue: project has no outputRoot')
    projectPath = os.path.abspath(args.project)
    statusLog = os.path.join(outputRoot, 'logs', 'exitStatus.log')
    parRunDir = resolveParRun(args)
    tasks = parRunTasks(parRunDir)
    existing, nPrinted, nSkipped, nQueued = {}, 0, 0, 0
    for track, asfFrame, burstFrame, first, second in rows:
        if track not in existing:
            existing[track] = existingFrames(outputRoot, track)
        if alreadyQueued(tasks, first['orbit'], second['orbit'],
                         burstFrame if burstFrame else asfFrame,
                         args.burstTolerance):
            nQueued += 1
            continue
        # Skip pairs already produced. Impossible when the burst frame is
        # unknown, so those are emitted rather than silently dropped.
        if burstFrame is not None and alreadyDone(existing[track],
                                                  first['orbit'], burstFrame,
                                                  args.burstTolerance):
            nSkipped += 1
            continue
        if args.maxPairs is not None and nPrinted >= args.maxPairs:
            continue
        # Both urls, even for a granule already on disk: it fast-paths in
        # milliseconds and makes the line self-repairing if that copy is moved.
        urls = ' '.join(m['url'] for m in (first, second))
        paths = [m['path'] if m['have']
                 else predictedPath(archiveDir, m['stem'])
                 for m in (first, second)]
        print(f"pullASF {urls} --archiveDir {archiveDir} && "
              f"runS1interferogram {paths[0]} {paths[1]} "
              f"--project {projectPath}; set rc = $status; "
              f'echo "`date +%FT%T` rc=$rc track={track} '
              f'frame={burstFrame if burstFrame else asfFrame} '
              f"orbits={first['orbit']}-{second['orbit']}\" >>! {statusLog}")
        nPrinted += 1
    if nSkipped:
        print(f'# skipped {nSkipped} pairs whose product already exists',
              file=sys.stderr)
    if nQueued:
        print(f'# skipped {nQueued} pairs already in {parRunDir} '
              f'(queue/running/done)', file=sys.stderr)


def main():
    # Exit quietly when the reader closes the pipe (`... | head`, quitting
    # less) instead of raising BrokenPipeError out of the print loop.
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    args = cmdLineParse()
    project = uwp.loadProject(args.project)
    archiveDir = project['archiveDir']
    baseline = int(project['temporalBaseline'])
    excludeMonths = set() if args.noExcludeMonths \
        else set(project.get('excludeMonths', []))
    firstDate = parseDate(args.firstDate)
    lastDate = parseDate(args.lastDate)
    # Precise orbits (POEORB) lag the acquisition by ~3 weeks, so an image that
    # is too young cannot be processed yet. Tighten lastDate rather than
    # filtering separately - both images of a pair are then covered.
    if args.minAge > 0:
        cutoff = (datetime.datetime.now()
                  - datetime.timedelta(days=args.minAge)).date()
        lastDate = min(lastDate, cutoff) if lastDate is not None else cutoff
    allTracks = any(t.lower() == 'all' for t in args.tracks)
    trackSet = set() if allTracks else {int(t) for t in args.tracks}

    # Works from a searchASF metadata file, not the local archive, so it runs
    # before the archive scan below (which opens zips for ascending node times
    # this mode does not need).
    if args.onlinePairs:
        onlinePairs(args, project, firstDate, lastDate, excludeMonths,
                    trackSet)
        return

    # Gather SAFE zips from the YYYY-MM subdirs
    zips = glob.glob(f'{archiveDir}/[0-9][0-9][0-9][0-9]-[0-9][0-9]/'
                     f'S1*_IW_SLC__*.zip*')

    # Parse names, filter by track / date / month / polarisation, dedup
    # .zip vs .zip.1 by stem
    polarization = projectPolarization(project)
    scenes, nWrongPol = {}, 0
    for path in zips:
        try:
            sat, orb, start, stem = parseName(path)
        except (ValueError, IndexError):
            continue
        track = computeTrack(orb, sat, start)
        if not allTracks and track not in trackSet:
            continue
        if firstDate and start.date() < firstDate:
            continue
        if lastDate and start.date() > lastDate:
            continue
        if start.month in excludeMonths:
            continue
        if not hasPolarization(stem, polarization):
            nWrongPol += 1
            continue
        scenes.setdefault(stem, {'sat': sat, 'orb': orb, 'start': start,
                                 'track': track, 'path': path})

    if nWrongPol:
        print(f'# skipped {nWrongPol} SAFEs carrying no {polarization} '
              f'channel (topsApp would fail obscurely on them)',
              file=sys.stderr)

    # Frame number (ascending-node time cached per absolute orbit)
    ascCache = {}
    for s in scenes.values():
        if s['orb'] not in ascCache:
            ascCache[s['orb']] = ascNodeTime(s['path'])
        s['frame'] = frameNumber(s['start'], ascCache[s['orb']])

    # Canonicalize frames per track: cluster first-burst frames whose
    # consecutive values differ by <= burstTolerance into one frame (labeled by
    # the cluster's smallest first burst), so acquisitions offset by a burst
    # (or a chained few) are treated as the same frame.
    tol = args.burstTolerance
    perTrack = {}
    for s in scenes.values():
        perTrack.setdefault(s['track'], []).append(s)
    for sc in perTrack.values():
        f0s = sorted(set(s['frame'] for s in sc))
        canon, cluster = {}, [f0s[0]]
        for f in f0s[1:]:
            if f - cluster[-1] <= tol:
                cluster.append(f)
            else:
                canon.update({m: cluster[0] for m in cluster})
                cluster = [f]
        canon.update({m: cluster[0] for m in cluster})
        for s in sc:
            s['cframe'] = canon[s['frame']]

    # Group by (track, canonical frame) and pull out baseline-separated pairs
    groups = {}
    for s in scenes.values():
        groups.setdefault((s['track'], s['cframe']), []).append(s)
    pairs = []
    for members in groups.values():
        members.sort(key=lambda x: x['start'])
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                dd = (members[j]['start'].date()
                      - members[i]['start'].date()).days
                if dd == baseline:
                    # same mission only (S1A-S1A, S1C-S1C, ...): a 12-day
                    # cross-mission pair is not a valid interferogram
                    if members[i]['sat'] == members[j]['sat']:
                        pairs.append((members[i], members[j]))
                elif dd > baseline:
                    break

    # Keep pairs whose frame falls in the track's allowed ranges (or, with
    # --excludedPairs, only those that fall outside).
    ranges = loadPhasePairRanges(project)
    if ranges is not None:
        def inRange(track, frame):
            return any(lo <= frame <= hi
                       for lo, hi in ranges.get(track, []))
        pairs = [p for p in pairs
                 if inRange(p[0]['track'],
                            p[0]['cframe']) != args.excludedPairs]

    # Drop frames that only fit at a shorter baseline than this project uses;
    # they come back on their own once temporalBaseline is lowered.
    limits = loadMaximumDT(project)
    if limits is not None:
        kept = [p for p in pairs
                if (maximumDTfor(limits, p[0]['track'], p[0]['cframe']) or
                    baseline) >= baseline]
        if len(kept) < len(pairs):
            print(f'# dropped {len(pairs) - len(kept)} pairs at '
                  f'{baseline} days exceeding their frame maximumDT',
                  file=sys.stderr)
        pairs = kept

    # Output, sorted by track / frame / date.
    asfLookup = loadAsfFrameLookup(project)
    pairs.sort(key=lambda p: (p[0]['track'], p[0]['cframe'], p[0]['start']))

    # Trim after sorting, so "the first N" is a stable, meaningful slice rather
    # than whatever the archive glob happened to return first. The note goes to
    # stderr so stdout stays pipeable and nobody mistakes a truncated list for
    # the whole scan. --queue applies its own limit further down, counting only
    # the pairs that still need running -- an N spent on already-done pairs
    # would hand back an empty command list.
    if args.maxPairs is not None and not args.queue \
            and len(pairs) > args.maxPairs:
        print(f'# showing first {args.maxPairs} of {len(pairs)} pairs '
              f'(--maxPairs)', file=sys.stderr)
        pairs = pairs[:args.maxPairs]

    def asfOf(track, frame):
        return (asfLookup.get(track, {}).get(frame, '-')
                if asfLookup is not None else '-')

    if args.queue:
        # runS1interferogram commands for the pairs with no product yet
        outputRoot = project.get('outputRoot')
        if not outputRoot:
            u.myerror('--queue: project has no outputRoot to check against')
        projectPath = os.path.abspath(args.project)
        # Each command records its own exit status. runS1interferogram cannot
        # do this for itself: a run killed by SIGKILL (137, OOM) or SIGSEGV
        # (139, a crash in GDAL/numpy) dies without running any of its own
        # code, which is exactly the case that has to be told apart.
        #
        # csh syntax, since that is what the queue is run under: $status not
        # $?, backticks not $(), and `>>!` so the append still works under
        # `set noclobber` when the file does not exist yet. $status is copied
        # into rc first - the date substitution would otherwise overwrite it.
        statusLog = os.path.join(outputRoot, 'logs', 'exitStatus.log')
        parRunDir = resolveParRun(args)
        tasks = parRunTasks(parRunDir)
        existing, nSkipped, nToRun, nPrinted, nQueued = {}, 0, 0, 0, 0
        for a, b in pairs:
            if a['track'] not in existing:
                existing[a['track']] = existingFrames(outputRoot, a['track'])
            if alreadyDone(existing[a['track']], a['orb'], a['frame'],
                           args.burstTolerance):
                nSkipped += 1
                continue
            if alreadyQueued(tasks, a['orb'], b['orb'], a['frame'],
                             args.burstTolerance):
                nQueued += 1
                continue
            # Keep scanning past the --maxPairs cut so the totals below still
            # describe the whole queue, not just the part that was printed.
            nToRun += 1
            if args.maxPairs is not None and nPrinted >= args.maxPairs:
                continue
            print(f"runS1interferogram {a['path']} {b['path']} "
                  f"--project {projectPath}; set rc = $status; "
                  f'echo "`date +%FT%T` rc=$rc track={a["track"]} '
                  f'frame={a["cframe"]} orbits={a["orb"]}-{b["orb"]}" '
                  f">>! {statusLog}")
            nPrinted += 1
        # to stderr so the command list on stdout stays pipeable
        if nPrinted < nToRun:
            print(f'# showing first {nPrinted} of {nToRun} to run '
                  f'(--maxPairs)', file=sys.stderr)
        print(f'# {nToRun} to run, {nSkipped} already done '
              f'under {outputRoot}'
              + (f', {nQueued} already in {parRunDir} '
                 f'(queue/running/done of {len(tasks)} tasks)'
                 if parRunDir else '; par_run not checked'), file=sys.stderr)
    elif args.summary:
        # per frame, count of pairs by year, in fixed per-year columns:
        #   track burstFrame asfFrame  count (year)  count (year) ...
        # a frame with no pairs in a year leaves that year's column blank.
        counts = {}
        for a, b in pairs:
            key = (a['track'], a['cframe'])
            byYear = counts.setdefault(key, {})
            y = a['start'].year
            byYear[y] = byYear.get(y, 0) + 1
        years = sorted({y for byYear in counts.values() for y in byYear})
        cw = len(f'{0:>3} ({years[0]})') if years else 0
        for track, frame in sorted(counts):
            byYear = counts[(track, frame)]
            cells = [f'{byYear[y]:>3} ({y})' if y in byYear else ' ' * cw
                     for y in years]
            print(f'{track:>3} {frame:>5} {str(asfOf(track, frame)):>5}  '
                  + '  '.join(cells))
    else:
        # track orbit1 orbit2 burstFrame asfFrame image1 image2
        for a, b in pairs:
            print(f"{a['track']} {a['orb']} {b['orb']} {a['cframe']} "
                  f"{asfOf(a['track'], a['cframe'])} "
                  f"{os.path.basename(a['path'])} {os.path.basename(b['path'])}")


if __name__ == "__main__":
    main()
