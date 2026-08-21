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
their precise orbit will not have been downloaded yet.

With --queue the pairs are instead printed as ready-to-run runS1interferogram
commands, skipping any pair whose product already exists under
outputRoot/track-N/<orbit1>_<frame1>.

Config keys (uwproject.yaml): archiveDir, temporalBaseline, excludeMonths,
outputRoot (--queue only).
"""
import argparse
import datetime
import glob
import os
import re
import signal
import sys
import zipfile

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
    parser.add_argument('--noExcludeMonths', action='store_true',
                        default=False,
                        help='ignore the project excludeMonths list')
    parser.add_argument('--minAge', type=int, default=25,
                        help='skip acquisitions younger than this many days, '
                        'so their precise orbit has been downloaded [25]')
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


def loadAsfFrameLookup(project):
    ''' {track: {burstFrame: ASFframe}} from the asfFrameLookup file, or None.
    Maps our burst frame to ASF's frame number (their internal scheme). '''
    f = project.get('asfFrameLookup')
    if not f or not os.path.exists(f):
        return None
    with open(f) as fp:
        return yaml.safe_load(fp)


def parseName(path):
    ''' sat, absolute orbit, acquisition-start datetime, and stem from a SAFE
    zip name (.zip or .zip.1). '''
    stem = os.path.basename(path).split('.zip')[0]
    f = stem.split('_')
    sat = f[0]
    absOrbit = int(f[7])
    start = datetime.datetime.strptime(f[5], '%Y%m%dT%H%M%S')
    return sat, absOrbit, start, stem


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


def alreadyDone(frames, orbit, frame, tolerance):
    ''' True if a product exists for this reference orbit at this frame. The
    product frame comes from the merged sensingStart while ours comes from the
    SAFE-name start, so allow the same burst slop used to canonicalize frames.
    '''
    return any(abs(f - frame) <= tolerance
               for f in frames.get(orbit, []))


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

    # Gather SAFE zips from the YYYY-MM subdirs
    zips = glob.glob(f'{archiveDir}/[0-9][0-9][0-9][0-9]-[0-9][0-9]/'
                     f'S1*_IW_SLC__*.zip*')

    # Parse names, filter by track / date / month, dedup .zip vs .zip.1 by stem
    scenes = {}
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
        scenes.setdefault(stem, {'sat': sat, 'orb': orb, 'start': start,
                                 'track': track, 'path': path})

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

    # Output, sorted by track / frame / date.
    asfLookup = loadAsfFrameLookup(project)
    pairs.sort(key=lambda p: (p[0]['track'], p[0]['cframe'], p[0]['start']))

    def asfOf(track, frame):
        return (asfLookup.get(track, {}).get(frame, '-')
                if asfLookup is not None else '-')

    if args.queue:
        # runS1interferogram commands for the pairs with no product yet
        outputRoot = project.get('outputRoot')
        if not outputRoot:
            u.myerror('--queue: project has no outputRoot to check against')
        projectPath = os.path.abspath(args.project)
        existing, nSkipped = {}, 0
        for a, b in pairs:
            if a['track'] not in existing:
                existing[a['track']] = existingFrames(outputRoot, a['track'])
            if alreadyDone(existing[a['track']], a['orb'], a['frame'],
                           args.burstTolerance):
                nSkipped += 1
                continue
            print(f"runS1interferogram {a['path']} {b['path']} "
                  f"--project {projectPath}")
        # to stderr so the command list on stdout stays pipeable
        print(f'# {len(pairs) - nSkipped} to run, {nSkipped} already done '
              f'under {outputRoot}', file=sys.stderr)
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
