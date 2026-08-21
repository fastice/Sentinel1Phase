#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
setupS1PhaseTracks - instantiate per-track tie_plan_header files for a
Sentinel-1 unwrapped-interferogram (phase) project.

Phase analogue of s1setup.setupS1Tracks, pared down to the one thing the phase
workflow needs: --copyFiles, which drops a tie_plan_header (with <TRACK>
substituted) into each track's tiepoints/ directory from templates/. The phase
project has no velocityStats / vel_thumb_header machinery, so none of that is
carried over here.

Run from the project root (the directory holding uwproject.yaml, templates/, and
the track-* subdirectories), or point --project at a uwproject.yaml elsewhere --
its directory is taken as the project root.

    setupS1PhaseTracks --copyFiles
    setupS1PhaseTracks --copyFiles --tracks track-25 track-90
    setupS1PhaseTracks --copyFiles --overWrite
"""
import argparse
import datetime
import glob
import os
import re
import subprocess

import yaml
import utilities as u
import uwproject as uwp


def cmdLineParse():
    parser = argparse.ArgumentParser(
        description='Instantiate per-track tie_plan_header files from '
                    'templates/ for a Sentinel-1 phase project.',
        epilog='Part of the Sentinel1Phase package.')
    parser.add_argument('--project', default='uwproject.yaml',
                        help='project yaml [uwproject.yaml]; its directory is '
                             'the project root (templates/, track-*)')
    parser.add_argument('--tracks', nargs='+', metavar='track-N',
                        help='restrict to these track dirs (e.g. --tracks '
                             'track-25 track-90); default: all track-* dirs')
    parser.add_argument('--copyFiles', action='store_true',
                        help='copy tie_plan_header (substituting <TRACK>, <DEM>, '
                             '<TIEFILE>) into each track tiepoints/, skipping '
                             'files that already exist')
    parser.add_argument('--runRefresh', action='store_true',
                        help='run refreshties.py -phase (maketies + makeframetie '
                             '-> tie_script) across the selected tracks and '
                             'return; honors --year, --nThreads, --overWrite')
    parser.add_argument('--year', type=int, nargs='+', metavar='YYYY',
                        help='years to refresh (default: 2015 through the '
                             'current year)')
    parser.add_argument('--nThreads', type=int, default=4, metavar='N',
                        help='concurrent tracks passed to refreshties -nThreads '
                             '[4]')
    parser.add_argument('--overWrite', action='store_true',
                        help='with --copyFiles, overwrite an existing '
                             'tie_plan_header; with --runRefresh, pass '
                             '--overWrite to refreshties (rerun existing tie products)')
    return parser.parse_args()


def defaultYears():
    return list(range(2015, datetime.datetime.now().year + 1))


def runRefreshTies(projectDir, trackDirs, years, overWrite=False, nThreads=4):
    ''' Run refreshties.py in -phase mode across the selected tracks, from the
    project root. Phase (unwrapped-interferogram) flavor of
    s1setup.setupS1Tracks.runRefreshTies -- the -phase flag routes maketies to
    the phase products. Aborts if refreshties fails. '''
    tracksStr = str([os.path.basename(d) for d in trackDirs]).replace(' ', '')
    overWriteFlag = '--overWrite ' if overWrite else ''
    yearsStr = ' '.join(str(y) for y in years)
    cmd = (f'refreshties.py -phase {overWriteFlag}-nThreads {nThreads} '
           f'-toRun="{tracksStr}" {yearsStr} -noPrompt')
    print(f'Running: {cmd}  (in {projectDir})')
    result = subprocess.run(['csh', '-c', cmd], cwd=projectDir)
    if result.returncode != 0:
        u.myerror(f'refreshties failed (exit {result.returncode}) in {projectDir}')


def getTrackDirs(projectDir, tracks=None):
    ''' sorted track-* dirs under projectDir, numeric on the track integer; an
    explicit --tracks list is validated to exist. '''
    if tracks:
        dirs = [os.path.join(projectDir, t) for t in tracks]
        missing = [d for d in dirs if not os.path.isdir(d)]
        if missing:
            u.myerror(f'Track directories not found: {missing}')
    else:
        dirs = glob.glob(os.path.join(projectDir, 'track-*'))
    return sorted(dirs,
                  key=lambda p: int(re.search(r'track-(\d+)', p).group(1)))


def loadTemplate(projectDir):
    ''' contents of templates/tie_plan_header, or None if absent. '''
    path = os.path.join(projectDir, 'templates', 'tie_plan_header')
    if not os.path.exists(path):
        return None
    with open(path) as fp:
        return fp.read()


def applySubstitutions(content, trackNum, dem, tiepointFile):
    ''' same <TRACK>/<DEM>/<TIEFILE> convention as s1setup.setupS1Tracks. '''
    return (content.replace('<TRACK>', trackNum)
                   .replace('<DEM>', dem)
                   .replace('<TIEFILE>', tiepointFile))


def regionDem(project):
    ''' DEM for the <DEM> template substitution: the `dem` from the project's
    region yaml (region / regionFile key), falling back to the project demTiff.
    Only the DEM is taken from the region file -- runS1interferogram keeps using
    the self-contained demTiff/velMap/icemask keys. '''
    regionPath = project.get('region') or project.get('regionFile')
    if regionPath and os.path.exists(regionPath):
        with open(regionPath) as fp:
            region = yaml.safe_load(fp) or {}
        if region.get('dem'):
            return region['dem']
    return project.get('demTiff', '')


def setupTrackDirs(projectDir, trackDirs, project, overwrite=False):
    ''' Ensure each track's tiepoints/ exists and holds a tie_plan_header
    instantiated from templates/tie_plan_header. Copy-if-absent unless
    overwrite. <DEM> is filled from the region yaml's `dem` (region/regionFile
    key; falls back to demTiff); <TIEFILE> from tiepointFile. '''
    template = loadTemplate(projectDir)
    if template is None:
        u.myerror(f'no templates/tie_plan_header under {projectDir}')
    dem = regionDem(project)
    tiepointFile = project.get('tiepointFile', '')

    nCreated, nSkipped = 0, 0
    for trackDir in trackDirs:
        trackNum = os.path.basename(trackDir).split('-')[1]
        tpdir = os.path.join(trackDir, 'tiepoints')
        if not os.path.isdir(tpdir):
            os.makedirs(tpdir)
            print(f'Created {tpdir}')

        dest = os.path.join(tpdir, 'tie_plan_header')
        exists = os.path.exists(dest)
        if exists and not overwrite:
            nSkipped += 1
            continue
        with open(dest, 'w') as fp:
            fp.write(applySubstitutions(template, trackNum, dem, tiepointFile))
        print(f'{"Overwrote" if exists else "Created"} {dest}')
        nCreated += 1

    if nSkipped:
        print(f'Skipped {nSkipped} existing (use --overWrite to replace), '
              f'wrote {nCreated}')


def main():
    args = cmdLineParse()
    projectDir = os.path.dirname(os.path.abspath(args.project))
    project = uwp.loadProject(args.project)
    trackDirs = getTrackDirs(projectDir, args.tracks)

    if args.copyFiles:
        setupTrackDirs(projectDir, trackDirs, project, overwrite=args.overWrite)
        return

    if args.runRefresh:
        years = args.year if args.year else defaultYears()
        runRefreshTies(projectDir, trackDirs, years, overWrite=args.overWrite,
                       nThreads=args.nThreads)
        return

    u.myerror('nothing to do: pass --copyFiles or --runRefresh')


if __name__ == '__main__':
    main()
