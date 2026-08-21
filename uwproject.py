#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Project-configuration loader for the Sentinel-1 unwrapped-interferogram
workflow (runS1interferogram).

A uwproject.yaml collects everything that used to be scattered across
isce2grimp/data/template.yml, sarfunc region yamls, and hardcoded strings:

  scratchDir  - per-machine local scratch (hostname -> path)
  threads     - optional per-machine topsApp cpu counts
  dataDir     - directory holding the input SAFE zips
  orbitDir    - S1 precise-orbit (OPOD) directory
  outputRoot  - where final low-volume GrIMP products are written
  dem         - ISCE-format DEM (.dem.wgs84) used by topsApp
  isce        - topsApp parameters (looks, polarization, ionosphere, ...)
  epsg/demTiff/velMap/icemask/icemaskTiff - region data for phase simulation
                and masking (self-contained), OR set regionFile: to reuse an
                existing sarfunc region yaml.
"""
import os
import socket
import yaml
import utilities as u
import sarfunc as sf


def loadProject(projectFile):
    ''' Load a uwproject.yaml and return it as a dict. '''
    if not os.path.exists(projectFile):
        u.myerror(f'loadProject: project file does not exist: {projectFile}')
    with open(projectFile) as fp:
        project = yaml.safe_load(fp)
    return project


def scratchForHost(project, host=None):
    ''' Return the local scratch dir for this host from the scratchDir map. '''
    if host is None:
        host = socket.gethostname()
    scratchMap = project.get('scratchDir', {})
    if host not in scratchMap:
        u.myerror(f'scratchForHost: no scratchDir entry for host "{host}"; '
                  f'add one to the project file (have: '
                  f'{list(scratchMap.keys())})')
    scratch = scratchMap[host]
    if not os.path.exists(scratch):
        u.myerror(f'scratchForHost: scratch dir does not exist: {scratch}')
    return scratch


def threadsForHost(project, host=None):
    ''' cpu count for this host: threads[host] else isce.cpus else 8. '''
    if host is None:
        host = socket.gethostname()
    default = project.get('isce', {}).get('cpus', 8)
    return project.get('threads', {}).get(host, default)


def regionDefsFromProject(project):
    ''' Build a sarfunc.defaultRegionDefs for phase-sim / masking. Uses
    regionFile if given, otherwise the self-contained keys in the project.
    Note: `region` may be a region-yaml path (setupS1PhaseTracks reads its DEM);
    here it is only used for the region name, so reduce a path to its stem. '''
    if project.get('regionFile'):
        return sf.defaultRegionDefs(None, regionFile=project['regionFile'])
    regionName = os.path.splitext(
        os.path.basename(str(project.get('region') or 'greenland')))[0]
    regionDef = {
        'name': regionName,
        'epsg': project['epsg'],
        'dem': project['demTiff'],
        'velMap': project['velMap'],
        'icemask': project['icemask'],
        'icemaskTiff': project.get('icemaskTiff'),
    }
    return sf.defaultRegionDefs(None, regionDef=regionDef)


def icemaskTiff(regionDefs):
    ''' tiff ice mask used for velocity masking. Falls back to the region dem
    directory convention if not explicitly set. '''
    return regionDefs.region.get('icemaskTiff')
