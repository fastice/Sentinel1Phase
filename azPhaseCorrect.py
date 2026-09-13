#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:12:15 2021

@author: ian
"""

import argparse
import time
import numpy as np
from bs4 import BeautifulSoup
import os
import utilities as u
import nisarfunc as nisar
import threading
import pyproj
import glob
import shutil
from datetime import datetime
import sarfunc as sf
from subprocess import call
from isceenv import isceShell
import rasterio
import rioxarray
from scipy.ndimage import gaussian_filter as gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import warnings
warnings.filterwarnings("ignore",
                        category=rasterio.errors.NotGeoreferencedWarning)


icemaskMap = {'greenland':
              '/Volumes/insar7/ian/gimp/mask/GimpIceMask_250m.tif',
              'antarctica': None}


def fileS1Args():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(description='\033[1m Apply burst by burst'
                                     ' phase correction to compensate for '
                                     'azimuth motion \033[0m',
                                     epilog='Notes:  ', allow_abbrev=False)
    parser.add_argument('isceDir', metavar='isceDir', type=str, nargs=1,
                        help='ISCE directory with product')
    parser.add_argument('--region', type=str,
                        choices=['greenland', 'amundsen', 'taku'],
                        default='greenland',
                        help='region code')
    parser.add_argument('--outputSuffix', type=str, default='',
                        help='suffix for gimpSuffix and intermediateSuffix '
                        'directories')
    parser.add_argument('--noAzCorrect', action='store_true', default=False,
                        help='Do not apply azimuth correction [False]')
    parser.add_argument('--noPhaseRemove', action='store_true', default=False,
                        help='Do not remove phase before unwrapping [False]')
    parser.add_argument('--gimpConvertOnly', action='store_true',
                        default=False,
                        help='Only do conversion from isce to gimp [False]')
    parser.add_argument('--unwrapOnly', action='store_true',
                        default=False,
                        help='Use prior merge and only unwrap & mv to gimp '
                        '[False]')
    parser.add_argument('--noRemap', action='store_true', default=False,
                        help='Do not remap to gimp formats [False]')
    args = parser.parse_args()
    if not os.path.exists(args.isceDir[0]):
        u.myerror(f'fileS1Args: basedir ({args.isceDir[0]}) '
                  'does not exist')
    return args.isceDir[0].rstrip('/'), args.region, args.noAzCorrect, \
        args.noPhaseRemove, not args.noRemap, args.gimpConvertOnly,\
        args.unwrapOnly, args.outputSuffix

# ---- IO Routines: Parameters


def logFail(message):
    with open(f'fail.{os.getpid()}', 'w') as fp:
        print(message, file=fp)


def burstSize(xmlFile):
    ''' Get burst size from xml'''
    xmlData = readXML(xmlFile)
    coord1 = xmlData.find('component', {'name': 'coordinate1'})
    nRange = int(coord1.find('property', {'name': 'size'}).value.text)
    coord2 = xmlData.find('component', {'name': 'coordinate2'})
    nAzimuth = int(coord2.find('property', {'name': 'size'}).value.text)
    return {'nRange': nRange, 'nAzimuth': nAzimuth}


def getDT(dataPath):
    ''' Get time interval between images as Secondary - Reference.
    Returns +dT for reference earlier images'''
    xmlData = readXML(f'{dataPath}/topsProc.xml')
    tr = datetime.strptime(
        xmlData.find('reference').find('ascendingnodetime').text.split()[0],
        '%Y-%m-%d')
    ts = datetime.strptime(
        xmlData.find('secondary').find('ascendingnodetime').text.split()[0],
        '%Y-%m-%d')
    return (ts-tr).days


def getBeamParams(dataPath, beams=['IW1', 'IW2', 'IW3']):
    ''' Extract the required params from the reference file beam xmls'''
    beamParams = {}
    params = ['pulserepetitionfrequency', 'azimuthsteeringrate',
              'numberoflines', 'numberofvalidsamples', 'radarwavelength']
    paramTypes = [float, float, int, int, float]
    for beam in beams:
        xmlData = readXML(f'{dataPath}/referencedir/{beam}.xml')
        results = {}
        for param, paramType in zip(params, paramTypes):
            results[param] = paramType(
                xmlData.find('property', {'name': param}).value.text)
        beamParams[beam] = results
    return beamParams


def readXML(xmlFile):
    ''' Open and existence check xml file and return a beautiful soup object'''
    if not os.path.exists(xmlFile):
        logFail(f'readXML: File does not exist: {xmlFile}')
        u.myerror(f'readXML: File does not exist: {xmlFile}')
        return None
    with open(xmlFile) as fpXML:
        xmlData = BeautifulSoup(fpXML, 'xml')
        return xmlData


# ---- IO Routines: Data


def getLatLonXY(baseName, regionDefs, nrlks=30, nalks=6):
    ''' Get lat/lon and convert to xy
    baseName should contain wildcard * in place of lat and lon'''
    latLonXY = {}
    for coord in ['lat', 'lon']:
        fileName = baseName.replace('*', coord)
        x = rasterio.open(fileName)
        # Downscale if needd
        lookedShape = int(x.shape[0]/nalks), int(x.shape[1]/nrlks)
        latLonXY[coord] = x.read(out_shape=lookedShape).reshape(lookedShape[0],
                                                                lookedShape[1])
    midPt = tuple((np.array(latLonXY['lat'].shape)/2).astype('int'))
    latLonXY['latC'] = latLonXY['lat'][midPt]
    latLonXY['epsg'] = regionDefs.epsg()
    if latLonXY['epsg'] is None:
        logFail('epsg is None, need to implement wkt')
        u.myerror('epsg is None, need to implement wkt')
    latLonXY['x'], latLonXY['y'] = lltoxy(latLonXY['lat'], latLonXY['lon'],
                                          epsg=latLonXY['epsg'])
    return latLonXY


def getLOS(fileName, nrlks=30, nalks=6):
    ''' get los from files and down scale if needed '''
    x = rasterio.open(fileName)
    lookedShape = (int(x.shape[0]/nalks), int(x.shape[1]/nrlks))
    los = x.read(out_shape=lookedShape)
    incidence = los[0]
    # ISCE is CCW angle from ground to sat
    # - to go to CW and 180 to look from sat to ground
    azAngle = -1 * los[1] + 180.
    return {'incidence': incidence, 'azAngle': azAngle}


_velCache = {}
_velCacheLock = threading.Lock()


def readVelCrop(velFile, llxy, pad=3000):
    ''' Crop every band of the velocity map to the llxy footprint in a single
    read, and return (data, x, y).

    The dataset handle is opened once and cached: with a handle per call the
    correction opened the map (and, for a VRT, its component tiffs) 27 times
    per pair, which is slow over NFS on any host that does not hold the map
    locally. The read is done under a lock because rasterio/GDAL dataset
    objects are not safe for concurrent reads, and because serializing keeps
    the 5 correction threads from hitting a remote mount at once. Only the
    read is locked; filtering and interpolation stay parallel. '''
    valid = np.abs(llxy['lat']) > 0
    with _velCacheLock:
        if velFile not in _velCache:
            _velCache[velFile] = rioxarray.open_rasterio(velFile)
        vel = _velCache[velFile]
        box = {}
        for coord in ['x', 'y']:
            box[f'min{coord}'] = max(np.round(np.min(llxy[coord][valid]) - pad,
                                              decimals=-3),
                                     vel[coord].min().values.item())
            box[f'max{coord}'] = min(np.round(np.max(llxy[coord][valid]) + pad,
                                              decimals=-3),
                                     vel[coord].max().values.item())
        velCrop = vel.rio.clip_box(**box)
        return velCrop.data, velCrop.x.values, velCrop.y.values


def interpVelBand(velData, band, velX, velY, llxy, sigma, noData=-2.0e9):
    ''' Interpolate one band (1=vx, 2=vy) of a cropped multi-band velocity map
    (from readVelCrop, in the same polar-stereo grid as llxy) onto the
    radar-grid points in llxy. '''
    valid = np.abs(llxy['lat']) > 0
    data = np.flipud(velData[band - 1]).astype(float)
    x = velX
    y = np.flipud(velY)
    data[data <= noData + 1] = np.nan
    data = gaussian_filter(np.nan_to_num(data), sigma=sigma)
    f = RegularGridInterpolator((y, x), data, method='linear',
                                bounds_error=False, fill_value=0)
    out = np.zeros(llxy['lat'].shape)
    out[valid] = f((llxy['y'][valid], llxy['x'][valid]))
    return out


_maskCache = {}
_maskCacheLock = threading.Lock()


def getIceMask(maskFile):
    ''' Load the tiff ice mask and cache it, keyed by file name.

    The mask covers the whole region, so it is the same for every burst and
    every frame, but readData reads the full map uncropped (~65 MB decoded for
    the 250 m Greenland mask). Without the cache runCorrections re-reads it
    once per burst -- 27 times per pair, 5 threads at a time -- which is the
    largest repeated read in the correction and the one that has hung on a
    stalled mount. The lock also stops the first read from being started by
    several threads at once. interpGeo and RegularGridInterpolator only read
    the object, so sharing one across threads is safe. '''
    with _maskCacheLock:
        if maskFile not in _maskCache:
            maskData = u.geoimage(geoType='scalar')
            maskData.readData(maskFile, tiff=True, dType='u1')
            # RegularGridInterpolator casts any integer input to float64, so
            # the uint8 mask would otherwise inflate 8x on the way in (65 MB
            # -> 517 MB for the 250 m Greenland mask). float32 is already far
            # more than a 0/1 mask needs and halves that.
            maskData.x = maskData.x.astype(np.float32)
            maskData.setupInterp(method='nearest')
            _maskCache[maskFile] = maskData
        return _maskCache[maskFile]


def getVel(llxy, regionDefs, sigma=1):
    ''' Read velocity and interpolate to the radar grid using xy from llxy'''
    # Ice mask (tiff mask used to zero out non-ice before correction)
    maskFile = regionDefs.region.get('icemaskTiff')
    maskData = getIceMask(maskFile)
    mask = maskData.interpGeo(llxy['x'] * 0.001, llxy['y'] * 0.001)  # needs km
    mask[np.isnan(mask)] = 0
    # Velocity: multi-band VRT/tiff (band1 vx, band2 vy) or legacy base name
    # with separate .vx.tif/.vy.tif companions.
    velFile = regionDefs.velMap()
    if velFile.endswith('.vrt') or velFile.endswith('.tif'):
        velData, velX, velY = readVelCrop(velFile, llxy)
        vx = interpVelBand(velData, 1, velX, velY, llxy, sigma)
        vy = interpVelBand(velData, 2, velX, velY, llxy, sigma)
    else:
        v = nisar.nisarVel(verbose=False)
        v.readDataFromTiff(velFile)
        v.vx = gaussian_filter(v.vx, sigma=sigma)  # Apply a light filter
        v.vy = gaussian_filter(v.vy, sigma=sigma)
        vx, vy = v.interp(llxy['x'], llxy['y'], noSpeed=True)
        noData = np.isnan(vx)
        vx[noData] = 0
        vy[noData] = 0
    return vx * mask, vy * mask


def readMergedOriginal(fileName):
    ''' Make a copy of the inteferogram as the original if it does not exist,
    then return the original '''
    if not os.path.exists(fileName):
        logFail(f'readMergedOriginal: file does not exist: {fileName}')
        u.myerror(f'readMergedOriginal: file does not exist: {fileName}')
    #
    copyISCEProduct(fileName, f'{fileName}.orig')
    # read and return data
    print(f'reading {fileName}.orig.vrt')
    x = rasterio.open(f'{fileName}.orig.vrt').read()[0]
    return x


def readComplexBurst(burstFile):
    ''' Read a complex burst '''
    if not os.path.exists(burstFile):
        logFail(f'readComplexBurst: file does not exist: {burstFile}')
        u.myerror(f'readComplexBurst: file does not exist: {burstFile}')
        return None
    return rasterio.open(f'{burstFile}.vrt').read()[0]


def readUW(fileName):
    ''' get uw from file '''
    x = rasterio.open(fileName)
    uw = x.read()
    return uw

# ---- Copy Routines


def copyISCEMeta(sourceFile, destFile, overwrite=False):
    ''' Make a copy of the meta file and change the file name internally'''
    if not os.path.exists(sourceFile):
        print(f'copyISCEMeta, file not found {sourceFile}')
        return
    # Strip suffix
    sourceName = '.'.join(os.path.basename(sourceFile).split('.')[0:-1])
    destName = '.'.join(os.path.basename(destFile).split('.')[0:-1])
    if os.path.exists(destFile) and not overwrite:
        return  # Copy already exists so don't copy unless overwrite is True
    # Read file and replace src baseName with dest baseName
    with open(sourceFile, 'r') as fpIn, open(destFile, 'w') as fpOut:
        # print(os.path.basename(sourceName), os.path.basename(destName))
        for line in fpIn:
            print(line.replace(sourceName, destName), file=fpOut, end='')


def copyISCEProduct(srcFile, destFile, overwrite=False):
    ''' make a copy of an isce product, including vrt and xmls with filename
    updates'''
    if not os.path.exists(srcFile):
        logFail(f'copyISCEProduct: file does not exist: {srcFile}')
        u.myerror(f'copyISCEProduct: file does not exist: {srcFile}')
    # Only copy file if it doesn't exist or overwrite True
    if os.path.exists(destFile) and not overwrite:
        return
    # Copy data
    shutil.copyfile(srcFile, destFile)
    # Copy meta
    for suffix in ['.xml', '.vrt']:
        copyISCEMeta(f'{srcFile}{suffix}', f'{destFile}{suffix}',
                     overwrite=overwrite)


# ---- Geocode


def lltoxy(lat, lon, epsg=3413):
    ''' ll to xy coordinates. Return values in metres (PS map units) '''
    toxy = pyproj.Transformer.from_crs('epsg:4326', f'epsg:{epsg}')
    x, y = toxy.transform(lat, lon)
    return x, y


# ---- Burst Routines


def setupBursts(dataPath, nalks=6, nrlks=30):
    ''' Get burst listing and make copies if they do not already exist '''
    beams = {}
    # Look dependent part of file names
    looks = f'.{nalks}alks_{nrlks}rlks'
    if nalks == 1 and nrlks == 1:
        looks = ''
    # Cycle through bursts and beams
    for beam in ['IW1', 'IW2', 'IW3']:
        beamBursts = sorted(glob.glob(
            f'{dataPath}/fine_interferogram/{beam}/burst_??{looks}.int'))
        # Compute names for originals
        origBursts = [x.replace('.int', '.int.orig') for x in beamBursts]
        for burst, origBurst in zip(beamBursts, origBursts):
            # Only copy if copy doesn't exists so subsequent corrections don't
            # overwrite
            copyISCEProduct(burst, origBurst, overwrite=False)
        beams[beam] = origBursts
    return beams


def azVel(vel, los, ll):
    ''' Rotate velocity to get azimuth directed component '''
    xyAngle = np.arctan2(-ll['y'], -ll['x'])
    rotAngle = np.deg2rad(los['azAngle']) - xyAngle
    va = vel['vx'] * np.sin(rotAngle) + vel['vy'] * np.cos(rotAngle)
    return va


def s1Squint(beamParams, beam, shape, nrlks=30, nalks=6):
    ''' Compute S1 beam dependent squint angle '''
    # print(beamParams[beam])
    nA, nR = shape
    azSquint = beamParams[beam]['azimuthsteeringrate'] / \
        beamParams[beam]['pulserepetitionfrequency'] * \
        ((np.arange(0, nA) - nA/2) * nalks)
    squint = np.repeat(azSquint, nR).reshape(nA, nR)
    return squint


def applyBurstCorrection(burstFile, phaseCorrection):
    ''' From apply the corretion to the "orig" burst and save the result in
    the standard isce file '''
    print('burst correction')
    # print(burstFile)
    burstData = readComplexBurst(burstFile)
    expCorrection = np.exp(phaseCorrection * (-1j))
    correctedBurst = burstData * expCorrection
    # print(burstData.shape)
    outputFile = burstFile.replace('.orig', '')
    correctedBurst.astype('complex64').tofile(outputFile)
    # print(outputFile)


def correctIndividualBurst(dataPath, beam, burstFile, beamParams,  nDays,
                           nrlks, nalks, regionDefs):
    ''' Fix a single burst '''
    # print('-', dataPath)
    burstNum = int(burstFile.split('burst')[-1].split('_')[1].split('.')[0])
    scaleFactor = nDays/365. * 4 * np.pi / beamParams[beam]['radarwavelength']
    llTemplate = f'{dataPath}/geom_reference/{beam}/*_{burstNum:02d}.rdr.vrt'
    llxy = getLatLonXY(llTemplate, regionDefs, nrlks=nrlks, nalks=nalks)
    losFile = f'{dataPath}/geom_reference/{beam}/los_{burstNum:02d}.rdr'
    los = getLOS(losFile, nrlks=nrlks, nalks=nalks)
    vel = dict(zip(['vx', 'vy', 'vv'], getVel(llxy, regionDefs)))
    va = azVel(vel, los, llxy)
    squint = s1Squint(beamParams, beam, vel['vx'].shape, nrlks=nrlks,
                      nalks=nalks)
    phaseCorrection = va * np.sin(squint) * scaleFactor
    applyBurstCorrection(burstFile, phaseCorrection)
    # print('+')


def duplicateUW():
    ''' Make a copy of the unwrapped file if not already present. In the
    single-unwrap flow (topsApp run only through burstifg) there is no prior
    unwrapped file to back up, so this is a no-op. '''
    if not os.path.exists('merged/filt_topophase.unw'):
        return
    if not os.path.exists('merged/filt_topophase.unw.orig'):
        copyISCEProduct('merged/filt_topophase.unw',
                        'merged/filt_topophase.unw.orig', overwrite=False)
        copyISCEProduct('merged/filt_topophase.unw.conncomp',
                        'merged/filt_topophase.unw.conncomp.unw',
                        overwrite=False)


def runCorrections(dataPath, beamFiles, beamParams,  nDays, regionDefs,
                   noAzCorrect, cpus=8, azThreads=3):
    ''' Run the azimuth corrections.

    azThreads is deliberately below the old hardcoded 5. Each thread holds
    full-resolution lat/lon/x/y for its burst (~1.1 GB), so the threads
    starting together produce an allocation spike -- and they start just as
    topsApp finishes, having left RAM full of page cache. Sustained reclaim
    over the following seconds is the PSI signature systemd-oomd kills on;
    all three known kills landed in exactly this window. At ~84 s per burst,
    dropping 5 -> 3 costs about 5 min per pair. '''
    if not noAzCorrect:
        threads = []
        for beam in beamFiles:
            for burstFile in beamFiles[beam]:
                myArgs = [dataPath, beam, burstFile, beamParams,  nDays, 1, 1,
                          regionDefs]
                thread = threading.Thread(target=correctIndividualBurst,
                                          args=myArgs)
                threads.append(thread)
        # Run threads
        u.runMyThreads(threads, azThreads, 'azCorrections', prompt=False)
    #
    duplicateUW()  # Copy previous if not already copied
    #
    isceShell('topsApp.py --start=mergebursts --end=filter', ompThreads=cpus)

# ---- Merged Routines


def getSlopes(llxy, demFile):
    ''' Get slopes from dem'''
    pad = 3000  # Pad area to compute slopes to avoid edge issues
    centerPt = (int(llxy['lat'].shape[0]/2), int(llxy['lat'].shape[1]/2))
    latC, lonC = llxy['lat'][centerPt], llxy['lon'][centerPt]
    valid = np.abs(llxy['lat']) > 0  # Only do points where there are data
    # Read dem
    dem = rioxarray.open_rasterio(demFile)
    # Setup a crop box
    box = {}
    for coord in ['x', 'y']:
        # make sure crop not outsize of dem (outer min/max ops)
        box[f'min{coord}'] = max(np.round(np.min(llxy[coord][valid]) - pad,
                                          decimals=-3),
                                 dem[coord].min().values.item())
        box[f'max{coord}'] = min(np.round(np.max(llxy[coord][valid]) + pad,
                                          decimals=-3),
                                 dem[coord].max().values.item())
    # Crop dem
    demCrop = dem.rio.clip_box(**box)
    # Grab data and compute gradients
    z = np.flipud(demCrop.data[0])
    x, y = demCrop.x, np.flipud(demCrop.y)
    # PS scale correction
    lengthScale = 1.0 / \
        pyproj.Proj(demCrop.rio.crs).get_factors(lonC, latC).parallel_scale
    dzdy, dzdx = np.gradient(z, y * lengthScale, x * lengthScale)
    # Setup interpolators and interpolate to radar coords for all valid points.
    fx = RegularGridInterpolator((y, x), dzdx, method='linear',
                                 bounds_error=False, fill_value=0)
    fy = RegularGridInterpolator((y, x), dzdy, method='linear',
                                 bounds_error=False, fill_value=0)
    dzdxR, dzdyR = np.zeros(llxy['lat'].shape), np.zeros(llxy['lat'].shape)
    dzdxR[valid] = fx((llxy['y'][valid], llxy['x'][valid]))
    dzdyR[valid] = fy((llxy['y'][valid], llxy['x'][valid]))
    return dict(zip(['dzdy', 'dzdx'], [dzdyR, dzdxR]))


def simPhase(vel, los, ll, dem, nDays, waveLength):
    ''' Compute phase due to motion using a velocity map '''
    xyAngle = np.arctan2(-ll['y'], -ll['x'])
    rotAngle = np.deg2rad(los['azAngle']) - xyAngle
    # Ground-range directed velocity
    vr = vel['vx'] * np.cos(rotAngle) - vel['vy'] * np.sin(rotAngle)
    # Compute vz for surface parallel flow
    grad = getSlopes(ll, dem)
    vz = grad['dzdx'] * vel['vx'] + grad['dzdy'] * vel['vy']
    # Compute horizontal and vertical corrections
    phaseScale = 4 * np.pi/waveLength / 365.25 * nDays
    p1 = vr * np.sin(np.deg2rad(los['incidence'])) * phaseScale
    p2 = - vz * np.cos(np.deg2rad(los['incidence'])) * phaseScale
    p = p1 + p2
    p[np.isnan(p)] = 0
    return p


def removePhase(interferogram, phaseSim):
    ''' Remove phase given as fp array from complex interferogram '''
    nonZero = np.abs(interferogram) > 0.0001  # Only correct nonZero values
    interferogram[nonZero] = interferogram[nonZero] * \
        np.exp(phaseSim[nonZero] * (-1j))
    interferogram[np.isnan(interferogram)] = 0 + 0j
    return interferogram, nonZero


def removeSimulated(dataPath, regionDefs, dem, nDays, waveLength):
    ''' Simulate phase from velocity and remove from complex inteferogram'''
    fileName = f'{dataPath}/merged/filt_topophase.flat'
    llxy = getLatLonXY(f'{dataPath}/merged/*.rdr.full.vrt', regionDefs,
                       nrlks=30, nalks=6)
    los = getLOS(f'{dataPath}/merged/los.rdr.vrt', nrlks=1, nalks=1)
    vel = dict(zip(['vx', 'vy', 'vv'], getVel(llxy, regionDefs, sigma=1)))
    phaseSim = simPhase(vel, los, llxy, dem, nDays, waveLength)
    print(phaseSim.shape, type(phaseSim))
    phaseSim.astype('float32').tofile(f'{dataPath}/merged/simPhase')
    # Use the ionosphere file to crate sim phase
    copyISCEMeta(f'{dataPath}/merged/topophase.ion.vrt',
                 f'{dataPath}/merged/simPhase.vrt', overwrite=False)
    #
    interferogram = readMergedOriginal(fileName)  # Read the interferogram
    interferogram, nonZero = removePhase(interferogram, phaseSim)
    # Overwrite prior result, which should now be backed up as orig
    interferogram.astype('complex64').tofile(fileName)
    return phaseSim, nonZero


def restoreSimulated(fileName, phaseSim, nonZero):
    ''' add the phase back to the unwrapped result'''
    uw = readUW(fileName)
    uw[1][nonZero] = uw[1][nonZero] + phaseSim[nonZero]
    uw = np.transpose(uw, axes=[1, 0, 2]).reshape(uw.shape[1], uw.shape[2]*2)
    uw.tofile(fileName)


def unWrapWithSimPhaseRemoved(dataPath, regionDefs, nDays, waveLength,
                              noPhaseRemove, cpus=8):
    ''' Simulate phase subtract, unwrap, add simulated phase back'''
    if not noPhaseRemove:
        duplicateUW()  # Make a copy first if one doesn't already exist
        demFile = regionDefs.dem(tiff=True)
        phaseSim, nonZero = removeSimulated(dataPath, regionDefs, demFile,
                                            nDays, waveLength)
    # Run unwrapper
    isceShell('topsApp.py --start=unwrap', ompThreads=cpus)
    uwFileName = f'{dataPath}/merged/filt_topophase.unw'
    cmpxFileName = f'{dataPath}/merged/filt_topophase.flat'
    # Now add back the phase
    if not noPhaseRemove:
        restoreSimulated(uwFileName, phaseSim, nonZero)
        # Restore the phase
        interferogram = rasterio.open(f'{cmpxFileName}.vrt').read()[0]
        #  interferogram = readMerged(cmpxFileName)  # Read the interferogram
        interferogram, nonZero = removePhase(interferogram, -phaseSim)
        # Overwrite prior result, which should now be backed up as orig
        interferogram.astype('complex64').tofile(cmpxFileName)


def otherDir(startDir, dirName):
    ''' Create new path, check dir exists, and if not create. Return path '''
    # in .../processingDir/track, so extract processingDir to replace
    baseName = os.path.basename(os.path.dirname(startDir))
    print(baseName, startDir)
    newDir = startDir.replace(baseName, dirName)
    if not os.path.exists(newDir):
        print(f'**** Creating {newDir}')
        os.mkdir(newDir)
    return newDir

# ---- Remap outputs


def remapOutput(dataPath, regionDefs, gimpDir, intermediatePath,
                resetMask=False):
    '''
    Create the intermediate ISCE product and the final GrIMP product.

    gimpDir          - track-level dir where the GrIMP product is written
    intermediatePath - full path to the intermediate product dir; its basename
                       must start "orbit1-orbit2" for downstream orbit parsing
    '''
    # Import here to avoid any import-order side effects
    from setupisceuw import setupUW
    print(f'gimpDir: {gimpDir}')
    print(f'intermediatePath: {intermediatePath}')
    for myDir in [gimpDir, os.path.dirname(intermediatePath)]:
        if not os.path.exists(myDir):
            os.makedirs(myDir)
    # Create intermediate product (convert_isce lives in the isce2grimp env)
    isceShell(f'convert_isce -i {dataPath} -o {intermediatePath}')
    # Copy stragglers
    simFileSource = f'{dataPath}/merged/simPhase'
    shutil.copyfile(simFileSource, f'{intermediatePath}/simPhase')
    corrSource = f'{dataPath}/merged/topophase.cor'
    shutil.copyfile(corrSource, f'{intermediatePath}/topophase.cor')
    shutil.copyfile(f'{corrSource}.vrt',
                    f'{intermediatePath}/topophase.cor.vrt')
    # Ionosphere estimate (if computed) - carried through for evaluation,
    # exported by setupUW as a separate tiff (not applied to the phase)
    ionSource = f'{dataPath}/merged/topophase.ion'
    if os.path.exists(ionSource):
        shutil.copyfile(ionSource, f'{intermediatePath}/topophase.ion')
        if os.path.exists(f'{ionSource}.vrt'):
            shutil.copyfile(f'{ionSource}.vrt',
                            f'{intermediatePath}/topophase.ion.vrt')
    # SE Tide correction (pass the ISCE dir so SETide finds merged/los.rdr
    # instead of reconstructing the old processing-tree path)
    command = f'SETide.py {intermediatePath} --isceDir {dataPath}'
    print(command)
    call(command, shell=True, executable='/bin/csh')
    # Actual conversion to GrIMP formats; returns the final product dir
    return setupUW(intermediatePath, gimpDir, regionDefs, resetMask=resetMask)


# ---- Run


def runAzPhaseCorrect(dataPath, regionDefs, gimpDir, intermediatePath,
                      noAzCorrect=False, noPhaseRemove=False, reMap=True,
                      gimpConvertOnly=False, unwrapOnly=False, resetMask=False,
                      cpus=8, azThreads=3):
    ''' Apply the burst-by-burst azimuth correction, unwrap with the simulated
    phase removed, and remap the result to intermediate and final GrIMP
    products.

    dataPath          - ISCE processing directory (topsApp run through burstifg)
    regionDefs        - sarfunc.defaultRegionDefs providing DEM/velMap/masks
    gimpDir           - track-level dir where the GrIMP product is written
    intermediatePath  - full path to the intermediate product dir
    '''
    processingDir = os.getcwd()
    if not gimpConvertOnly:
        os.chdir(dataPath)
        workPath = '.'
        nrlks, nalks = 1, 1
        # Make a copy of the bursts
        print('setting up bursts')
        beamFiles = setupBursts(workPath, nrlks=nrlks, nalks=nalks)
        print('getting burst params')
        beamParams = getBeamParams('.')
        nDays = getDT(workPath)
        waveLength = beamParams['IW1']['radarwavelength']
        # Loop through beam files
        if not unwrapOnly:
            print('run az corrections')
            t0 = time.time()
            runCorrections(workPath, beamFiles, beamParams,  nDays, regionDefs,
                           noAzCorrect, cpus=cpus, azThreads=azThreads)
            print(f'[TIME] azCorrect + mergebursts/filter: '
                  f'{time.time() - t0:.1f} s')
        #
        print('run unwrap')
        t0 = time.time()
        unWrapWithSimPhaseRemoved(workPath, regionDefs, nDays, waveLength,
                                  noPhaseRemove, cpus=cpus)
        print(f'[TIME] simPhase-remove + unwrap: {time.time() - t0:.1f} s')
        os.chdir(processingDir)
    #
    if reMap:
        print('convert to gimp')
        t0 = time.time()
        product = remapOutput(dataPath, regionDefs, gimpDir, intermediatePath,
                              resetMask=resetMask)
        print(f'[TIME] convert_isce + SETide + setupUW: '
              f'{time.time() - t0:.1f} s')
        return product
    return None


# ---- Main (standalone CLI)


def buildRegionDefs(region):
    ''' Legacy CLI helper: build regionDefs from a named region and inject the
    tiff ice mask that used to be hardcoded in icemaskMap. '''
    regionDefs = sf.defaultRegionDefs(region)
    if regionDefs.region.get('icemaskTiff') is None and region in icemaskMap:
        regionDefs.setRegionField('icemaskTiff', icemaskMap[region])
    return regionDefs


def main():
    ''' Standalone CLI entry. Preserves prior behavior: derive gimp/
    intermediate dirs as siblings of the processing dir and parse orbits from
    the ISCE dir name (track-frame-orbit1-orbit2). '''
    dataPath, region, noAzCorrect, noPhaseRemove, reMap, gimpConvertOnly, \
        unwrapOnly, outputSuffix = fileS1Args()
    processingDir = os.getcwd()
    regionDefs = buildRegionDefs(region)
    # Legacy directory derivation
    gimpDir = otherDir(processingDir, 'gimp' + outputSuffix)
    intermediateDir = otherDir(processingDir, 'intermediate' + outputSuffix)
    info = dict(zip(['track', 'frame', 'orbit1', 'orbit2'],
                    os.path.basename(dataPath).split('-')))
    if len(info.keys()) != 4:
        logFail(f'main: Could not parse path {dataPath}')
        u.myerror(f'main: Could not parse path {dataPath}')
    intermediatePath = \
        os.path.join(intermediateDir,
                     f'{info["orbit1"]}-{info["orbit2"]}-{info["frame"]}')
    runAzPhaseCorrect(dataPath, regionDefs, gimpDir, intermediatePath,
                      noAzCorrect=noAzCorrect, noPhaseRemove=noPhaseRemove,
                      reMap=reMap, gimpConvertOnly=gimpConvertOnly,
                      unwrapOnly=unwrapOnly)


if __name__ == "__main__":
    main()
