#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 12:14:01 2019

@author: ian
"""
import utilities as u
import shutil
import argparse
import os
import numpy as np
import glob
from bs4 import BeautifulSoup
from datetime import datetime
import sarfunc as s
import rioxarray
import rasterio
from osgeo import gdal
from subprocess import call

# file name and path info
dataFiles = {'unw': 'filt_topophase.unw',
             'unwcc': 'filt_topophase.unw.conncomp',
             'icemaskMap': {'greenland':
                            '/Volumes/insar7/ian/gimp/mask/GimpIceMask_90m',
                            'antarctica': None},
             'icemask': 'icemask',
             'SECorrection': 'SECorrection.vrt',
             'region': None,
             'regionData': None}


def setupISCEArgs():
    ''' Handle command line args'''
    parser = argparse.ArgumentParser(description='\n\n\033[1mMap ISCE '
                                     'directory int-orbit1-orbit2 to '
                                     'orbit1_frame\033[0m\n\n')
    parser.add_argument('--region', type=str, default='greenland',
                        choices={'greenland', 'antarctica'})

    parser.add_argument('--resetMask', action='store_true', default=False,
                        help='If an icemask exists, generate a new version'
                        ' [keep old mask]')
    parser.add_argument('isceProductPath', metavar='iscePath', type=str,
                        nargs=1, help='Path to reduced ISCE unwrapped product')
    parser.add_argument('gimpDir', metavar='gimpPath', type=str,
                        nargs=1,
                        help='Track level dir where gimp products are saved')
    args = parser.parse_args()
    #
    print(args.region)
    regionDefs = s.defaultRegionDefs(args.region)
    #
    if not os.path.exists(args.isceProductPath[0]):
        u.myerror(f'Input file - {args.isceProductPath[0]} - does not exist')

    orbit1, orbit2 = parseOrbits(args.isceProductPath[0].rstrip('/'))
    return args.isceProductPath[0].rstrip('/'), args.gimpDir[0].rstrip('/'), \
        orbit1, orbit2, args.resetMask, regionDefs


def returnProperty(xmlData, myProperty):
    myProps = []
    for x in xmlData.find_all('property'):
        if x['name'] == myProperty:
            myProps.append(x.getText())
    return myProps


def simIceMask(outputDir, geodat, regionDefs, resetMask):
    #
    maskFile = f"{outputDir}/{dataFiles['icemask']}"
    # return if not resetMask and there is already a mask
    if not resetMask and os.path.exists(maskFile):
        return
    # generate mask
    args = f"-mask {regionDefs.dem()} " \
        f"{regionDefs.icemask()} {geodat} {maskFile}"
    command = 'siminsar'
    u.callMyProg(command, myArgs=args.split(), screen=True)
    


def readXML(iscePath, fileName, myName):
    myXmlFile = f'{iscePath}/{fileName}'
    fp = checkAndOpenMyFile(myXmlFile)
    xmlfp = fp.read()
    xmlData = BeautifulSoup(xmlfp, 'xml')
    fp.close()
    return returnProperty(xmlData, myName)


def parseOrbDateFromSafeName(safe):
    # split on S1 in case _ in path name
    pieces = safe.split('S1')[1].split('_')
    orbit = int(pieces[7])
    date = datetime.strptime(pieces[5].split('T')[0], "%Y%m%d")
    return {orbit: date}


def parseDates(iscePath):
    # read safe names from xml and use to generate dict for orb > date
    safe1, safe2 = readXML(iscePath, 'topsApp.xml', 'safe')
    orbDate = parseOrbDateFromSafeName(safe1)
    orbDate.update(parseOrbDateFromSafeName(safe2))
    return orbDate


def parseOrbits(iscePath):
    pieces = os.path.basename(iscePath).split('-')
    try:
        orbit1, orbit2 = int(pieces[0]), int(pieces[1])
    except Exception:
        u.myerror(f'Problem parsing orbit numbers from path: {iscePath}')
    return orbit1, orbit2


def getFrames(iscePath):
    frameFile = glob.glob(f'{iscePath}/frames.*.*')
    if len(frameFile) > 0:
        pieces = frameFile[0].split('frames')[1].split('.')
        return int(pieces[1]), int(pieces[2])
    else:
        u.myerror('Could not find frames.*.* file ')


def nLooksFromGeodat(geodat):
    lString = os.path.basename(geodat).split('.')[0].split('geodat')[1]
    pieces = lString.split('x')
    return int(pieces[0]), int(pieces[1])


def checkMyFile(myFile, function=''):
    if not os.path.exists(myFile):
        u.myerror(f'{function}: {myFile} does not exist')


def checkAndOpenMyFile(myFile):
    #
    checkMyFile(myFile, function='checkAndOpenMyFile')
    return open(myFile, 'r')


def setupOutputDir(iscePath, outputPath, orbit1):
    # get frames
    frame1, frame2 = getFrames(iscePath)
    #
    outDir = f'{outputPath}/{orbit1}_{frame1}'
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    # link from gimp back to intermediate to establish change of custody
    linkIntermediate = f'{outDir}/{os.path.basename(iscePath)}'
    try:
        if not os.path.exists(linkIntermediate):
            os.symlink(iscePath, linkIntermediate)
    except Exception:
        print('Could not create linke')  # Shouldn't happen, but avoid fail
    #
    # Copy the geodats (geojson primary + secondary, plus legacy .in) into the
    # product; use the primary geojson downstream.
    geodats = glob.glob(f'{iscePath}/geodat*.geojson') + \
        glob.glob(f'{iscePath}/geodat*.in')
    if len(geodats) == 0:
        u.myerror(f'{iscePath}/geodat*.geojson not found')
    for g in geodats:
        shutil.copyfile(g, f'{outDir}/{os.path.basename(g)}')
    primary = [g for g in geodats
               if g.endswith('.geojson') and '.secondary.' not in g]
    geodat = primary[0] if primary else geodats[0]
    newGeodat = f'{outDir}/{os.path.basename(geodat)}'
    #
    nlr, nla = nLooksFromGeodat(geodat)
    return outDir, newGeodat, frame1, frame2, nlr, nla


def applyConnectedComponents(ccFile, uw):
    ''' Keep only largest cc'''
    cc = readCC(ccFile)
    labels = np.delete(np.unique(cc), [0])
    counts = dict(zip(labels,
                      [np.count_nonzero(cc == x) for x in labels if x != 0]))
    maxKey = max(counts, key=counts.get)
    uw[cc != maxKey] = -2.0e9


def applyMask(maskFile, georxa, uw):
    #
    if os.path.exists(maskFile):
        print(maskFile)
        iceMask = u.readImage(maskFile, georxa.nr, georxa.na, 'u1')
        uw[iceMask != 1] = -2.0e9


def applySETide(SETideFile, georxa, uw):
    if os.path.exists(SETideFile):
        print('applying solid earth')
        SETide = np.squeeze(rasterio.open(SETideFile).read())
        SETide *= 4. * np.pi / georxa.wavelength
        print(uw.shape, SETide.shape)
        valid = uw > -1.99e9
        print(np.nanstd(uw[valid]))
        print(np.sum(valid), np.nanmean(SETide[valid]))
        # tide caculation is ~ (... -Vz) , which -> uw = uw -  SEtide
        uw[valid] = uw[valid] - SETide[valid]
        print(np.nanstd(uw[valid]))
    else:
        u.mywarning('No Solid Earth Correction')


def readUW(fileName):
    ''' get uw from file '''
    x = rasterio.open(fileName)
    uw = x.read()
    return uw


def readCC(ccFile):
    ''' Read connected components file '''
    checkMyFile(ccFile, function='readCC')
    return rasterio.open(ccFile).read()[0]


def writeRawFloatVrt(rawFile, nr, na, description='', noData=-2.0e9):
    ''' Write a GDAL VRT beside a raw big-endian float32 (>f4) radar-grid file
    so it can be viewed with showImage/gdal. Mirrors the icemask.vrt layout
    (identity geotransform, VRTRawRasterBand) but Float32/MSB. '''
    if not os.path.exists(rawFile):
        return
    vrt = f'''<VRTDataset rasterXSize="{nr}" rasterYSize="{na}">
  <GeoTransform> -5.0000000000000000e-01,  1.0000000000000000e+00,  0.0000000000000000e+00, -5.0000000000000000e-01,  0.0000000000000000e+00,  1.0000000000000000e+00</GeoTransform>
  <VRTRasterBand dataType="Float32" band="1" blockXSize="{nr}" blockYSize="1" subClass="VRTRawRasterBand">
    <Metadata>
      <MDI key="Description">{description}</MDI>
    </Metadata>
    <NoDataValue>{int(noData)}</NoDataValue>
    <SourceFilename relativeToVRT="1">{os.path.basename(rawFile)}</SourceFilename>
    <ImageOffset>0</ImageOffset>
    <PixelOffset>4</PixelOffset>
    <LineOffset>{nr * 4}</LineOffset>
    <ByteOrder>MSB</ByteOrder>
  </VRTRasterBand>
</VRTDataset>
'''
    with open(f'{rawFile}.vrt', 'w') as fp:
        fp.write(vrt)


def writeRadarTiff(data, filename, description='', noData=-2.0e9):
    ''' Write a numpy array as a GeoTIFF on the radar grid (identity transform,
    LZW) plus a sidecar VRT wrapping it. The .tif is the viewable product; the
    .vrt is what the GrIMP tie/mosaic tools (tiepoints -yaml, intfloat -inputVRT,
    mosaic3d) read -- the tie workflow globs by extension and locates phase via
    a fixed-name phase.uw.vrt symlink (see mapUW). Geometry for mosaic3d comes
    from the geojson geodat, not this transform. '''
    na, nr = data.shape
    transform = rasterio.transform.Affine(1.0, 0.0, -0.5, 0.0, 1.0, -0.5)
    with rasterio.open(filename, 'w', driver='GTiff', height=na, width=nr,
                       count=1, dtype='float32', nodata=noData,
                       compress='lzw', transform=transform) as dst:
        dst.write(data.astype('float32'), 1)
        if description:
            dst.update_tags(1, Description=description)
    # Sidecar VRT wrapping the tif (NISAR x.tif + x.vrt convention). Translate,
    # not BuildVRT: gdalbuildvrt rejects this radar grid's positive NS
    # resolution ('does not support positive NS resolution'), skips its only
    # source and returns None. Written beside the tif, so the source reference
    # stays relative and the product survives being copied out of scratch.
    vrt = gdal.Translate(f'{filename[:-4]}.vrt', filename, format='VRT')
    vrt.FlushCache()
    vrt = None


def mapUW(iscePath, outDir, orbit1, orbit2, frame1, geodat, haveIceMask):
    ''' get the sign based on convention that orbit1 is earliest in time -
    negate if otherwise '''
    orbDates = parseDates(iscePath)
    print(orbDates, orbit1, orbit2)
    print(iscePath, outDir, orbit1, orbit2, frame1, geodat)
    mySign = np.sign((orbDates[orbit2] - orbDates[orbit1]).days)
    #
    georxa = u.geodatrxa(file=geodat)
    unw = readUW(f'{iscePath}/{dataFiles["unw"]}')
    uw = mySign * unw[1]
    p = unw[0]
    # con comp mask
    applyConnectedComponents(f'{iscePath}/{dataFiles["unwcc"]}', uw)
    # ice mask
    if haveIceMask:
        applyMask(f"{outDir}/{dataFiles['icemask']}", georxa, uw)
    # SE Correction
    if dataFiles['SECorrection'] is not None:
        applySETide(f'{iscePath}/{dataFiles["SECorrection"]}', georxa, uw)
    #
    nr, na = georxa.nr, georxa.na
    uwName = f'{orbit1}_{frame1}.{orbit2}_{frame1}.{georxa.nlr}x' \
        f'{georxa.nla}.isce.uw'
    powName = f'{orbit1}_{frame1}.{georxa.nlr}x{georxa.nla}.pow'
    corName = f'{orbit1}_{frame1}.{orbit2}_{frame1}.{georxa.nlr}x' \
        f'{georxa.nla}.isce.cor'
    # Unwrapped phase + power as GeoTIFF (+ sidecar VRT) — final product format
    writeRadarTiff(uw, f'{outDir}/{uwName}.tif', 'unwrapped phase')
    writeRadarTiff(p, f'{outDir}/{powName}.tif', 'power')
    # Fixed-name symlink tieScript.process_phase_yaml() locates the phase by
    # (phase.uw.vrt -> <uw>.vrt), regardless of the S1 product naming.
    phaseLink = f'{outDir}/phase.uw.vrt'
    if os.path.lexists(phaseLink):
        os.remove(phaseLink)
    os.symlink(f'{uwName}.vrt', phaseLink)
    # Interp: intfloat reads/writes big-endian (>f4); write a temp, interp,
    # read back as >f4 (byte order must match intfloat's MSB), tiff, clean up.
    tmpUw = f'{outDir}/{uwName}'
    u.writeImage(tmpUw, uw, '>f4')
    call(f'intfloat -wdist -nr {nr} -na {na} -ratThresh 1 -thresh 100 '
         f'-islandThresh 300 {tmpUw} > {tmpUw}.interp',
         shell=True, executable='/bin/csh')
    interp = u.readImage(f'{tmpUw}.interp', nr, na, '>f4')
    writeRadarTiff(interp, f'{outDir}/{uwName}.interp.tif',
                   'unwrapped phase (interp)')
    os.remove(tmpUw)
    os.remove(f'{tmpUw}.interp')
    # Simulated (topo+motion) phase
    if os.path.exists(f'{iscePath}/simPhase'):
        phase = np.fromfile(f'{iscePath}/simPhase',
                            dtype='float32').reshape(na, nr)
        writeRadarTiff(phase, f'{outDir}/simPhase.tif', 'simulated phase')
    # Correlation
    if os.path.exists(f'{iscePath}/topophase.cor.vrt'):
        corr = readUW(f'{iscePath}/topophase.cor.vrt')
        writeRadarTiff(corr[1], f'{outDir}/{corName}.tif', 'correlation')
    # Ionosphere estimate - EXPORTED for evaluation, NOT applied to the phase
    if os.path.exists(f'{iscePath}/topophase.ion.vrt'):
        ion = readUW(f'{iscePath}/topophase.ion.vrt')
        writeRadarTiff(ion[-1], f'{outDir}/ionosphere.tif',
                       'ionosphere phase (not applied)')
    return orbDates, georxa


def makePairInfo(outDir, orbit1, orbit2, orbDates, georxa):
    fp = open(os.path.join(outDir, f'{orbit1}.{orbit2}.pairinfo'), 'w')
    print(f'{orbit1}  {orbit2}  {orbDates[orbit1].strftime("%Y-%m-%d")} '
          f'{orbDates[orbit2].strftime("%Y-%m-%d")} {georxa.nlr} {georxa.nla}',
          file=fp)
    fp.close


def setupMotion(outDir, orbit1, orbit2, frame1):
    motionDir = f'{outDir}/motion'
    if not os.path.exists(motionDir):
        os.mkdir(motionDir)
    baselineFile = f'{outDir}/motion/baselines.orig'
    fpBase = open(baselineFile, 'w')
    print('0. 0. 0. 0.\n0. 0. 0. 0.\n', file=fpBase)
    fpBase.close()


def setupUW(iscePath, gimpDir, regionDefs, resetMask=False):
    ''' Convert an intermediate ISCE product (int-orbit1-orbit2) into a final
    GrIMP product under gimpDir/orbit1_frame1. Importable entry point used by
    azPhaseCorrect.remapOutput and runS1interferogram. '''
    iscePath = iscePath.rstrip('/')
    gimpDir = gimpDir.rstrip('/')
    orbit1, orbit2 = parseOrbits(iscePath)
    #
    outDir, geodat, frame1, frame2, nlr, nla = setupOutputDir(iscePath,
                                                              gimpDir, orbit1)
    #
    haveIceMask = regionDefs.icemask() is not None
    if haveIceMask:
        simIceMask(outDir, geodat, regionDefs, resetMask)
    #
    orbDates, georxa = mapUW(iscePath, outDir, orbit1, orbit2, frame1, geodat,
                             haveIceMask)
    #
    makePairInfo(outDir, orbit1, orbit2, orbDates, georxa)
    #
    setupMotion(outDir, orbit1, orbit2, frame1)
    #
    print(outDir, geodat, frame1, frame2)
    print(orbit1, orbit2)
    return outDir


def main():
    ''' Convert isce intermediate product to gimp formats'''
    iscePath, gimpDir, orbit1, orbit2, resetMask, regionDefs = setupISCEArgs()
    setupUW(iscePath, gimpDir, regionDefs, resetMask=resetMask)


if __name__ == "__main__":
    main()
