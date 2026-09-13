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
import re
from bs4 import BeautifulSoup
from datetime import datetime
import sarfunc as s
import rioxarray
import rasterio
import yaml
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

# A product with less than this much of its ice unwrapped is automatically
# given a hard Exclude by summaryQA.
EXCLUDEPERCENTVALID = 5.0
# Marks an Exclude this code wrote, so a later re-run can clear its own
# verdict without ever touching an Exclude placed by hand.
AUTOEXCLUDETAG = 'AUTO-EXCLUDED by setupisceuw'


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


def applyConnectedComponents(ccFile, uw, iceMask=None):
    ''' Keep a single connected component: the one covering the most ice, or -
    with no ice mask - the largest overall. Label 0 is snaphu's unreliable
    class and is never kept. Returns the component array for summaryQA.

    Selecting on ice rather than on raw size matters: a large component lying
    mostly on ocean/rock can otherwise beat a smaller one sitting entirely on
    ice, discarding most of the usable phase. '''
    cc = readCC(ccFile)
    labels = [int(label) for label in np.unique(cc) if label != 0]
    if not labels:
        # snaphu found nothing reliable - the whole frame is invalid
        u.mywarning(f'applyConnectedComponents: no components in {ccFile}')
        uw[:] = -2.0e9
        return cc
    if iceMask is not None:
        counts = {label: int(((cc == label) & iceMask).sum())
                  for label in labels}
    else:
        counts = {label: int(np.count_nonzero(cc == label))
                  for label in labels}
    uw[cc != max(counts, key=counts.get)] = -2.0e9
    return cc


def readIceMask(maskFile, georxa):
    ''' Boolean ice mask, or None when there is no mask file. Read before the
    connected-component step, which now selects on it. '''
    if not os.path.exists(maskFile):
        return None
    print(maskFile)
    return u.readImage(maskFile, georxa.nr, georxa.na, 'u1') == 1


def applyMask(iceMask, uw):
    ''' Blank everything off the ice mask. '''
    if iceMask is not None:
        uw[~iceMask] = -2.0e9


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


def autoExclude(outDir, percentValid, measure):
    ''' Write a hard Exclude when too little of the frame unwrapped, and
    return whether it is now excluded.

    A stale auto-Exclude is cleared when a re-run comes out above threshold,
    but only if it carries AUTOEXCLUDETAG - a hand-placed Exclude is a human
    decision and is never removed or overwritten here. '''
    excludeFile = f'{outDir}/Exclude'
    excluded = (percentValid is not None
                and percentValid < EXCLUDEPERCENTVALID)
    existing = ''
    if os.path.exists(excludeFile):
        with open(excludeFile) as fp:
            existing = fp.read()
    if existing and AUTOEXCLUDETAG not in existing:
        u.mywarning(f'autoExclude: leaving hand-placed Exclude in {outDir} '
                    f'untouched ({measure} = {percentValid}%)')
        return excluded
    if excluded:
        with open(excludeFile, 'w') as fp:
            print(f'{datetime.now()}: {AUTOEXCLUDETAG}: {measure} = '
                  f'{percentValid}% (below {EXCLUDEPERCENTVALID}%) - too '
                  f'little of this frame unwrapped to be useful for ties',
                  file=fp)
        print(f'wrote {excludeFile}: {measure} = {percentValid}%')
    elif existing:
        os.remove(excludeFile)
        print(f'cleared stale auto-Exclude in {outDir} '
              f'({measure} = {percentValid}%)')
    return excluded


def summaryQAfromProduct(productDir):
    ''' Rebuild summaryQA.yaml from a finished product, for products whose ISCE
    scratch is gone. Everything is read back from the product's own tiffs, so
    the connected-component metrics cannot be recovered and are written None.

    Note the coverage numbers reflect whatever component selection produced the
    product - reprocessing with -remapOnly, where the scratch survives, gives
    both the current selection and the full metrics. '''
    productDir = productDir.rstrip('/')

    def readIfPresent(stem):
        ''' Read <stem>.tif, or fall back to <stem>.vrt for the older products
        written as raw MSB rasters with a VRT sidecar. None if neither. '''
        for pattern in (f'{stem}.tif', f'{stem}.vrt'):
            hits = glob.glob(f'{productDir}/{pattern}')
            if hits:
                return rasterio.open(hits[0]).read()[0]
        return None

    uw = readIfPresent('*.isce.uw')
    if uw is None:
        u.mywarning(f'summaryQAfromProduct: no unwrapped phase in '
                    f'{productDir}')
        return None
    sim = readIfPresent('simPhase')
    corr = readIfPresent('*.isce.cor')
    # Ice mask is a raw byte image; take its dimensions from the phase grid
    iceMask = None
    maskFile = f"{productDir}/{dataFiles['icemask']}"
    if os.path.exists(maskFile):
        na, nr = uw.shape
        iceMask = u.readImage(maskFile, nr, na, 'u1') == 1
    # Provenance from the product dir name and the pairinfo file
    orbit1, frame = [int(v) for v in
                     os.path.basename(productDir).split('_')[:2]]
    meta = {'orbit1': orbit1, 'orbit2': None, 'frame': frame,
            'date1': None, 'date2': None, 'temporalBaseline': None,
            'looks': None, 'ionosphereEstimated':
                os.path.exists(f'{productDir}/ionosphere.tif')}
    pairFiles = glob.glob(f'{productDir}/*.pairinfo')
    if pairFiles:
        with open(pairFiles[0]) as fp:
            pieces = fp.read().split()
        if len(pieces) >= 4:
            meta['orbit2'] = int(pieces[1])
            meta['date1'], meta['date2'] = pieces[2], pieces[3]
            meta['temporalBaseline'] = abs(
                (datetime.strptime(pieces[3], '%Y-%m-%d')
                 - datetime.strptime(pieces[2], '%Y-%m-%d')).days)
    looks = re.search(r'\.(\d+x\d+)\.isce\.uw',
                      ' '.join(os.listdir(productDir)))
    if looks:
        meta['looks'] = looks.group(1)
    meta['rebuiltFromProduct'] = True
    return summaryQA(productDir, uw, None, iceMask, sim, corr, meta)


def summaryQA(outDir, uw, cc, iceMask, sim, corr, meta):
    ''' Write summaryQA.yaml beside the product: coverage, scatter about the
    simulated phase, correlation, and connected-component statistics. Written
    into the product dir so it travels with the product when it is copied out
    of scratch. All arrays are on the radar grid; phase values are radians.

    Every key is always written, None where the input needed for it was not
    available - QA rebuilt from a finished product (summaryQAfromProduct) has
    no connected-component file, so those entries come out None. '''
    def scalar(value, digits=4):
        ''' Plain float (not np.float32) so yaml writes a number, not a blob;
        None stays None. '''
        return None if value is None else round(float(value), digits)

    valid = uw > -1.99e9
    anyValid = bool(valid.any())
    qa = dict(meta)
    qa['nr'], qa['na'] = int(uw.shape[1]), int(uw.shape[0])
    qa['validPixels'] = int(valid.sum())
    qa['percentValid'] = scalar(100.0 * valid.mean(), 2)
    # Coverage over ice is what actually feeds the tie points; percentValid
    # over the whole frame mixes that up with how much of it is ocean/rock.
    qa['icePixels'] = int(iceMask.sum()) if iceMask is not None else None
    qa['percentValidOnIce'] = scalar(
        100.0 * (valid & iceMask).sum() / max(int(iceMask.sum()), 1), 2) \
        if iceMask is not None else None
    # Scatter about the simulated phase: how far the unwrapped result departs
    # from the velocity/topography model. Only the standard deviation is
    # meaningful - the unwrapped phase carries an arbitrary constant, so the
    # mean of the residual says nothing.
    qa['sigmaRelativeToSim'] = scalar(np.nanstd((uw - sim)[valid]), 3) \
        if sim is not None and anyValid else None
    qa['meanCorrelationAll'] = scalar(np.nanmean(corr), 4) \
        if corr is not None else None
    qa['meanCorrelationValid'] = scalar(np.nanmean(corr[valid]), 4) \
        if corr is not None and anyValid else None
    # Connected components. keptComponent is the one applyConnectedComponents
    # kept; largest* vs bestOnIce* show what selecting on raw size instead of
    # on ice would have given, so the difference stays visible.
    for key in ('numberOfConnectedComponents', 'keptComponent',
                'largestComponent', 'largestComponentPixels',
                'largestComponentOnIcePixels', 'bestOnIceComponent',
                'bestOnIceComponentPixels'):
        qa[key] = None
    if cc is not None:
        sizes = {int(label): int(np.count_nonzero(cc == label))
                 for label in np.unique(cc) if label != 0}
        qa['numberOfConnectedComponents'] = len(sizes)
        if sizes:
            largest = max(sizes, key=sizes.get)
            qa['largestComponent'] = largest
            qa['largestComponentPixels'] = sizes[largest]
            qa['keptComponent'] = largest
            if iceMask is not None:
                onIce = {label: int(((cc == label) & iceMask).sum())
                         for label in sizes}
                best = max(onIce, key=onIce.get)
                qa['largestComponentOnIcePixels'] = onIce[largest]
                qa['bestOnIceComponent'] = best
                qa['bestOnIceComponentPixels'] = onIce[best]
                qa['keptComponent'] = best
    # Exclude on the on-ice coverage where there is a mask; without one the
    # whole-frame number is all we have.
    if iceMask is not None:
        measure, percent = 'percentValidOnIce', qa['percentValidOnIce']
    else:
        measure, percent = 'percentValid', qa['percentValid']
    qa['automaticallyExcluded'] = autoExclude(outDir, percent, measure)
    qaFile = f'{outDir}/summaryQA.yaml'
    with open(qaFile, 'w') as fp:
        fp.write('# summaryQA.yaml - per-product quality metrics written by '
                 'setupisceuw.summaryQA\n')
        fp.write('# phase values are radians; sigmaRelativeToSim is the '
                 'scatter of (unwrapped - simulated)\n')
        yaml.safe_dump(qa, fp, default_flow_style=False, sort_keys=False)
    print(f'wrote {qaFile}')
    return qa


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
    # ice mask is read first: the connected-component choice is made on it
    iceMask = None
    if haveIceMask:
        iceMask = readIceMask(f"{outDir}/{dataFiles['icemask']}", georxa)
    # con comp mask - keeps the component covering the most ice
    cc = applyConnectedComponents(f'{iscePath}/{dataFiles["unwcc"]}', uw,
                                  iceMask)
    # ice mask
    applyMask(iceMask, uw)
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
    phase, correlation = None, None
    if os.path.exists(f'{iscePath}/simPhase'):
        phase = np.fromfile(f'{iscePath}/simPhase',
                            dtype='float32').reshape(na, nr)
        writeRadarTiff(phase, f'{outDir}/simPhase.tif', 'simulated phase')
    # Correlation
    if os.path.exists(f'{iscePath}/topophase.cor.vrt'):
        corr = readUW(f'{iscePath}/topophase.cor.vrt')
        correlation = corr[1]
        writeRadarTiff(correlation, f'{outDir}/{corName}.tif', 'correlation')
    # Ionosphere estimate - EXPORTED for evaluation, NOT applied to the phase
    haveIon = os.path.exists(f'{iscePath}/topophase.ion.vrt')
    if haveIon:
        ion = readUW(f'{iscePath}/topophase.ion.vrt')
        writeRadarTiff(ion[-1], f'{outDir}/ionosphere.tif',
                       'ionosphere phase (not applied)')
    # Quality metrics, from the arrays already in hand
    summaryQA(outDir, uw, cc, iceMask, phase, correlation,
              {'orbit1': orbit1, 'orbit2': orbit2, 'frame': frame1,
               'date1': orbDates[orbit1].date(),
               'date2': orbDates[orbit2].date(),
               'temporalBaseline':
                   abs((orbDates[orbit2] - orbDates[orbit1]).days),
               'looks': f'{georxa.nlr}x{georxa.nla}',
               'ionosphereEstimated': haveIon})
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
