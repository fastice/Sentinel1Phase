# setupisceuw.py

Final "remap to GrIMP" stage of the Sentinel-1 phase pipeline. Converts one ISCE
unwrapped intermediate product into a GrIMP-format frame directory — applying the
time-sign convention, connected-component and ice masks, and the solid-Earth-tide
correction, then writing big-endian GrIMP files and gap-filling the phase.

Called by [`azPhaseCorrect.remapOutput`](azPhaseCorrect.md) as the last step:

```
setupisceuw.py --region <region> <iscePath> <gimpDir>
```

```
setupisceuw.py [--region {greenland,antarctica}] [--resetMask]
               <isceProductPath> <gimpDir>
```

| Arg | Meaning |
|---|---|
| `isceProductPath` | ISCE intermediate product dir, basename `orbit1-orbit2-frame` (from `convert_isce`, after `SETide.py`). |
| `gimpDir` | Track-level directory under which the `orbit1_frame1/` output dir is created. |
| `--region` | `greenland` (default) or `antarctica`; selects DEM, velocity, mask, EPSG via `sarfunc.defaultRegionDefs`. |
| `--resetMask` | Regenerate the radar-geometry ice mask even if one already exists (default: keep existing). |

Input `orbit1-orbit2-frame` is remapped to an output frame directory
`gimpDir/orbit1_frame1/`. **`main()` runs on import** (no `__name__` guard), same
as `azPhaseCorrect.py`.

---

## What it produces

For a pair `orbit1`/`orbit2` at frame `frame1`, in `gimpDir/orbit1_frame1/`:

| File | Written by | Contents |
|---|---|---|
| `orbit1_frame1.orbit2_frame1.NlrxNla.isce.uw` | `mapUW` → `u.writeImage('>f4')` | unwrapped phase, sign-corrected + masked + tide-corrected (big-endian float32) |
| `…isce.uw.interp` | `intfloat` | gap-filled version of the `.uw` (distance-weighted fill) |
| `orbit1_frame1.NlrxNla.pow` | `mapUW` | power/amplitude band (band 0 of the ISCE `.unw`) |
| `orbit1_frame1.orbit2_frame1.NlrxNla.isce.cor` | `mapUW` | correlation (band 1 of `topophase.cor`) |
| `simPhase` | `mapUW` | motion-simulated phase copied from the ISCE product, byteswapped to big-endian |
| `icemask` | `simIceMask` (`siminsar`) | ice/no-ice mask in radar geometry |
| `geodatNlrxNla.in` | `setupOutputDir` (copied) | radar-geometry metadata for the frame |
| `orbit1.orbit2.pairinfo` | `makePairInfo` | `orbit1 orbit2 date1 date2 Nlr Nla` |
| `motion/baselines.orig` | `setupMotion` | zeroed baseline placeholder |
| `<iscePath basename>` (symlink) | `setupOutputDir` | link back to the intermediate product ("chain of custody") |

`Nlr`/`Nla` are the range/azimuth look counts, parsed from the geodat filename
(`geodatNlrxNla.in`) by `nLooksFromGeodat`.

---

## Pipeline (main → helpers)

```
setupISCEArgs      parse args; region → defaultRegionDefs; parseOrbits(basename)
setupOutputDir     getFrames (frames.F1.F2) → outDir=gimpDir/orbit1_frame1
                   mkdir, symlink back to iscePath, copy geodat*.in,
                   nLooksFromGeodat → nlr,nla
simIceMask         (greenland) siminsar -mask <dem> <icemask90m> <geodat> <icemask>
mapUW              the core conversion (below)
makePairInfo       write orbit1.orbit2.pairinfo
setupMotion        mkdir motion/, write zeroed baselines.orig
```

---

## Core conversion — `mapUW` (lines 220–261)

### 1. Time-sign convention (222–230)

GrIMP requires the unwrapped phase referenced so that **`orbit1` is the earlier
acquisition** (positive phase = motion over `date2 − date1`). ISCE does not
guarantee that ordering, so the sign is set from the acquisition dates:

```
orbDates = parseDates(iscePath)          # safe names in topsApp.xml → {orbit: date}
mySign   = np.sign((orbDates[orbit2] − orbDates[orbit1]).days)
uw       = mySign * unw[1]               # band 1 = unwrapped phase; band 0 = power
```

`parseDates` reads the two `safe` product names from `topsApp.xml` and
`parseOrbDateFromSafeName` extracts orbit number (field 7) and date (field 5) from
each S1 SAFE name.

### 2. Connected-component mask — `applyConnectedComponents` (173–180)

Keeps **only the largest connected component** of the unwrap; everything else is
set to the GrIMP noData sentinel `−2.0e9`:

```
cc      = readCC(filt_topophase.unw.conncomp)   # ISCE conn-comp labels, 0 = invalid
labels  = unique(cc) minus 0
maxKey  = label with the largest pixel count
uw[cc != maxKey] = −2.0e9
```

This discards isolated unwrapped islands whose integer-cycle offset relative to
the main region is unknown.

### 3. Ice mask — `applyMask` (183–188)

Reads the radar-geometry ice mask (`icemask`, produced by `simIceMask`) as `u1`
and blanks non-ice pixels:

```
uw[iceMask != 1] = −2.0e9
```

### 4. Solid-Earth-tide correction — `applySETide` (191–204)

Removes the solid-Earth tide contribution. The tide file (`SECorrection.vrt`,
written by the upstream `SETide.py` step) is a **displacement**; it is scaled to
phase and **subtracted** on valid pixels only:

```
SETide *= 4π / wavelength                 # metres → radians (one-way scale as stored)
valid   = uw > −1.99e9                     # exclude the −2e9 noData sentinel
uw[valid] = uw[valid] − SETide[valid]
```

The in-code note records the sign rationale: the tide enters like a `−Vz` term, so
the corrected phase is `uw − SETide`.

### 5. Output + gap fill (241–260)

- Writes the corrected `uw` and the power band as big-endian float32
  (`u.writeImage(…, '>f4')`) — the GrIMP MSB convention.
- Runs `intfloat` to distance-weight-fill gaps:
  `intfloat -wdist -thresh 100 -nr <nr> -na <na> -islandAreaThresh 300` →
  `…isce.uw.interp`.
- Copies `simPhase` (ISCE little-endian float32 → **byteswapped** to big-endian).
- Copies the correlation band (`topophase.cor` band 1) to `…isce.cor`.

`georxa = u.geodatrxa(geodat)` supplies `nr, na, nlr, nla, wavelength`.

---

## Conventions and gotchas

- **noData sentinel** is `−2.0e9`, matching `intfloat`/GrIMP `minValue`; the valid
  test is `> −1.99e9`.
- **Byte order**: all GrIMP outputs are big-endian (`'>f4'`); `simPhase` is
  explicitly byteswapped because ISCE stores it little-endian.
- **`orbit1` = earliest** is enforced by `mySign`, not assumed from the path order.
- **Frame numbering**: the output dir and filenames use `frame1` for both ends
  (`orbit1_frame1.orbit2_frame1`), i.e. a single frame per pair.
- **Region config** (`defaultRegionDefs`): `.dem()` and the hardcoded
  `icemaskMap` 90 m GIMP mask feed `siminsar`; `antarctica` has `icemaskMap` =
  `None`.
- **`main()` runs on import** (no `__name__` guard).

## Known issues (flagged, not fixed)

- **`makePairInfo` never closes its file** — line 269 is `fp.close` (missing
  `()`), so `orbit1.orbit2.pairinfo` is only flushed on interpreter exit.
- **The `None` guards check the wrong object**: `main` and `mapUW` test
  `dataFiles['icemaskMap'] is not None` / `dataFiles['SECorrection'] is not None`,
  but those are a dict and a string literal — **always truthy**. The intended test
  is the per-region entry `icemaskMap[region]`. For `--region antarctica`
  (`icemaskMap['antarctica'] = None`) `simIceMask` therefore passes the literal
  string `None` to `siminsar`; the Greenland (default) path is unaffected.

## External programs invoked

`siminsar` (radar-geometry ice-mask simulation) and `intfloat` (gap fill) — both
GrIMP C binaries — plus `sarfunc`/`utilities` Python helpers
(`geodatrxa`, `readImage`/`writeImage`, `callMyProg`, `defaultRegionDefs`).
