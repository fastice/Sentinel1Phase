# azPhaseCorrect.py

Burst-by-burst and merged-scene InSAR phase correction for Sentinel-1 IW (TOPS)
interferograms, tuned for glacier/ice-sheet velocity. It sits between an ISCE
`topsApp.py` interferogram and the GrIMP unwrapped-phase products, applying two
physically distinct corrections and then remapping to GrIMP formats.

```
python azPhaseCorrect.py <isceDir> [--region {greenland,amundsen,taku}]
                         [--noAzCorrect] [--noPhaseRemove]
                         [--gimpConvertOnly] [--unwrapOnly] [--noRemap]
                         [--outputSuffix STR]
```

`<isceDir>` is an ISCE TOPS product directory named `track-frame-orbit1-orbit2`
(the basename is parsed on `-` in the remap stage).

---

## Why two corrections are needed

A repeat-pass interferogram over ice contains phase from surface motion. Two
Sentinel-1/TOPS-specific problems make the raw ISCE interferogram hard to use:

1. **Azimuth-motion × TOPS squint** produces a phase ramp that *jumps at every
   burst boundary*. This is a systematic, burst-periodic artifact unique to
   TOPS acquisitions and must be removed **before the bursts are merged**.
2. **Fast, steep flow** produces steep fringe gradients that the unwrapper
   trips over. Subtracting a *simulated* motion phase (from an external velocity
   map) flattens the interferogram so the unwrapper only has to resolve the
   residual, after which the simulated phase is added back.

The two corrections are independent, each with its own opt-out flag
(`--noAzCorrect`, `--noPhaseRemove`).

---

## Algorithm 1 — Azimuth burst phase correction

Functions: `correctIndividualBurst`, `s1Squint`, `azVel`, `applyBurstCorrection`,
driven by `runCorrections`.

### The physics

Sentinel-1 IW uses **TOPS** (Terrain Observation by Progressive Scans): within
each burst the antenna is electronically steered in azimuth at a constant rate
`azimuthsteeringrate` (rad/s). Consequently a ground target is imaged at a
**squint angle that depends on its azimuth position within the burst** — zero at
burst centre, growing (with opposite sign) toward the two burst edges.

With a non-zero squint θ, the radar line-of-sight has an along-track (azimuth)
component `sin θ`. Surface motion in the azimuth direction therefore leaks into
the measured range/phase, scaled by `sin θ`. Because θ varies linearly across the
burst and **resets each burst**, the resulting phase is a sawtooth that steps
discontinuously at burst seams in the merged product. The correction removes it
in the per-burst geometry, where it is well defined.

### The squint field — `s1Squint` (lines 301–309)

For a burst of azimuth length `nA` lines:

```
azSquint(n) = (azimuthsteeringrate / PRF) · (n − nA/2) · nalks       [rad]
```

- `azimuthsteeringrate / PRF` = radians of beam steering per azimuth line
  (`PRF` = `pulserepetitionfrequency`, lines/s).
- `(n − nA/2)·nalks` = azimuth line offset from burst centre (in single-look
  lines; `nalks` = azimuth looks, = 1 in the burst-level call).
- Result is broadcast across all range columns (`np.repeat … reshape(nA, nR)`):
  the squint is range-independent within a burst.

So the squint is **zero at burst centre, linear and antisymmetric to the edges**.

### The azimuth velocity component — `azVel` (lines 293–298)

The external velocity map is in polar-stereographic (PS) grid components
`(vx, vy)` (m/yr). It must be projected onto the satellite **along-track
(azimuth)** direction, accounting for PS grid convergence:

```
xyAngle  = arctan2(−y, −x)                    # PS grid convergence toward pole
rotAngle = deg2rad(azAngle) − xyAngle         # azimuth heading in the vx/vy frame
va       = vx·sin(rotAngle) + vy·cos(rotAngle)
```

- `azAngle` is the satellite azimuth (sat→ground, clockwise) derived from the
  ISCE LOS band: `azAngle = −los[1] + 180` in `getLOS` (ISCE stores a CCW
  ground→sat angle; the negation + 180° flips it to CW sat→ground). See
  **LOS conventions** below.
- `xyAngle = arctan2(−y,−x)` is the direction from the pixel toward the PS
  origin (the pole); it captures the meridian convergence so that a map-frame
  azimuth is expressed correctly against the local `x`/`y` grid axes.

### Assembling and applying the correction — `correctIndividualBurst` (326–341)

```
scaleFactor     = nDays/365. · 4π / radarwavelength           # m(azimuth disp.) → 2-way rad
phaseCorrection = va · sin(squint) · scaleFactor              [rad], shape (nA, nR)
```

Reading it right to left:

| factor | meaning |
|---|---|
| `va` | along-track velocity (m/yr) |
| `· nDays/365.` | → along-track displacement over the pair interval (m) |
| `· sin(squint)` | project azimuth displacement onto the squinted LOS |
| `· 4π/λ` | LOS displacement (m) → two-way interferometric phase (rad) |

`applyBurstCorrection` (312–323) then multiplies the **original** complex burst
by `exp(−1j · phaseCorrection)` and writes the corrected burst under the standard
ISCE name:

```
correctedBurst = burstData · exp(−1j · phaseCorrection)
```

`burstData` is read from the `.orig` copy (see **Backup convention**), so the
operation is idempotent across re-runs.

### Orchestration — `runCorrections` (355–374)

For every burst of every beam (`IW1/IW2/IW3`), a `correctIndividualBurst` thread
is queued and run 5-at-a-time via `u.runMyThreads`. Per burst it:
`getLatLonXY` (lat/lon → PS x/y) → `getLOS` → `getVel` (velocity interpolated to
the radar grid, ice-masked) → `azVel` → `s1Squint` → apply. The burst number is
parsed from the filename (`burst_NN`). After all bursts:

```
topsApp.py --start=mergebursts --end=filter
```

merges the corrected bursts and produces the filtered merged interferogram.

---

## Algorithm 2 — Simulated-phase removal for robust unwrapping

Functions: `simPhase`, `getSlopes`, `removePhase`, `removeSimulated`,
`restoreSimulated`, driven by `unWrapWithSimPhaseRemoved`.

### The idea

Simulate the LOS motion phase `phaseSim` from an external velocity map + DEM,
**subtract** it from the merged complex interferogram, unwrap the (now nearly
flat) residual, then **add** `phaseSim` back to the unwrapped result. This keeps
the unwrapper on well-behaved gradients where fast/steep flow would otherwise
alias.

### The simulated LOS phase — `simPhase` (417–432)

Two motion components are projected into the LOS and summed:

```
rotAngle = deg2rad(azAngle) − xyAngle                 # same rotation as azVel
vr = vx·cos(rotAngle) − vy·sin(rotAngle)              # ground-range-directed velocity
vz = dzdx·vx + dzdy·vy                                 # vertical vel., surface-parallel flow
phaseScale = 4π/λ / 365.25 · nDays                     # m/yr → 2-way rad over the pair

p1 =  vr · sin(incidence) · phaseScale                 # horizontal (range) LOS term
p2 = −vz · cos(incidence) · phaseScale                 # vertical LOS term
phaseSim = p1 + p2   (NaN → 0)
```

- `vr` is the **cross-track / ground-range** component of horizontal velocity —
  the orthogonal projection to `va` above, using the same `rotAngle`.
- **Surface-parallel flow**: the vertical velocity `vz` is inferred as the
  horizontal velocity dotted with the surface slope, `vz = v · ∇z`. Slopes
  `dzdx, dzdy` come from `getSlopes`.
- LOS projection: LOS displacement = `vr·sin(inc) − vz·cos(inc)`, the standard
  right-looking geometry (range-directed motion scaled by `sin(incidence)`,
  vertical by `cos(incidence)`), converted to two-way phase by `4π/λ`.

### DEM slopes — `getSlopes` (379–414)

Reads the region DEM (`sf.defaultRegionDefs(region).dem`), crops a padded box
around the scene footprint, and computes `dzdx, dzdy` with `np.gradient`. A
**polar-stereographic scale correction** is applied so gradients are in true
ground metres:

```
lengthScale = 1 / pyproj.Proj(crs).get_factors(lonC, latC).parallel_scale
dzdy, dzdx  = np.gradient(z, y·lengthScale, x·lengthScale)
```

Gradients are then interpolated (`RegularGridInterpolator`, `fill_value=0`) onto
the radar-grid `(x, y)` of the valid pixels.

### Subtract / unwrap / restore — `unWrapWithSimPhaseRemoved` (472–493)

```
removeSimulated(...):                                   # 444–461
    phaseSim = simPhase(vel, los, llxy, dem, nDays, λ)
    write merged/simPhase (+ .vrt cloned from topophase.ion.vrt)
    interferogram = readMergedOriginal(filt_topophase.flat)   # backs up → .orig
    interferogram, nonZero = removePhase(interferogram, phaseSim)   # × exp(−1j·phaseSim)
    overwrite filt_topophase.flat

topsApp.py --start=unwrap

restoreSimulated(unw, phaseSim, nonZero):               # 464–469
    uw[1][nonZero] += phaseSim[nonZero]                 # add motion phase back to band 1
    write filt_topophase.unw

# finally re-flatten the complex file with the opposite sign
interferogram, _ = removePhase(interferogram, −phaseSim)
overwrite filt_topophase.flat
```

`removePhase` (435–441) only touches pixels with `|interferogram| > 1e−4`
(`nonZero`) and zeroes NaNs, so masked/void pixels are left untouched and the
same `nonZero` mask is reused on restore.

> **Grid note**: in `removeSimulated`, lat/lon is read from `*.rdr.full.vrt`
> downsampled `30×6` (`getLatLonXY` default) while `los.rdr` is read at `1×1`
> — the multilooked `los.rdr` already sits on the `30×6` grid, so shapes match.
> Changing one look factor without the other breaks this.

---

## Supporting routines

| Function | Lines | Role |
|---|---|---|
| `fileS1Args` | 36–70 | argparse; returns the 8-tuple of options. |
| `burstSize` | 80–87 | burst `nRange`/`nAzimuth` from a coordinate XML. |
| `getDT` | 90–100 | `nDays` = secondary − reference ascending-node dates (`topsProc.xml`). |
| `getBeamParams` | 103–116 | PRF, steering rate, wavelength, line counts per beam. |
| `readXML` | 119–127 | BeautifulSoup XML reader with existence check. |
| `getLatLonXY` | 133–154 | lat/lon rasters → PS `x,y` (metres), epsg from region. |
| `getLOS` | 157–166 | incidence + azimuth angle from ISCE `los.rdr` (`azAngle=−los[1]+180`). |
| `getVel` | 169–190 | external velocity map (`nisar.nisarVel`) + ice mask → radar grid. |
| `lltoxy` | 262–266 | pyproj EPSG:4326 → PS transform. |
| `setupBursts` | 272–290 | glob bursts per beam, back up each to `.orig`. |
| `copyISCEProduct` / `copyISCEMeta` | 225–256 | copy data + rewrite basename inside `.xml`/`.vrt`. |
| `duplicateUW` | 345–352 | back up `filt_topophase.unw`(.conncomp) once. |
| `otherDir` | 496–505 | derive/create sibling `gimp*`/`intermediate*` output dirs. |
| `remapOutput` | 510–546 | `convert_isce` → copy stragglers → `SETide.py` → `setupisceuw.py`. |
| `logFail` | 75–77 | drop a `fail.<pid>` marker on error. |

---

## Conventions and gotchas

- **Backup / idempotency**: first run copies each correctable product to
  `<name>.orig`; corrections read `.orig` and write the standard name. Re-runs
  never double-correct. Delete the `.orig` files to force a clean re-baseline.
- **`main()` runs on import** — no `if __name__` guard (bottom of file).
- **Look factors**: burst-level az correction runs at full res (`nrlks,nalks=1,1`
  set in `main`); merged sim-phase runs on the `30×6` multilooked grid.
- **Year length**: `nDays/365.` (az correction) vs `nDays/365.25` (sim phase) —
  a minor inconsistency, harmless at current precision.
- **LOS conventions**: `getLOS` returns `incidence` (band 0) and
  `azAngle = −los[1] + 180` (band 1), converting ISCE's CCW ground→sat angle to a
  CW sat→ground azimuth used by both `azVel` and `simPhase`.
- **Sign conventions**: interferogram corrections are applied as
  `× exp(−1j·phase)`; sim phase is subtracted before unwrap and added back after
  (`restoreSimulated`), then the flat file is re-flattened with `−phaseSim`.
- **Region config** comes from `sf.defaultRegionDefs(region)` (epsg, `velMap`,
  `dem`); the Greenland ice mask path is hardcoded in `icemaskMap`.
- **Path parsing**: `remapOutput` requires the product-dir basename to split into
  exactly `track-frame-orbit1-orbit2`.

## External programs invoked

`topsApp.py` (ISCE; `mergebursts`/`filter`/`unwrap` stages, run via `csh` with
`OMP_NUM_THREADS` set), and the GrIMP remap CLIs `convert_isce`, `SETide.py`,
`setupisceuw.py`.
