# runS1interferogram

Single command that turns two Sentinel-1 SLC SAFE zips into one GrIMP
unwrapped-interferogram product. It consolidates the old
`prep_pair` → `run_isce` → `azPhaseCorrect` → `convert_isce`/`SETide`/`setupisceuw`
chain into one step, driven by a `uwproject.yaml`.

Run it from your **base** conda environment. The Python workflow runs in base;
the ISCE-dependent shell commands (`topsApp.py`, `convert_isce`, `gdalwarp`,
`gdal2isce_xml.py`) are dispatched into the `isce2grimp` env with `conda run` by
[`isceenv.isceShell`](../isceenv.py), so the ISCE env never has to be activated.

Pair *selection* is not done here — a valid pair is assumed to be given on the
command line. Use [`findS1InsarPairs`](../findS1InsarPairs.py) to produce pairs.

---

## Usage

```
runS1interferogram REFERENCE.zip SECONDARY.zip [--project uwproject.yaml]
```

```
# basenames resolved against dataDir / archiveDir
runS1interferogram S1C_IW_SLC__1SDH_20251108T095750_..._B629.zip \
                   S1C_IW_SLC__1SDH_20251120T095750_..._A4F1.zip

# explicit track, more threads, keep the ISCE files for debugging
runS1interferogram ref.zip sec.zip --track 25 --cpus 16 --debug
```

Feeding it from the pair lister: `findS1InsarPairs --queue` emits ready-to-run
commands for the pairs with no product yet (a count of queued vs already-done
goes to stderr, so stdout stays pipeable):

```
findS1InsarPairs 25 --firstDate 2025-11-01 --queue > runQueue
csh runQueue
```

**The queue lines are csh, not sh** — they use `set rc = $status`, backticks,
and `>>!`. Run them with `csh` (or paste them at a tcsh prompt); `sh runQueue`
will not work.

To also process pairs whose images are not downloaded yet, add `--onlinePairs`:
the lines then begin with a [`pullASF`](../../asfSearchAndDownload/Documents/pullASF.md)
that fetches whichever images the pair is missing before running it.

```
findS1InsarPairs all --onlinePairs --queue --firstDate 2024-01-01 --lastDate 2024-02-15 > runQueue
csh runQueue
```

```
pullASF <url1> <url2> --archiveDir <A> && runS1interferogram <f1> <f2> --project <p>; set rc = ...
```

One `pullASF` per pair, not one per granule: downloads are throttled to one at a
time archive-wide, so fetching a whole pair under a single acquisition of that
slot lets each job start computing while the next one downloads. stderr reports
how many pairs need a download and **how much data that is** — a wide window can
queue terabytes, and `--maxPairs` caps it.

Each line appends its own exit status to `<outputRoot>/logs/exitStatus.log`:

```
2026-08-21T13:00:21 rc=0 track=2 frame=426 orbits=63424-63599
```

That record exists because `runS1interferogram` cannot report on its own death.
A run killed by SIGKILL (`rc=137`, typically the OOM killer) or by a segfault in
a native library (`rc=139`, e.g. GDAL/numpy) never executes any of its own code,
so it writes no `Fail.*` record — the two cases are indistinguishable from the
product tree alone, and `rc` is what separates them. See Failure reporting below.

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `reference` | — | Reference (earlier) SLC SAFE zip. Path, or a basename resolved in `dataDir` then `archiveDir/YYYY-MM/`. Both `.zip` and `.zip.1` are tried at each location, so the name resolves whichever suffix the file currently has — `.1` is the marker `fileS1` writes once a granule is unpacked into the assembly tree, and a queue line generated hours earlier must still resolve if filing renamed it since. |
| `secondary` | — | Secondary SLC SAFE zip, same resolution rules. |
| `--project PATH` | `uwproject.yaml` | Project yaml (see keys below). |
| `--track N` | from manifest | Relative-orbit/track number. Read from `relativeOrbitNumber` in the reference `manifest.safe` if not given; aborts if neither is available. |
| `--cpus N` | `threads[host]`, else `isce.cpus`, else 8 | `OMP_NUM_THREADS` for `topsApp` and the correction stage. |
| `--debug` | False | Keep the ISCE and intermediate files in scratch. |
| `--keepScratch` | False | Same effect on cleanup as `--debug` (do not delete the scratch working dir). |
| `--noAzCorrect` | False | Skip the azimuth burst correction (passed through to `azPhaseCorrect`). |
| `--noPhaseRemove` | False | Skip simulated-phase removal before unwrapping. |

Scratch is **always** kept when the run raises — the failure path logs
`scratch kept for debugging: <path>` and re-raises.

---

## What it does

1. **Config / host setup.** `uwproject.loadProject` reads the yaml;
   `isceenv.setEnv` picks the ISCE conda env; `scratchForHost` /`threadsForHost`
   select the machine-local scratch dir and cpu count by hostname;
   `regionDefsFromProject` builds the `sarfunc.defaultRegionDefs` used for phase
   simulation and masking.
2. **Resolve the pair.** `resolveSafe` finds each SAFE (existing path → `dataDir`
   → `archiveDir/YYYY-MM/`), `parseSafe` pulls mission / absolute orbit /
   acquisition start and stop from the filename, and the reference
   `manifest.safe` is read once for the track (`trackFromManifest`) and the
   burst frame (`frameFromManifest`).
3. **Verify precise orbits.** `checkOrbits` requires a POEORB `.EOF` in
   `orbitDir` whose validity window spans each acquisition, and aborts before
   any processing if one is missing (see below).
4. **topsApp inputs.** `linkSafe` symlinks each zip into the scratch ISCE dir
   under a plain `.zip` name (ISCE's reader dispatches on the extension, so a
   `.zip.1` must be linked). `prepareIsceDem` builds the lat/lon DEM,
   `writeTopsAppXML` writes `topsApp.xml`.
5. **Run topsApp through burst interferograms only** — `topsApp.py --end=burstifg`,
   in the scratch ISCE dir, with `OMP_NUM_THREADS=<cpus>` and
   `OMP_PLACES=sockets(1)`. Merging, filtering and unwrapping happen later,
   inside `azPhaseCorrect`, because the burst correction must be applied to
   `fine_interferogram/` *before* `mergebursts`.
6. **Correct, unwrap, remap** — `azPhaseCorrect.runAzPhaseCorrect` applies the
   azimuth burst correction, removes the simulated phase, unwraps once, restores
   the simulated phase, then remaps via `convert_isce` → `SETide.py` →
   `setupisceuw.setupUW`. See [azPhaseCorrect.md](azPhaseCorrect.md) and
   [setupisceuw.md](setupisceuw.md).
7. **Copy the product out and clean up.** `copyProduct` copies the final
   (low-volume) product to `outputRoot/track-N/`, then the whole per-pair scratch
   tree is deleted unless `--debug`/`--keepScratch`.

Each stage is wrapped in `timed()`, so the log carries `[TIME] <stage>: <s>`
lines plus a `[TIME] TOTAL pipeline (excl. DEM)`.

---

## Scratch layout and disk strategy

The heavy ISCE files (~2 GB/pair) live **only** on the machine-local scratch disk
to avoid NFS contention; only the final ~150 MB product is written to the shared
output volume.

```
<scratchDir[host]>/<track>-<orbit1>-<orbit2>-<refStart>/
├── isce/                       topsApp.xml, dem.wgs84, SAFE symlinks, merged/, ...
├── gimp/track-<N>/             GrIMP product built by setupUW
└── intermediate/<o1>-<o2>/     convert_isce intermediate + simPhase, topophase.cor[.ion]
```

The scratch tag includes the reference **acquisition start time** (`frameTag`),
not just the orbits: frames of one datatake share absolute orbit numbers, so
without it two frames of the same pass would collide in scratch.

Outputs that survive cleanup:

```
<outputRoot>/track-<N>/<orbit1>_<frame1>/                 final GrIMP product
<outputRoot>/logs/<YYYY-MM-DD>/<track>.<orbit1>.<frame>.<runStamp>.log
<outputRoot>/logs/<YYYY-MM-DD>/Fail.<track>.<orbit1>.<frame>.<runStamp>
<outputRoot>/logs/exitStatus.log        one line per --queue run, all machines
```

Logs are grouped by the **run** date, so one `logs/YYYY-MM-DD/` directory holds
everything processed that day across all tracks and frames. `<runStamp>` is the
run start time (`YYYYMMDDTHHMMSS`), so repeat runs of the same pair never
collide and the fail records accumulate as a history rather than overwriting.

`<frame>` is the burst frame computed up front from the manifest's ascending-node
time — the same calculation as `findS1InsarPairs` and
`isce2grimp.convert_isce.get_frame_number` — so log names line up with the
`<orbit1>_<frame1>` product dirs. If the manifest can't be read it falls back to
`0` for naming only (non-fatal).

`copyProduct` drops symlinks (the back-link to the intermediate dir) so the
copied product is self-contained, and stages through `<dest>.new` /`<dest>.old`
renames so a re-run over an existing product works even on NFS dirs holding
`.nfs*` leftovers.

The log is opened fresh (`w`) per run and written under
`outputRoot/logs/<runDate>/`, so it survives the scratch cleanup; ISCE command
output is teed into the same file.

---

## Precise-orbit check

`checkOrbits` runs before anything is created, and requires that `orbitDir`
contains a POEORB `.EOF` for the right mission whose validity window spans each
acquisition:

```
S1A_OPER_AUX_POEORB_OPOD_<produced>_V<validStart>_<validStop>.EOF
```

A pair with a missing orbit fails immediately rather than an hour into
`topsApp`. The matched file is logged for each acquisition.

POEORB production lags the acquisition by about 20 days, which is why
`findS1InsarPairs --minAge` (default 25 days) keeps too-young acquisitions out
of the queue in the first place. That flag is an *age heuristic*; `checkOrbits`
is the actual check, and catches gaps in the archive that age alone wouldn't.

---

## Failure reporting

On any failure the program logs the traceback, keeps the scratch tree, writes a
fail record beside the run log, and **exits with a meaningful code**:

```
exitCode: 2
error: isceShell: command failed (exit 2): gdalwarp -q -overwrite ...
log: <outputRoot>/logs/2026-08-20/25.8348.571.20260820T165821.log
```

The exit code is the failing ISCE command's own code where there is one —
`isceenv.isceShell` raises `IsceCommandError` (a `RuntimeError` subclass
carrying `returncode`) — otherwise 1.

Three subtleties this handles:

- `azPhaseCorrect` and `setupisceuw` abort via `u.myerror`, which calls
  `sys.exit()` and raises `SystemExit` — **not** an `Exception`. The handler
  catches `(Exception, SystemExit, KeyboardInterrupt)` so those aborts still
  produce a fail record.
- That bare `SystemExit` carries no message (myerror prints to the console), so
  `failMessage` falls back to the source line that aborted, skipping `myerror`
  itself:
  `SystemExit at runS1interferogram.py:271: u.myerror('prepareIsceDem: need ...`
- `installSignalHandlers` turns SIGHUP and SIGTERM into `KilledBySignal`, whose
  `returncode` is `128 + signal`, so a terminated run still records why:
  `exitCode: 143 / error: killed by SIGTERM (signal 15)`. Ctrl-C keeps its
  normal behaviour and is caught as `KeyboardInterrupt` (130).

### What a *missing* record means

SIGKILL and SIGSEGV cannot be intercepted — the process runs none of its own
code — so those leave **no** `Fail.*` record at all. With the handlers above in
place that absence is now diagnostic rather than ambiguous, and the queue's
`exitStatus.log` says which one it was:

| symptom | cause |
|---|---|
| `Fail.*` naming a signal | SIGHUP / SIGTERM (terminal closed, killed) |
| no `Fail.*`, `rc=137` | SIGKILL — usually the OOM killer |
| no `Fail.*`, `rc=139` | SIGSEGV — crash in a native library (GDAL, numpy) |
| no `Fail.*`, no `rc` | run was not launched through a `--queue` line |
| `rc` 64–70, no log at all | the `pullASF` ahead of it failed, so the run never started — see [pullASF](../../asfSearchAndDownload/Documents/pullASF.md) for the code |

A failure during setup — an unresolvable SAFE, or a track that cannot be
determined — now exits **1**. It previously exited 0, because `u.myerror` calls a
bare `sys.exit()`; those checks run before the log dir exists, so they report to
the console and exit non-zero rather than writing a `Fail.*` record whose name
would be built from the very values that failed to resolve.

The handlers are installed just before the main `try`, *after* module import.
Import takes ~5 s (rasterio, rioxarray, sarfunc, GDAL), so a signal inside that
first few seconds still goes unrecorded.

> Not to be confused with `azPhaseCorrect.logFail`, which drops an unrelated
> `fail.<pid>` marker into the *current directory* (i.e. inside scratch, removed
> with it). The `Fail.*` records described here are the durable ones.

---

## The ISCE DEM

`topsApp` needs a lat/lon (EPSG:4326) DEM, but the project's `demTiff` is the
polar-stereo GrIMP DEM. `prepareIsceDem`:

- uses `dem` from the project if it is set, not `TBD`, and exists; otherwise
- computes the frame footprint from the reference `manifest.safe`
  `gml:coordinates` (`frameBbox`, +0.3° margin) and crops/reprojects `demTiff`
  into `isce/dem.wgs84` with `gdalwarp -t_srs EPSG:4326 -r bilinear -of ISCE` at
  `isce.demPosting` (default 0.000277°), then runs `gdal2isce_xml.py` on it.

Deriving the DEM per frame is the normal path — it keeps the project
self-contained and avoids maintaining a stitched `.dem.wgs84` per region.

> **Gotcha:** each `gdalwarp`/`gdal2isce_xml.py` is issued as its own
> `isceenv.isceShell` call. A single `conda run` only dispatches the *first*
> command of a `;`-chain into the env, so chaining them silently runs the second
> one in base.

---

## `topsApp.xml`

`writeTopsAppXML` starts from the `isce2grimp` template (`data/template.yml`,
read with `dinosar.read_yaml_template`) and overrides:

| topsinsar key | Source |
|---|---|
| `demfilename` | DEM from `prepareIsceDem` |
| `azimuthlooks` / `rangelooks` | `isce.azimuthlooks` / `isce.rangelooks` |
| `unwrappername` | `isce.unwrappername` |
| `doionospherecorrection` | `isce.doionospherecorrection` |
| `reference|secondary.safe` | linked zip names in the scratch ISCE dir |
| `reference|secondary.orbit directory` | `orbitDir` |
| `reference|secondary.polarization` | `isce.polarization` (default `hh`) |
| `reference|secondary.output directory` | `referencedir` / `secondarydir` |

Anything not listed keeps the template's value.

---

## `uwproject.yaml` keys

Read by this program (loader: [`uwproject.py`](../uwproject.py)):

| Key | Required | Description |
|---|---|---|
| `scratchDir` | yes | `{hostname: path}`. Aborts if the current host has no entry, or the path does not exist. |
| `threads` | no | `{hostname: cpus}` override for `--cpus`. |
| `dataDir` | yes | Directory holding the input SAFE zips. |
| `archiveDir` | no | Archive root with `YYYY-MM/` subdirs; searched second, so cross-month pairs resolve. |
| `orbitDir` | yes | Sentinel-1 precise-orbit (OPOD) directory, e.g. `/Volumes/insar9/ian/Data/SentinelGreenland/OPOD`. Checked by `checkOrbits` and written into `topsApp.xml` as `orbit directory`. |
| `outputRoot` | yes | Where final products and `logs/` are written. |
| `demTiff` | yes* | Polar-stereo DEM used to derive the ISCE DEM (and for phase simulation). |
| `dem` | no | Explicit ISCE-format `.dem.wgs84`; skips the per-frame DEM build. *Required only if `demTiff` is absent.* |
| `isce` | no | topsApp parameters: `azimuthlooks`, `rangelooks`, `unwrappername`, `doionospherecorrection`, `polarization`, `demPosting`, `cpus`. |
| `isceEnv` | no | Conda env name for ISCE commands (default `isce2grimp`). |
| `epsg`, `velMap`, `icemask`, `icemaskTiff`, `region` | yes† | Self-contained region data for phase simulation and masking. |
| `regionFile` | no | † Path to an existing `sarfunc` region yaml, used instead of the self-contained keys. |

Also present in the same file but read by the companion tools:
`archiveDir`, `temporalBaseline`, `excludeMonths`, `phasePairRanges`,
`asfFrameLookup` (`findS1InsarPairs`); `tiepointFile`, `region`/`regionFile`
(`setupS1PhaseTracks`).

---

## Failure modes worth knowing

- **No `scratchDir` entry for the host** — hard error naming the hosts that do
  have entries. Add the machine to the map rather than pointing it at NFS.
- **Missing precise orbit** — aborts up front; see the orbit-check section.
- **`isceShell` non-zero exit raises** (`set -o pipefail` keeps the `tee` pipe
  from masking it), so a failed `topsApp` aborts the pipeline instead of
  silently continuing into `azPhaseCorrect` with a half-built merge.
- **SAFE not found** — `resolveSafe` reports all three locations it tried.
- **Track undeterminable** — a corrupt/truncated zip makes `trackFromManifest`
  return `None`; pass `--track`.
- **Any failure** leaves the scratch tree in place and logs its path; the run
  log and `Fail.*` record under `outputRoot/logs/<runDate>/` are what survive
  cleanup.

---

## Related

- [`findS1InsarPairs`](../findS1InsarPairs.py) — lists candidate pairs at the
  project temporal baseline; its `image1 image2` columns are this program's two
  positional arguments. Shares `uwproject.yaml` and the `--project` convention.
- [`setupS1PhaseTracks`](setupS1PhaseTracks.md) — instantiates the per-track
  `tie_plan_header` for the phase products this program writes.
- [`azPhaseCorrect`](azPhaseCorrect.md) — the burst azimuth correction and
  simulated-phase unwrap stage (step 5); also usable standalone on an existing
  ISCE directory.
- [`setupisceuw`](setupisceuw.md) — the intermediate → GrIMP frame conversion at
  the end of that stage.
- [`isceenv.py`](../isceenv.py) — the `conda run` dispatch into the ISCE env.
- [`uwproject.py`](../uwproject.py) — project-yaml loader and host/region helpers.
