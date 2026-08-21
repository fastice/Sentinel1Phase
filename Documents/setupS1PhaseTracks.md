# setupS1PhaseTracks

Instantiates the per-track `tie_plan_header` file for a Sentinel-1
unwrapped-interferogram (phase) project. Pared-down phase analogue of
[`s1setup.setupS1Tracks`](../../s1setup/s1setup/setupS1Tracks.py): the only mode
is `--copyFiles`, which drops a `tie_plan_header` (with `<TRACK>` substituted)
into each track's `tiepoints/` directory from the project's `templates/`
directory. The phase project has no velocity mosaics, so the `vel_thumb_plan` /
`vel_thumb_header_<range>` machinery, the `velocityStats` run modes, and the
`secondaryDirectories` cascade of `setupS1Tracks` are all omitted.

Run from the project root (the directory holding `uwproject.yaml`, `templates/`,
and the `track-*` subdirectories), or point `--project` at a `uwproject.yaml`
elsewhere — its directory is taken as the project root.

---

## Usage

```
setupS1PhaseTracks --copyFiles
setupS1PhaseTracks --copyFiles --tracks track-25 track-90
setupS1PhaseTracks --copyFiles --overWrite
setupS1PhaseTracks --copyFiles --project /path/to/uwproject.yaml
```

---

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--project PATH` | `uwproject.yaml` | Project yaml; **its directory is the project root** (`templates/`, `track-*`). |
| `--tracks track-N …` | all `track-*` | Restrict to these track dirs (validated to exist). Numeric-sorted on the track integer, same as `setupS1Tracks`. |
| `--copyFiles` | — | Instantiate `tie_plan_header` into each track's `tiepoints/`, skipping files that already exist. |
| `--runRefresh` | — | Run `refreshties.py -phase` (maketies + makeframetie → `tie_script`) across the selected tracks and return. Honors `--year`, `--nThreads`, `--overWrite`. |
| `--year YYYY …` | 2015…current | Years to refresh (passed as positional years to `refreshties.py`). |
| `--nThreads N` | 4 | Concurrent tracks passed to `refreshties -nThreads`. |
| `--overWrite` | False | With `--copyFiles`, overwrite an existing `tie_plan_header`; with `--runRefresh`, pass `--overWrite` to `refreshties` (rerun existing tie products). |

One of `--copyFiles` / `--runRefresh` is required; without either the program
errors ("nothing to do").

---

## What it does

For each selected track directory:

1. Ensures `<track>/tiepoints/` exists (created if absent).
2. Reads `templates/tie_plan_header` and substitutes placeholders (see below).
3. Writes `<track>/tiepoints/tie_plan_header` — copy-if-absent unless
   `--overWrite`. A summary line reports how many were skipped vs written.

The run is idempotent: a second `--copyFiles` with no `--overWrite` skips every
existing header and writes nothing.

---

## `--runRefresh`

Runs `refreshties.py` in `-phase` mode from the project root, across the
selected tracks:

```
refreshties.py -phase [--overWrite ]-nThreads <N> -toRun="[track-25,track-90,...]" <years> -noPrompt
```

The `-phase` flag is what distinguishes this from
[`setupS1Tracks --runRefresh`](../../s1setup/Documents/setupS1Tracks.md) (the
legacy-baseline velocity flavor): it routes `maketies` to the phase
(unwrapped-interferogram) products. The run aborts if `refreshties` fails.

---

## Template substitutions

Same `<TRACK>` / `<DEM>` / `<TIEFILE>` convention as
`setupS1Tracks.applySubstitutions`. In practice the phase `tie_plan_header`
template only carries `<TRACK>` (the track dir stub); `<DEM>` and `<TIEFILE>` are
no-ops when the template doesn't use them.

| Placeholder | Filled from | Notes |
|---|---|---|
| `<TRACK>` | track dir integer (`track-25` → `25`) | e.g. `tie_dir = .../track-<TRACK>/tiepoints`. |
| `<DEM>` | region yaml `dem` (`uwproject.yaml` `region`/`regionFile`; falls back to `demTiff`) | Only the DEM is taken from the region file — runS1interferogram keeps its self-contained `demTiff`/`velMap`/`icemask`. |
| `<TIEFILE>` | `uwproject.yaml` → `tiepointFile` | Empty string if the key is absent. |

The template itself lives at `<project root>/templates/tie_plan_header`; see it
for the baked `DEM`, `default_nDays`, `extraties`, and `base_nlooks` values.

---

## Related

- [`s1setup.setupS1Tracks`](../../s1setup/s1setup/setupS1Tracks.py) — the fuller
  velocity-project orchestrator this is modeled on (`--copyFiles` plus
  `--runVelstatsregions` / `--runVelStats` and the secondary-directory cascade).
- [`findS1InsarPairs`](../findS1InsarPairs.py) — lists the interferometric pairs
  for a track; shares the `uwproject.yaml` config and `--project` convention.
- [`runS1interferogram`](../runS1interferogram.py) — produces the per-pair GrIMP
  phase product the tie plan operates on.
