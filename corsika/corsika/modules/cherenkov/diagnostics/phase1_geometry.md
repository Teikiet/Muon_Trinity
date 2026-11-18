# Phase 1 — tracer geometry: stage-local kappa and Hermite closure

Status: **shipped** (both tracers). Absorption, projection modes, routing and the
table's contents are unchanged. This note records what was measured, so the numbers
can be re-derived rather than trusted.

## What was wrong

`rayTraceSphericalToPlane()` and `rayTraceSpherical()` were already classical RK4 on
`(position, direction)`, but the ray equation's coefficient

```
du/ds = (n'/n)(r_hat - (u.r_hat)u)          kappa = (dn/dr)/n
```

was computed **once per step, at the step's start radius**, and the four RK4 stages
reused it. That is a left-endpoint rule for `int kappa ds`, i.e. first order in
`ds/H`, not the 4th order the RK4 promises. A second, independent error was the plane
crossing: the fraction was exact for the (affine) plane distance, but the position was
interpolated **linearly along the chord** between the two step ends, leaving the
chord-to-arc sagitta (`c^2/(8*rho)`, `rho ~ 5e6 m`) in the impact point.

Measured on the frozen telescope acceptance set (48 pairs, truth = the same RK4 at
ds = 1 m), before the fix:

```
max |dL|/L            = 4.117e-07
max impact error      = 3.695e-01 m      <-- 370x the 1 mm target
max step-to-step jump = 4.819e-08
```

## What was done

1. **Stage-local kappa** (both tracers). `kappaAt(r)` is evaluated inside `rayDeriv`,
   at the radius of the stage that needs it. With a separable index it is two O(1)
   lookups in `RefractiveIndexTable` (`N`, `dN/dz`) and **no** `refractiveIndexFunc_`
   call per step; without one it is two function calls per stage via
   `(ln n(h+1) - ln n(h-1)) / 2`, which is `d(ln n)/dr` exactly.
2. **Cubic Hermite closure** (`HermiteStep`). The step's endpoint tangents are already
   in the state, so the position, the arrival direction and the path fractions inside
   the step are taken from the cubic through `(p0, ds*u0, p1, ds*u1)`. The crossing is
   then solved by Newton (started from the chord fraction, 6 iterations, clamped, with
   the chord value as fallback). This also removes the end-of-step direction error of
   `~(1-alpha)*ds*kappa` (~8e-5 rad at 400 m) from the recorded arrival direction and
   from the `xi` used by the absorption table.
3. **4th-order group path.** `Sigma = int N ds` is carried as an extra RK4 quadrature
   on the same stages and the optical path is `L + f_group*Sigma`, instead of the
   left-endpoint `n_group(step start)*ds` (bias ~0.4 m at 400 m steps, ~1.2 ns).
4. **`RefractiveIndexTable` caches its cell log-slopes** (`B_`, 96 kB for the
   12 000-node production grid), so `N`, `A` and `dN/dz` no longer recompute
   `ln(N_{i+1}/N_i)` per call. This is what makes (1) cheaper than the frozen kappa it
   replaced, in the production (table) path.

The finite-difference branch of (1) is kept verbatim for non-separable indices and for
runs without a table, but it now costs 2 function calls per stage instead of 2 per
step. Measured on the no-table path (`testCherenkovCurvedTracerGrid`):
86.34 us/call vertical (was 83 us), 264.29 us/call grazing.


## Root cause, isolated

2-D vertical section, `n = 1 + 3.27e-4 exp(-h/8 km)`, h_em = 20 km, d = 200 km
(the worst row of the frozen set), reference = the same scheme at ds = 0.25 m:

```
scheme          dz_frozen[mm]   dz_stage_local[mm]
fixed 1200 m        768.6               7.2
fixed  400 m        254.5               0.041
fixed  100 m         63.7               0.043
controller          255.4               0.255
```

The frozen-kappa column reproduces the 0.37 m seen in production; the stage-local
column is three orders inside the 1 mm target.

## Validation (all after the fix)

`testCherenkovPlaneTracerAcceptance`, frozen set, truth ds = 1 m:

```
max |dL|/L            = 1.701e-12      (gate 1e-9)
max impact error      = 6.184e-07 m    (gate 1e-3 m)
max step-to-step jump = 3.231e-11      (gate 1e-9)
in-FoV 27: pass 17, flagged 0, fail 0  (unchanged)
```

`CHER_TRUTH_STEP_M=<m>` refines the truth (added to the test) so its own convergence
can be shown; the default 1 m keeps archived numbers reproducible.

`testCherenkovTracerAcceptance` (sphere tracer): **now passes** (11 assertions). Before,
it failed at `worstOne <= 1e-9` with 0.968, and the n = 1 probe read 9.682e-01; it now
reads 3.792e-11, and the tau step-halving check is 1.156e-03.

Second test case in the same file (curved table, telescope mode):
```
max |dL|/L      2.833e-08 -> 9.742e-09
max impact      0.20 m    -> 0.01 m
maxTauRel       5.613e-02 -> 5.598e-02   (still above its 5e-2 gate; see below)
```

Other suites: `testCherenkovAbsorptionChord` 21 assertions, `testCherenkovChordRoute`
17, `testCherenkovCurvedTracerGrid` 5, `testCherenkovE2E` 62 -- all pass.

## Production, seed 600 (10^4 GeV muon, 15 m dish, 2944 m, 43 km emission)

Straight line (A) is untouched, as expected: 2079 hits, Cherenkov 9.07 s
(production 5.46 + bunch loop 3.10), run 46.0 s -- the same numbers as the archived
`m_A_s600` (2079, 9.18 s, 47.75 s).

Curved (C, legacy tracer; E, table):

| | C | E | archived C (old binary) |
|---|---|---|---|
| hits | 2075 | 2073 | 2074 |
| gated (traced) photons | 199 315 (1.51 % of 13.18 M bunches) | 199 315 | -- |
| table hits | -- | 97 248 | table not used (Telescope was bypassed) |
| table rejections | -- | 102 067 (upgoing 65 033, grazing 5 256, belowTable 802) | -- |
| legacy fallback | 66 223 ok / 133 092 fail | 61 078 ok / 40 989 fail | -- |
| bunch loop | 101.8 s | 80.8 s | 3.4 s |

The table route confirms the legacy tracer to +-3 hits, which is the first time the
table has been exercised in Telescope mode (the archived six-configuration comparison
ran a binary whose eligibility predicate excluded Telescope).

## Not fixed here

* `maxTauRel = 5.598e-02` against its 5e-2 gate in the table/telescope test case. This
  is the absorption-tau disagreement documented in `b_table/README.md` (25.4b); it moved
  by 0.015 pp, i.e. it is not a geometry effect.
* `CurvedRayTable::buildColumn` still freezes its log-slope per **substep**. Substeps
  are 4-32 per 10-100 m cell, so the effect is far below the table's mm-level budget,
  but it is the same first-order-in-dh/H term and is left as is deliberately.
* `testCherenkovTable3km` fails its own criteria (11 % of cells within 1 %). That test
  evaluates the rejected coarse-theta absorption-table design of `b_table` sections
  9.5/12; it is not a tracer test and its failure is the measurement it was written to
  make.

## Phase 2 -- cost: exact upgoing cut, log rate-limiting, and a corrected diagnosis

### 2.1 The cut

In a spherically symmetric medium whose `n*r` increases monotonically (the condition
`CurvedRayTable`'s constructor verifies), a photon emitted with a non-negative radial
direction cosine has a monotonically increasing radius.  Its crossing `X` of the
(infinite) dish plane then satisfies

```
|X| >= r_emit   =>   |X - C| >= r_emit - r_obs   =>   sqrt(x^2+y^2) >= r_emit - r_obs
```

because the dish coordinates are the plane components of `X - C`.  So an emission more
than one aperture above the dish cannot produce a hit, and the photon is refused before
tracing.  This is a proof, not a tunable gate.  `rayTraceSpherical` gets the equivalent
sphere-target cut.  `CHER_NO_UPGOING_CUT=1` disables it for A/B (same pattern as
`CHER_AB_PLANE_BACK`).

### 2.2 Log rate-limiting

The tracers emitted one `WARN` per escaping photon -- **40 032 lines on one event**.
They now always increment `escapedAtmosphere_`/`maxStepsReached_` and log the first
three; the totals appear once as
`Cherenkov tracer cost: upgoing pre-cut N (of M gated), escaped atmosphere E, max-steps S`.
Same events: **241 log lines** instead of 40 032+, with the same 40 028 escapes counted.

### 2.3 A/B, same binary, same seed (yields must be identical)

Moving dish (`--telescope-position-correction 1`, the benchmark of the sections above):

```
              hits   bunchloop   upgoing pre-cut   escaped
cut ON        2075   101.882 s            1663     40028
cut OFF       2075   100.527 s               0     40028
```

Static dish (`--telescope-position-correction 0`):

```
              hits   bunchloop   upgoing pre-cut   escaped
cut ON          45    73.528 s               2         0
cut OFF         45    74.031 s               0         2
```

Hits are identical in both pairs: the cut loses nothing, as proved.  The measured
delta is -0.503 s (static) and +1.355 s (moving, i.e. the cut run was slower there),
which brackets zero: both are inside the +-1.4 % run-to-run spread of these 75-100 s
runs.  The cut removes 1 663 (moving) / 2 (static) traces of work, which is too
little to resolve against that spread -- it was never going to be a speed-up, it is
a proof-based refusal that costs nothing.

### 2.4 The diagnosis in the Phase-1 section was WRONG, and this is the correction

The Phase-1 note predicted that refusing the 65 033 `Upgoing` table rejections would
recover a large part of the 80-100 s.  It does not, for two reasons that the A/B makes
visible:

1. **The table judges `Upgoing` against the fixed `groundZ_m_`, but the hit test uses
the telescope centre**, and `--telescope-position-correction 1` moves that centre onto
the shower.  With the dish moved to 43 km, `r_emit - r_obs ~ 15 m`, so the proof does
not apply and only 1663 photons are cut.  With a static dish the proof applies, but the
gated population is then descending by construction (45 hits, 2 cut).
2. **The rejected rays are not waste.**  Of the 102 067 table rejections in E, 61 078
succeed in the legacy plane tracer -- they really do cross the plane, tens of km above
the dish, where the aperture then discards them.  Their trace is what establishes that.

Where the time actually goes: 199 315 traces at ~0.5 ms each.  A trace costs
`steps x ~0.25-0.5 us` and this benchmark's paths are 200-700 km (a horizontal 43 km
muon), i.e. ~1000-2000 steps.  For a representative vertical shower the same tracer is
20-80 us per photon (`testCherenkovCurvedTracerGrid`: 86 us vertical, 264 us grazing),
so this configuration is 10-20x more expensive *per trace* than the typical case, and
the cost is path length, not routing.

### 2.5 What Phase 2 therefore recommends for cost

* Keep both changes: the cut is exact and free; the log fix removes a 40 000-line flood.
* Do **not** expect a big saving from routing: the table cannot serve these rays (its
  contract is the sphere crossing), and the exact trace is what proves they miss.
* The remaining levers are (a) a plane-parametrised table for Telescope mode -- design D
  of `b_table` section 26, whose phi axis was refuted in 26.4, so this is research, not
  a fix; (b) accepting the cost for horizontal configurations, which are 10-20x more
expensive per photon than vertical ones.

## Phase 3 -- adaptive DOPRI5 5(4): measured, adopted in the application

### 3.1 What was added

* `corsika/modules/cherenkov/Dopri5.hpp` -- allocation-free Dormand-Prince 5(4) with FSAL,
  the quartic dense output (Hairer's P matrix, the one scipy's RK45 uses), per-component
  scaled error control, and the bounded controller `h_new = h*clip(0.2, 5, 0.9 e^-1/5)`
  with `e == 0` handled explicitly.  Scratch storage is fixed (`k[7][N]`), so the
  stepping loop allocates nothing.
* `TracerIntegrator { RK4, DOPRI5 }`, `setDopriTolerances(atol_pos[m], atol_dir,
  atol_opt[m])`, `setDopriStepLimits(min[m], max[m])`.  **RK4 stays the default, the
  reference and the fallback**; DOPRI5 is taken only when it is enabled AND the index is
  separable AND the refractive-index table is valid, because its Sigma component is
  scaled by f_group.
* `rayTracePlaneDopri` -- the same ODE, frame and tau ownership as the RK4 tracer.  The
  plane crossing is located on the DENSE OUTPUT with a bracketed solve, sampling the
  plane distance at t = 0, 1/4, 1/2, 3/4, 1 so a dip that returns to the same side is
  still seen; position, direction, group path and tau are all taken at that same t.
  The tau quadrature keeps its altitude-fixed 10 m / 50 m rule, now driven by the dense
  output, so a large step does not coarsen absorption.
* Counters: accepted/rejected steps, derivative evaluations, min/max step, capped steps
  (clipped by the maximum), crossing iterations, min-step failures, and which tracer
  served the ray.  Reported by the application once per run.
* `tests/modules/testCherenkovTracerBenchmark.cpp` -- 12 frozen rays (200/400/700 km
  horizontal, 10-40 km downward, 100/200 km near-tangent, two layer-crossing, one steep,
  one escaping) against a CONVERGED reference (DOPRI5 atol_pos = 1e-6 m, no cap).

### 3.2 Reference convergence (the criterion the first attempt got wrong)

```
reference self-check (100x tighter tolerance): max 0.000339 mm, inconsistent 0
```

A radial cap is **not** a step cap: at 0.5 m/step a steep ray needs ~3e4 steps and hits
RK4's 1e4 step limit, so the first "self-check" (RK4 at a 0.5 m radial cap) was
ill-posed and was replaced by a tolerance-convergence check on the reference itself.

### 3.3 Accuracy and cost on the frozen set (per-call microseconds)

```
RK4 (production caps): max |dP| 0.000341 mm, |dU| 1.15e-11, |dt| 0.0003 ns, 1790 us/set

maxstep   dopri_us  speedup  horiz_x  max_dP_mm     max_dU  max_dt_ns   steps  rej
   5000       445.7     4.02     4.39    0.002520   8.43e-11    0.0000    9300  160
  20000       196.5     9.11    12.73    0.118117   4.98e-09    0.0004    3640  300
  50000       162.1    11.05    16.73    0.123695   4.98e-09    0.0004    2700  640
```

All three caps are far inside the 1 mm impact, 1e-7 direction and 0.05 ns timing
limits.  Reproducibility, three consecutive runs of the same binary: 9.19 / 9.23 / 9.22x
at 20 km and 11.18 / 11.21 / 11.13x at 50 km (spread +-0.5 %), with the accuracy columns
bit-identical because the tracer is deterministic.  20 km was chosen (9.1x overall, 12.7x on the horizontal band): the 20 -> 50 km
step adds 20 % speed but takes in whole-atmosphere steps, and the crossing sampler is a
5-point scan, so the extra risk is not worth 20 %.  With the RK4 caps retained the same
integrator is 0.38x -- *slower* than RK4, because 7 evaluations per step at the same step
size is pure overhead.  That is the measured answer to "would adaptivity alone help".

### 3.4 Production, seed 600 (43 km muon, 15 m dish, 2944 m, moving dish)

Configuration C (curved legacy tracer, `--disable-cherenkov-curved-table`), same binary,
same seed, alternating runs:

```
                    hits   bunch loop   Cherenkov total   run time
tracer rk4          2075     102.378 s        108.467 s      145.0 s
tracer dopri5       2075      18.888 s         24.858 s       58.9 s
speedup                        5.4x             4.4x          2.5x
```

197 652 rays traced, 11.9 M accepted / 7.5 M rejected steps, 116.7 M derivative
evaluations, h in [1, 2e4] m, 1.6 M capped steps, 2.7 M crossing iterations,
**0 min-step failures**.

Configuration E (curved + CurvedRayTable, so 179 226 rays reach the exact tracer):
86.3 s -> 32.5 s (2.7x), bunch loop 80.3 s -> 26.5 s (3.0x), **hits 2073 identical**.

### 3.5 The recorded quantities, not just the counts (CSV A/B, 2075 rows)

```
rows                     2075 vs 2075  (1:1 pairing, no hit/miss change)
x_m       max |d|        0.01 mm
y_m       max |d|        3.58 mm      (2 rows > 1 mm of 2075; median 0, p99 0.01 mm)
z_m, wavelength, production height, all direction cosines, time, weight, n-1: equal to
                         the CSV's 6-digit print precision
optical path             max |d| 1 mm
```

Of the two rows above 1 mm, one is a **measured RK4 error**: re-running DOPRI5 at
atol_pos = 1e-5 m reproduces the default run exactly (0.000 mm) while RK4 sits 2.07 mm
away.  The other is ill-conditioned: RK4, DOPRI5-default and DOPRI5-tight all differ
(the tight run moves that impact by 10 mm), i.e. the coordinate has a large derivative
with respect to the trajectory there -- the same conditioning cliff established in
Phase 1, not an integrator defect.  Tightening the tolerance from 1e-3 to 1e-5 m cost
2.3 % of the bunch loop (18.888 -> 19.313 s).

### 3.6 Rejections: the next lever, measured

39 % of attempted steps are rejected (7.5 M of 19.4 M).  The cause is visible in the
model: `RefractiveIndexTable` stores a PIECEWISE-CONSTANT log-slope, so kappa jumps at
every 10 m node.  A 20 km step spans 2000 of them and the embedded error estimate sees
the node noise, not the smooth curvature.  Splitting the step at every node would cap h
at 10 m and destroy the win, so the honest options are (a) leave it (already 5.4x),
(b) evaluate kappa from a C1 local fit over the neighbouring nodes -- a model change
that needs its own acceptance study, or (c) tighten the atmosphere grid.  Each costs
7 evaluations per accepted step today (9.8 including the rejections) against RK4's 4.

### 3.7 Adoption decision (against the agreed criteria)

| criterion | result |
|---|---|
| no unexplained hit/miss changes | 2075 vs 2075 (C), 2073 vs 2073 (E); no mismatches on the frozen set |
| accuracy vs converged references | 0.12 mm impact (limit 1 mm), 5e-9 direction (1e-7), 0.0004 ns (0.05 ns) |
| no new failures / silent min-step exits | 0 min-step failures in 11.9 M accepted steps |
| >= 20 % faster on the expensive horizontal benchmark | 12.7x on that band, 9.1x on the whole frozen set |
| no material regression on typical showers | every ray in the set improved; the steep 5 km ray is the smallest gain |

**Adopted**: the application's `--cherenkov-tracer` default is now `dopri5` with
`--cherenkov-tracer-max-step 20000` and `atol pos/dir/optical = 1e-3 m / 1e-9 / 1e-4 m`;
`--cherenkov-tracer rk4` restores the previous behaviour exactly.  The MODULE default
stays `RK4`, so every acceptance test and every archived number in this note still
measure the RK4 reference -- the switch is explicit at the application level.

Scope, stated rather than implied: the DOPRI5 path serves the PLANE tracer (Telescope
mode and its legacy fallback).  `rayTraceSpherical` (Ground mode) still uses RK4, so
Ground-mode curved runs are unchanged in cost and accuracy.

Deviation from the specification, with its measurement: explicit splitting at
atmospheric-layer boundaries was NOT implemented, because the refractive-index model's
nodes are 10 m apart and splitting there would cap the step at 10 m and erase the whole
6-9x.  The embedded error estimate is what reacts to those jumps instead; the measured
cost of that choice is the 39 % rejection rate of section 3.6.
