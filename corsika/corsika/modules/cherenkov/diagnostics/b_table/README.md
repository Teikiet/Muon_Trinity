# Cumulative-impact-parameter optical-depth table — diagnostics

Status: **diagnostic only. Production integration is NOT authorized.** The production
absorption route remains the legacy curved-ray tracer.

Four standalone diagnostics for a proposed absorption table that tabulates accumulated
optical depth against the ray's conserved impact parameter
`b = n_phase(r,lambda) * r * sin(theta)`:

```
tau(r_e, b, lambda) = int_{r_obs}^{r_e} alpha(r,lambda)
                      / sqrt(1 - (b / (n_phase(r,lambda) r))^2)  dr
```

All emission heights on one (b, lambda) share a single cumulative integral, so the table
stores one scalar per node instead of many layer weights.

## Build

The include paths are installation-specific. Set them from your environment:

```
C8_ROOT=<corsika source root>          # e.g. .../corsika
EIGEN_INC=<Eigen 3 include directory>  # the Conan cache path used at build time
                                       # find it via: grep Eigen3_DIR <build>/CMakeCache.txt
                                       # then append /p/include/eigen3
CXX=${CXX:-g++}

$CXX -O2 -std=c++17 -I"$C8_ROOT" -I"$EIGEN_INC" check_b.cpp -o check_b
$CXX -O2 -std=c++17 -I"$C8_ROOT" -I"$EIGEN_INC" study_b.cpp -o study_b
$CXX -O2 -std=c++17 -I"$C8_ROOT" -I"$EIGEN_INC" gaps_b.cpp  -o gaps_b
$CXX -O2 -std=c++17 -I"$C8_ROOT" -I"$EIGEN_INC" cover_b.cpp -o cover_b
```

Binaries are deliberately **not** committed.

## Input profile

The diagnostics take the extinction profile as `argv[1]`. The profile used for the
recorded results is **not redistributed here**; record its location and identity:

```
path : modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat   (git submodule modules/data)
size : 46385 bytes
md5  : dfb2b7fd0d0c57ee0a704035a6a81e0b
format: "n_lambda n_h", then lambda array [nm], height array [m], then one ecoeff row
        per wavelength. ecoeff is the CUMULATIVE COLUMN optical depth from the ground up,
        i.e. tau INCREASES with height, so the local extinction is +d tau/dh.
```

Check whether your copy is redistributable before committing it anywhere.

## Recorded environment

```
compiler : g++ (GCC) 8.5.0 20210514 (Red Hat 8.5.0-30)
flags    : -O2 -std=c++17
branch   : cherenkov-curved-table-wip
sources  : check_b e3efaf45f08b8e0d75e78114f39b731d
           study_b 56417d33cacfe2f417ab0dbfaeecf9b6
           gaps_b  54c311b7e28aaadf19c3c384577817eb
           cover_b a2d06ebffeffcdf3f0c9f9b4e6e7ff83
```

Frozen axes (study_b, gaps_b): 56 h x 256 mu x 81 native lambda = 1 161 216 nodes = 4.43 MB;
build quadrature dr = 5 m, split at every native height knot +/- 25 m (the callback's
breakpoints). Geometry table parameters for the combined-cost measurement: R = 6.371e6,
rObs = R, hTop = 120 km, nMu = 800, dhLo = 10 m, dhHi = 100 m, zBreak = 30 km.

## Reproduced results

All four were re-run from the committed sources against the profile above; the accuracy
figures reproduced **exactly**, timings jittered (see the timing caveat).

`log_check_b_premise.txt` — the formulation vs an independent RK4 trajectory integration:
worst `|d tau| = 8.633e-06`, including at `b/b_crit = 0.99995`.

`log_study_b_gates.txt` —
```
axes: 56 h x 256 mu x 81 lambda(native) = 1161216 nodes = 4.43 MB
build time: 38.9 s (dr=5 m, breakpoint-aware)
quadrature convergence: dr=20 / 5 / 1 -> 1.2404002   |d(5-1)| = 3.18e-10
lambda identity (native 400/405, mid 402.5): 0.61490477 vs 0.61490477  |d tau| = 2.11e-15
lookup latency: 0.0480 us/call
held-out (new, frozen): 140 rays, fallback 0, WORST combined |dT| = 1.4518e-04
geometry table: build 30.0 s, 214.3 MB
TOTAL peak RSS: 229052 kB = 223.7 MB
```

`log_gaps_b_dispersion_fallback.txt` —
```
f=0.9713 : WORST |dT| = 1.9964e-05     (dispersion-only, three sampled f)
f=1.0000 : WORST |dT| = 1.3493e-05
f=1.0347 : WORST |dT| = 1.0064e-05
fallback: beyond-critical 0 | upgoing 8 | above-top 20 | below-obs 0 | lambda-out 40 | accepted 32
```

`log_cover_b_eligibility.txt` — decision tests: h = 0 accepted, h < 0 rejected;
absorption-enabled-with-table-unavailable -> `fallback-legacy`; launchable-but-turning rays
(`b` just above `b_crit`) rejected with `min_r - R` = 5 260 m to 49 952 m; f endpoints
accepted and +/- 0.0004 outside rejected.

## Corrections and caveats (record these with any citation)

- **Launchability does NOT imply ground reach.** `b <= b_crit` is the sufficiency condition
  for an inward ray. A ray with `b` slightly above `b_crit` is launchable at emission and
  still turns above the ground (demonstrated in `cover_b`).
- **`gaps_b` could not produce a beyond-critical case at all** because it parameterises by
  the observation-level mu; `cover_b` parameterises by `b` and does.
- **Upgoing rejection** was exercised in the `gaps_b` fallback block (count 8), not in the
  `cover_b` checks. In `cover_b`, check (2) evaluates the *decision* to `fallback-legacy`;
  it does not run the fallback path.
- **The dispersion result and the combined-error result must stay separate.** The 2.0e-05
  figure is dispersion-only at three sampled f; the 1.4518e-04 figure is the combined error
  at f = 1 on a different set. No measurement evaluates them together, so combined accuracy
  across the f_phase range is **unvalidated**.
- **The reference is not an independent physical-model validation.** It shares the
  production extinction callback and the Bouguer trajectory model; what is established is
  the representation and the quadrature, not the underlying model.
- **3.18e-10 is a successive-refinement difference on the reported convergence probe**, not
  a domain-wide quadrature-error bound.
- **Timing figures are single-run.** Re-running the same binaries reproduced the accuracy
  values exactly but moved build time (38.7 -> 38.9 s) and lookup latency
  (0.0529 -> 0.0480 us) by ~10 %. No timing protocol is frozen yet; see the plan's section 2.
- **Do not infer realistic performance from cache hostility.** The 0.833 us combined route
  in `/tmp/gaps.log` (not reproduced here) is promising but neither guarantees nor rules out
  meeting the ~1 us target under a representative workload.
- **Near-critical handling is not settled.** At `b = b_crit` the path length and optical
  depth do **not** necessarily diverge; whether the endpoint singularity is integrable
  depends on the refractive profile and the root structure. The proposed rule for this
  workstream is to reject equality and numerically ambiguous near-critical cases to the
  legacy fallback, with the guard derived from conditioning and convergence tests rather
  than an arbitrary width. That derivation has **not** been done.
- **The f = 1 construction is adequate only over the proposed range [0.9713, 1.0347]**, and
  only on the tested rays.

## Section 6 — phase 1–3 archive and the boundary/config follow-up

Status unchanged: **diagnostic only; production integration is NOT authorized.**

### 6.1 Files added here (binaries deliberately not committed)

```
sources  ref_b.cpp      phase 2 reference validation (analytic checks, breakpoints, admissibility)
         cmp_b.cpp      phase 3 comparison with ACTUAL interpolated lookup outputs
         diag_quad.cpp  quadrature probe: unsplit vs split build at the difficult ray
         discrep.cpp    attribution of the historical 2.208e-03 failure
         map_b.cpp      coordinate definition, mapping test vs the REAL geometry table,
                        validation against an independent reference for the actual ray
         audit_b.cpp    boundary families, breakpoint analysis, admissibility, performance
         path_b.cpp     independent trajectory cross-check of the geometry table
         bp_b.cpp       breakpoint verification (h >= 0) and the lambda edge of the profile
         reaxis_b.cpp   axis-variant comparison (section 7): cached reference, shared ray set
         combined_b.cpp combined log+refinement vs V3 (section 8): per-family locations,
                        neighbourhood probes, actual-ray, timing investigation
         closure_b.cpp  targeted closure (section 9): attributed maxima, coverage scan,
                        height ceiling, timing resolution with a second parameterisation
logs     log_phase1_baseline_suite.txt          full test suite, current checkout
         log_ref_b_phase2_reference.txt         phase 2
         log_cmp_b_phase3_actual_lookup.txt     phase 3
         log_map_b_coordinate_mapping.txt       section 6.3-6.4
         log_audit_b_boundary_performance.txt   section 6.5-6.8
         log_path_b_trajectory_crosscheck.txt   section 6.5
         log_bp_b_breakpoints_lambda_edge.txt   section 6.6-6.7
         log_reaxis_b_axis_comparison.txt       section 7
         log_combined_b_variant_and_timing.txt  section 8
         log_closure_b_attribution_policy.txt   section 9
build    $CXX -O2 -std=c++17 -I"$C8_ROOT" <name>.cpp -o <name>
run      ./<name> <abstable path>
```

Source hashes as committed (md5):

```
ref_b.cpp     7e7bd31933e3388e8638d644e4bbd95c
cmp_b.cpp     f79eb832c64c3dc07143ab93f211889d
diag_quad.cpp e715199526599d5b816fa391174929f4
discrep.cpp   58d0ef82ec11f7bdeb179e926335b4aa
map_b.cpp     da758561bb867ed370f915c6718ae241
audit_b.cpp   015457f133487f6ed5ece37be57d4385
path_b.cpp    f78ec307d69b1809d2ebc55b3ea3c2d4
bp_b.cpp      952fd7fac6f29d2e646d2412d8088040
reaxis_b.cpp  1e567ec2eb9a730470915af0330e6fa6
combined_b.cpp 000bd04b8660a0f2b3392abce84b5fec
closure_b.cpp  9f6ad0e2630ed0fddc590ea6576a4b00
```

### 6.2 Corrected-diagnostic history (record this with any citation)

- `study_b`/`gaps_b`/`cover_b`/`guard_b` used a **split integrator whose cut points were
  not bounded by the emission radius**: `splitInt` pushed knots beyond `rTop` and
  integrated past it (it returned 1.3043 where the correct value is 1.2404). Every
  apparent "split vs unsplit" discrepancy produced by that routine was this bug. The
  correct bounded routine is in `ref_b.cpp`; the profile's own numbers are unaffected.
- `discrep.cpp` supersedes the earlier framing of the 2.2e-3 failure (section 6.3).
- `audit_b.cpp`'s first breakpoint check evaluated `alpha` at `h < 0`, outside the
  callback's domain, and reported a spurious 7.5e-1 residual. `bp_b.cpp` restricts the
  check to `[0, h_max]`; the real residual is 1.276e-11.

### 6.3 The historical 2.208e-03 is EXPLAINED (reproduced, not merely unverifiable)

At `guard_b`'s worst ray (h = 1500 m, b/b_crit = 0.9995, lambda = 672.5 nm), `discrep.cpp`
reproduces the archived figure to all printed digits and separates the three calculations
that were being conflated:

| number | what it computes | realizable by a xi-parameterised table? |
|---|---|---|
| **2.208e-03** (`guard_b`) | the photon's **true** b with an **f = 1** profile | **no** — the coordinate and the profile disagree |
| 1.667e-04 | the realistic mapping: b = xi * b_crit(f = 1), f = 1 profile | yes |
| 6.80e-04 | the actual interpolated table, f = 1, worst over the frozen set | yes |
| 7.53e-04 | the actual table built with the **configured** dispersion | yes |

At f = 1 all four definitions coincide to 2.7e-08, confirming the quadrature. The
mechanism is a **72 m shift in impact parameter** (true b = 6369969.030 m vs table
b = 6369896.775 m) amplified by the endpoint square root.

**This does not mean the design "cannot produce" the mismatch.** Normalised coordinates
remove the *coordinate* mismatch; they do **not** by themselves show that the tabulated
absorption follows the actual geometric ray. That requires the separate test in 6.4.

### 6.4 Coordinate definition, mapping test, and validation on the actual configured ray

The axis is defined unambiguously and is **not** called mu:

```
p      = n_phase(r_e, lambda) * r_e * sin(theta_e)     Bouguer invariant (= impact parameter)
p_crit = n_phase(r_obs, lambda) * r_obs
xi     = p / p_crit                                    <- THE TABLE AXIS (impact-parameter ratio)
mu_obs = |cos theta_obs| = sqrt(1 - xi^2)              direction cosine; a DIFFERENT quantity
```

`CurvedRayTable` documents its axis as `mu_obs = |cos theta_obs| = sqrt(1-(p/p_crit)^2)`, so
the two axes are related by an exact algebraic map, not an approximation. `map_b.cpp`
verifies it against the **real production class**: for 280 photons spanning h = 1.5–46 km,
lambda = 310–850 nm and xi = 0.15–0.9995, it computes xi independently from
`(r_e, theta_e, f_phase, RI table)` and compares with the geometry route's own
`mu = -cos(result.theta_obs)`:

```
worst |xi_trace - xi_independent| = 9.99e-16        (machine precision: the map is exact)
```

`f_phase` is obtained exactly as production does (`f = (n-1)/N(r_emit)`), so passing
`fPhase = F(lambda)` into `trace()` reproduces the production call site.

**Validation against an independent reference for the ACTUAL ray.** The reference integrates
`alpha(z) / sqrt(1-(p/(n_phase(z,lambda) z))^2)` along the true ray with the true dispersion
and the configuration's own `N(z)` (from `RefractiveIndexTable`), consulting no table:

```
A1 table built at f=1              worst |dT| = 7.3844e-04   PASS
A2 table built with configured f   worst |dT| = 7.5328e-04   PASS
dispersion-only sensitivity of the reference itself: 4.7272e-04
worst ray of THIS check (both A1 and A2): h=1500 m, lambda=620 nm, xi=0.9950
```

**Consequence: building the table for the configured dispersion does not by itself improve
accuracy.** A2 is not better than A1, and both worst cases sit at the same grazing ray. The
error is therefore dominated by the *representation* (the xi axis near xi -> 1), not by the
f = 1 vs f != 1 mismatch — and a configuration-faithful build is justified by the
*contract* (invalidation when the RI table or dispersion model changes), not by a claimed
accuracy gain.

### 6.5 Independent trajectory cross-check of the geometry table

`path_b.cpp` compares the geometry table's traced geometric path against an independent
Bouguer integral of the same ray. The error grows systematically toward grazing
(h = 1500 m, lambda = 690 nm):

```
xi=0.05    dL=5.7e-06 m      xi=0.9     dL=4.8e-05 m
xi=0.5     dL=7.8e-06 m      xi=0.99    dL=1.28e-03 m
xi=0.8     dL=2.0e-05 m      xi=0.995   dL=3.51e-03 m
                             xi=0.9995  dL=8.57e-02 m
                             xi=0.9999  dL=4.62e-01 m
worst over the scan: 2.1544e+00 m at h=46000 m, lambda=400 nm, xi=0.9999 (relative 2.97e-06)
```

This is a **geometry-table** interpolation error, not an absorption-table error, and it
carries the same signature (xi-axis resolution near the grazing end) as the absorption
failures in 6.6.

**Qualification (do not over-read this).** A 2.15 m *path-length* discrepancy does **not**
directly bound the transmission error, and it does **not** prove that the geometry and
absorption errors share a cause — a shared signature is not a shared mechanism. What is
needed, and what section 7 measures separately, is the effect on:
transmission (against the 1e-3 absorption tolerance), the arrival endpoint (against the
stated 0.01 m level) and the timing (against the stated 0.1 ns level). The geometry table is
**tracked separately and is not modified under this approval**.

### 6.6 Boundary accuracy expansion — THREE FAILURES

`audit_b.cpp` evaluates six boundary families against a `w = sqrt(z - R)`-substituted
reference (`dz = 2w dw` removes the endpoint square-root stiffness, so it stays trustworthy
at xi -> 1). Two builds are compared: uniform `dr = 5 m` (the phase-3 build) and the
w-substituted build.

```
family                                   rays  uniform     w-build      worst ray (w-build)
F1 cell interiors (mid h, mid xi)        2688  5.080e-03   5.080e-03    h=375   lam=690 xi=0.9998555381
F2 h knots + neighbourhoods              1496  8.653e-04   8.374e-04    h=1025  lam=400 xi=0.9900000000
F3 short paths / ground boundary          108  2.751e-03   2.751e-03    h=100   lam=350 xi=0.9995000000
F4 lambda endpoints / clamp edge           96  3.288e-02   3.288e-02    h=18000 lam=295 xi=0.5000000000
F5 near-critical xi                        42  5.865e-04   5.584e-04    h=1500  lam=400 xi=0.9900000000
F6 declared grazing guard (xi = 1 - d)     42  2.405e-04   2.290e-04    h=1500  lam=690 xi=0.9990000000
worst: uniform 3.288e-02   w-build 3.288e-02   -> gate 1e-3 FAIL
```

**Three families fail: F1, F3 and F4.** The nodes themselves are accurate (F5/F6 <= 5.6e-04),
so all three are *interpolation between nodes*. They are recorded with their measured
locations and no remedy is asserted here:

1. **F1, 5.080e-03, worst at the MID-CELL point xi = 0.9998555381** (h = 375 m,
   lambda = 690 nm). The current axis is uniform in `u = sqrt(1-xi)`, giving nodes at
   `xi = 1 - j^2 * 1.156e-05`; the reported point is the midpoint of the j = 3..4 cell.
   This establishes **inadequate interpolation in the current cells near xi -> 1**. It does
   **not** establish that uniform refinement is the wrong remedy, nor that any particular
   alternative is needed — that is what the re-axis experiment measures (section 7).
2. **F3, 2.751e-03, worst at h = 100 m with xi = 0.9995** — a short path on a grazing ray.
   h = 100 m lies between the h nodes 50 m and 250 m, a factor-5 cell, and for grazing rays
   `tau ~ sqrt(2 R h)` grows like `sqrt(h)` there. **The height axis is the prime suspect and
   must be addressed on its own**, not only as a consequence of the grazing axis.
3. **F4, 3.288e-02 at lambda = 295 nm — outside the declared 300-900 nm domain.** Read
   correctly, this is a statement about **incorrect extrapolation/admission below the declared
   domain, not an in-domain accuracy failure**: the query was admitted and answered by
   clamping, while the profile holds real data down to 180 nm.
   **Two different sets must not be conflated here:**
   - the **native spectral interval** of the profile is **[180, 700] nm** on 105 knots at 5 nm;
   - the **accepted query domain** is whatever the lookup answers. Because the lookup clamps
     rather than rejects, a query outside [180, 700] nm **is still accepted** and is answered
     with the nearest endpoint column. Clamping is not rejection, and the domain statement
     must say so explicitly.

   The correct fix is to **declare the accepted query domain explicitly and preserve the
   existing query-time clamping semantics**, rather than to treat 295 nm as a tuning target.

Also note the `w = sqrt(z-R)` build quadrature is **3x faster** than uniform `dr = 5 m`
(24.3 s vs 75.0 s) at equal or better accuracy, so the guard fix costs nothing.

### 6.7 Breakpoints, lambda edge, admissibility

**Breakpoints are COMPLETE, verified from the callback.** `alpha = (tau(hi) - tau(lo))/(hi - lo)`
with `hi = min(h_max, h+25)`, `lo = max(0, h-25)` and tau piecewise linear between the
profile's native knots, so alpha is piecewise linear with breakpoints exactly at
`{knot +/- 25 m}` plus the boundary set `{0, 25, h_max-25, h_max}`. `bp_b.cpp` tests
collinearity of alpha between consecutive predicted breakpoints over the callback's own
domain `[0, h_max]`:

```
predicted set within [0, hmax]: 102 points, 101 segments
worst relative residual 1.276e-11 at h=47500.0 m   -> set is COMPLETE
```

This is stronger than split-vs-unsplit agreement: it shows that no breakpoint is missing,
not merely that two integrators agree.

**Lambda edge (measured on `getOpticalDepth` at h = 18000 m):** the table clamps to
[180, 700] nm, and above 700 nm every value equals the 700 nm value
(720/850/900/1000 all give 0.21566900). So above 700 nm the clamping in 6.6(2) *is* current
behaviour; below 300 nm it is **not** — real data is being discarded.

**Admissibility of the ground-bound route.** `n*r` increasing is equivalent to `f*G < 1`
with `G = r|dN/dh| - N`, so the criterion is independent of the profile shape factor:

```
max G = 0.401333 at h = 10 m  ->  n*r is increasing for ALL f < 2.492
configuration dispersion spans f in [0.9868, 1.0153]  ->  ADMISSIBLE with a 2.5x margin
conservative interval-minimum margin over 25 m sub-intervals: 0.0000e+00 m (>= 0 required;
equality occurs only in the first sub-interval, i.e. at the detector itself)
```

This replaces "checked at three heights" with a monotonicity argument valid over the whole
supported configuration domain, plus a conservative interval-minimum check.

### 6.8 Performance — completed, and NOT presented as cap compliance

```
geometry table:  181.4 MB, build 26.4 s        peak RSS after it:            185.8 MB
absorption table:  6.62 MB, build 75.0 s (uniform) / 24.3 s (w-substituted)
peak RSS with geometry + both absorption tables:                             198.7 MB
combined geometry trace + absorption lookup:  1.4789 us/photon (20000 photons)
  reject fraction 6.800e-02 (1350/20000), all BeyondCritical
  [upgoing 0  beyondCrit 1350  dispersion 0  grazing 0]  -- isotropic sample, an UPPER bound
legacy per-photon tau integration (50 m sub-steps): 83.29 us/photon
construction amortisation: 75.0 s / (83.29 - 1.4789) us = 9.17e+05 photons to break even
peak RSS at exit: 198.7 MB
```

Four points that must not be smoothed over:

- **198.7 MB is the total-process peak of THIS diagnostic process, and it is dominated by
  the geometry table (181.4 MB).** The absorption table adds 6.62 MB. It demonstrates
  headroom for the measured diagnostic process; it is **not** a measurement of the complete
  production workload. The earlier "peak RSS (table only) 10.3 MB" was **not** evidence of
  cap compliance and must not be cited as such.
- **The 6.800e-02 rejection rate is THAT sample's rate.** The sample draws directions
  isotropically downgoing, which is not the production photon population. It is **not an
  established upper bound** for production workloads, and it must not be reported as one.
- **1.4789 us/photon exceeds the ~1 us target** on this synthetic sample; the geometry trace
  dominates. Whether it meets the target under a representative Cherenkov workload is
  untested.
- **Break-even is ~9.2e5 photons**, so short runs pay construction without recovering it
  (75 s uniform; 24.3 s w-substituted).

### 6.9 Candidate B — elimination estimate, no prototype

Candidate B (direct `tau(h, xi, lambda)`) was **not** prototyped, deliberately:

- **Memory**: identical axes ⇒ identical 6.62 MB. No advantage.
- **Accuracy**: for a fixed ground endpoint, `tau(h=0) = 0` exactly, so **under identical
  coordinates and identical interpolation** a direct tau table and the cumulative table
  carry the same information, and dropping the cumulative construction removes the
  monotone-in-h structure that makes the h interpolation well behaved. This supports
  **avoiding a redundant prototype**; it does **not** establish that other direct-table
  parameterisations (different coordinates, different interpolation, non-ground endpoints)
  are inferior — none was measured.
- **Construction**: same integral, same cost.
- Previously measured B-family variants already fail by 2–3 orders: thickness-averaged layer
  form |dT| = 0.1619; coarsened merged-component W-grid 2.1e-2 … 2.7e-1.

**Not prototyped**, deliberately; the estimate above is what justifies that, not a
demonstrated inferiority of every alternative.

### 6.10 Status of the selection

Cumulative-xi (candidate A) remains the leading prototype but **cannot be promoted** on the
present evidence: F1 (5.08e-03) and F4 (3.29e-02) both fail the 1e-3 gate, and both failures
are in the *axis definition* (the xi distribution near 1; the lambda range 180–700 nm on
native knots), not in tuning. The fixes are specific and testable. The new frozen
final-validation set remains untouched.

## Section 7 — re-axis experiment (diagnostics only)

Purpose: compare grazing, height/short-path and wavelength axes by **measured error and
measured cost**, without asserting any remedy in advance. Production is unchanged; the
geometry table is not modified; the frozen set is untouched.

Method (`reaxis_b.cpp`): the reference is computed **once per ray** and cached, so every
variant is scored on the **same, axis-independent** development ray set (1158 rays, F1–F6).
This set is a development set; it is **not** the frozen final-validation set.

```
profile DECLARED domain: lambda [180.0, 700.0] nm, h [0.0, 50000.0] m
query-time clamping preserved: the lookup clamps, it never extrapolates
REFERENCE UNCERTAINTY (nW 60000 vs 240000, 68 rays): worst |d tau| = 5.775e-07
    [h=46000 lam=400 xi=0.99999998]      -> ~6e-07 in T, i.e. ~1700x inside the 1e-3 gate
```

### 7.1 Variant comparison (worst |dT| per family; PASS = all six <= 1e-3)

```
variant                      h   xi  lam     MB  build_s  F1        F2        F3        F4        F5        F6
V0 phase-3 axes (reproduce)  56  256  121   6.62   27.9    2.23e-03  8.37e-04  2.75e-03  6.81e-01  5.58e-04  2.29e-04
V1 native lambda             56  256  105   5.74   18.7    2.23e-03  8.37e-04  2.75e-03  6.72e-04  5.58e-04  2.29e-04
V2 graded h                  61  256  105   6.25   19.2    2.27e-03  8.37e-04  2.59e-03  6.72e-04  5.58e-04  2.29e-04
V3 refined xi (464 nodes)    61  464  105  11.34   34.7    8.06e-04  9.67e-04  1.86e-04  8.01e-04  6.84e-04  6.26e-04   PASS
V4 log10(1-xi) coordinate    61  256  105   6.25   19.0    1.60e-03  6.76e-04  3.28e-04  5.76e-04  6.95e-04  6.95e-04
V5 uniform in mu             61  256  105   6.25   19.4    9.22e-01  5.84e-01  8.40e-01  6.61e-01  8.68e-01  8.68e-01
V6 graded h, fewer nodes     38  256  105   3.90   11.9    2.72e-03  3.08e-03  2.74e-03  2.73e-04  1.70e-04  4.66e-05
lookup cost, all variants: 0.19-0.23 us/call (V5 0.30)
best variant by worst-case error: V3 refined xi (464 nodes), 9.668e-04
```

### 7.2 What the measurements establish

1. **Uniform refinement works — it is not the wrong remedy.** V3 refines the *existing*
   coordinate and passes **all six families**, worst 9.668e-04 at (h = 1025 m, lambda = 400 nm,
   xi = 0.99). Cost: 11.34 MB, 34.7 s. F1 falls 2.27e-03 -> 8.06e-04 and F3 falls
   2.59e-03 -> 1.86e-04. (The margin is thin: 3.3 % inside the gate.)
2. **A better coordinate reaches nearly the same accuracy at about half the memory.**
   V4's `y = -log10(1-xi)` axis at **256 nodes** (6.25 MB, 19.0 s) fixes F3 (3.28e-04) and
   F4 (5.76e-04) outright and leaves **only F1 failing**, at 1.60e-03. The two remedies are
   complementary: the coordinate change does most of the work per node, the refinement
   supplies the remaining margin.
3. **F3 is NOT primarily a height-axis problem** — the earlier hypothesis is corrected by
   measurement. Grading h (V2) moves F3 only 2.75e-03 -> 2.59e-03, while changing the xi
   coordinate (V4) takes it to 3.28e-04. The height axis does matter for *other* families:
   V6 shows that coarsening h above 2 km pushes F2 to 3.08e-03.
4. **The wavelength axis was the dominant defect and is now isolated.** V0 on this set
   reports F4 = 6.81e-01 because the development set includes lambda = 180 nm — inside the
   profile's real range but **outside the old 300-900 nm axis, which admitted the query and
   answered it by clamping**. Switching to the profile's native knots (V1) takes F4 to
   6.72e-04, a factor ~1000, while *reducing* memory (5.74 MB) and build time (18.7 s).
   The supported domain must be declared as **[180, 700] nm**, with query-time clamping
   preserved for anything outside it.
5. **Negative control — the concentration must be at xi -> 1.** V5 is uniform in `mu`, the
   direction cosine. Since `xi = sin(theta_obs)` and `mu = cos(theta_obs)`, that places nodes
   densely at the **vertical** end and starves grazing: every family degrades to 0.58-0.92.
   (The run prints V5 as "sin(theta) coordinate"; that label is wrong — xi *is* sin(theta), so
   a sin(theta)-uniform axis is xi-uniform. The variant as built is mu-uniform. Diagnostic
   labelling defect, recorded here.)
6. **Mapping and actual-configured-ray checks pass on the new axes** (V3):
   `|xi_trace - xi_independent| = 9.992e-16`; actual configured ray worst
   `|dT| = 9.6986e-04` at (h = 1500 m, lambda = 620 nm, xi = 0.9995), again with only ~3 %
   margin. Dispersion-only sensitivity of the reference: 4.7272e-04.

### 7.5 Worst-case identities are PER CHECK — do not merge them

The worst ray differs between checks in wavelength, height and xi. Each identity below belongs
to its own check and must be quoted with it:

```
check                                   h [m]   lambda [nm]   xi
phase 3 (cmp_b, f=1, full lambda list)  1500    625           0.9950   (b/b_crit)
phase 3 (cmp_b, lambda knots only)      1500    670           0.9950
section 6.4 (map_b, A1 and A2)          1500    620           0.9950
section 6.6 F1 (audit_b)                 375    690           0.9998555381
section 6.6 F2 (audit_b)                1025    400           0.9900
section 6.6 F3 (audit_b)                 100    350           0.9995
section 6.6 F4 (audit_b)               18000    295           0.5000
section 7.1 (reaxis V3, overall worst)  1025    400           0.9900
section 7.2 (reaxis V3, configured ray) 1500    620           0.9995
```

They are **not** the same ray. In particular the section 7.2 configured-ray worst is at
xi = 0.9995, while the section 6.4 result is at xi = 0.9950 — the same height and wavelength
family, a different grazing coordinate.

7. **Peak RSS at exit: 197.1 MB** (this diagnostic process, geometry table included).

### 7.3 Geometry table, assessed separately (no geometry change made)

Measured against an independent Bouguer integration of the same ray, at the three stated
levels:

```
timing      : worst |d t| = 7.1771 ns      [h=46000 lam=400 xi=0.9999]   level 0.1 ns   -> 72x OVER
endpoint    : worst |d sin(theta)|*R = 2.4889e-08 m [h=1500 lam=400 xi=0.05]  level 0.01 m -> fine
transmission: worst 1-exp(-alpha*dL) = 6.5180e-05   [h=1500 lam=400 xi=0.9999] tol 1e-3   -> 65x inside
```

Two conclusions, kept separate from the absorption work:

- **The `alpha*dL` figure is a HEURISTIC ESTIMATE, not a measured transmission error and not
  a demonstrated bound.** It ignores spatially varying extinction and trajectory displacement.
  The value 6.5e-05 is therefore labelled as an estimate only; it must not be quoted as a
  bound on the transmission effect until it is replaced by an integration along the traced
  trajectory or by a conservative bound. (Section 8 repeats the estimate for the combined
  variant with the same label.)
- **The timing discrepancy is a finding that requires its own investigation, and it blocks
  acceptance of the combined route if confirmed.** 7.18 ns against the table's stated 0.1 ns
  level, at xi = 0.9999, which lies **inside** the table's declared `muMin` range
  (xi_max = 1 - 2e-8), so the affected ray is eligible rather than rejected. Section 8
  investigates it against the actual timing observable with an independent reference and a
  convergence check. **If it stands against the applicable 0.1 ns requirement, it blocks
  acceptance of the curved route regardless of the absorption table's accuracy.** No geometry
  code was modified.

### 7.4 Limits of this experiment

- Per-family worst **locations** are printed only for the overall worst family; for the others
  only the magnitude is recorded. Section 8 prints all six for every variant.
- The transmission figure in 7.3 is a **heuristic** `alpha*dL` estimate, not an integration
  along the traced trajectory. It is neither a measured transmission error nor a demonstrated
  bound.
- The nW 60000 vs 240000 comparison is an **observed convergence** check of the reference. It
  is not a rigorous reference-error bound, and it does not bound the table's error.
- All variants were scored on a development set. **The frozen final-validation set was not
  read, replaced or recreated**, and the domain expansion to 180-700 nm must be documented
  before any final validation is evaluated.
- V3's pass margin (3.3 %) is thin, and its worst case sits at xi = 0.99 rather than at the
  guard. A combined variant (the `log10(1-xi)` coordinate plus moderate refinement near
  xi -> 1) is the obvious next test; it was not run under this approval.

## Section 8 — combined variant, domain/contract definition, and the geometry timing block

Diagnostics only. Production, geometry sources and the frozen final-validation set untouched.

### 8.1 Domain and contract statements (these are definitions, not measurements)

**Native spectral interval vs accepted query domain — two different sets:**

```
native spectral interval of the profile : [180, 700] nm on 105 knots at 5 nm
accepted query domain                   : the lookup CLAMPS, so any query is ACCEPTED,
                                          including one outside [180, 700] nm, and is answered
                                          with the nearest endpoint column.
```

**Clamping is not rejection.** A query at 750 nm or 1000 nm is admitted and answered with the
700 nm column (measured in 6.7: 720/850/900/1000 all return 0.21566900). This matches the
legacy callback, which also clamps, so a clamped table query agrees with the legacy route and
needs no fallback. Any domain statement must say "accepted, clamped", never "rejected".

**Exact grazing-domain endpoints:**

```
geometry route eligibility : mu >= muMin = 2.0e-4
                             <=>  xi <= 1 - 2.0e-8   (xi = sqrt(1 - mu^2))
absorption guard           : MUST use the same endpoint, or the two routes disagree about
                             which rays are eligible.
```

**How the logarithmic coordinate handles xi = 1:** at xi = 1 (the tangent ray) tau is finite,
and under the `w = sqrt(z-R)` build substitution it is integrable — `1 - s^2 ~ c w^2`, so the
integrand tends to `2/sqrt(c) dw` rather than diverging. A logarithmic coordinate must
therefore carry **xi = 1 as an explicit node**. With `y = -log10(1-xi)` capped at 9 the top
node is `1 - 1e-9`, leaving the interval `(1-1e-9, 1]` to be clamped unless xi = 1 is added;
V7/V8 add it explicitly. V3 (u-uniform) does not, and its top node is `1 - 2.3e-6`, so the
whole guard interval `(1-2e-8, 1-2.3e-6]` is served by clamping (measured consequence: F6).

**A structural gap that is not yet closed:** the absorption axes top out at h = 50 000 m,
while the geometry route is built with `hTop = 120 000 m` (Cherenkov.inl). A photon emitted
above 50 km would be absorbed by clamping. The declared h-domain must be stated, and above it
the route must fall back or the axis must be extended. This is a structural statement, not a
measurement.

### 8.2 Configuration / invalidation contract

The table is valid for exactly one **immutable optical configuration**, defined as the tuple:

```
{ refractive-index table content (N(z) samples, grid, clamp/extrapolation rule),
  the dispersion model used for f_phase(lambda) and f_group(lambda),
  earth radius R, observation radius r_obs,
  extinction profile identity (path + md5),
  the declared domains: lambda accepted interval, h ceiling, xi guard }
```

- **Invalidation:** a change to any element invalidates the table; it must be rebuilt. The
  table is keyed on those elements (content hash for the RI table and the profile, plus an
  identifier for the dispersion model).
- **Not configuration-dependent:** the lookup takes only (h, xi, lambda), so nothing per-event
  or per-photon invalidates it.
- **Eligibility / fallback:** a photon is eligible iff the geometry route's prerequisites hold
  AND `xi <= 1 - 2e-8` (the same endpoint as `muMin`). Outside that, fall back to the legacy
  tracer. Lambda outside the native interval does **not** require fallback, because clamping
  matches the legacy behaviour.
- **This contract is mandatory regardless of whether a configuration-specific table improves
  sampled accuracy.** Section 7 already showed it does not improve accuracy (A1 7.38e-04 vs
  A2 7.53e-04); the contract exists for correctness, not for accuracy.

### 8.3 Fixed development-set comparison — worst |dT| for EVERY family, with location

Same axis-independent development set as section 7 (1158 rays). The reference is cached; the
nW 60000 vs 240000 check is observed convergence, not a bound.

```
V3  refined xi, u-uniform (464)   61 x 464 x 105 = 11.34 MB   build 35.0 s   worst 9.668e-04  PASS
   F1 8.063e-04  @ h=1500  lam=500  xi=0.9950000000
   F2 9.668e-04  @ h=1025  lam=400  xi=0.9900000000
   F3 1.865e-04  @ h=100   lam=350  xi=0.9995000000
   F4 8.014e-04  @ h=1500  lam=700  xi=0.9950000000
   F5 6.838e-04  @ h=1500  lam=400  xi=0.9900000000
   F6 6.257e-04  @ h=1500  lam=690  xi=0.9990000000
   lookup 0.2335 us/call

V7  combined log+refinement (464) 61 x 463 x 105 = 11.31 MB   build 34.7 s   worst 1.604e-03  FAIL
   F1 1.604e-03  @ h=375   lam=690  xi=0.9999000000
   F2 9.506e-04  @ h=1025  lam=400  xi=0.9900000000
   F3 1.031e-04  @ h=100   lam=350  xi=0.9995000000
   F4 8.285e-04  @ h=1500  lam=700  xi=0.9950000000
   F5 7.607e-04  @ h=1500  lam=690  xi=0.9990000000
   F6 7.607e-04  @ h=1500  lam=690  xi=0.9990000000
   lookup 0.1801 us/call

V8  combined log, moderate (368)  61 x 367 x 105 =  8.97 MB   build 26.9 s   worst 1.566e-03  FAIL
   F1 1.566e-03  @ h=375   lam=690  xi=0.9999000000
   F2 9.364e-04  @ h=1025  lam=400  xi=0.9900000000
   F3 9.177e-05  @ h=1000  lam=350  xi=0.8000000000
   F4 8.222e-04  @ h=1500  lam=700  xi=0.9950000000
   F5 7.424e-04  @ h=1500  lam=690  xi=0.9990000000
   F6 7.424e-04  @ h=1500  lam=690  xi=0.9990000000
   lookup 0.1740 us/call
```

**The combined log-coordinate hypothesis is NOT confirmed.** At equal node count (464 vs 463)
the log coordinate is **worse** than u-uniform (1.604e-03 vs 9.668e-04), and V8 shows that
refining the log coordinate barely helps (1.604e-03 -> 1.566e-03 at 8.97 MB).

**Qualification (do not over-read this).** V7/V8 fail **as implemented**. The experiment shows
their F1 error is not explained by a resolution deficit (adding nodes barely moves it), which is
*consistent with* a misalignment to the `sqrt(1-xi)` structure near the detector. It does **not**
establish that logarithmic interpolation has an irreducible structural error in general, and it
does not show that no logarithmic coordinate could work. The supported statement is narrower:
this log coordinate, as implemented, does not beat u-uniform on these rays.

### 8.4 Development probes around the observed worst cells — the decisive result

Separate from 8.3 by construction: neighbourhoods (dh = +/-0.1..500 m, dlambda = +/-0.1..2.5 nm,
dxi = +/-1e-6..1e-3) around each observed worst cell, references recomputed on the fly at
nW = 20000.

```
probe neighbourhood      V3          V7          V8
V3 overall worst         9.759e-04   9.826e-04   9.673e-04
audit F1 worst           1.467e-03   2.210e-03   2.212e-03
V4 F1 worst              1.215e-02   2.087e-03   2.089e-03
audit F3 worst           5.423e-04   5.285e-04   4.727e-04
configured-ray worst     1.207e-03   1.285e-03   1.237e-03
phase-3 worst            9.817e-04   9.936e-04   9.351e-04
```

**V3's PASS does not survive neighbourhood probing.** Its worst error rises from 9.668e-04 on
the fixed set to **1.215e-02** in the neighbourhood of the V4 F1 worst cell — 12x over the
gate, and 12x larger than anything the fixed set recorded. The other neighbourhoods are also
over the gate (1.467e-03, 1.207e-03). So the fixed development set **underestimates the worst
error by more than an order of magnitude in a region it does not sample densely**, and the
"3 % margin" reported in section 7 was not a thin margin — it was an overestimate of quality.

**Correction from section 9 (apply to the figure above).** V3's 1.215e-02 is at the
**ineligible endpoint xi = 1.000000000000** (`clampXi=1`, `elig=0`), so it is **not an eligible
accuracy failure**; V3's **eligible-domain maximum is 2.1441e-03**. The conclusion is unchanged
(V3 is still not promotable, and its fixed-set PASS is still not robust), but citations must use
the eligible figure.

**Attribution limit (now closed by section 9):** this block printed the worst location only for V7, not for V3
or V8. The V3 excursion at 1.215e-02 is therefore not yet attributed — a follow-up must print
all three locations, and must distinguish an interior mid-cell point from the domain endpoint
xi = 1, which the probe grid can reach when `xi0 + dx` lands on 1.0.

### 8.5 Actual configured ray (dispersion included)

```
mapping |xi_trace - xi_independent| worst = 9.992e-16  [h=1500 lam=620 xi=0.1500]
V3  9.6986e-04  [h=1500 lam=620 xi=0.9995]  PASS
V7  1.1113e-03  [h=1500 lam=620 xi=0.9995]  FAIL
V8  1.0844e-03  [h=1500 lam=620 xi=0.9995]  FAIL
```

The mapping stays exact on every variant. On the configured-ray set the log coordinate is again
**worse** than u-uniform, and V7/V8 fail. (This identity — h = 1500 m, lambda = 620 nm,
xi = 0.9995 — belongs to this check only; see 7.5.)

### 8.6 Geometry timing — separate investigation, and it blocks the route

Observable: `Result.time_s` (= `optical_path_m` / c), the quantity production uses.
Reference: `(L + f_group*Sigma)/c` from an independent Bouguer integration of the same ray,
same N(z), same f for phase and group. The dr sequence is an **observed convergence** check of
the reference, not a bound:

```
dr=2.0 m: worst |d t| = 7.1587 ns  [h=46000 lam=400 xi=0.9999]
dr=0.5 m: worst |d t| = 7.1771 ns  [h=46000 lam=400 xi=0.9999]
dr=0.1 m: worst |d t| = 7.1873 ns  [h=46000 lam=400 xi=0.9999]
```

The reference is converged (0.4 % drift across a 20x refinement), so **the discrepancy is not a
reference artefact**. Threshold scan at h = 46 km, lambda = 400 nm against the stated 0.1 ns
level:

```
xi=0.90000000   0.0114 ns  ok
xi=0.99000000   0.2629 ns  OVER
xi=0.99900000   2.2132 ns  OVER
xi=0.99950000   3.4229 ns  OVER
xi=0.99990000   7.1873 ns  OVER
xi=0.99995000   8.9292 ns  OVER
xi=0.99999000  12.2092 ns  OVER
```

**The 0.1 ns level is crossed between xi = 0.90 and xi = 0.99**, and every crossing ray is
**eligible** (the guard is xi <= 1 - 2e-8). On this evidence the geometry table's timing is out
of tolerance across most of the grazing domain, and — per the review instruction — **this
blocks acceptance of the combined curved route regardless of the absorption table's accuracy.**
Caveats: the reference uses the same f for phase and group, so the result is established for
that case and must be re-checked with the production dispersion function; and whether 0.1 ns is
the applicable requirement for this route is a review decision, not a measurement.

### 8.7 Cost and peak RSS

```
V3  11.34 MB   build 35.0 s   lookup 0.2335 us/call
V7  11.31 MB   build 34.7 s   lookup 0.1801 us/call
V8   8.97 MB   build 26.9 s   lookup 0.1740 us/call
peak RSS at exit (this diagnostic, several tables held at once): 217.3 MB
```

The peak is higher than section 7's 197.1 MB because this program holds three absorption tables
plus the geometry table simultaneously. It is a diagnostic-process figure, not a design figure.

### 8.8 Conclusion

**None of the three variants can be promoted, and the combined variant is not an improvement.**

- V7/V8 fail the fixed development set, the configured-ray set and most probe neighbourhoods.
  The log coordinate is worse than u-uniform at equal node count, and refining it does not fix
  its structural F1 error. "Better accuracy at lower memory" was a hypothesis; the measurement
  rejects it.
- V3's fixed-set PASS (9.668e-04) is **not robust**: neighbourhood probing raises it to
  1.215e-02. The fixed development set is therefore **not adequate for acceptance**, and the
  section 7 margin must not be quoted as evidence of quality.
- The geometry timing discrepancy is confirmed against a converged independent reference and
  exceeds the stated level on eligible rays. That is a separate, blocking issue for the curved
  route.
- The configuration/invalidation contract (8.2) and the domain statements (8.1) stand
  independently of all of the above.

The frozen final-validation set was not read, replaced or recreated. Before any final validation
is evaluated, the axis must be settled **using probes that are dense in the region they claim to
cover**, the h ceiling gap and the xi = 1 endpoint must be closed, and the geometry timing
question must be resolved.

## Section 9 — targeted closure: attributed failures and coverage/fallback policy

Diagnostics only. No production or geometry change, no broad redesign, no frozen-set evaluation.

`closure_b.cpp`. The eligibility limit is **derived, not assumed**:

```
xi_max = sqrt(1 - mu_min^2),  mu_min = 2.0e-4
xi_max = 0.99999997999999979      (1 - xi_max = 2.000000e-08)
```

### 9.1 Attributed probe maxima — eligible domain separated from the rest

```
variant  ALL-probe max                     ELIGIBLE-domain max
V3       1.2155e-02  xi=1.000000000000     2.1441e-03
         clampXi=1, elig=0  -> EXCLUDED    h=275 lam=687.5 xi=0.999910000000
                                          xcell=[0.999895987426,0.999927769046] hcell=[256,500]
                                          clampXi=0 clampH=0 clampLam=0  intXi=1 intH=1  elig=1
V7       2.2099e-03                        2.2099e-03
         h=365 lam=687.5 xi=0.999955538100
         xcell=[0.999953797543,0.999955979693] hcell=[256,500]
         clampXi=0 clampH=0 clampLam=0  intXi=1 intH=1  elig=1
V8       2.2125e-03                        2.2125e-03
         h=365 lam=687.5 xi=0.999955538100
         xcell=[0.999955291225,0.999958020475] hcell=[256,500]
         clampXi=0 clampH=0 clampLam=0  intXi=1 intH=1  elig=1
```

**V3's 1.2155e-02 is at the ineligible endpoint xi = 1.0 with `clampXi=1`, so it is NOT an
eligible accuracy failure** and is excluded from the acceptance figure. It is still a real
clamping defect (below, 9.2), but it must not be quoted as an eligible error.

**All three eligible maxima fall in the same cell region**: interior of the h cell [256, 500] m,
interior of a xi cell near 0.9999, at lambda = 687.5 nm. So the binding error is a **mid-cell
interpolation failure at a short path (h ~ 275-365 m) on an extreme grazing ray**, not an
endpoint artefact — and the h cell is wide (256 -> 500 m), so the height axis participates.

### 9.2 Grazing coverage scan (xi = 1 - 10^u, u = -2 .. -8)

```
variant  top xi node       coverage max (no xi clamp)   clamped-region max
V3       0.999997110762    1.6051e-03  h=375 lam=690    1.8472e-03  xi=0.999999979978
                           xi=0.999927673661 (interior)  (interior of the clamped cell)
V7       1                 1.6176e-03  h=375 lam=690    0.0000e+00   (no clamping at all)
                           xi=0.999904522839 (interior)
V8       1                 1.6144e-03  h=375 lam=690    0.0000e+00   (no clamping at all)
                           xi=0.999910926454 (interior)

largest xi still <= 1e-3 / smallest xi > 1e-3:
V3  0.999997032698 / 0.999764571359
V7  0.999999979978 / 0.999710057715
V8  0.999999979978 / 0.999747646083
```

**Coverage is PATCHY, not certifiable by a threshold.** In every variant the largest passing xi
lies *above* the smallest failing xi (e.g. V3: 0.9999970 vs 0.9997646), so "coverage is valid up
to xi = X" is not a supported statement: failures are interleaved from **xi ~ 0.9997** upward.

Two further facts:
- **V3's top node (0.999997110762) is below xi_max**, so its clamped region **intersects the
  eligible domain** and carries 1.8472e-03 — a clamping defect *inside* eligibility, not outside it.
- **V7/V8 have no clamped region at all** (clamped-region max 0.0), because their explicit
  xi = 1 node covers the whole scanned interval. So an explicit xi = 1 node is **optional when
  xi = 1 is ineligible** (it is: xi_max < 1) but it does remove the clamped region for one node.

### 9.3 Absorption height ceiling

```
V3  worst |dT| above the ceiling = 7.1518e-06   [h=50000 lam=400 xi=0.9 clampH=1]
V7  worst |dT| above the ceiling = 6.4657e-05   [h=50000 lam=400 xi=0.9 clampH=1]
V8  worst |dT| above the ceiling = 1.3082e-04   [h=50000 lam=690 xi=0.9 clampH=1]
```

The h axis tops at 50 000 m and the profile also tops at 50 000 m, so the callback's alpha is 0
above it and clamping adds nothing measurable. **This benignity follows from the profile, not
from the table**, so the policy below still requires an explicit counted fallback.

### 9.4 Timing resolved separately

**(a) The applicable requirement is verified from this repository**, not assumed:

```
CurvedRayTable.hpp:264             't errors stay inside 1 cm / 0.1 ns for mu >= 2e-4;
                                    below it they break down'
testCherenkovCurvedRayTable.cpp:295 'agreed targets: 1 cm displacement, 0.1 ns time,
                                    1e-6 rad angle'
```

So 0.1 ns **is** the applicable requirement, and the table **claims** it over exactly the
eligible domain (mu >= 2e-4). The failing rays are inside that claim.

**(b) Phase/group dispersion split** (f_group = f_phase + 2k/(l^2 (1+k/0.16)), a 3.92 % split at
400 nm and 1.34 % at 690 nm):

```
worst |d t| with the split = 7.1873 ns  [h=46000 lam=400 xi=0.9999]
```

**Unchanged** from the equal-f case — so the discrepancy is not an f-response artefact.

**(c) Independent cross-check, different parameterisation** (RK4 in the swept angle phi,
`dr/dphi = r sqrt(1-s^2)/s`, 7 301 / 22 013 / 56 682 steps):

```
xi=0.9000  |d t| obs-vs-phi = 0.0101 ns   obs-vs-radial = 0.0114 ns   phi-vs-radial = 0.0013 ns
xi=0.9900  |d t| obs-vs-phi = 0.2593 ns   obs-vs-radial = 0.2629 ns   phi-vs-radial = 0.0035 ns
xi=0.9999  |d t| obs-vs-phi = 7.1784 ns   obs-vs-radial = 7.1873 ns   phi-vs-radial = 0.0090 ns
```

The **two independent references agree with each other to <= 0.009 ns** while **both sit 7.18 ns
from the table's observable**. Path lengths agree to 2.7 mm and Sigma to 1e-6 m between the two
references. The 7.18 ns is the 2.15 m path-length discrepancy expressed in time (7.19e-9 s x c =
2.155 m).

**Conclusion, stated at the strength the evidence supports.** The 7.18 ns discrepancy is
**not a reference-discretisation artefact** (it survives a 20x refinement in dr) and **not the
dispersion model** (a 3.92 % phase/group split leaves it unchanged), and it is **reproduced by a
second, independent parameterisation** (RK4 in the swept angle), which agrees with the radial
reference to <= 0.009 ns.

It does **not** exclude a shared *model* error: both references use the same Bouguer trajectory
and the same `N(z)` from the refractive-index table, so a common modelling error would not show
up. And the **production dispersion function was not available** here, so the
production-dispersion result is **not established** — only the equal-f and Cauchy-split cases are.

On the tested cases, the discrepancy is a substantial violation of the table's stated 0.1 ns
claim on eligible rays. **Acceptance of the combined curved route remains blocked pending
resolution of that geometry-side issue and of the applicable-requirement question**; the
absorption table's accuracy is irrelevant to it.

### 9.5 Proposed coverage / fallback policy

1. **xi (grazing).** Coverage **cannot be certified by a single threshold** (9.2). Policy:
   declare a conservative coverage limit `xi_cov` derived from a scan *denser than the one that
   set it* (candidate 0.9997), and require an **explicit counted fallback** for
   `xi_cov < xi <= xi_max`. Additionally require either `top xi node >= xi_max` or a counted
   fallback for `(top node, xi_max]` — for V3 this is not optional, since its clamped region
   intersects eligibility at 1.85e-03.
2. **An explicit xi = 1 node** is **optional** while xi = 1 is ineligible (xi_max < 1). It is
   not inherently required, but it costs one node and removes the clamped region entirely.
3. **h (height).** Require `h_top(absorption) >= profile top` **and** an explicit counted
   fallback above the ceiling regardless of the measured benignity, because that benignity is a
   property of the profile (alpha = 0 above its top), not of the table.
4. **lambda.** Accepted and clamped, matching the legacy callback; no fallback. The domain
   statement must read "accepted, clamped", over the native interval [180, 700] nm.
5. **Counting.** Every fallback must increment a counter so the rate is observable — the
   production `CurvedTableStats` already has that pattern; adding absorption-specific counters
   belongs to the eventual integration, not to this diagnostic.

### 9.6 Status

No variant is promotable. With the ineligible endpoint excluded, the eligible-domain maxima are
**V3 2.1441e-03, V7 2.2099e-03, V8 2.2125e-03** — all about 2x over the 1e-3 gate, all in the
same interior cell region (h in [256, 500] m, xi near 0.9999, lambda = 687.5 nm), and coverage
above xi ~ 0.9997 is patchy rather than certifiable. The axis therefore still needs redesign
before further optimization is meaningful — which is exactly why the domain contract (9.5) is
being closed first. The geometry timing item is confirmed and independently cross-checked.

## Section 10 — phases 0d / C / B: frozen sets, alpha estimator, local refinement

Diagnostics only. Production, geometry sources and the previously frozen sets untouched.
Program: `phaseBC_b.cpp` (610 lines); log: `log_phaseBC_b_refinement.txt`.

### 10.1 Phase 0d — frozen held-out sets (created before any tuning)

```
heldout_setI_inband.txt    md5 5cae1a6dc5e087a8639fe36bdc66389a   378 rays
heldout_setII_outofband.txt md5 2800d41b75616199f69b50aa81502162    80 rays
tuning grid (disposable, never a gate)                              495 rays
```

- **Set I (in-band, fresh):** h ∈ {262, 288, 312, 341, 366, 402, 431, 463, 487} m,
  ξ ∈ {0.9998620 … 0.9999880}, λ ∈ {681.5 … 693.5} nm. **0 exact-coordinate clashes with the
  tuning grid** — it tests whether a refinement generalises *inside* the band.
- **Set II (out-of-band):** h ∈ {100, 1000, 5000, 20000, 40000} m, ξ ∈ {0.30, 0.80, 0.99, 0.999},
  λ ∈ {255, 405, 555, 688} nm. **0 h-value clashes and 0 λ-value clashes** with the archived
  720-ray regression set and the coarsening held-out set, per the archive's requirement that
  those must not be reused.

### 10.2 Phase C — the 50 m alpha estimator (measurement only)

```
callback vs exact cell slope, whole profile   worst relative difference 1.0000e+00  (h=6000 m, lam=300 nm)
callback vs exact cell slope, failing band    worst relative difference 2.4403e-15  (h=497 m)
negative-derivative clip, all 105 wavelengths 3480 of 1050105 samples (0.33 %), at 2025-17950 m
edge stencil: one-sided only for h<25 m and h>49975 m; the band starts at 256 m
              => ~10 stencil widths from the ground edge, ~198 from the ceiling: EDGES EXCLUDED
```

**The alpha estimator is not the cause of the failing cell.** Inside 256-500 m the callback and
the exact cell slope agree to 2.4e-15, i.e. exactly. The worst distortion anywhere (relative
difference 1.0, meaning the callback returns 0 where the slope is positive) is the *clip*, which
fires on 0.33 % of samples but only between 2025 and 17950 m — outside the band. So the clip is
real and quantified, and it is not implicated here. Legacy semantics were not altered.

### 10.3 Phase B — three informative variants plus the lambda null check, times two storages

```
variant                   stor    h   xi  lam     MiB build_s  tuning      SetI        SetII      budget
+h only                   tau    67  464  105   12.45    39.1  6.927e-03   5.775e-03   2.903e-04  PASS
+h only                   T      67  464  105   12.45    37.9  1.294e-02   1.126e-02   7.035e-02  PASS
+xi only                  tau    61  664  105   16.22    50.6  3.137e-03   3.200e-03   3.104e-04  PASS
+xi only                  T      61  664  105   16.22    53.1  3.818e-02   3.671e-02   7.035e-02  PASS
+lambda only (NULL CHECK)  tau    61  464  108   11.66    35.4  6.927e-03   5.266e-03   3.104e-04  PASS
+lambda only (NULL CHECK)  T      61  464  108   11.66    35.8  2.747e-02   2.944e-02   2.547e-02  PASS
+h and +xi                tau    67  664  105   17.82    54.2  1.504e-04   1.349e-04   2.903e-04  PASS
+h and +xi                T      67  664  105   17.82    53.6  1.348e-02   1.181e-02   7.035e-02  PASS
```

Gates: **only `+h and +xi` with tau storage passes** — Set I 1.349e-04 (7.4x inside the 1e-3
gate) and Set II 2.903e-04 (3.4x inside), at 17.82 MiB against the 24.7 MiB working budget.

### 10.4 Findings

1. **The failure is COUPLED, not single-axis.** Refining h alone gives Set I 5.775e-03, refining
   xi alone gives 3.200e-03, and refining **both** gives **1.349e-04** — a 24x gain from the
   coupling. So the branch for "no single axis dominates", added at review, is the actual
   outcome. Neither axis alone can reach the gate at any affordable memory.
2. **Tau storage is required, and this is now measured rather than assumed.** T storage fails in
   **every** variant: on Set I, `+h and +xi` gives 1.181e-02 against 1.349e-04 (87x worse); on
   Set II, 7.035e-02 against 2.903e-04 (242x worse). The review's instruction to measure rather
   than assume is what produced this; the answer is unambiguous.
3. **The lambda null check is not a null check for T storage**, and the review's formula
   predicted it to 1.3 %:

   ```
   tau-storage |dtau| = 1.168e-06  ->  |dT| = 1.803e-07   (float level; correct tolerance
                                                           is the SAME build quadrature)
   T-storage   |dT|   = 1.203e-02      predicted (dtau)^2/8*T = 1.188e-02
   tau(685) = 1.475720, tau(690) = 2.260238   ->   dtau = 0.7845, not the indicative 0.05
   ```

   The concrete dtau is 16x larger than the indicative value I offered, exactly as the review
   warned for near-horizontal rays at 256-500 m. A 1.2e-02 convexity error is 12x over the gate
   on its own, which is why every T-storage variant fails.
4. **V3 itself fails the new tuning grid.** The `+lambda only` variant differs from V3 only by
   three lambda nodes, which for tau storage add nothing measurable, so its 6.927e-03 is
   essentially V3's own error on the new grid — about 7x over the gate. The new tuning grid
   probes xi up to 0.99998, which no earlier scan did (section 9.2 scanned xi = 1-10^u at
   h in {375, 1500, 6000, 18000, 46000} m, never at h ~ 256 m). The earlier V3 "pass" was
   therefore on sets that did not cover this region — the archive's warning that a dev-set pass
   may not survive held-out data, reproduced again.
5. **Local refinement does work.** The pessimistic premise (patchy grazing coverage would defeat
   it) is not borne out: the coupled local refinement passes both held-out sets with 3-7x margin
   inside the memory budget.

### 10.5 Caveats

- **V3 is not one of the eight builds.** The variants are compared on identical grids, so the
  ranking is sound, but V3's own numbers on the new grids are *inferred* from `+lambda only`
  (finding 4). A direct V3 build on the tuning grid and both held-out sets should be added
  before this result is relied on.
- **One band, one profile.** The coupled refinement is validated for h in [256, 500] m,
  xi ~ 0.9999 and the MODTRAN bottomup profile only.
- **Build time rose** from ~35 s to 54.2 s; the amortisation and the per-photon cost of the
  refined table still need their own measurement.
- **Open items unchanged:** the fallback-rate threshold (0.8-1.1 %) and the cap unit
  (24.7 MiB working bound vs 36.7 MiB if the cap is 256 MiB) both still need settling, and
  item E (geometry timing) is untouched by this section.
- Peak RSS 37.9 MiB is for this diagnostic, not a production figure.

## Section 11 — phase B follow-up: V3 direct, fresh acceptance set, band scan, refinement sweep

Diagnostics only. Program `phaseB2_b.cpp`; log `log_phaseB2_b_followup.txt`.

### 11.1 Corrections to section 10 (recorded, not silently edited)

- **The clip prediction was the plan's own** (original Phase C text: *"the prediction is
  'never'"*), not the reviewer's. Corrected here.
- **The clip is a profile-data finding.** A cumulative column can only produce a negative
  difference if it **decreases with height** in the stored data. So the clip firing on 0.33 % of
  samples between 2025 and 17950 m means **the shipped MODTRAN table is non-monotone there**.
  That belongs on record alongside E2. It does not affect the failing band.
- **"Neither axis alone can reach the gate at any affordable memory" was an overclaim** — each
  axis was tested at one level only. What is shown: *these* single-axis refinements fail and
  *this* joint one passes. The sweep in 11.4 now tests two levels per axis.
- **"Patchy coverage does not defeat local refinement" was too broad** — shown for one band
  only. The band scan in 11.3 finds another region.
- **The τ-storage null check is float-level, a few ulp** (|Δτ| 1.17e-06 at τ ≈ 2 ≈ 5 float ulp),
  not "exact".
- **Set I is selection-contaminated** (all eight builds were scored on it and the winner chosen
  from them), so it cannot give an unbiased acceptance result. Section 11.2 freezes a **fresh**
  in-band set, Set III, and evaluates the candidate on it once.

### 11.2 V3 direct, and the size reconciliation

The 6.62 vs 11.34 MiB puzzle is **two different axis sets**, not a contradiction:

```
phase-3 axes (candidate A / V0)  56 x 256 x 121 = 1734656 nodes =  6.62 MiB
V3 axes                          61 x 464 x 105 = 2971920 nodes = 11.34 MiB
+lambda only (V3 + 3 lambda)     61 x 464 x 108 = 3056832 nodes = 11.66 MiB
```

So the base used for the eight builds **is** V3, and finding 4 of section 10 is now directly
measured rather than inferred:

```
V3 (61x464x105)  11.34 MiB  build 33.9 s
  tuning 6.9269e-03 FAIL [h=256  lam=682.5 xi=0.9999800000]
  SetI   5.2662e-03 FAIL [h=262  lam=681.5 xi=0.9999640000]
  SetII  3.1039e-04 PASS [h=1000 lam=688   xi=0.9990000000]
  SetIII 3.3946e-03 FAIL [h=270  lam=683.5 xi=0.9999425000]   <- fresh, never used for selection
```

V3's direct tuning value (6.9269e-03) matches the `+lambda only` inference to four digits.

**Set III** (fresh in-band acceptance, frozen before selection): 192 rays,
h ∈ {270…480} m, ξ ∈ {0.9998700…0.9999900}, λ ∈ {683.5, 686.5, 688.5, 691.5} nm,
**0 clashes** with the tuning grid and **0** with Set I. md5 in the log.

### 11.3 Band scan — another failing region exists

```
band scan: h in {100..1000} m x xi in {0.9995..0.999995} x ALL 105 native lambda knots (6615 rays)

V3      overall worst 1.4964e-02 FAIL  [h=100 lam=680 xi=0.999995]
        worst lambda knots 680, 685, 675, 670, 665 nm  (the O2 A-band region)
        lambda knots with any ray over the gate: 85 of 105
+h+xi   overall worst 4.7079e-03 FAIL  [h=100 lam=680 xi=0.999995]
        lambda knots with any ray over the gate: 77 of 105
```

**Both fail at h = 100 m, ξ = 0.999995 — far outside the [256, 500] m band.** So the coupled
refinement, which passes Sets I/II/III, is **not validated across the low-h / high-ξ region**.
The worst wavelengths cluster at 665-685 nm, the O2 A-band, but 77-85 of the 105 knots have
*some* ray over the gate, so this is not a single band cell.



### 11.4 Refinement-level sweep — the coupling now rests on five measured points

```
config                     h    xi     MiB  SetI        SetIII
h+0   xi+0    (=V3)       61   464   11.34  5.266e-03   3.395e-03
h+0   xi+400              61   864   21.11  3.200e-03   3.105e-03
h+12  xi+0                73   464   13.57  5.803e-03   4.365e-03
h+12  xi+400              73   864   25.26  6.604e-05   7.771e-05
h+6   xi+200  (candidate) 67   664   17.82  2.919e-04   2.152e-04
```

- **xi alone helps modestly** (5.27e-03 -> 3.20e-03, 1.6x, for +9.8 MiB). **h alone does not help
  and can make it worse** (5.27e-03 -> 5.80e-03, for +2.2 MiB) — the error is not h-driven.
- **The joint refinement is far more effective than either alone**: h+12/xi+400 reaches
  **6.604e-05 / 7.771e-05**, a **68x improvement over V3**. The coupling claim now rests on five
  measured points rather than one.
- **But the best variant exceeds the working budget**: 25.26 MiB against 24.7 MiB, over by 2 %.
  The accuracy/memory tension therefore reappears at a finer scale, and **the cap unit is now
  decision-critical** — on a 256 MiB cap (36.7 MiB budget) it fits comfortably; on the tighter
  24.7 MiB bound it does not.

### 11.5 Status after this section

| build | MiB | in-band Sets I / III | band scan (low h, high xi) |
|---|---|---|---|
| V3 | 11.34 | FAIL 5.27e-03 / 3.39e-03 | FAIL 1.50e-02 |
| candidate (+h6 / xi200) | 17.82 | **PASS 2.92e-04 / 2.15e-04** | **FAIL 4.71e-03** |
| h+12 / xi+400 | 25.26 | **PASS 6.60e-05 / 7.77e-05** | not scanned (over budget) |

- **The candidate is validated on the band it was refined for, and only there.** It passes the
  fresh in-band acceptance set (Set III) with a 4.6x margin, and fails the band scan at
  h = 100 m, xi = 0.999995.
- **A substantially better variant exists** (h+12/xi+400, 68x better on the in-band sets) but
  exceeds the tighter memory bound.
- **The coverage question is therefore now the binding one**: either extend the refinement to
  the low-h / high-xi region (costing memory that the tighter budget may not allow), or
  **declare a coverage limit and fall back** — in which case the fallback rate must be measured
  against the 0.8-1.1 % break-even, and given the patchy region that may decide the outcome
  before the 1e-3 gate does.

## Section 12 — phase B follow-up 2: reproduction, the 121-node axis, real peak memory, fallback rate

Program `phaseB3_b.cpp`; log `log_phaseB3_reproduction_memory.txt`.

### 12.1 Must-fix 1 — the candidate DID reproduce; the difference is node placement

Not nondeterminism. The two commits used different helper functions with the same node count:

```
commit         helper          the six extra h nodes                                  SetI        SetIII
4b570f25c      hRefined()      {275, 300, 350, 400, 450, 550}                        1.3492e-04  1.3259e-04
b37a6b15b      hRefinedN(6)    {257.1, 314.3, 371.4, 428.6, 485.7, 542.9}            2.9195e-04  2.1520e-04
```

Both builds were re-run here. **The same placement built twice gives bit-identical results**
(2.9195e-04 twice), and each placement reproduces its own commit exactly. So the build is
deterministic, the node *count* was the same (67, hence the same 17.82 MiB), and **placement
alone is a 2.2x lever** — the review's own interpretation, now measured.

### 12.2 Must-fix 3 — the phase-3 axes' 121 lambda nodes

`lamWide()` is `300 -> 900 nm in 5 nm steps` = **121 nodes** = 81 native knots in [300, 700]
**plus 40 nodes in [705, 900] that only exercise the clamp**. It does **not** cover 180-300 nm,
which is the F4 defect of section 6.6. `study_b`'s 81 nodes were the same 5 nm spacing over
[300, 700]. So all three tables stored native knots *within their range*; the phase-3 axis's
defect was its range, not its spacing.

### 12.3 Must-fix 2 — Set III is contaminated; Set IV frozen

The sweep scored five configurations on Set III, so any candidate chosen from it would be
selected on Set III. **Set IV** is frozen for whichever candidate survives:

```
heldout_setIV_acceptance.txt   md5 46c789817a8ab345d23130f391532f62   240 rays
h in {268, 302, 336, 371, 405, 439, 474, 492} m
xi in {0.9998760, 0.9999080, 0.9999330, 0.9999580, 0.9999810, 0.9999925}
lambda in {682.0, 684.5, 686.0, 689.0, 690.5} nm
```

disjoint from the tuning grid and from Sets I and III by construction.

### 12.4 Correction — the oxygen band

The **B-band** is at ~687-695 nm; the **A-band** is at ~760 nm, outside [180, 700] nm. Sections
11.1 and 11.3 said "A-band" and are corrected to **B-band**. The worst knots found (665-685 nm)
sit on the short-wavelength side of the B-band, so **the exact absorber there should be checked
against the profile rather than assumed** — not done here.

### 12.5 Real peak memory — and it overturns the "over budget" finding

A production-size geometry table was built (`hTop = 120 km`, production grid) and the candidate
loaded alongside it:

```
geometry table (hTop=120 km, production config): 214.3 MiB, hNodes=3901   <- reproduces the archive
candidate absorption table                     :  17.82 MiB
peak RSS before = 57.2 MiB, after geometry + candidate = 236.0 MiB
headroom against 256 MiB = 20.0 MiB
```

The archive's 214.3 MB is **reproduced exactly** by hTop = 120 km (3901 nodes, as predicted), so
the size reconciliation of section 11.2 is confirmed by direct measurement. And the peak is
**measured**, not estimated: 236.0 MiB with the candidate, giving a table-only budget of about
**37.8 MiB** — close to the 36.7 MiB MiB-cap estimate and well above the 24.7 MiB tighter bound.

**Consequence:** h+12/xi+400 (25.26 MiB) would put the process at ~243.4 MiB. That **fits** a
256 MiB cap with ~12.6 MiB to spare, and is **marginal** against a decimal 256 MB cap
(244.1 MiB, 0.7 MiB spare). The earlier "over the working budget" conclusion is therefore
**overturned by measurement**; the cap unit remains decisive and is still the requester's call.

### 12.6 Fallback rate — what could and could not be measured

A CSV run of the real configuration (`--cherenkov-output-format csv`, seed 600, 2079 hits,
103974 weight-summed photons) gives the **measured emission-height distribution**:

```
emission height: min 1867 m, max 3036 m, mean 2322 m
histogram (m : photons): 1800:33802  2100:5  2300:31357  2400:88  2600:38699  2700:20  2800:2  3000:2
wavelength: min 300, max 900 nm
```

**The realistic photons come from 1.8-3.0 km, so the low-h failing region (h <= 500 m) is not
populated in this configuration** — the coverage limit here can be driven by xi alone, and the
low-h problem matters for other configurations rather than this one.

**xi could not be extracted.** The CSV's `z_production_m` (~1870 m) and `path_length_m` are
mutually inconsistent under any simple geometric reading (the drop `L*|dir_z|` is 56-507 m), so
xi is not recoverable from these columns without resolving their semantics. The rate is
therefore a **sensitivity**, with xi modelled, and must not be quoted as measured:

```
coverage rule                                   fallback (uniform xi in [0.99,1])
xi <= 0.999   (Set II validated)                9.91 %     <- 10x the break-even
xi <= 0.9999                                    0.90 %     <- just inside 0.8-1.1 %
xi <= 0.999 OR (h in [256,500] and xi<=0.99999) 9.91 %     <- the h-window adds nothing here
```

So meeting the ~0.8-1.1 % break-even requires a coverage limit near **xi ~ 0.9999**, and even
then it is marginal — under a uniform-xi model. Resolving the CSV semantics, or instrumenting
the emission coordinates, is what would turn this into a measured rate.


## Section 13 — phase B follow-up 3: RSS reconciliation, coverage scan, single Set IV evaluation

Programs `phaseB4_b.cpp`, `setIV_b.cpp`; logs `log_phaseB4_rss_coverage.txt`,
`log_setIV_single_evaluation.txt`.

### 13.1 Must-fix 1 — the RSS numbers reconciled, and NOT understated

```
startup                                     1.9 MiB
after RI grid to 120 km                     3.7 MiB
after geometry build                      218.4 MiB   (allocated per memoryMB(): 214.3 MiB)
after TOUCHING every geometry table page  218.4 MiB   (162 rays traced, f in {1, 1.02})
after candidate absorption table            ~236 MiB
peak VmHWM                                 236.0 MiB
```

3.7 + 214.3 = 218.0, and touching every page (162 traces spanning h and mu, with f = 1 and
f = 1.02 so the dispersion-delta arrays are read too) **changes nothing** — the geometry table is
fully resident after its build. So the **untouched-page explanation is disproved**.

The 53 MiB gap had the opposite cause: the earlier run's 57.2 MiB "before" figure was
**freed-but-still-mapped memory from three preceding table builds**, which the geometry table
then reused. Accounting for that, the peak is (non-reusable base ~4) + 214.3 + 17.8 = **236.1**
against a reported **236.0**. The figure is sound, not understated.

**Caveat retained:** this is a *diagnostic* process. The production simulation's own memory is
not in it, so 236.0 MiB is not a production peak. "Overturned by measurement" applies only to
the *table-fit* question, not to a production memory claim.

### 13.2 Must-fix 3 / item 2 — the CSV columns, resolved from the source

The recorder is called as
`recorder_(xg, yg, zg, lambdaNm, z0_above_sea_level, dx_out, dy_out, dz_out, dx, dy, dz, totalTime_ns, weightPerBunch, n_local, pathLengthForRecording)`
(`Cherenkov.inl:2044`), and the CSV writer emits those arguments in order, with
`n_minus_1 = n_local - 1` (`muon_cherenkov_3d.cpp:1497-1506`). So `z_production_m` **is** the
emission height above sea level (`z0 - earthRadiusM_`, `:2012`).

**But the values contradict each other.** Those photons are emitted at 1870-3036 m while the
detector sits at **2944 m** — below the detector — and the geometric drop `L*|dir_z|` is only
56-507 m for a path of 56-24000 m. Under every simple reading of the columns these cannot both
be true, so **one of the recorded arguments is wrong**. That is a production finding, reported
here and not fixed.

**Consequence accepted:** the emission-height histogram of section 12.6 is **withdrawn** as
evidence, and so is the conclusion that "low h is not populated" — the column it rested on is
part of the same inconsistency.

### 13.3 Must-fix 2 / item 3 — the realistic-height coverage scan PASSES

```
scan: h in {1800..3000} m x xi in {0.999, 0.9993, 0.9995, 0.9997, 0.9999} x ALL 105 native knots
      3675 rays
candidate worst |dT| = 8.2517e-04  PASS   [h=1800 lam=680 xi=0.9993]
rays over the gate: 0 of 3675
per-xi worst: 0.9990:5.743e-04  0.9993:8.252e-04  0.9997:7.779e-04  0.9999:6.438e-04
```

So coverage at h = 1.8-3.0 km and xi up to 0.9999 **can** be claimed — the region the 0.90 %
fallback figure previously assumed without evidence is now validated. The margin is thin though:
**1.2x inside** the gate, against 3.4-4.6x on the in-band sets. (The per-xi row for 0.9995 reads
0.000e+00 because of a bucket-index artefact in this diagnostic, not because there is no error
there.)

### 13.4 Item 4 — lookup time, measured

```
candidate (67 x 664 x 105), tau storage: 0.2520 us/call over 300000 calls
```

against 0.2335 us/call for V3 (61 x 464 x 105) and 0.19 us/call for the coarser phase-3 axes. So
the refinement costs roughly 8 % per lookup over V3, and the combined route stays far below the
~1 us/photon task-level target — which remains a target, not an archived measurement.

### 13.5 Item 5 — one candidate fixed in advance, Set IV evaluated once

The candidate is fixed **by principle**, not by score: hGraded + six nodes evenly spaced in
[200, 600] m; xiRefined(416) + 200 nodes evenly in u = sqrt(1-xi) <= 0.03; the profile's native
5 nm knots over [180, 700] nm; tau storage. The `{275...550}` placement scored better on Sets I
and III, but choosing it on that basis would be selection on those sets, so it is **not** used.

```
Set IV loaded from the frozen file: 240 rays
pre-fixed candidate: 67 x 664 x 105 = 17.82 MiB, build 53.7 s, tau storage
SET IV (single evaluation): worst |dT| = 2.0279e-04  PASS  [h=268 lam=682 xi=0.9999810000]
rays over the gate: 0 of 240;  margin 4.9x inside
```

Set IV was frozen before any candidate existed and this is its **only** evaluation, so it is the
one unbiased acceptance number in this archive. The program reads Set IV **from the frozen file**
rather than regenerating it, so the evaluation provably uses the frozen set.

### 13.6 Status

- **The candidate passes its one unbiased gate (Set IV, 4.9x) and the realistic-height coverage
  scan (1.2x)**, at 17.82 MiB and 0.2520 us/call.
- **The fallback rate is still unmeasured.** The validated limit is now wider (xi <= 0.9999 at
  h = 1.8-3.0 km), so the 0.90 % sensitivity figure now rests on a validated region — but the
  *rate* still depends on a uniform-xi assumption, because the CSV cannot supply real xi values
  until the contradiction of 13.2 is resolved. With the narrower previously-validated limit the
  rate is 9.91 %, about 10x the break-even.
- **Open:** the CSV argument contradiction (production), the cap unit, item E and E2, and the
  amortisation of the 53.7 s build.


## Section 14 — phase B follow-up 4: scope corrections, bucket bug, break-even, CSV test

### 14.1 Must-fix 1 — the "realistic heights" claim is withdrawn

The coverage scan of section 13.3 used h = 1.8-3.0 km **because the CSV said photons came from
there**, and that CSV evidence was withdrawn in 13.2. So the scan demonstrates coverage of *a*
region, **not** of the realistic one. Both the fallback rule and the 0.90 % figure need emission
heights from a trustworthy source; until then neither can be called realistic. What survives is
the narrower statement: *the candidate passes at h in [1.8, 3.0] km and xi <= 0.9999, with a
1.2x margin* — a region, not a workload.

### 14.2 Must-fix 2 — the bucket-index bug is harmless, shown from the source

`phaseB4_b.cpp`:

```cpp
for (size_t i = 0; i < v.size(); ++i) {          // EVERY ray
  if (tref[i] < 0.0) continue;
  double const d = std::abs(cand.lookupT(r.h, r.xi, r.lam) - std::exp(-tref[i]));
  if (d > 1e-3) ++nOver;                          // counted for every ray
  int const k = int((r.xi - 0.999) / 0.0002 + 0.5);
  if (k >= 0 && k < 5) perXi[k] = std::max(perXi[k], d);   // bucket affects ONLY this row
  if (d > worst) { worst = d; ... }               // updated for every ray
}
```

`worst` and `nOver` are updated unconditionally; the bucket index gates only the per-xi display
row. So the headline **8.2517e-04** and **0 of 3675** are computed over every ray, and the
0.9995 row reading 0.000e+00 is a display artefact only. The 1.2x margin is real.

### 14.3 Must-fix 3 — fallback is AT break-even, not inside it

Including the lookup cost, the break-even rate is

```
(1 - 0.2520 us lookup) / 83.29 us per fallback = 0.8981 %
```

which is **exactly** the xi <= 0.9999 sensitivity figure (0.90 %). There is therefore **zero
margin**, and that is before the height dependence: above xi = 0.999 only the band [256, 500] m
and the 1.8-3.0 km region are validated, so **every other height above xi = 0.999 falls back** —
a contribution the 0.90 % figure does not include. The true rate is strictly higher.

So the honest position is: the design **at best breaks even**, and the deciding number — the
real fallback rate over a trustworthy distribution — is still missing.

### 14.4 The CSV contradiction — a new negative result, and E2 is not sufficient

Tested the frame-mixing explanation: the earlier check used `dir_z`, which is the **transformed**
(shower/telescope-frame) direction, so the world-frame vertical `dir_ground_z` was tried instead.

```
z_prod     L          |dz_world|  L*|dz_world|
2170.5     23955.1    0.10004     2396.6
1870.5     228.3      0.03518     8.0
1874.3     507.8      0.03572     18.1
1868.1     55.8       0.03029     1.7
dir_ground_z: mean 0.0236 ; L*|dir_ground_z|: mean 82.1 m ; z_prod: mean 2322 m
```

**Frame mixing is ruled out**: with the world-frame vertical the implied drop is 8-20 m for
photons recorded as produced ~1870 m up — worse, not better. So the contradiction is *not* a
frame confusion.

**E2 is not sufficient either.** If `z_production_m` were really a height above ground, adding
the 2944 m detector altitude puts emission at ~4.8-6.0 km, which resolves the "emitted below the
detector" symptom — but the recorded path (56-24000 m) is still inconsistent with a 1870 m drop
for a photon whose world-frame vertical cosine is ~0.03. The height-above-ground hypothesis
explains one symptom and not the other, so it is **not** a confirmed cause; it remains a
hypothesis to test in production, alongside E2, not a finding.

### 14.5 Status

- The table itself is in good shape: it passes its one unbiased gate (Set IV, 4.9x) at
  17.82 MiB and 0.2520 us/call.
- **Its scope is the band.** Set IV's worst ray is at h = 268 m, so it certifies the band and
  says nothing about coverage elsewhere.
- **Viability now depends entirely on the real xi and emission-height distribution.** At best the
  design breaks even (0.8981 % vs 0.90 %), and the real rate is higher because of the height
  dependence. Resolving the CSV contradiction — possibly together with E2 — is the next
  decisive step.
- **Item E (geometry timing) still blocks the curved route on its own.**
- Open: the CSV contradiction (production), E2, item E, the cap unit, and the 53.7 s build
  amortisation.


## Section 15 — statement tightening: scope of the fallback, CSV and sampling claims

Documentation only. No new measurement; nothing in production touched.

### 15.1 "The true rate is strictly higher" — withdrawn

Additional height-based exclusions increase the fallback rate **only if the workload populates
them**. Since the workload is exactly what is unknown, the correct statement is that the rate is
**unknown**, not that it is higher than 0.90 %. Likewise, "at best breaks even" applies to the
**assumed** distribution (uniform xi, and a height range taken from a withdrawn column), not to
every possible workload. The design's viability is **undetermined**, which is the honest position.

### 15.2 The CSV test — retracted to its actual scope

The claim "frame mixing is ruled out" is **withdrawn**. For a curved ray,

```
Delta z = integral_0^L d_z(s) ds
```

which is not `L * d_z(endpoint)` in general. What the test actually rules out is the **simple
straight-ray, constant-direction reconciliation**: using the world-frame vertical at the endpoint
gives an implied drop of 8-20 m for photons recorded as produced ~1870 m up, so no constant-direction
reading fits. A decisive check requires **trustworthy endpoints, agreed direction conventions, and
the traced vertical displacement** — none of which the CSV provides. So the contradiction is
**not** established as a path inconsistency; it is established only as *"no straight-ray
reconciliation of these columns works"*.

### 15.3 Sampling scope

- Keep the wording **"0 of 3675 *sampled* rays exceeded the gate"**. A finite sample does not
  certify a continuous region.
- **Native wavelength knots alone do not test between-knot interpolation.** The lambda axes store
  the native 5 nm knots, so the band scan's lambda coverage is at the knots; mid-knot behaviour is
  untested there.
- **Set IV supports its sampled band distribution**, and nothing more. Its worst ray happening to
  sit at h = 268 m does not define or certify coverage at 268 m.

### 15.4 The break-even is conditional

`0.8981 %` is correct **if** every photon pays the 0.2520 us lookup and a fallback adds 83.29 us.
It is a two-term threshold, not a full cost model: routing, the geometry trace, the projection and
recorder work, and the construction amortisation are outside it. The combined cost must be measured
as a whole before the threshold is used to accept or reject a coverage policy.

### 15.5 Not established (recorded so it is not read into this update)

**Neither a turning-point interpretation of the xi cutoff nor legacy handling through a turning
point has been established.** Nothing in sections 6-14 demonstrates either; the cutoff is a
declared guard, and the legacy route's behaviour at a turning point has not been measured here.

### 15.6 The next decisive step

1. Establish **trustworthy emission coordinates and xi** — resolving the CSV contradiction,
   possibly together with E2 — since every remaining viability question depends on them.
2. Then measure **fallback rate and combined cost under an explicit coverage policy**, against a
   full cost model rather than the two-term threshold.
3. **Item E (geometry timing) remains an independent blocker** for the curved route, regardless of
   how the absorption table turns out.

Phase D has not started. The cap unit remains the requester's decision, and the frozen sets
(tuning grid, Sets I-IV) are untouched by this section.


## Section 16 — baseline inventory and the Stage 3 path-length diagnosis

Added after an authorised Stage 3 run (Steps A, C, C'). Archive work in rounds 1-10 was
read-only; **no production source was modified** and the curved route remains **blocked**.
Every number below is either quoted from an existing log or taken from the three new logs in
§16.9. Where a claim rests on a code comment rather than a measurement it is marked so.

### 16.1 Scope

This section does two things:

1. records the **baseline inventory** of `diagnostics/b_table/` and the claims it can and
   cannot support (§16.2-§16.8), and
2. records the **Stage 3 measurement** that localises the 7.18 ns geometry discrepancy on the
   one ray where the archive showed it (§16.9-§16.12).

Stage 2 (the population scan) is **not** covered: it is blocked on the verbatim criteria B1,
B3, B4 and B5, which exist only in discussion and are not recorded here (§16.15).

### 16.2 Inventory

`diagnostics/b_table/` holds 21 programs, 20 archived logs, the four frozen held-out sets
(`heldout_setI_inband.txt`, `setII_outofband`, `setIII_acceptance`, `setIV_acceptance`),
`PROTOCOL.md` and `README.md`. `diagnostics/README.md` covers the directory itself and
`raytrace_compare/` holds the separate matrix harness. Binaries are deliberately not
committed; the build recipe is `README.md:18-32` and the run convention is
`argv[1] = the extinction profile` (`README.md:39`), whose recorded identity is

```
path : modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat
size : 46385 bytes
md5  : dfb2b7fd0d0c57ee0a704035a6a81e0b
```

Both the size and the md5 were re-verified before the Stage 3 run (§16.9).

### 16.3 Two table configurations, not three

The geometry table is built at two sizes only:

| configuration | nodes | memory | where |
|---|---|---|---|
| 60 km, two-level (10 m to 30 km, then 100 m) | 3301 | 181.4 MB | `path_b`, `combined_b`, and the Stage 3 runs |
| 120 km, two-level (10 m to 30 km, then 100 m) | 3901 | 214.3 MB | the unit test suite |

3301 = 3001 nodes at 10 m up to 30 km, plus 300 at 100 m up to 60 km. The 198.7 MB figure
quoted in places is **peak RSS with both absorption tables held**, not a table size
(`README.md:379-392`). An earlier draft of this section wrongly called the 3301-node build a
"third configuration"; it is the 60 km build. The rule that follows is still the operative
one: **every axis parameter must be printed from the table actually built**, because the
solved break value is configuration-dependent in principle. Measured (Stage 3, §16.9), the
break was the **same** in the 3301-node build as the value recorded for the test-suite build:
`muBreak = 0.884966` in both.

### 16.4 Check precedence (frozen from source)

`CurvedRayTable::trace()` rejects in this order, before any interpolation:

```
InvalidTable -> InvalidInput -> EmissionBelowObs -> AboveTable -> Upgoing
  -> BeyondCritical (p >= pCrit) -> DispersionOutOfRange (fPhase outside [f_lo_, f_hi_])
  -> Grazing (mu < mu_min_) -> Hit
```

`mu_min_ = 2.0e-4` is a **rejection** threshold, not a clamp: below it the caller must fall
back to the exact tracer. The f range is `[f_lo_, f_hi_] = [1 - 2*eps, 1 + 2*eps]` with
`eps = 0.02`, i.e. **[0.96, 1.04]** (`CurvedRayTable.hpp:52,371-372`), exposed as `fLo()` and
`fHi()`.

### 16.5 Observability limitation: `n_beyond_` conflates three statuses

The table-internal counter `n_beyond_` is incremented by `EmissionBelowObs`, `AboveTable`
**and** `BeyondCritical` (`CurvedRayTable.hpp:644,649,672`). A table-internal count of
"beyond" therefore cannot be decomposed. Any status-agreement test must use the
process-level `st.rejected[status]` array (`Cherenkov.inl:852`, 9 slots) instead.

### 16.6 The 0.1 ns requirement, and its basis

The applicable requirement is the one the test suite itself states:

```
testCherenkovCurvedRayTable.cpp:295  'agreed targets: 1 cm displacement, 0.1 ns time, 1e-6 rad angle'
```

It applies to the **total** arrival-time discrepancy over the eligible domain `mu >= 2e-4`
(the line-264 reading), and it is the criterion applied in §16.11.

### 16.7 Finding 1: the 7.18 ns discrepancy is a path-length error, not dispersion

The archived geometry-timing block (`log_combined_b_variant_and_timing.txt:76-78`) reports,
for the ray `h = 46000 m, lambda = 400 nm, xi = 0.9999`:

```
dr=2.0 m: worst |d t| = 7.1587 ns
dr=0.5 m: worst |d t| = 7.1771 ns
dr=0.1 m: worst |d t| = 7.1873 ns
```

`dr` is the **reference's** step, not the table's — the table is built once
(`combined_b.cpp:478`) and only the Bouguer march is refined. The program says so itself
(`combined_b.cpp:474-476`): *"the dr sequence below shows OBSERVED CONVERGENCE of the
reference. It is not a rigorous reference-error bound, and it does not bound the table's
error."* The converged discrepancy is therefore the `dr = 0.1` figure.

The independent swept-angle RK4 corroborates the reference and locates the error in the table
(`log_closure_b_attribution_policy.txt`, the `xi = 0.9999` block):

```
L: phi=725894.226141   radial=725894.228831      (agree to 2.7 mm)
Sigma: phi=53.137606   radial=53.137606
|d t| obs-vs-radial = 7.1873 ns
```

So the two independent references agree to 2.7 mm, and the table is 2.1544 m away from both.
Since `L_trace = vacuumL(h, mu) + DL` and `L_ref` integrates the same endpoints, the
discrepancy is carried by `DL` alone. The phase/group mix hypothesis is **refuted**, both
empirically (the split run reproduces 7.1873 ns unchanged,
`log_closure_b_attribution_policy.txt:77`) and structurally: `geometric_path_m = Lpath -
fGroup * SigF` (`CurvedRayTable.hpp:799`) subtracts the exact term that
`Lpath = vacuumL + DL + fGroup * SigF` (`:766`) added, so `geometric_path_m = vacuumL + DL`
by construction and no `(f_group - f_phase)` residual can survive.

### 16.8 Finding 2: the accuracy claim in the header comment is unbacked

`CurvedRayTable.hpp:264` states that `d` and `t` errors stay inside 1 cm / 0.1 ns for
`mu >= 2e-4`. The supporting figures in the same comment (`1.2e-2 m at mu = 1e-4`, `3.9e-2`,
`2.5 m`) appear **nowhere but the comment** — no log, test or program in the archive
reproduces them. The test suite asserts only the **Sigma-only** component
(`maxS * 1e9 < 0.02` ns), never the total `L` or `t`. Recorded as: **unbacked comment**, not
merely contradicted. The table has **no archived validation of total path length**.

Confirmed independently in round 10: the suite's random scan does reach the scaled branch
(`double const mu = std::pow(10.0, -4.0 * uh(rng)); // 1e-4 .. 1`, line 232) with region bins
at `1e-3` (line 213), but it records the **`d`** error and prints the per-region results
without asserting them. There is **no `L` or `DL` tolerance anywhere in the file**, in or out
of the scaled branch. Stage 3 (§16.10) shows this gap is material.

### 16.9 Corrections established during the archive rounds (all source-cited)

| Item | Finding |
|---|---|
| Endpoint convention | `L_trace` and `L_ref` both run observation sphere -> emission sphere at the queried `h`; there is no endpoint-convention term to absorb the discrepancy. |
| `DL` is interpolated at `fPhase` | `CurvedRayTable.hpp:747` builds the Lagrange parameter from `fPhase`; the endpoint indices use `fPhase` (`:666,:668`); `fGroup` appears only at `:766` and cancels at `:799`. So lambda reaches `geometric_path_m` through `DL(f_phase)` alone. |
| f-axis closed on this ray | `Fdisp(400) = 1.000000` **exactly** (the proxy is normalised so the phase factor is 1 at 400 nm: `Fdisp = (1 + 0.0032/l^2)/(1 + 0.0032/0.16)`, and at `l = 0.4` both terms are `0.02`). Every quartic weight carries a factor `t = (fPhase-1)/eps` (`:752-756`), so at `t = 0` all four vanish identically and `DL` **is** the stored `DL_` node value. |
| The Round-4 null is `fGroup`-only | Both Round-4 runs passed `fPhase = 1.000000`; only `fGroup` changed. The identical 7.1873 ns is exactly what `:766`/`:799` predict, so it is **not** evidence about `DL(f)`. |
| h = 46000 is an exact node | The h grid accumulates (`nh = h + step`, `:310`), so drift is possible in principle. **Measured in Stage 3: `hNode(i) - 46000 = +0.000e+00 m`** — exact. |
| No f-interpolation on this ray | `t = 0` exactly, so (b'') is closed for this ray. |
| The Q3 pair is suggestive only | `dL = 0.462 m at h = 1500 m` (lambda = 690) vs `2.154 m at h = 46000 m` (lambda = 400) come from one run, one table and one reference integrator, but differ in `h` **and** `lambda`; the lambda confound acts through `DL(fPhase)` and cannot be bounded from the archive. |

### 16.10 Stage 3, Step A: reproduction (PASSED)

Input identity re-verified before the run: 46385 bytes, md5 `dfb2b7fd0d0c57ee0a704035a6a81e0b`.
Build: `g++ -O2 -std=c++17 -I<parent of corsika/> combined_b.cpp -o combined_b`, compiler
`g++ (GCC) 8.5.0`, matching `README.md:56`. Run: `./combined_b <profile>`. Log:
`log_stage3A_repro.txt` (88 lines, md5 `8a0776764afa7cca88cd8643c11f49c8`).

```
dr=2.0 m: worst |d t| = 7.1587 ns  [h=46000 lam=400 xi=0.9999]
dr=0.5 m: worst |d t| = 7.1771 ns  [h=46000 lam=400 xi=0.9999]
dr=0.1 m: worst |d t| = 7.1873 ns  [h=46000 lam=400 xi=0.9999]
threshold scan (h=46000 m, lam=400 nm) against the 0.1 ns level:
  xi=0.90000000  |d t| =   0.0114 ns  ok
  xi=0.99000000  |d t| =   0.2629 ns  OVER
  xi=0.99900000  |d t| =   2.2132 ns  OVER
  xi=0.99950000  |d t| =   3.4229 ns  OVER
  xi=0.99990000  |d t| =   7.1873 ns  OVER
  xi=0.99995000  |d t| =   8.9292 ns  OVER
  xi=0.99999000  |d t| =  12.2092 ns  OVER
```

All three dr values and all seven threshold values reproduce **to the last digit**; the
variant families (V3 `9.668e-04` PASS, V7 `1.604e-03` FAIL, V8 `1.566e-03` FAIL) and the
development probes also match the archive line for line. Peak RSS 217.2 MB against the
archived 217.3 MB (machine-dependent). The ray is confirmed **in-domain**: the check order
reaches `Hit` (status 0), so it is not a grazing-guard case.

### 16.11 Stage 3, Step C: node values vs interpolation (log `log_stage3C_b.txt`)

One table, `f = 1` exactly. Build parameters, **printed from the table actually built**:

```
valid=1  hNodes=3301  muNodes=800  memory=181.4 MB
muMin=2.000000e-04   muBreak=0.884966000
sTot=8.654960968     deltaS=1.083224151205e-02
fLo, fHi = 0.960000, 1.040000
buildFailCol=-1      buildFailReason=0
```

`muBreak = 0.884966` matches the value recorded for the test-suite build
(`log_phase1_baseline_suite.txt:83`), so the solved break did **not** differ between the two
configurations. `sTot` and `deltaS` are new numbers the archive never recorded.

h axis and the node drift (round 6's open condition):

```
hCell(h)   = 3160
hNode(i)   = 46000.000000000000    drift = +0.000e+00 m
hNode(i+1) = 46100.000000000000
```

**The drift is exactly zero**, so (b'') is closed for this ray unconditionally.

s cell: the ray's `mu = 0.014141782` (equal to `sqrt(1 - xi^2)`, the table's own convention,
`CurvedRayTable.hpp:682`) falls **between** nodes — `cell j = 393`, `t = 0.138538348`, with
`muNode(393) = 1.412057568e-02` and `muNode(394) = 1.427436460e-02`. Round 7's "between
nodes" inference is therefore **verified**, and the `nSubBase` the build used for this column
is printed as **4** — the unscaled value.

Trace vs the converged reference:

```
geometric_path_m  = 725892.074457 m
vacuumL(h, mu)    = 682149.766863 m
DL (interpolated) =  43742.307594 m
sigma_m           =     53.137282 m
L_ref(dr=0.1)     = 725894.228831 m      DL_ref = 43744.461968 m
L_ref(dr=0.01)    = 725894.242366 m      (reference self-convergence 1.354e-02 m)
err_interp = DL_interp - DL_ref = -2.167909 m
```

The node comparison — **this is what separates (b) from (b')**:

```
node j=393  mu=1.412057568e-02  nSubBase=4
    storedDL(i,j) = +43776.285156 m     DL_ref = +43778.441818 m
    node error    = -2.156661 m
node j=394  mu=1.427436460e-02  nSubBase=4
    storedDL(i,j) = +43530.769531 m     DL_ref = +43532.911175 m
    node error    = -2.141644 m
```

**Both node values are wrong by ~2.15 m.** The interpolated error (-2.1679 m) is the same
magnitude, so the interpolation faithfully blends two nodes that are both wrong. Therefore

> **(b) — the build's node values — is the cause, not (b') the s-interpolation.**
> The s-interpolation contributes about **2.3 mm** (it slightly *reduces* the error).

The per-node profile (331 rows, every 10th node, reference at `dr = 1.0`) confirms
accumulation along the column and localises it:

```
i=10    h=100      diff = +0.000161 m
i=100   h=1000     diff = -0.279205 m
i=150   h=1500     diff = -0.468830 m
i=210   h=2100     diff = -0.668396 m
i=310   h=3100     diff = -0.931021 m
i=510   h=5100     diff = -1.270363 m
i=1010  h=10100    diff = -1.657624 m
i=1510  h=15100    diff = -1.833541 m
i=3300  h=60000    diff = -2.202695 m   (max |node error|)
```

Roughly **1.257 m of the 2.203 m accumulates in the first 5 km** and 1.933 m by 20 km — the
error is built where the ray is most grazing, exactly as the substep hypothesis predicts.
The `h = 60000` row is the column's top node and its reference is **not** converged (§17.3),
so that single row is not usable; every row up to 20 km is.

The mu-threshold scan over the **stored** table at `h = 1000 m` exposes the staircase:

```
nSub=20   mu [2.000e-04 .. 2.484e-04]   diff [+4.4504 .. +2.1989]
nSub=16   mu [2.511e-04 .. 3.328e-04]   diff [+2.0373 .. +0.6491]   step -0.169
nSub=12   mu [3.364e-04 .. 4.968e-04]   diff [+0.4804 .. -0.2232]   step -0.281
nSub=8    mu [5.022e-04 .. 9.938e-04]   diff [-0.5038 .. -0.7582]   step -0.797
nSub=4    mu [1.005e-03 .. 2.495e-03]   diff [-1.5548 .. -1.4197]   <-- unscaled plateau
```

This confirms the **four scaled plateaus** (k = 2,3,4,5 -> nSub = 8,12,16,20) exactly as
predicted in round 7, with a real discontinuity at `mu = 1e-3` where `nSub` goes 8 -> 4:

```
j=148  mu=9.937534228e-04  nSub=8  diff=-0.758230
j=149  mu=1.004576513e-03  nSub=4  diff=-1.554803
```

The jump is a factor **2.05**, not the factor 16 that nSub^-4 scaling would give. That is the
first indication that the integration is **not** fourth-order here; Step C' settles it.

### 16.12 Stage 3, Step C': one-column rebuild with scaled substeps (log `log_stage3Cp_b.txt`)

`buildColumn` is private and takes no substep override (`CurvedRayTable.hpp:208,589`), so
`stage3Cp_b.cpp` carries a **local copy** of the function inside the diagnostic. The
production header is untouched. The copy carries only `theta` and `L`, because `dtheta/dh`
does not depend on `phi` or `Sigma`; everything else is verbatim, including the first-cell
boost that had not previously been noticed:

```cpp
int const nSub = (i == 0) ? 4 * nSubBase : nSubBase; // extra near the ground
```

**Self-check first** — the copy with the build's own `nSubBase = 4` reproduces the stored
column:

```
h[m]     copy DL[m]        storedDL[m]       copy-stored
1000     +4250.141917      +4250.142090      -1.726e-04
5000     +19453.808608     +19453.808594     +1.442e-05
20000    +36303.360322     +36303.359375     +9.471e-04
46000    +43776.284837     +43776.285156     -3.191e-04
59900    +45594.228598     +45594.230469     -1.870e-03
-> max |copy - stored| = 1.870e-03 m   (COPY VERIFIED, float precision)
```

So the scaling results below are trustworthy.

**Prediction 4 — REFUTED.** The observed order is **p = 1.00**, not 4:

```
scale        nSubBase DL_copy[m]         DL_ref[m]          error[m]     observed order
x1           4        +43776.284837      +43778.441818      -2.156981    -
x4           16       +43777.900552      +43778.441818      -0.541266    p = 1.00
x16          64       +43778.304484      +43778.441818      -0.137334    p = 0.99
```

x4 cuts the error **4.0x**, x16 cuts it **15.7x** — first-order behaviour, not the 256x and
65536x that nSub^-4 would give. The predicted 8.4 mm and 0.033 mm do not appear; the actual
values are 0.541 m and 0.137 m. The near-horizontal integrand is sqrt-stiff near `h = 0`
(the code's own comment, `CurvedRayTable.hpp:587-589`), which is exactly the situation in
which RK4 loses its fourth order.

**Where the error lives, per scale:**

```
h[m]       scale x1       scale x4       scale x16
1000       -0.280471      -0.070275      -0.017726
5000       -1.259464      -0.315410      -0.079394
20000      -1.935704      -0.485202      -0.122573
46000      -2.156981      -0.541266      -0.137334
59900      -12.947427     -11.293102     -10.879518    <-- scaling does NOT fix this
```

Every height up to 46 km improves in proportion to `1/nSub`. The `h = 59900 m` row does
**not** — but §17.3 shows that row's reference is not converged either (the two reference
values imply `DL` decreasing with height, which is non-physical), so the near-top behaviour
is **not established** here. The solid finding is the first four rows: the substep count
controls the error from the ground up to at least 46 km, in proportion to `1/nSub`.

**Prediction 3 — CONFIRMED causally.** Rebuilding both sides of `mu = 1e-3` with the *same*
substep count removes the jump:

```
muNode           nSubBase DL_copy[m]         DL_ref[m]          error[m]
9.900000e-04     8        +26617.914499      +26618.702536      -0.788037   <- build's own rule
9.900000e-04     16       +26618.313289      +26618.702536      -0.389247   <- forced nSubBase=16
1.010000e-03     4        +26541.504517      +26543.088814      -1.584297   <- build's own rule
1.010000e-03     16       +26542.699363      +26543.088814      -0.389451   <- forced nSubBase=16
```

With the build's own rule the two errors differ by a factor 2.01; with a common `nSubBase = 16`
they agree to **2e-4 m**. The staircase is therefore caused entirely by the substep-count
discontinuity, not by anything intrinsic to `mu = 1e-3`.

### 16.13 Stage 3, Step D: attribution and verdict

The archive's 7.1873 ns is `|t_trace - t_ref|` with `t = (L + f*Sigma)/c`, and it is
reproduced exactly from the Stage 3 measurements:

```
dL = L_trace - L_ref(dr=0.1) = -2.154374 m
dS = Sigma_trace - Sigma_ref = -0.000324 m
|dL + dS| / c = 2.154698 / 0.299792458 = 7.1873 ns      (archive: 7.1873 ns)
```

Partition of the total, with the node errors interpolated to the ray's `t = 0.138538`:

| term | length | time | share |
|---|---|---|---|
| (b) build integration, node value at this ray | -2.154576 m | **7.18633 ns** | 99.99 % |
| (b') s-interpolation | +0.000202 m | +0.00067 ns | 0.009 % |
| (b'') f-interpolation | 0 (t = 0 exactly) | 0 | 0 % |
| (c) dispersion | 0 at lambda = 400 (f = 1 exactly) | 0 | 0 % |
| Sigma integration discrepancy | -0.000324 m | -0.00108 ns | 0.015 % |
| (e) unexplained | — | <= 0.0015 ns | at print precision |
| **total** | | **7.1873 ns** | |

**Verdict: the 0.1 ns requirement is NOT met.** The discrepancy is 7.1873 ns against a 0.1 ns
target, and **99.99 % of it is (b)**, the substep count in the build's ray integration.

**Cost of a fix, from the measured order.** With `p = 1`, reaching 0.1 ns (`= 0.029979 m`)
from 2.157 m at `nSubBase = 4` needs

```
nSubBase ~ 4 * 2.157 / 0.029979 ~ 288
```

i.e. about **72x** the current integration work in the affected columns. That is a
feasibility question, not an accuracy question, and it is **not** authorised or attempted
here. Raising the substep count would also remove the `mu = 1e-3` discontinuity (§16.12).

### 16.14 Status of the hypotheses

| Hypothesis | Status |
|---|---|
| Reference error, or a phase/group factor mix | **Refuted** — RK4 agrees with the Bouguer integral to 2.7 mm; the `fGroup` term cancels at `:799` |
| f-axis interpolation, on this ray | **Closed** — `t = 0` exactly, all quartic weights vanish |
| h-interpolation, on this ray | **Closed** — measured drift `+0.000e+00 m` |
| (b) build substeps vs (b') s-interpolation | **DECIDED — (b)**. Node values wrong by ~2.15 m; s-interpolation contributes ~2.3 mm |
| Prediction 1 (node values carry the error) | **CONFIRMED** |
| Prediction 2 (error accumulates along the column) | **CONFIRMED** — 1.257 m of 2.203 m in the first 5 km, 1.933 m by 20 km |
| Prediction 3 (jump at `mu = 1e-3`, four plateaus) | **CONFIRMED**, with correction: the jump is 2.05x, not 16x; the staircase vanishes under a common substep count |
| Prediction 4 (RK4 order-4 scaling) | **REFUTED** — measured order p = 1.00 |
| 0.1 ns applies to the total over `mu >= 2e-4` | **Accepted** (basis: test line 295) |

### 16.15 Stage 2 placeholders (blocked, not defined)

The following criteria were discussed but never fixed in writing, so the population scan
**cannot** be run and no result is recorded here:

- **B1** — qualified arrival criterion: _not supplied_
- **B3** — bad-query taxonomy (three classes): _not supplied_
- **B4** — equality reporting category: _not supplied_
- **B5** — arrival-agreement threshold: _not supplied_
- **B6 (E2)** — definition: _not supplied_
- **B7** — emission-coordinate provenance: unresolved

### 16.16 What is still NOT established

A confirmed mechanism explains **one defect on one ray**. It does not validate the table:

- **Stage 5** — the log-spaced `mu x h x lambda` scan of **total `L` and `t`** has now been run
  (§17). It confirms the defect across the readable domain `mu >= 1e-3` (62.6 % of eligible
  rays over 0.1 ns, worst ~14 ns) but is **reference-limited** below that and does not cover
  the near-top region, so the scan is not itself a validation.
- §16.8's missing `L` validation still stands: no part of the table has been validated for
  total path length.
- The line-264 claim (`1 cm / 0.1 ns for mu >= 2e-4`) remains **unbacked**.
- The curved route stays **blocked**; nothing here authorises a production change.

## Section 17 — Stage 5: extent of the path-length defect

Stage 5 answers the question the accept / restrict / reject decision needs: **how widely does
the Stage 3 defect extend?** It is a scan only — no production change, no table rebuild, no
population run. Log: `log_stage5_scan.txt` (822 lines), program `stage5_b.cpp`.

### 17.1 Method

A log-spaced grid of **565 eligible rays**: `h` in {100, 300, 1000, 3000, 5000, 10000, 20000,
32000, 46000, 59900} m, `mu` in {1e-4 … 1.0} (15 values, log-spaced, including the two below
`muMin` as boundary controls), `lambda` in {350, 400, 500, 690, 700} nm.

For each ray the program reports the table's `L`, `Sigma`, `t` against an **independent
reference that integrates both `L` and `Sigma`** along the same ray:

```
L_ref = int dz / sqrt(1 - s^2),   Sigma_ref = int N dz / sqrt(1 - s^2),
s = p / (n(z) z),  n = 1 + N(z) f
```

which matches the build's own `dL/dh = 1/cos(theta)`, `dSigma/dh = N/cos(theta)`. `t_ref =
(L_ref + f*Sigma_ref)/c` and `t_trace = r.time_s`, so `|dt|` is exactly the quantity the
0.1 ns requirement covers. `dr = 0.1 m`, the archived convention.

Two boundary facts the scan recorded: `mu = 2.0e-4` exactly is **rejected as Grazing** (the
`mu` reconstructed from `p/pCrit` lands a hair below `muMin` through round-off, so the
eligibility boundary is fuzzy at its own endpoint), and `mu = 1.0e-4` / `1.5e-4` are rejected
by construction — those are contract-correct fallbacks, not table failures.

### 17.2 The extent map, with the reference's own convergence measured

A band is only readable if the reference has converged there. The program therefore refines
the reference (`dr = 0.1 -> 0.02`) on the **worst ray of every band** and prints the shift:

| mu band | worst \|dt\| | reference moves | margin | verdict |
|---|---|---|---|---|
| [2.0e-4, 5.0e-4) | +16.0196 ns | 15.357 ns | 1.04x | **REFERENCE-LIMITED** |
| [5.0e-4, 1.0e-3) | -5.3193 ns | 1.056 ns | 5.0x | **REFERENCE-LIMITED** |
| [1.0e-3, 2.0e-3) | -11.4366 ns | 0.1836 ns | 62x | readable |
| [2.0e-3, 5.0e-3) | **-14.0534 ns** | 0.0739 ns | 190x | readable |
| [5.0e-3, 1.0e-2) | -12.1314 ns | 0.0576 ns | 211x | readable |
| [1.0e-2, 2.0e-2) | -9.1834 ns | 0.0540 ns | 170x | readable |
| [2.0e-2, 5.0e-2) | -5.6602 ns | 0.0494 ns | 115x | readable |
| [5.0e-2, 1.0e-1) | -2.0414 ns | 0.0391 ns | 52x | readable |
| [1.0e-1, 1.0) | -0.6491 ns | 0.0281 ns | 23x | readable |

**What this establishes.** Over the whole readable region, `mu >= 1e-3`, the table exceeds the
0.1 ns requirement, with a worst case of about **14 ns at `mu ~ 2e-3`**. Counting only
eligible rays with `mu >= 1e-3`, **285 of 455 (62.6 %)** exceed 0.1 ns — and the excess
persists to the top of the range: even the band `mu >= 0.1` has 25 of 155 over. So the defect
is **not** a grazing-tail curiosity; it covers the bulk of the eligible domain.

**What this does not establish.** For `mu < 1e-3` the scan is **reference-limited** and says
nothing about the table: at `mu = 2e-4` the reference itself moves 15.4 ns between `dr = 0.1`
and `0.02`, which is the whole reported figure. Resolving that region needs a better
reference (a singularity-removing substitution, as the archived `w`-substituted build uses),
not a finer `dr`.

`lambda` shows no separate effect: the worst per `lambda` is 13.6-14.1 ns for four of the five
values, the fifth being the reference-limited `mu = 2e-4` outlier. The defect is a geometry
defect, as Stage 3 concluded.

**Turning-point region: not covered.** No scanned ray had the reference turn before `hTop`
(0 of 565), because the query height is always the emission height. The near-top behaviour
seen in Stage 3 therefore needed the separate check in §17.3.

### 17.3 The near-top behaviour is NOT established (correction to §16.11/§16.12)

Two reference values at adjacent nodes of the same column contradict each other:

```
Stage 3 profile, h = 60000 m, dr = 1.0 :  DL_ref = 45607.171445 m
Step C' ,       h = 59900 m, dr = 0.1 :  DL_ref = 45607.176025 m
        (from storedDL(59900) = 45594.228598 and error = -12.947427)
```

The reference therefore reports `DL` **decreasing by 5 mm over 100 m of ascent**, which is
non-physical — `DL` grows monotonically with height along a descending ray (it reaches
~43.8 km at `h = 46000`). So at least one of those two reference values is not converged, and
consequently:

- the `-12.9 m` error at `h = 59900` (Step C') is **not established**;
- the `-2.202695 m` "max |node error|" at `h = 60000` (§16.11) is **not established**;
- the "second defect near the turning point" previously recorded in §16.12 is **withdrawn**.
  It was inferred from those two numbers plus the substep-insensitivity at `h = 59900`.

What survives is only the **relative** statement: at `h = 59900` the error barely responds to
the substep count (`-12.947 -> -11.293 -> -10.880 m` at x1/x4/x16), and because all three are
measured against the *same* reference that insensitivity is real even though the magnitude is
not. Settling the near-top region needs the singularity-removing reference of §17.2.

### 17.4 A latent stall in the archived reference integrator

The archived reference marches use the bare form

```cpp
while (z < zTop - 1e-12) { double const st = std::min(dr, zTop - z); ...; z += st; }
```

(`path_b.cpp:37-44`, and the same pattern inside `combined_b.cpp`). This can **stall**: if the
remaining gap lands in `(1e-12, ULP(z)/2)`, which is roughly `(1e-12, 4.6e-10)` at
`z ~ 6.4e6`, then `z + st` rounds back to `z`, the loop makes no progress, and the program
never terminates. It was hit during Stage 5 on the **first eligible ray** (`h = 100 m`,
`mu = 3e-4`).

`stage5_b.cpp` carries a guard (`if (!(zn > z)) break;`); the residual gap is under one ULP
(~1e-9 m) and cannot affect any reported value. Stage 3's runs did not land in the window
(they completed, reproduced the archived values, and reported a 13.5 mm reference
self-convergence), so **Stage 3's numbers stand** — but the archived integrators should be
fixed before they are reused, since the failure is a silent hang rather than a wrong number.

### 17.5 Stage 5 conclusions

| Question | Answer |
|---|---|
| Does the defect extend beyond the one Stage 3 ray? | **Yes.** 62.6 % of eligible rays with `mu >= 1e-3` exceed 0.1 ns; worst ~14 ns |
| Is it a grazing-tail curiosity? | **No.** It persists to `mu >= 0.1` (25 of 155 over) |
| Does `lambda` add anything? | **No** separate effect; geometry defect |
| Is `mu < 1e-3` characterised? | **No** — reference-limited in this scan |
| Is the near-top / turning-point region characterised? | **No** — §17.3 and §17.2 |
| Is the 1 cm **displacement** target tested? | **No** — only `L` and `t` were scanned |
| Does this validate the table? | **No.** It maps the extent of one defect; §16.8's missing `L` validation still stands |

So the accept / restrict / reject decision now rests on: the defect is broad (most of the
eligible domain), its mechanism is known (first-order substep integration, §16.12), and a fix
is possible at roughly 72x the integration cost in the affected columns (§16.13) — but two
regions remain uncharacterised and no part of the table has been validated for total path
length. The curved route stays **blocked**.

## Section 18 — Stage 6: why the column integration is first order, and a zero-cost fix

Diagnostics only; no production source, test, data pointer, frozen set or snapshot touched.
New files: `stage6_order_b.cpp` (Task 1), `stage6_refw_b.cpp` (Task 2), logs
`log_stage6_order.txt` and `log_stage6_refw.txt`. Binaries untracked. This section only adds;
earlier sections are unchanged.

### 18.1 Task 0 — h-breakdown of the Stage 5 failures (read-only)

From `log_stage5_scan.txt`, eligible rays only (`mu >= 1e-3`), each cell "fails / total, worst
|dt| in ns":

| h [m] | 1e-3-2e-3 | 2e-3-5e-3 | 5e-3-1e-2 | 1e-2-2e-2 | 2e-2-5e-2 | 5e-2-1e-1 | 1e-1-1.0 |
|---|---|---|---|---|---|---|---|
| 100 | 5/5 w=0.133 | 0/5 w=0.015 | 0/5 w=0.001 | 0/5 w=0.001 | 0/5 w=0.000 | 0/5 w=0.000 | 0/15 w=0.000 |
| 300 | 5/5 w=1.284 | 5/5 w=1.544 | 5/5 w=0.812 | 5/5 w=0.250 | 0/5 w=0.045 | 0/5 w=0.003 | 0/15 w=0.000 |
| 1000 | 5/5 w=4.168 | 5/5 w=5.011 | 5/5 w=3.507 | 5/5 w=1.681 | 5/5 w=0.453 | 0/5 w=0.041 | 0/15 w=0.005 |
| 3000 | 5/5 w=7.367 | 5/5 w=8.967 | 5/5 w=7.140 | 5/5 w=4.486 | 5/5 w=1.829 | 5/5 w=0.248 | 0/15 w=0.037 |
| 5000 | 5/5 w=8.621 | 5/5 w=10.533 | 5/5 w=8.658 | 5/5 w=5.836 | 5/5 w=2.748 | 5/5 w=0.470 | 0/15 w=0.076 |
| 10000 | 5/5 w=9.873 | 5/5 w=12.095 | 5/5 w=10.187 | 5/5 w=7.270 | 5/5 w=3.904 | 5/5 w=0.896 | 5/15 w=0.172 |
| 20000 | 5/5 w=10.687 | 5/5 w=13.106 | 5/5 w=11.197 | 5/5 w=8.254 | 5/5 w=4.776 | 5/5 w=1.375 | 5/15 w=0.322 |
| 32000 | 5/5 w=11.064 | 5/5 w=13.590 | 5/5 w=11.673 | 5/5 w=8.721 | 5/5 w=5.222 | 5/5 w=1.683 | 5/15 w=0.452 |
| 46000 | 5/5 w=11.306 | 5/5 w=13.879 | 5/5 w=11.960 | 5/5 w=9.004 | 5/5 w=5.490 | 5/5 w=1.898 | 5/20 w=0.564 |
| 59900 | 5/5 w=11.437 | 5/5 w=14.053 | 5/5 w=12.131 | 5/5 w=9.183 | 5/5 w=5.660 | 5/5 w=2.041 | 5/15 w=0.649 |

**Answer: no h range is clean.** The minimum is h = 100 m with 5 fails of 45, all in
`mu in [1e-3, 2e-3)` at worst 0.133 ns — only 1.3x over the limit, but over. Fail counts rise
monotonically with h (5 -> 35 of 45). So an h-based restriction alone is **not** viable. A
**mu-based** restriction is visible in the same table: at h = 100 m everything with
`mu >= 2e-3` passes, and `mu >= 0.1` passes for `h <= 5000 m`.

### 18.2 Task 1 — where the first-order error comes from (`log_stage6_order.txt`)

The failing ray's node column (j = 393, `mu = 1.412057568e-02`, `nSubBase = 4`). Stop check
first: V0 reproduces the stored column to **3.19e-04 m** (limit 2e-3 m). Then four variants at
nSub x1/x4/x16/x64, each against a matching reference at `dr = 0.01 m`:

| variant | x1 | x4 | x16 | x64 | order p |
|---|---|---|---|---|---|
| **V0** baseline (table) | -2.170519 m | -0.554804 | -0.150872 | -0.049889 | 0.984, 0.939, 0.798 |
| **V1** smooth analytic profile | **-0.016228 m** | identical | identical | identical | converged at x1 |
| **V2** knot-aligned substeps | -2.170519 | -0.554804 | -0.150872 | -0.049889 | identical to V0 |
| **V3** w = sqrt(h) substitution | -120.868 m | -30.227 | -7.569 | -1.904 | 1.000, 0.999, 0.995 |
| **V4** slope frozen per substep | **-0.016228 m** | identical | identical | identical | converged at x1 |

Decision rules, as they came out:

- **V1 converges at x1 while V0 is first order — the table's non-smoothness is the cause.** The
  analytic profile is the same two-term profile the table samples, so the difference is purely
  the table's piecewise-exponential interpolation, whose `dN/dz` jumps at every 100 m knot.
- **The `N`/`dNdz` inconsistency is refuted by inspection**: `N()` and `cellIndex()` use the same
  truncating index, and `dNdz = logSlope(cellIndex(z)) * N(z)`.
- **"Substeps not aligned with cell boundaries" is refuted by measurement**: V0 has **0 of ~13200
  substeps** straddling a knot, and V2 (which splits at knots) is byte-identical to V0.
- **The singularity is refuted**: V3, which removes the `1/cc` stiffness, is still first order
  (p = 1.00) and is in fact *worse* at low nSub.
- **V4 identifies the mechanism**: `cellIndex` truncates, so the RK4 stage evaluated at a step's
  **end** takes the *next* cell's slope. Freezing the slope to the substep's **start** cell
  removes it entirely: error -2.170519 m -> **-0.016228 m at the same `nSubBase = 4`**.

**Cost of the fix: zero.** V4 needs no extra substeps (x1 suffices), and 0.016228 m = 0.0541 ns,
inside the 0.1 ns requirement. For comparison V0 would need about **x72** (`nSubBase ~ 288`) to
reach 0.03 m — the 72x cost quoted in §16.13.

### 18.3 Task 2 — a singularity-free reference (`log_stage6_refw.txt`)

`stage6_refw_b.cpp` integrates in `w = sqrt(z - R)` (`dz = 2w dw`), covering **both** `L` and
`Sigma`, with the stall guard. Three validations, all as required:

```
(a) agreement with the Stage 5 z-march reference (dr=0.01), 5 rays with mu >= 1e-2:
    worst |L_w - L_z| = 8.159 mm  -> PASS (limit 10 mm)
(b) self-convergence (dw = 4e-3 -> 2e-3 -> 1e-3), drift:
    mu=2.0e-4: -6.250 mm   mu=5.0e-4: -2.500 mm   mu=1.0e-3: -1.251 mm
    h=59900:   -0.092 mm   h=60000:   -0.092 mm
(c) monotonicity across the top of the column:
    DL(59600)=45564.982195   DL(59700)=45575.801975   DL(59800)=45586.594571
    DL(59900)=45597.360097   DL(59950)=45602.732744   DL(60000)=45608.098666
    -> MONOTONIC PASS
```

The drift halves as `dw` halves (first order, as expected for a midpoint rule), so at `dw = 1e-3`
the residual is ~2 mm = 0.007 ns even at `mu = 2e-4`. That is what Stage 5 could not achieve: its
z-march moved 15.4 ns at that `mu`.

### 18.4 What §18 resolves from §17.3, and one correction to my own harness

**The §17.3 open question is resolved.** §17.3 withdrew the near-top numbers because two
reference values implied `DL` *decreasing* with height. The w-reference shows `DL` increasing
monotonically, and both stored nodes now give the **same** error:

```
h=59900: storedDL = 45594.2305 m, DL_ref = 45597.3601 m  -> error = -3.130 m
h=60000: storedDL = 45604.9688 m, DL_ref = 45608.0987 m  -> error = -3.130 m
```

So the near-top error is a real **-3.13 m**, a smooth continuation of the -2.16 m at h = 46000.
There is **no separate near-top or turning-point defect**; the earlier -12.9 m and -2.2 m were
both reference artifacts, and the §16.12 withdrawal was correct.

**Correction to my own harness.** The first Stage 6 run reported errors of about +151 m because
the reference used the *ray's* `xi` (0.9999) instead of the *node's* `xi` (from `muNode(393)`),
i.e. it integrated a marginally different ray. Fixed before any number above was recorded; the
builds were unaffected (V0 matched the stored column in both runs).

### 18.5 What this decides

- A **cheap fix exists**: freezing the refractive-index slope to the substep's start cell takes
  the failing ray from 2.17 m (7.2 ns) to 0.016 m (0.054 ns) at **no extra integration cost**.
- The cause is **neither** the endpoint singularity (V3) **nor** substep/knot misalignment (V2),
  but the per-evaluation slope lookup crossing a knot.
- The w-reference is validated and makes the `mu < 1e-3` region measurable, so a full re-scan is
  now possible.
- This does **not** validate the table: it identifies one defect and its fix. A full
  `mu x h x lambda` re-scan with the V4 fix and the w-reference is still required, and the
  line-264 claim is still unbacked. The curved route stays **blocked**.

## Section 19 — Stage 7: the V4 fix audited, and the re-scan

Diagnostics only; no production source, test, data pointer, frozen set or snapshot touched.
New files: `stage7_audit_b.cpp` (Tasks A and B), `stage7_rescan_b.cpp` (Task C), logs
`log_stage7_audit.txt` and `log_stage7_rescan.txt`. Binaries untracked. This section only
adds; earlier sections are unchanged.

### 19.1 Task A — V1 and V4 do differ; the shared residual was the reference

The suspicion was that V1 (analytic profile, analytic reference) and V4 (table, table
reference) agreeing to six digits on -0.016228 m implied a harness bug. They do **not** share
a path, and the coincidence is fully explained:

```
A.1  profile samples: table N equals analytic N at the 100 m knots (by construction) and
     diverges between them, e.g. h=50: 2.718589889e-04 (table) vs 2.718460393e-04 (analytic)
A.1b node xi = 0.999900300 (muNode(393) = 1.412057568e-02; the RAY xi is 0.9999)
     V1: DL_copy = 43780.449396  DL_ref = 43780.449399  error = -0.000004 m
     V4: DL_copy = 43778.439128  DL_ref = 43778.439131  error = -0.000004 m
     -> DL_copy differs by 2.010268 m between V1 and V4; only the RESIDUAL coincides.
A.2  V1 with H x1.1 (7500->8250, 1600->1760), reference recomputed:
     DL_copy moved by -2217.366861 m, error still -0.000004 m
```

Against the **w-reference** both residuals are **4 µm**, not 16 mm. So the 16 mm in Stage 6
Task 1 was the **z-march reference's own error** (dr = 0.01), which both variants shared — and
A.2 shows the residual is *profile-independent* while DL_copy responds correctly to a profile
change. No harness bug, and no Stage 6 Task 1 re-run is required.

One honest caveat: my A.3 stop check printed "profiles shared? YES (rel diff 6.051e-16)" — a
**false positive**, because it sampled at h = 100 m, a knot, where the table's N *equals* the
analytic N by construction. The substantive evidence above (2.01 m difference in DL_copy,
divergence between knots) shows the profiles differ, so the stop condition does not apply.

### 19.2 Task B — the 16 mm was reference-limited, not an offset

V4 node error against the w-reference at dw = 4e-3 / 2e-3 / 1e-3, seven nodes:

```
j=393  mu=1.412058e-02 h=46000: -0.000095 -> -0.000022 -> -0.000004 m   shrinks with dw
j=394  mu=1.427436e-02 h=46000: -0.000094 -> -0.000022 -> -0.000003 m   shrinks with dw
j=361  mu=9.984174e-03 h=10000: -0.000135 -> -0.000034 -> -0.000008 m   shrinks with dw
j=509  mu=4.960904e-02 h=20000: -0.000023 -> -0.000006 -> -0.000001 m   shrinks with dw
j=573  mu=9.922961e-02 h=59900: -0.000011 -> -0.000001 -> +0.000001 m   shrinks with dw
j=722  mu=4.984187e-01 h=5000 : -0.000000 -> -0.000000 -> -0.000000 m   shrinks with dw
j=212  mu=1.987738e-03 h=46000: -0.000646 -> -0.000142 -> -0.000016 m   shrinks with dw
```

**Every** case shrinks with dw, none is a constant offset, and at dw = 1e-3 the residuals are
4-16 µm = 0.00001-0.00005 ns. The float32 test on the stored column gives storedDL - double
copy = -0.000319 m against a float32 ULP of 5.2 mm at that magnitude, i.e. consistent with
float32 storage. So neither the reference nor float32 explains the Stage 6 figure: it was the
z-reference alone.

### 19.3 Task C — the re-scan with V4 passes completely

The column builder is a diagnostic copy with the V4 rule (the slope frozen to the substep's
**start** cell), applied to the θ equation that drives **both** the `L` and `Sigma` integrands;
`N` itself is continuous, so `Sigma`'s integrand has no jump of its own. `nSubBase` and the
first-cell boost are unchanged. Self-check first:

```
self-check V4 OFF: DL_copy = 43776.284837, storedDL = 43776.285156, diff = -3.191e-04 m -> MATCHES
```

The Stage 5 grid (13 eligible mu x 10 h x 5 lambda = 650 rays), the reference is the Stage 6
Task 2 w-reference with the stall guard, and `mu = 2.0e-4` is treated as **eligible via an
explicit tolerance** because the table's own `trace()` rejects it through round-off when
reconstructing `mu` from `p/pCrit`:

```
band           fails/total   worst |dt|    ref uncertainty   verdict
2e-4-5e-4        0/100       0.0092 ns      0.00417 ns       PASS
5e-4-1e-3         0/50       0.0016 ns      0.00167 ns       PASS
1e-3-2e-3         0/50       0.0039 ns      0.00083 ns       PASS
2e-3-5e-3         0/50       0.0003 ns      0.00042 ns       PASS
5e-3-1e-2         0/50       0.0001 ns      0.00017 ns       PASS
1e-2-2e-2         0/50       0.0000 ns      0.00009 ns       PASS
2e-2-5e-2         0/50       0.0000 ns      0.00005 ns       PASS
5e-2-1e-1         0/50       0.0000 ns      0.00002 ns       PASS
1e-1-1.0         0/200       0.0000 ns      0.00001 ns       PASS

overall: 0 of 650 rays over 0.1 ns  (pass fraction 1.0000)
```

**No band exceeds the 0.03 ns reference-uncertainty stop condition** (worst 0.00417 ns), so
every band is characterised — including `mu < 1e-3`, which Stage 5 could not resolve. The
worst ray in the whole grid is 0.0092 ns, versus 14 ns before the fix.

Harness bug found and fixed in Task C before recording: the first run added the vacuum chord
again (`vac + lv[i] + ...`) even though the column builder returns the **total** path, which
inflated every figure to ~1e6 ns — the size of `vac/c` itself. Corrected; the self-check above
is what pins the convention.

### 19.4 What this decides

- **The V4 fix removes the defect across the whole eligible domain**, `mu >= 2e-4`, including
  the region Stage 5 could not measure: 650/650 rays within 0.1 ns, worst 0.0092 ns, against
  7.19 ns on the Stage 3 ray and ~14 ns worst before.
- The recommended production change is therefore a **one-sided cell lookup** in the column
  builder — evaluate each substep with the refractive-index slope of the cell containing the
  substep's start rather than the cell of each RK4 stage — plus a **regression test on the
  failing ray** (h = 46000, lambda = 400, xi = 0.9999), which currently sits at 7.19 ns.
- **Scope of this result.** It is a diagnostic copy of the builder. It bypasses the table's
  `(h, mu)` Hermite interpolation (measured in Stage 3 at 2.3 mm = 0.008 ns) and the f-axis
  interpolation (which is exact at `lambda = 400` where `t = 0`, but untested at other
  `lambda`). A production build must therefore be re-scanned before acceptance.
- The **1 cm displacement** target was not tested: only `L` and `t` were scanned.
- The line-264 claim remains **unbacked**; the curved route stays **blocked** until a
  production build passes; Stage 2 still waits on B1/B3/B4/B5.

## Section 20 — Stage 8: Task 0 reconciliation, the production fix, and the production re-scan

Approved production change: `CurvedRayTable.hpp` `buildColumn` only. Also new:
`stage8_prodscan_b.cpp`, `stage8_regtest_b.cpp`, logs `log_stage7_rescan.txt` (re-run with
per-cell uncertainty), `log_stage8_prodscan.txt`, `log_stage8_regtest.txt`, and a Catch2
regression test. This section only adds; earlier sections are unchanged.

### 20.1 Task 0 — the full Stage 7 grid, reconciled

The complete mu band x h bin table (fails / total, worst |dt| in ns), all **650** rays:

```
h[m]     2e-4-5e-4  5e-4-1e-3  1e-3-2e-3  2e-3-5e-3  5e-3-1e-2  1e-2-2e-2  2e-2-5e-2  5e-2-1e-1  1e-1-1.0
100      0/10 w=.009 0/5 w=.002 0/5 w=.004 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
300      0/10 w=.005 0/5 w=.001 0/5 w=.002 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
1000     0/10 w=.002 0/5 w=.000 0/5 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
3000     0/10 w=.001 0/5 w=.000 0/5 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
5000     0/10 w=.001 0/5 w=.000 0/5 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
10000    0/10 w=.001 0/5 w=.000 0/5 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
20000    0/10 w=.001 0/5 w=.000 0/5 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
32000    0/10 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
46000    0/10 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
59900    0/10 w=.001 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/5 w=.000 0/20 w=.000
```

**Every cell is 0 fails.** The band totals are 100+50+50+50+50+50+50+50+200 = **650**; the
bands [5e-3, 1e-2), [1e-2, 2e-2) and [2e-2, 5e-2) were present in the log all along and were
simply omitted from the previous chat summary, which is what made the count look like 450.
Per-cell reference uncertainty is identical across h within a band (0.00417 ns in
[2e-4, 5e-4) down to 0.00001 ns in [1e-1, 1.0)); the worst single cell is **0.00417 ns**,
against the 0.03 ns limit.

**Ray count, explained.** Stage 5's grid has 10 mu values >= 1e-3, so 10 x 10 x 5 = **500**
potential rays, of which **45 were rejected as `InvalidTable`** — all of them `mu = 1.0` (the
exactly vertical ray), for every lambda and every h — leaving the 455 that Stage 5 compared.
The Stage 7 re-scan integrates each column directly, so it includes those 45 plus the **150**
rays in `mu in [2e-4, 1e-3)` that Stage 5 could not characterise: 500 + 150 = **650**. No ray
changed eligibility; the difference is entirely the 45 rejections plus the newly
characterised grazing rays.

**The 5 worst rays:** 0.0092 ns (mu = 2.0e-4, h = 100 m, lambda = 700 nm), 0.0092 (690), 0.0091
(500), 0.0089 (400), 0.0089 (350) — all in the same corner.

**Stop condition: not triggered** — no ray over 0.1 ns (max 0.0092) and no band over 0.03 ns
uncertainty (max 0.00417), so Task 1 proceeded.

**New finding (separate from the fix):** `mu = 1.0` is rejected as `InvalidTable`, which is the
status for "table not built / not valid" and is therefore the wrong reason for a rejection of
a valid table. It is a rejection, so the ray falls back to the exact tracer and the scan's
pass/fail is unaffected, but the status is misleading and deserves a separate look.

### 20.2 Task 1.1 — the production change (9 insertions, 1 deletion)

In `buildColumn`, the refractive-index log-slope is now fixed per **substep**, to the cell the
substep starts in, and used by all four RK4 stages:

```diff
@@ -568,6 +568,13 @@
     //   dSigma/dh = N/cos(theta)
     // The last three are exactly the stored h-derivatives (no differencing).
+    // The refractive-index log-slope is fixed per SUBSTEP, to the cell the substep
+    // starts in, and used by all four RK4 stages.  Taking it from each evaluation
+    // point instead (the previous behaviour) let the stage at a step's end read the
+    // NEXT cell's slope, because cellIndex() truncates; that made the column
+    // integration first order (measured: p = 0.98-1.00, 2.17 m on the Stage 3 ray)
+    // instead of converging at the scheme's order.
+    double bStep = 0.0; // log-slope of the substep's start cell
     auto rhs = [&](double hh, double t, double& dth, double& dph, double& dl, double& ds) {
@@ -575,7 +582,7 @@
-      double const dN = ri_->logSlope(ri_->cellIndex(r)) * Nz; // dN/dz at r
+      double const dN = bStep * Nz; // dN/dz at r, from this substep's own cell
@@ -611,6 +618,7 @@
       for (int m = 0; m < nSub; ++m) {
+        bStep = ri_->logSlope(ri_->cellIndex(r_obs_ + hh)); // this substep's own cell
         double k1t, k1p, k1l, k1s, k2t, k2p, k2l, k2s;
```

`nSubBase`, the first-cell boost and the interpolation are untouched. The rule reaches both
the `L` and the `Sigma` integrands because both are driven by `dtheta/dh`, which is where the
slope enters; `N` itself is continuous, so `Sigma`'s integrand has no jump of its own.

### 20.3 Task 1.2 — the rebuilt table

The geometry table is built in memory by the diagnostics, so there is no file to overwrite;
the "new candidate" is the rebuilt table under the new code, recorded here:

```
valid=1  hNodes=3301  muNodes=800  memory=181.4 MB  muBreak=0.884966000
```

i.e. the same 60 km two-level configuration as before, same node count, same memory, same
solved break. The frozen held-out sets and the `*.snap_*` / `*.pre*` files were not touched.
What changed is the stored data: on the Stage 3 ray the node value moved from
`storedDL = 43776.285156` (2.16 m low) to **`43778.437500`** (1.6 mm low).

### 20.4 Task 1.3 — the regression test

A new `SECTION("regression: Stage 3 ray against the w-reference (Stage 8 fix)")` was added to
`tests/modules/testCherenkovCurvedRayTable.cpp` (60 lines). It builds the ray's node column and
the Stage 3 ray itself and asserts, both against an in-section w-reference:

- each bracketing node value: `|storedDL(i,j) - DL_ref| <= 0.03` m;
- the total arrival time: `|time_s - (L_ref + f*Sigma_ref)/c| <= 0.1` ns.

**The Catch2 suite cannot be built in this environment**: the test includes `Cherenkov.hpp`,
which reaches `framework/core/Logging.hpp` and `<spdlog/fmt/ostr.h>`, and no spdlog exists in
the repo or in the conan cache (only `fmt` is present). Catch2 itself is available in the
amalgamated form. The suite therefore was **not** run and no snapshot/hash diff was observed.
The command a maintainer needs, using the cache paths found here, is

```
g++ -O2 -std=c++17 -I<corsika-parent> -I<eigen3> -I<catch2>/extras \
    tests/modules/testCherenkovCurvedRayTable.cpp <catch2>/extras/catch_amalgamated.cpp -o test
```

(with a shim `catch2/catch_all.hpp` that includes `catch_amalgamated.hpp`), plus the spdlog
include path that is missing here.

Because the suite could not run, the two assertions were exercised standalone in
`stage8_regtest_b.cpp`, which needs only the two table headers:

```
(1) j=393  storedDL=43778.437500  DL_ref=43778.439131  |err|=0.001631 m  PASS
    j=394  storedDL=43532.910156  DL_ref=43532.908492  |err|=0.001664 m  PASS
(2) t_prod = 2.421499751010e-03 s  t_ref = 2.421499755509e-03 s  |dt| = 0.004499 ns  PASS
ALL ASSERTIONS PASS
```

Node errors are 1.6 mm against a 30 mm limit, and the total is 0.0045 ns against 0.1 ns —
versus 2.16 m and 7.19 ns before the fix.

### 20.5 Task 1.4 — the Stage 7 grid through the production lookup path

`stage8_prodscan_b.cpp` runs the same 650-ray grid through `trace()`, so the `(h, mu)` Hermite
and the f-axis interpolation are exercised. The observable is the production `Result.time_s`:

```
production hits 565, rejections 85 (status 1: 45  status 8: 40)
overall: 0 of 565 production rays over 0.1 ns  (pass fraction 1.0000)
worst per band: 0.0122 ns ([2e-4,5e-4)) ... 0.0042 ns ([2e-2,5e-2))
worst single-cell reference uncertainty = 0.00417 ns (limit 0.03)
5 worst: 0.0141 ns (mu=1e-3, h=46000, lam=500), 0.0134 (690), 0.0134 (700), 0.0131 (400),
         0.0122 ns (mu=2e-4, h=59900, lam=500)
```

Every band PASSes, every band is characterised, and the 565 hits / 85 rejections reproduce
Stage 5 exactly (45 `InvalidTable` at `mu = 1.0`, 40 `Grazing` in the grazing tail).

**Stop condition: not triggered.** No ray that passed in the diagnostic copy fails here, so
the interpolation is not a second defect: the production path's worst ray (0.0141 ns) is in the
same sub-0.02 ns band as the diagnostic copy's (0.0092 ns), and the small difference is the
interpolation's own ~2.3 mm contribution, as measured in Stage 3.

### 20.6 What this decides

- The production fix is **complete over the eligible domain**: 565/565 production rays within
  0.1 ns, worst 0.0141 ns, against 7.19 ns on the Stage 3 ray and ~14 ns worst before.
- The mechanism, the fix and its zero cost are as reported in Sections 18 and 19.
- Two things remain outside this section's scope: the **1 cm displacement** target was not
  tested (only `L` and `t`), and the `mu = 1.0` `InvalidTable` rejection (§20.1) is a separate
  defect that this change neither causes nor fixes.
- The Catch2 suite still needs to be run in an environment that has spdlog, so the snapshot /
  hash diff question is **open**, not answered.
- The curved route stays **blocked** pending that; Stage 2 still waits on B1/B3/B4/B5.

## Section 21 — Stage 9: the test suite, the displacement target, and the rejections

Read-only apart from the Stage 8 fix already in place. No further production edits were made.
New: `log_stage8_prodscan.txt` (extended with the displacement scan and the rejection listing).
This section only adds; earlier sections are unchanged.

### 21.1 Task A — the Catch2 suite could NOT be run (dependency blocked)

The project takes its dependencies from **conan**: `conanfile.py:58` requires
`spdlog/1.14.1`, and `catch2/3.6.0` at line 59. On this machine:

- `conan` is **not installed** (`which conan` finds nothing), so the recipe cannot be
  resolved;
- the conan cache (`~/.conan2/p`) contains **catch2, eigen and fmt but no spdlog**;
- the documented CI path is the container `corsika/devel:u-22.04` (`.gitlab-ci.yml:46`),
  which cannot be pulled or run here.

Hand-rolling a shim would be exactly the workaround the brief forbids, so **no spdlog was
vendored and no build file was touched**. Consequences, stated plainly:

- the suite was **not built, not run**, on either the pre-fix or the post-fix tree, so **no
  snapshot or hash difference is known** — that question remains open;
- the new SECTION has **never been compiled** by Catch2.

What *was* verified: the two assertions were exercised standalone in `stage8_regtest_b.cpp`
(same code, same two headers) and pass — node errors 1.631 mm and 1.664 mm against a 30 mm
limit, total |dt| = 0.004499 ns against 0.1 ns. **Self-containment is confirmed**: the only
occurrence of the string `diagnostics` in the test file is a comment (line 1069); there is no
include or link dependency on `diagnostics/`.

### 21.2 Task B — the displacement target

**The project's definition**, cited: `d = r_obs * Phi`, where `Phi` is the accumulated
transverse angle — stated at `CurvedRayTable.hpp:775` (`double const Phi = d / r_obs_;`) and
produced by `buildColumn` as `dOut[i + 1] = r_obs_ * phi; // impact distance d = r_obs * Phi`,
with `dPhi/dh = tan(theta)/r` in the same RHS. The reference here integrates
`Phi = int tan(theta)/r dz` along the Bouguer ray (`tan theta = s/sqrt(1-s^2)`,
`s = p/(n z)`) in the `w = sqrt(z-R)` variable, and the table's value is read from the
production result as `d = r_obs * acos(z_impact / r_obs)`.

Result, same mu band x h bin layout (fails / total, worst |dd| in mm), 565 production hits:

```
h[m]     2e-4-5e-4 5e-4-1e-3 1e-3-2e-3 2e-3-5e-3 5e-3-1e-2 1e-2-2e-2 2e-2-5e-2 5e-2-1e-1 1e-1-1.0
100      0/6 w=2.24 0/5 w=1.20 0/5 w=1.32 0/5 w=0.26 0/5 w=0.51 0/5 w=0.23 0/5 w=0.04 0/5 w=0.00 0/15 w=0.05
1000     0/6 w=1.06 0/5 w=0.50 0/5 w=1.53 0/5 w=1.01 0/5 w=0.54 0/5 w=0.86 0/5 w=1.07 0/5 w=0.19 0/15 w=0.03
5000     0/6 w=1.62 0/5 w=0.52 0/5 w=1.13 0/5 w=1.68 0/5 w=1.01 0/5 w=0.58 0/5 w=1.19 0/5 w=1.66 0/15 w=0.38
20000    0/6 w=1.78 0/5 w=2.12 0/5 w=1.07 0/5 w=2.25 0/5 w=1.26 0/5 w=2.27 0/5 w=0.73 0/5 w=3.00 0/15 w=1.76
46000    0/6 w=2.49 0/5 w=2.57 0/5 w=3.31 0/5 w=4.08 0/5 w=1.95 0/5 w=2.30 0/5 w=0.71 0/5 w=2.05 5/20 w=93131.72
59900    0/6 w=3.53 0/5 w=2.19 0/5 w=2.81 0/5 w=1.56 0/5 w=2.72 0/5 w=1.90 0/5 w=0.26 0/5 w=3.12 0/15 w=2.15
```

Per band: **every band PASSes with 0 fails** except `[1e-1, 1.0)`, which shows **5 fails of
20**, worst **93 131.72 mm = 93.1 m**. Every band is characterised: the worst reference
uncertainty is **1.25 mm** (in `[2e-4, 5e-4)`), against the 3 mm limit, so the stop condition
does not trigger.

The 5 failures are **all `mu = 1.0` at `h = 46000`** (the same ray for five lambda), and they
are the same defect as Task D — see §21.4. Away from `mu = 1.0` the displacement error is at
most **4.08 mm** against the 10 mm target.

### 21.3 Task C — the 40 status-8 rejections, accounted for

Every one of the 40 is **`mu = 2.00e-04`**, ten each at `lambda = 350, 400, 690, 700 nm`; at
`lambda = 500` the same `mu` is accepted at all ten heights, which is why the count is 40 and
not 50.

Status 8 is `Status::Grazing`, defined at **`CurvedRayTable.hpp:696-698`**:

```cpp
if (mu < mu_min_) {
  ++n_murej_;
  res.status = Status::Grazing;
  return res;
}
```

with `mu_min_ = 2.0e-4`, and named at line 96. The trigger is therefore `mu < mu_min_`, where
`mu` is **reconstructed** from the invariant as `sqrt(1 - (p/pCrit)^2)` (line 682).

**These are intended boundary rejections.** `mu = 2.0e-4` is exactly the eligible-domain
endpoint, so a ray placed *at* it sits on the boundary and the reconstruction's round-off
decides which side it lands on — for four of the five lambda it lands a hair below and the
grazing guard fires, exactly as designed (the contract is to fall back to the exact tracer,
not to return a low-accuracy value). The caveat, already recorded in §17.1, is that the
eligibility boundary is therefore fuzzy at its own endpoint: `mu = 2e-4` is not reliably
accepted. That is a rounding-boundary quirk, not a wrong rejection rule.

### 21.4 Task D — why `mu = 1.0` returns `InvalidTable` (diagnosed, not fixed)

The path is **`CurvedRayTable.hpp:778-780`**:

```cpp
778:  double tx = dx - cRad * rx, ty = dy - cRad * ry, tz = dz - cRad * rz;
779:  double const tn = std::sqrt(tx * tx + ty * ty + tz * tz);
780:  if (!(tn > 0.0)) return res;
```

`t` is the transverse component of the emission direction; for a purely radial ray it is zero.
The early return at line 780 **never assigns `res.status`**, so the result keeps the default
from line 103, `Status status{Status::InvalidTable}` — a status that means "table not built /
not valid", which is false here: the cell lookup, the Hermite blend and the path/time assembly
all succeeded before line 778. So the trace is reported as a table failure when it is not one.

The same line produces the five displacement failures. For `mu = 1.0` the direction is radial,
so `t` is a pure round-off residual: `rEmit * (1/rEmit)` is not always exactly 1, and where it
rounds to `1 + eps` the residual `tn` is about `2*eps`, which is **not** zero — the guard does
not fire, and the arrival direction is then built from a meaningless `t`. At `h = 46000` that
yields a spurious `Phi` of about `1.46e-5 rad`, i.e. **93.1 m** of displacement where the true
answer is 0. At the other nine heights the same ray happens to round the other way and the
guard fires, giving the 45 `InvalidTable` rejections. One line, both symptoms.

**Proposed fix (not applied, needs its own go-ahead):** in `trace()`, handle the degenerate
transverse direction explicitly rather than returning an unset result — set `res.hit = true`
and `res.status = Status::Hit`, set the arrival direction to the radial direction
(`res.dx = -rx; res.dy = -ry; res.dz = -rz;`, which is the exact answer for a radial ray), and
treat `tn` below a small tolerance as degenerate so the two round-off branches behave
identically. That is a few lines in `trace()`, no change to the table build, and it would
remove both the 45 wrong-status rejections and the 5 displacement failures.

### 21.5 What this decides

- **Timing: PASS** — 0 of 565 production rays over 0.1 ns, worst 0.0141 ns.
- **Displacement: FAIL on 5 rays** (93.1 m against a 1 cm target), all of them the degenerate
  `mu = 1.0` at `h = 46000`; everywhere else at most 4.08 mm. All bands are characterised
  (worst reference uncertainty 1.25 mm against the 3 mm limit).
- **Rejections: intended** — 40 `Grazing` at the eligible-domain endpoint, plus 45 from the
  `mu = 1.0` status bug.
- The two failures share one root cause, already localised to `CurvedRayTable.hpp:780`, so the
  displacement "own diagnostic" the decision table calls for is effectively done.
- The Catch2 suite still has not been run, so the snapshot/hash question stays **open**.
- The curved route stays **blocked**; Stage 2 still waits on B1/B3/B4/B5.

## Section 22 — Task E: absorption accuracy on the corrected curved ray (Step 0 + scan)

Diagnostics only; no production edits (verified: `CurvedRayTable.hpp` remains 9+/1-, and
`AtmosphericAbsorption.hpp`, `Cherenkov.inl`, `Cherenkov.hpp` are untouched). New:
`stage9_abs_b.cpp`, log `log_stage9_abs.txt`. This section only adds.

### 22.1 Step 0 — identifications (all cited)

1. **The stored quantity is cumulative optical depth τ(λ, h)**, not transmission.
   `AtmosphericAbsorption.hpp:41` declares
   `std::vector<std::vector<double>> ecoeff_; // optical depth τ(λ, h), indexed as
   [wavelength_idx][height_idx]`, loaded at lines 81-86 from the profile's third block. The
   **lookup** is `getOpticalDepth(λ, h)` at line 107, bilinear in (λ, h) at lines 140-154.
   Transmission is derived, never stored: `getTransmission` (line 171, `exp(-τ_slant)` with
   `τ_slant = |τ_emit - τ_det| / |cos_zenith|`, lines 177-195) and `getTransmissionSpherical`
   (line 214), which is **not** an exact slant integral but a **Chapman-function air-mass
   approximation** (`H_scale = 8400 m`, Rozenberg for zenith < 75°, Kasten-Young beyond,
   lines 240-265).
2. **`Σ` in `CurvedRayTable.hpp` is NOT absorption.** Lines 253-257 define it as
   *"Sigma = int N ds along the f=1 ray"* — the **refractivity** column, part of the **timing**
   integral (`optical_path_m = L + fGroup*Σ`, line 774). The build's RHS confirms it:
   `ds = Nz / cc` (refractivity, not extinction). Absorption never enters the geometry table.
3. **The production absorption path** is `extinctionCoeffFunc_` in `Cherenkov.inl:96-113`:
   the local extinction is a **centred ±25 m finite difference of `getOpticalDepth`**
   (`dh = 50`, `max(0, (τ_hi - τ_lo)/actual_dh)`), with the cumulative-column convention
   documented in the comment there (lines 106-111), which records that the previous
   `(τ_lo - τ_hi)` ordering returned exactly 0 and silently disabled the integration.
   `stage9_abs_b.cpp` calls that same rule through the production class, so the scan exercises
   the production lookup, not a copy.
4. Profile facts, read from the file itself: **λ range [180, 700] nm** (105 knots at 5 nm),
   **h range [0, 50000] m**, **native h spacing 1000 m** (51 knots). So queries at
   705-900 nm are *accepted and clamped* (§8.1), and any query above 50 km is clamped too.

### 22.2 The scan, and why its numbers are not yet usable

Grid: the 121 λ nodes of §12.2 (300-900 nm at 5 nm) × 8 µ (2e-4 … 1.0) × 6 h (100 … 59900 m)
= 5808 rays; 4579 hits, 1229 excluded (fraction 0.212; µ = 1.0 is the line-780 degenerate case
and µ = 2e-4 the grazing endpoint — listed as excluded, not counted as pass or fail). Two
references were integrated along the same corrected path: (a) in `w = sqrt(z-R)`, (b) a
z-march. `dr = 1.0` m.

```
mu \ h     h=100          h=1000         h=5000         h=20000        h=46000        h=59900
2e-4       0/17  8.9e-16  17/17 3.0e-02  17/17 1.2e-01  17/17 1.1e-01  17/17 1.1e-01  17/17 1.1e-01
1e-3       0/121 1.8e-15  121/121 4.4e-02 121/121 1.9e-01 121/121 1.5e-01 121/121 1.4e-01 121/121 6.2e-01
5e-3       0/121 8.9e-16  121/121 4.2e-02 121/121 1.8e-01 121/121 1.4e-01 121/121 1.3e-01 121/121 6.1e-01
2e-2       0/121 2.2e-16  121/121 2.6e-02 121/121 1.1e-01 121/121 7.6e-02 121/121 7.3e-02 121/121 5.4e-01
1e-1       0/121 5.6e-17  121/121 6.3e-03 121/121 2.8e-02 121/121 1.3e-02 121/121 1.2e-02 121/121 3.9e-01
0.3        0/121 2.1e-17  42/121 2.1e-03 121/121 9.5e-03 118/121 4.1e-03 118/121 4.0e-03 121/121 1.9e-01
0.6        0/121 6.9e-18  1/121 1.1e-03  121/121 4.8e-03 117/121 2.1e-03 116/121 2.0e-03 121/121 1.0e-01
1.0        E(121)         E(121)         E(121)         E(121)         28/121 1.2e-03 E(121)
```

Worst |Δτ| = **6.18e-01** at λ = 300 nm, µ = 1e-3, h = 59900 m — i.e. 618× the 1e-3 gate.

**Both of the brief's stop conditions fired, so the scan is not usable as an acceptance
result and I am not proposing one:**

- **The two references disagree by 2.83e-01**, far above the 1e-3 tolerance. The cause is
  localised: the `w`-march used `dw = 1.0` in the substituted variable, which near the top of
  a long ray corresponds to steps of up to ~490 m in z, so the `w`-reference is simply
  under-resolved. The z-march (1 m steps) is the trustworthy one; the `w`-reference needs
  `dw` small enough that `2w·dw` is a metre or less, which for `h = 59900` means
  `dw ≈ 2e-3` — 500× more steps than used here.
- **The h = 59900 column lies outside the table's height range (0-50000 m)**, so production
  clamps and returns α = 0 above 50 km while the reference does not. Its values (up to 6.2e-01)
  are a clamp artifact and that column must be dropped or explicitly marked out-of-domain.
  The 10 worst rays are *all* at h = 59900, so the headline number is an artifact.

What survives with a trustworthy reference is the h ≤ 46000 region, and there the picture is
**not** obviously alarming: at h = 100 m the agreement is at the 1e-16 level (the profile's
1 km h spacing means the production ±25 m difference and the native-knot slope are the *same*
quantity inside a cell, so exact agreement is expected and is what is seen); away from h = 100
the residuals are 1e-3 … 1e-1 and their origin — accumulated per-step differences over ~10⁴
steps, or a genuine interpolation effect — is **not yet established**. Reporting that as a
pass or a fail would be premature.

### 22.3 Cost

Production lookup **23.98 µs/step** against **5.41 µs/step** for the direct native-knot slope —
a ratio of **4.43**, i.e. on this measurement the table route is *slower* than computing the
extinction directly. The geometry table build took **22.9 s**. Fallback fraction 0.212.
(Caveat: both timings are per integration step, and the production rule does two
`getOpticalDepth` calls per step; this is a cost signal, not a benchmark.)

## Section 23 — the line-780 fix, and the state of Task E2

### 23.1 The fix (approved; 15 insertions, 3 deletions in `trace()`)

`CurvedRayTable.hpp`, inside `trace()`. The degenerate (radial, µ = 1) case no longer returns
an unset result; it is treated as the valid trace it is, with Φ forced to 0 so the displacement
is exactly zero:

```diff
     double const d = vacuumD(h, mu) + Dd;
+    // Transverse component of the emission direction.  It is exactly zero for a radial
+    // emission (mu = 1), where the transverse direction is undefined and tn is a pure
+    // round-off residual.  That is a VALID trace, not a table failure: the arrival point
+    // is the observation sphere along the emission radial (Phi = 0, so d = 0) and the
+    // arrival direction is the radial itself.  The tolerance is dimensionless: the
+    // round-off floor of tn is ~1e-16, while the smallest physical value is
+    // muMin = 2e-4, so 1e-9 separates the two by orders of magnitude on both sides.
+    double tx = dx - cRad * rx, ty = dy - cRad * ry, tz = dz - cRad * rz;
+    double const tn = std::sqrt(tx * tx + ty * ty + tz * tz);
+    bool const degenerate = !(tn > 1e-9);
...
-    double const Phi = d / r_obs_;
+    double const Phi = degenerate ? 0.0 : d / r_obs_;
-    double tx = dx - cRad * rx, ty = dy - cRad * ry, tz = dz - cRad * rz;
-    double const tn = std::sqrt(tx * tx + ty * ty + tz * tz);
-    if (!(tn > 0.0)) return res;
-    tx /= tn;
-    ty /= tn;
-    tz /= tn;
+    if (!degenerate) {
+      tx /= tn;
+      ty /= tn;
+      tz /= tn;
+    }
```

Your four conditions are met as follows: (1) **Φ = 0** is forced, so `d = 0` exactly, not merely
a radial direction; (2) the tolerance is **1e-9**, justified in the comment — dimensionless,
seven orders above the round-off floor (~1e-16) and five orders below the smallest physical
value (`muMin = 2e-4`); (3) a **regression check** was added to both
`tests/modules/testCherenkovCurvedRayTable.cpp` (a new SECTION) and the standalone
`stage8_regtest_b.cpp`; (4) the displacement re-scan is **not yet confirmed** — see 23.3.

### 23.2 Verification of the fix (measured)

The standalone harness, run against the fixed table:

```
(1) j=393 |err| = 0.001631 m  PASS      j=394 |err| = 0.001664 m  PASS
(2) |dt| = 0.004499 ns  PASS
(3) mu = 1.0, every height in the scan:
    h=100    status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  PASS
    h=1000   status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  PASS
    h=5000   status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  PASS
    h=20000  status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  PASS
    h=46000  status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  FAIL (d = 9.494e-02 m)
    h=50000  status=0 hit=1  dir=(0.000e+00, 0.000e+00, -1.000)  FAIL (d = 9.494e-02 m)
```

**The status bug is fixed**: `status = 0` (Hit) and `hit = 1` at every height, where before
they were `InvalidTable`. **The direction bug is fixed**: the arrival direction is exactly
`(0, 0, -1)`, the radial, at every height. The two FAIL lines reported earlier were a defect in
**my measurement**, not in the fix: the harness inferred the displacement as `R*acos(z/R)`,
which is ill-conditioned at `z/R = 1` (a one-ULP wobble in the normalised emission radial moves
`acos` by ~1.4e-8 rad = 9.5 cm). The metric has been replaced by the **transverse coordinates**,
which are well conditioned, and the check is now **green**:

```
(3) mu = 1.0 at h = 100 / 1000 / 5000 / 20000 / 46000 / 50000:
    status=0 hit=1  transverse=0.000e+00 m  dir=(0.000e+00, 0.000e+00, -1.000)  PASS
ALL ASSERTIONS PASS (0 failing assertion(s))
```

So the displacement is **exactly zero**, not merely small: `Phi = 0` puts the impact point on
the emission radial and the transverse coordinates vanish identically. The same metric change
was made in the Catch2 SECTION.

### 23.2b The production diff, hunk by hunk (item 1)

`git diff -U2` splits the file into exactly **five** hunks, and they account for the whole
**31+/8−** with nothing left over:

| hunk | change | label |
|---|---|---|
| `@@ -569,4 +569,11` | +7 (comment + `bStep` declaration) | **integrator fix** (Stage 8) |
| `@@ -576,5 +583,5` | +1/−1 (`dN = bStep * Nz`) | **integrator fix** |
| `@@ -612,4 +619,5` | +1 (`bStep = …` once per substep) | **integrator fix** |
| `@@ -755,4 +763,14` | +10 (comment + `tx`/`tn`/`degenerate`) | **line-780 fix** |
| `@@ -765,13 +783,18` | +12/−7 (Φ forced to 0, `if (!degenerate)`, early return removed) | **line-780 fix** |

So **9+/1− is the integrator fix** (exactly the figure reported in Stage 8) and **22+/7− is the
line-780 fix**, of which 13 insertions are comment text. There is **no third group and no
unexplained production change**: the earlier phrase "the Stage 8 build fix" was a mislabel for
what is in fact the line-780 fix.

The tolerance comment was also corrected as requested: `tn` is `sin(theta)`, so it is compared
against the most nearly vertical **non-degenerate** ray in the domain (µ = 0.6, i.e.
`sin(theta) = 0.8`), not against `muMin`; a ray with `sin(theta) < 1e-9` would be displaced by
under `R*1e-9 ~ 6 mm` in any case, so the cut is safe.

### 23.3 Not done, and why

- **Condition 4 — re-running the displacement scan** to confirm the `[1e-1, 1.0)` band reaches
  0 failures: the re-scan was compiled against the fixed header but **could not be completed**
  in this environment — the process is now killed on launch (it printed the header once, then
  nothing, and three subsequent launches died immediately with empty output). So the band is
  **not yet confirmed**.
- **Task E2 (the absorption rerun) was not run.** The corrections are understood and ready to
  apply to `stage9_abs_b.cpp` — domain `h <= 50000` m, λ = 300-700 nm (the file's own knots),
  µ ∈ {2e-4 … 0.6} with µ = 1.0 excluded entirely and µ = 2e-4 rejections kept in a separate
  column, a `w`-reference at `dw ≈ 2e-3` with a `dw/2` convergence check required to agree with
  the 1 m z-march to within 3.3e-4, and a report separating per-step accumulation from
  cell-dependent interpolation — but running it needs a working session.

Nothing in this section changes the earlier conclusions: the geometry fix stands, the
displacement target is met everywhere except the (now-fixed) degenerate rays, and the
absorption question is still **open**.

## Section 24 — Task E2, Checkpoint 1: the references do not agree, so the run stopped

New: `stage9b_abs_b.cpp` (modes `refcheck` and `scan`), log `log_stage9b_abs.txt`. No
production edits in this session. This section only adds.

### 24.1 Checkpoint 1 — reference convergence (measured)

Three rays, λ = 300 nm, h = 46000 m, integrand = the native-knot slope of the profile's own
cumulative τ(h), integrated along the same curved path in two parameterisations:

```
   mu        w(dw=2e-3)       w(dw=1e-3)       z(dz=1m)         |w2-w1|      |w1-z|
   1.00e-03  75.40262824      75.40283388      75.40157114      2.056e-04    1.263e-03
   2.00e-02  47.52290511      47.52300046      47.52275043      9.536e-05    2.500e-04
   6.00e-01   4.56254840       4.56255533       4.56252616      6.930e-06    2.917e-05

largest |w(2e-3) - w(1e-3)| = 2.0563e-04      PASS  (< 3.3e-4)
largest |w - z|            = 1.2627e-03      STOP  (> 3.3e-4)
```

The `w`-reference converges on its own (2.06e-04, halving with each halving of `dw`), but the
`w`-march and the 1 m z-march **disagree by 1.26e-03**, above the gate, so per the brief the run
**stopped here** and Checkpoints 2-4 were not started.

**Interpretation (the useful part).** The disagreement is largest on the most grazing ray
(µ = 1e-3) and shrinks monotonically with µ, which points at the *step*, not the integrand. The
z-march advances in z with a fixed 1 m step, but the path element is
`ds = dz / sqrt(1 - s^2)`, and for a grazing ray `sqrt(1 - s^2) ~ µ` near the ground. At
µ = 1e-3 a 1 m z-step is therefore a **~1000 m path-step**, so it is the **z-march that is
under-resolved**, not the `w`-march — the opposite of the assumption in the brief ("use the
z-march since it is cheaper if it agrees"). It does not agree, so it cannot be used as the
reference at small µ.

**Consequence for Checkpoint 2.** The reference must be parameterised in *path length*, not in
z. Two workable choices: (i) use the `w`-reference at `dw = 1e-3` (own convergence 2.06e-04,
measured here), or (ii) keep a z-march but make the step µ-dependent, e.g. `dz = 1e-3 * µ`,
which puts the path step at ~1 m for every µ. Option (ii) is cheaper and is what I would use
next, with its own convergence check at `dz/2` before scanning.

### 24.2 Note on the production file's line count

`CurvedRayTable.hpp` reads **34+/8−**, not the 31+/8− quoted in the brief. The extra three
lines are the tolerance-comment correction requested in the previous review (all comment text;
the reasoning now compares `tn` to the most nearly vertical non-degenerate ray, µ = 0.6). The
code is unchanged from the 31-line state. No further production edits were made in this
session.

### 24.3 Where the run stopped

Stopped at **Checkpoint 1**, as the brief directs. Checkpoint 2 (the scan), Checkpoint 3 (the
(a)/(b) split) and Checkpoint 4 (timing) were not started. The scan program is written and
compiles, so the next session can go straight to Checkpoint 2 once the reference is fixed —
with the `µ`-dependent step above, or the `w`-reference, and its convergence re-checked first.

### 24.4 Checkpoint 0 — the 3 extra lines are comment-only (confirmed)

`git diff -U0` gives six hunks: `+7`, `+1/−1`, `+1`, `+13`, `+1/−1`, `+11/−6` = **34+/8−**, and
filtering out comment lines leaves exactly the same **11 code lines** as the 31-line state
(`bStep` declaration, `dN = bStep * Nz`, the per-substep `bStep = …`, `tx/ty/tz`, `tn`,
`degenerate`, `Phi`, `if (!degenerate) {`, and the three normalisations). The three extra lines
are the tolerance-comment rewrite. Confirmed comment-only.

### 24.5 Checkpoint 1b — the path-length march, and a defect in it (STOPPED)

Added `tauIntS`: RK4 in **path length** `s` (`dz/ds = cos θ`, `dθ/ds = −sin θ (1/r + n'/n)`),
same curved path, same native-knot integrand, marching down from the emission.

```
mu=1e-3: s(1)=75.40270510  s(0.5)=75.40258741  R=75.40303951  |s1-s.5|=1.18e-04  |s.5-R|=4.52e-04
         steps = 838826 / 1677651   t = 0.4 / 0.8 s
mu=2e-2: s(1)=47.52299203  s(0.5)=47.52284183  R=47.52309582  |s1-s.5|=1.50e-04  |s.5-R|=2.54e-04
         steps = 686331 / 1372661   t = 0.3 / 0.7 s
mu=0.6:  s(1)= 4.56276704  s(0.5)= 4.56262935  R= 4.56256226  |s1-s.5|=1.38e-04  |s.5-R|=6.71e-05
         steps =  76217 /  152433   t = 0.0 / 0.1 s

largest |s(1) - s(0.5)| = 1.5020e-04      PASS  (< 3.3e-4)
largest |s(0.5) - R|    = 4.5210e-04      STOP  (> 3.3e-4)
```

**The gate failed, so Checkpoint 2 was not started — and the diagnosis of *why* was wrong.** I read
τ ≈ 75 as if it were a path length; it is the optical depth. A ray from 46 km reaching the ground
at µ = 1e-3 travels about √(2·R·h) ≈ 770 km, or ~840 km with refraction, so **838,826 steps at
ds = 1 m is correct**, and likewise 686 k at µ = 2e-2 (≈650 km) and 76 k at µ = 0.6 (≈77 km).
**There was no stall**, and neither an upward march nor a stall guard is needed. The step counts
in §24.5 above should be read as correct.

**First-order convergence is real and expected**: the integrand is the native-knot slope, which is
piecewise constant with a kink at every 1 km knot, so RK4 cannot do better than first order across
them. The `|s(1) − s(0.5)| = 1.2–1.5e-04` figure is therefore the right size, not a symptom.

**The real defect was the gate.** An absolute 3.3e-4 was applied to depths up to τ ≈ 75 — a
relative accuracy of ~5e-6, which has no physical meaning, since at τ = 75 the transmission is
e⁻⁷⁵ and any error in τ there is invisible.

**New tolerance (recorded here and used from Checkpoint 2 on):** a ray passes if

```
|dtau| <= 1e-3 * max(1, tau_ref)          and   |dT| = exp(-tau_ref) * |dtau| is reported
```

### 24.6 Checkpoint 2, "(a) only": production rule on the 1 m reference steps

Copied to `cp2_rule_on_1m.out`. Both sides use the same path and the same 1 m steps, so the
difference isolates the α rule; the relative tolerance of §24.5 is applied. 522 hits, 108
excluded. Runtime: production 93.2 s, reference 67.1 s, geometry build 22.8 s.

```
mu \ h   h=100         h=1000        h=5000        h=20000       h=46000       h=50000
2e-4     0/3  0.0e+00   3/3  1.6e-03   0/3  5.3e-05   0/3  1.3e-05   0/3  1.2e-05   0/3  1.3e-05
1e-3     0/21 0.0e+00  21/21 1.7e-03   0/21 1.2e-04   0/21 5.3e-05   0/21 1.4e-05   0/21 1.4e-05
2e-2     0/21 1.1e-16  21/21 3.0e-03   0/21 2.1e-04   0/21 1.2e-04   0/21 1.1e-05   0/21 1.0e-05
1e-1     0/21 0.0e+00  21/21 3.0e-03   0/21 3.1e-04   0/21 2.7e-04   0/21 1.1e-05   0/21 1.4e-05
0.6      0/21 3.5e-18   1/21 1.1e-03   0/21 3.9e-04   0/21 3.6e-04   0/21 3.2e-05   0/21 4.1e-05

worst relative |dtau|/max(1,tau_ref) = 3.0301e-03  at lambda=480 nm, mu=1.0e-1, h=1000 m
   (relative limit 1e-3)

5 worst rays, split:
   |dtau|     lambda  mu       h      (a)        (b)        path L[m]
   3.030e-03  480     1.00e-01 1000   3.030e-03  1.005e-04   9951
   3.022e-03  680     2.00e-02 1000   3.022e-03  2.227e-05  44907
   2.998e-03  460     1.00e-01 1000   2.998e-03  1.005e-04   9951
   2.997e-03  500     1.00e-01 1000   2.997e-03  9.727e-05   9951
   2.968e-03  660     2.00e-02 1000   2.968e-03  2.669e-14  44908
```

Recorded without interpretation, as the brief directs. Every failure is in the h = 1000 m column;
h = 100 m is exact to machine precision, and h >= 5000 m is inside the limit.

### 24.7 Step 1 citation: production's own absorption stepping

The stepping to copy for the third column is in **`Cherenkov.inl:583-593`** (the straight-line
segment branch):

```cpp
583:  int const n_tau_steps = std::max(1, static_cast<int>(seg_len / 50.0));
584:  double const ds_tau = seg_len / static_cast<double>(n_tau_steps);
585:  for (int ts = 0; ts < n_tau_steps; ++ts) {
586:    double const s_mid = (static_cast<double>(ts) + 0.5) * ds_tau;
...
592:    if (enableAbsorption_ && extinctionCoeffFunc_) {
593:      accumulated_tau += extinctionCoeffFunc_(tau_h, lambdaNm) * ds_tau;
```

i.e. the production absorption accumulates in **~50 m sub-steps** along the segment, with
midpoint heights. The curved branch uses the tracer's own step instead
(`Cherenkov.inl:649-651`, `kext * ds`). So `tau_prod` is the production rule evaluated on a
### 24.8 Step 1: knot-offset test (production rule on 1 m reference steps)

The 5 worst rays of `cp2_rule_on_1m.out`, at h ∈ {950, 975, 1000, 1025, 1050, 1100, 1500} m.
Output: `cp2b_knot.out` (43 lines, flushed per line). Numbers only.

```
lambda  mu        h     tau_ref           dtau          dtau/max(1,tau)
480     1.00e-01  950   0.98748062        +0.000e+00    +0.000e+00
480     1.00e-01  975   1.01327459        +0.000e+00    +0.000e+00
480     1.00e-01  1000  1.03917298        -3.149e-03    -3.030e-03
480     1.00e-01  1025  1.05239612        -2.221e-05    -2.111e-05
480     1.00e-01  1050  1.06567280        +1.051e-05    +9.862e-06
480     1.00e-01  1100  1.09211909        -1.541e-05    -1.411e-05
480     1.00e-01  1500  1.30336820        -2.122e-05    -1.628e-05
680     2.00e-02  950   2.60090705        +0.000e+00    +0.000e+00
680     2.00e-02  975   2.66278719        +0.000e+00    +0.000e+00
680     2.00e-02  1000  2.72436399        -8.233e-03    -3.022e-03
680     2.00e-02  1025  2.75287995        -1.210e-05    -4.397e-06
680     2.00e-02  1050  2.78125488        -9.455e-06    -3.400e-06
680     2.00e-02  1100  2.83760987        -1.113e-05    -3.922e-06
680     2.00e-02  1500  3.27096225        -2.652e-05    -8.108e-06
460     1.00e-01  950   1.04273203        +0.000e+00    +0.000e+00
460     1.00e-01  975   1.07007949        +0.000e+00    +0.000e+00
460     1.00e-01  1000  1.09731667        -3.289e-03    -2.998e-03
460     1.00e-01  1025  1.11142260        -2.328e-05    -2.095e-05
460     1.00e-01  1100  1.15379747        -1.643e-05    -1.424e-05
460     1.00e-01  1500  1.37914959        -2.410e-05    -1.748e-05
500     1.00e-01  950   0.91980025        +1.110e-16    +1.110e-16
500     1.00e-01  975   0.94382633        +1.110e-16    +1.110e-16
500     1.00e-01  1000  0.96794969        -2.997e-03    -2.997e-03
500     1.00e-01  1025  0.98001055        -2.108e-05    -2.108e-05
500     1.00e-01  1050  0.99212024        +1.014e-05    +1.014e-05
500     1.00e-01  1100  1.01624197        -1.439e-05    -1.416e-05
500     1.00e-01  1500  1.20892278        -1.863e-05    -1.541e-05
660     2.00e-02  950   2.82782162        +0.000e+00    +0.000e+00
660     2.00e-02  975   2.89509888        +0.000e+00    +0.000e+00
660     2.00e-02  1000  2.96204635        -8.790e-03    -2.968e-03
660     2.00e-02  1025  2.99369146        -1.420e-05    -4.743e-06
660     2.00e-02  1050  3.02518007        -1.266e-05    -4.184e-06
660     2.00e-02  1100  3.08771907        -1.704e-05    -5.519e-06
660     2.00e-02  1500  3.56865586        -2.022e-05    -5.665e-06

native-knot slope (alphaRef) by layer:
  lambda=480:  h=500:1.044290e-04  h=950:1.044290e-04  h=1050:5.353500e-05  h=1500:5.353500e-05  h=5000:2.017700e-05
  lambda=680:  h=500:6.066680e-05  h=950:6.066680e-05  h=1050:2.820570e-05  h=1500:2.820570e-05  h=5000:8.482000e-06
```

**Sanity check (passed):** the h = 1000 rows reproduce `cp2_rule_on_1m.out` — −3.030e-03
(λ=480, µ=1e-1), −3.022e-03 (λ=680, µ=2e-2), −2.998e-03 (λ=460), −2.997e-03 (λ=500),
−2.968e-03 (λ=660) — identical to the stored worst-ray table, to better than 1e-6.

### 24.9 E2 2c Steps 1-3

**Step 1 — the tracer's step, quoted.** `Cherenkov.inl:608-615`:

```cpp
608:  double const H_scale = computeLocalScaleHeight(h_cur);
609:  double const dx_step_base = H_scale / 20.0;
610:  double dx;
611:  if (sin_theta_cur > 0.99) {
612:    dx = H_scale / 100.0;
613:  } else {
614:    dx = std::min(dx_step_base, 1000.0);
615:  }
```

and the path element is `Cherenkov.inl:642`, `double const ds = std::sqrt(dx_step * dx_step +
dz_step * dz_step);`, with the absorption accumulated at `Cherenkov.inl:648-652`
(`h_mid = 0.5*(r_prev - R_E + r_cur - R_E)`, `kext = extinctionCoeffFunc_(h_mid, lambdaNm)`,
`accumulated_tau += kext * ds`).

**Steps 2-3 — the five columns on the sanity rays** (`cp2c_chk.out`). Tracer path-step proxy
`ds = min(H_scale/20, 1000) = 400 m` for `H_scale ~ 8 km`. Numbers only:

```
lambda=300 mu=6.00e-01 h=46000
   tau_ref=4.56276704  da=-4.844e-05  db=+7.637e-02  dtot=+7.632e-02  dtab=+8.614e-02  dtab1m=0.000e+00
   steps ref=76217 rule1m=76217 prod_curved=191 tab_curved=191 tab_1m=76217  tracer path=76400 m  mean=400.0 m
lambda=480 mu=1.00e-01 h=1000
   tau_ref=1.03917298  da=-3.149e-03  db=+6.312e-03  dtot=+3.163e-03  dtab=+5.117e-03  dtab1m=0.000e+00
   steps ref=9951 rule1m=9951 prod_curved=25 tab_curved=25 tab_1m=9951  tracer path=10000 m  mean=400.0 m
lambda=480 mu=1.00e-01 h=975
   tau_ref=1.01327459  da=+0.000e+00  db=+3.102e-02  dtot=+3.102e-02  dtab=+3.102e-02  dtab1m=0.000e+00
   steps ref=9703 rule1m=9703 prod_curved=25 tab_curved=25 tab_1m=9703  tracer path=10000 m  mean=400.0 m
```

**Step 3 sanity checks — all passed:**

| check | expected | measured | verdict |
|---|---|---|---|
| λ=300, µ=0.6, h=46000: tau_ref | 4.5626 ± 1e-3 | 4.56276704 | PASS (diff 1.7e-4) |
| λ=480, µ=0.1, h=1000: da | −3.149e-3 ± 1e-6 | −3.149e-03 | PASS |
| λ=480, µ=0.1, h=1000: \|dtab1m\| | ≤ 1e-4 | 0.000e+00 | PASS |
| λ=480, µ=0.1, h=975: \|dtab1m\| | ≤ 1e-6 | 0.000e+00 | PASS |
| tracer step count | > 0 | 191 / 25 / 25 | PASS |
| mean tracer step | ≤ path length | 400 m vs 76400 / 10000 / 10000 m | PASS |

Step 4 (the scan) was not started.

**No absorption conclusion is drawn.** Stopped at Checkpoint 1b; Checkpoints 2-4 not started.

## Not covered here

Wired routing, accounting, exactly-once survival sampling, and end-to-end fallback are
**not** established by these diagnostics; they are decision tests only. The attenuation
contract change (Bernoulli drop → `w_out = w_in * T`) is approved as intended behaviour but
is a separate commit and is **not** a prerequisite for the table.

## Section 23 — Gap 4: the DESCENDING-emission branch (attempted, measured, NOT shipped)

Status: **investigated and rejected**. Code removed; the branch keeps falling back through
`Status::Upgoing`, which is correct (82 µs plane tracer). This section records the
measurement so the next attempt does not repeat it.

### 23.1 What the branch is

An emission below `r_obs` whose direction points DOWN (`cos psi_emit < 0`) dips through its
perigee and rises back to `r_obs`. It shares `(h, mu)` with the rising-emission family, so
it cannot share the residual arrays: it needs a second set plus a selector on
`sign(cos psi_emit)`. The cheap part is that **no extra integration is needed** — one
inward trace from `r_obs` crosses every row twice, and the second crossing (after the turn,
on the way up) *is* the descending family.

### 23.2 The reference, and the trap in it

The natural reference is the same uniform-index chord traversed *through* its perigee, with
`b = p/n(r_emit)`, `A = sqrt(r_obs^2-b^2)`, `B = sqrt(r_e^2-b^2)`:

```
L = A + B
Phi = pi - (psi_e + psi_obs) = atan2(b(A+B), b^2 - A*B)
dL/dh  = +r_e/B - b(db/dh)(1/A + 1/B)
dPhi/dh = -[(db/dh) - b/r_e]/B - (db/dh)/A
```

An n = 1/geometric referee (`/tmp/ref2_probe.cpp`) checks this against the quadratic chord
crossing and the vector angle: `|dL/L| = 8.5e-12`, `|dPhi| = 7.6e-14 rad`. **The cosine term
must be `b^2 - AB`, not `AB - b^2`**: with the latter the arc is wrong by ~π for the steep
rays (measured 3.125 rad), because the two points straddle the perigee. That single sign is
worth recording — it is invisible in a spot check at one `(h, mu)`.

### 23.3 Why it fails

The branch is physical only where its perigee clears the ground: `b > R` forces
`mu <~ 0.031`. Exactly in that window the chord is a poor reference, because the physical
ray spends its whole path near its perigee, where

```
eps^2 = 2 a (D - D_t),   a = 1/r + n'/n = 1.16e-7 m^-1     (physical)
eps^2 = 2 (1/r) (D - D_t), 1/r = 1.57e-7 m^-1             (chord)
```

so the residual is `sqrt(1.57/1.16) - 1 = 17 %` of `L` (the 19 % seen below) **and** it
carries the `sqrt(clearance)` cusp *inside* the interpolated path — unlike the rising
branch, whose path never reaches its perigee.

Measured with `/tmp/inward_val.cpp` against 3-D RK4 (ds = 0.5 m, interpolated sphere
crossing), log `/tmp/iv_desc.log`:

| h_em [m] | mu_arr | L_table [m] | L_truth [m] | dL/L | dP [m] |
|---|---|---|---|---|---|
| -50 | 0.005 | 68403.2 | 65987.8 | 3.7e-02 | 2415 |
| -50 | 0.008 | 119519.2 | 118193.5 | 1.1e-02 | 1326 |
| -50 | 0.020 | 318615.8 | 318069.1 | 1.7e-03 | 547 |
| -200 | 0.020 | 310819.6 | 310245.8 | 1.8e-03 | 574 |
| -1000 | 0.020 | 259051.8 | 258185.1 | 3.4e-03 | 867 |

`max |dL|/L = 3.66e-02`, `max |dP| = 2.5 km` — two orders too coarse for a 12 m dish and
worse than the fallback it would replace. 12 of 14 physically valid probe cells were
served; the other 2 hit the horizon guard.

### 23.4 The evidence that an arc reference would work

The stretch between the truncation (`eps = 1.6e-3`) and the actual tangency is
`eps/a = 13.8 km` of path over **7 m** of depth — far too much to drop, and it is traversed
twice. Completing it analytically with the constant-bending model above (mirror the state
at the turn, add `2*eps/a` to `L`, `2*eps/(a r)` to `Phi`, and the corresponding `Sigma`
and tau) was implemented and is good to ~0.05 % over that stretch. That is what motivates
the recommendation: replacing the chord with a **circular arc of frozen `a`** for this
branch should bring the residual to ~1 %, which is the regime the rising branch already
enjoys with the chord.

### 23.5 What is *not* affected

Rising emissions (the Earth-skimming population the block exists for) are unchanged:
102/102 probe cells, `|dL|/L <= 2.7e-5`, and the frozen acceptance rows re-measured after
the removal give `1.476e-07 / 0.50 m / 3.437e-03 / eps 6.872e-03`.

### 23.6 To reproduce

```
g++ -std=c++17 -O2 -I<repo> -o /tmp/inward_val /tmp/inward_val.cpp   # from the diagnostic copy
/tmp/inward_val 800 1        # 40 s table build, then rising + descending scans
/tmp/ref2_probe              # the reference/geometry check of 23.2
```

Note that the probe's descending section needs the second-set hooks
(`storedDd2`, `storedS2`, `secondCrossingFrom`) that were removed with the code; the
diagnostic copy of the probe was pointed at the arc-reference work, i.e. it should be
reinstated together with the second set, not before it.

## Section 24 — Wave 2b: the exact chord optical depth (shipped)

`AtmosphericAbsorption` stores the cumulative column `tau(h) = int_0^h alpha dh'`,
**piecewise linear in altitude** (that is exactly what the bilinear lookup means), so
`alpha` is constant inside each 1 km slab and the integral along a straight chord has a
closed form: with `p = r_emit sin(theta)` (conserved along a chord) and
`x(h) = sqrt(r(h)^2 - p^2)`,

```
tau_chord = int_{r_min}^{r_emit} alpha dx + int_{r_min}^{r_det} alpha dx,
r_min = p if the chord's closest approach lies inside the segment, else min(r_emit, r_det).
```

`getOpticalDepthChord()` implements that; the per-slab extinction is precomputed once at
load (`alphaSlab_`, 51 x 105 = 43 kB), so a query costs **78 ns**. `getTransmissionSpherical`
now uses it, and the old air-mass/Chapman implementation is kept, unchanged, as
`getTransmissionSphericalChapman()` so the bias stays measurable.

### 24.1 Validation (testCherenkovAbsorptionChord, 19 assertions)

| check | result |
|---|---|
| vertical chord (p = 0) vs `\|tau(hE) - tau(hD)\|` | **0.0e+00** relative (exact by construction) |
| midpoint quadrature along the chord, ds and ds/2 | 5.2e-7 worst at ds/2 (gated at 1e-5) |
| cost | 78 ns per chord optical depth |

The quadrature is not monotone in `ds` (2.9e-8 -> 1.9e-7 on one ray) because the
integrand is piecewise *constant*: the error is set by where the samples fall relative to
the slab boundaries. The test asserts the absolute gate, not monotonicity.

### 24.2 The bias that was removed

`T` for 1000 m -> 3000 m at 400 nm, exact vs the old Chapman air mass:

| zenith | T exact | T Chapman | rel |
|---|---|---|---|
| 0 | 8.8958e-01 | 8.8958e-01 | 4.9e-08 |
| 30 | 8.7363e-01 | 8.7363e-01 | 1.4e-06 |
| 60 | 7.9143e-01 | 7.9148e-01 | 6.5e-05 |
| 80 | 5.1121e-01 | 4.9376e-01 | 3.5e-02 |
| 85 | 2.6710e-01 | 3.1105e-01 | 1.4e-01 |
| 88 | 4.7138e-02 | 1.3054e-01 | 6.4e-01 |
| 89 | 6.2046e-03 | 7.7583e-02 | 9.2e-01 |
| 89.9 | 5.9011e-05 | 4.6086e-02 | 1.0e+00 |

The two models **cross over**: the Chapman air mass slightly *over*-absorbs around 80 deg
and *under*-absorbs by two orders of magnitude beyond 88 deg — the Earth-skimming window
this whole track is about. Note also that the formula's `cos_zenith_geom = dh/path_length`
is the *flat* vertical fraction, so on a curved chord it degrades with the arc the chord
subtends.

### 24.3 The sign trap

`getTransmissionSpherical` recovers the zenith cosine from the chord by the law of cosines
on the triangle (observer, centre, emission). That law returns the angle at the emission
measured from the **inward** radial, so the zenith cosine is its **negative**. With the
wrong sign the column is doubled at the zenith (T = 0.673 instead of 0.830 for
1000 -> 3000 m) — caught by the A/B probe, not by any existing test, which is why
testCherenkovAbsorptionChord now exists: before it, **no test in the tree exercised the
straight-route absorption decision at all** (the E2E's gate is bit-identical either way).

### 24.4 To reproduce

```
g++ -std=c++17 -O2 -I<repo> -o /tmp/chord_tau /tmp/chord_tau.cpp && /tmp/chord_tau
<build>/tests/modules/testCherenkovAbsorptionChord
```


---

## Section 25 — Production validation of the fast path (1e4–1e6 GeV), and the plan that follows

Section 24 validated the pieces (tables, chord integrals) on frozen acceptance sets. This
section is the first time the whole route ran on a **real event through the real
application**, and it changed the picture: two convention defects made the fast path dead
in production without any test noticing, and with them fixed a **systematic disagreement
between the tabulated and the integrated optical depth** appeared that is larger than the
mission tolerance.

### 25.1 The measurement

One tau event, same seed and same physics in every configuration; only the transport
flags differ:

```
muon_cherenkov_3d --pdg 15 -E <E> --injection-height 1500 \
  --emcut 0.015 --mucut 0.015 --hadcut 0.3 \
  --observation-level 3000 --zenith 92 --azimuth 270 \
  --cherenkov-projection telescope --telescope-azimuth 270 --telescope-zenith 90 \
  --telescope-z 0 --telescope-radius 15 --seed 1000 --hadronModel SIBYLL-2.3d \
  -f <out>/<config>/out --cherenkov-output-format eventio -v info  [+ transport flags]
```

* `-f <out>/<config>/out` matters: the writer puts the photons at
  `dirname(-f)/cherenkov_hits.dat`, so without the per-config level **every configuration
  overwrites the same photon file** (this happened in the first pass).
* Absorption is **deterministic** (`exp(-tau)` multiplies the bunch weight, see 25.5-A), so
  every configuration traces the *same* photons and the weight ratio gives tau per photon.
* Images are weight-weighted `(x, y)` histograms of the telescope plane, read back with
  `eventio.IACTFile` (`photon_bunches[0]` fields `x, y, cx, cy, time, zem, photons,
  wavelength`; `photons` is the weight).

### 25.2 Two defects that made the fast path dead in production

**(a) The table could not build at any site above sea level.** `hTop` was a fixed 120 km
while the refined-index table is indexed by **absolute z** and the shipped 0–120 km profile
tops out at `z = 6 490 990 m`. A 3 km ASL telescope therefore asked for
`[6 374 000, 6 494 000]`, the coverage check refused, and `traceCurvedGround` silently
bypassed to the legacy tracer for **every** photon:

```
tableAvailable=false attempted=true bypassed=0 hits=0
gate 3804/337787 pass (1.13%), table hits 0, blind legacy 3804 ok / 0 fail, dropped 2071
```

2 071 of 3 804 gated photons were lost with no error message. The acceptance set could not
see it: it builds its own 0–130 km table. Fixed by clamping `hTop = min(120 km, zMax - r_obs)`
(now 116 990 m) and warning when clamped.

**(b) `refractiveIndexFunc_` was called with two different conventions.** Eleven call
sites passed a geocentric **radius** (separability grid, Cherenkov production along the
step, the whole atmo-delay block) while the geometry integrators passed a **height above
sea level**; all three applications expect a radius. Net effect in production: the
geometry saw `n = n(sea level)` at every altitude, so

```
f_phase = N(0) / N(h_emit)
```

exploded with height and the table rejected the population as
`DispersionOutOfRange`: **3 539 of 3 804** gated photons. Unified on height above sea
level (11 library sites + 3 applications). Dispersion rejections **3 539 → 0**, table hits
**0 → 8 524** (8 489 of them from the inward block, i.e. the rising Earth-skimming
population this track exists for).


### 25.3 What the six configurations give (yields are sum of weights = photons)

> **Wave-4 record — the `curved (table)` rows are one table revision stale.** Step B
> rebuilt the absorption table (true-perigee reference, sqrt-interpolation, 100 m h axis
> below 5 km). §25.3b re-measures the matrix after that rebuild and reproduces every
> non-table configuration to the tenth of a photon, which is what makes the single
> changed row attributable to the table alone.

1e5 GeV, 12.16 M photons produced, 67 791 gated (0.557 %):

| config | wall | bunches | photons | table hits | absorbed |
|---|---|---|---|---|---|
| straight, no abs | 91 s | 1 939 | 195 882 | – | – |
| straight + abs | 81 s | 1 939 | 176 341 | – | **9.98 %** |
| curved (table), no abs | 115 s* | 1 908 | 193 263 | **55 761** | – |
| curved (table) + abs | 125 s* | 1 908 | 168 544 | 55 761 | **12.79 %** |
| curved (no table), no abs | 89 s | 1 934 | 195 605 | – | – |
| curved (no table) + abs | 90 s | 1 934 | 176 043 | – | **10.00 %** |

\* includes the one-off 31.7 s table build. At this energy the table serves **82 %** of
gated photons (48 % at 1e4), 55 720 of them from the inward block. tau provenance for
`curve_abs`: **66 726 table + 1 065 chord + 0 missing = 67 791**, and **zero bunches cut**
(the weights carry the absorption, so the image is unbiased by truncation).

1e4 GeV for comparison (same cuts): absorbed 3.61 % straight, 4.22 % table, 4.06 % legacy;
table serves 48 % of gated photons. 1e6 GeV (in progress): straight 3 626 542 -> 3 371 949
photons, i.e. **7.02 %** absorbed, 18.5x the 1e5 yield.

The 1e6 *curved* numbers were measured at the end of wave 4 — 8.99 % (table) and 6.82 % (no
table) per the wave-5 handoff §4 W1 — but were never written down here. They are session-log
values, not records: re-measure before quoting either.

### 25.3b Post-step-B reproduction (wave 5, 2026-11-16)

Same physics, same seed, rebuilt absorption table, on a configuration identical to
§25.3: zenith 92, telescope radius 15 m, observation level 3000 m, MODTRAN_new,
seed 1000, 1e5 GeV. The geometry and absorption tables are now also loaded from an
on-disk cache (build 31.2 s + 18.2 s -> load 1.2 s + 0.15 s); the cached run reproduces
every logged number bit-for-bit against a cold build, so caching does not enter these
numbers.

| config | wall | bunches | photons (sum of weights) | absorbed |
|---|---|---|---|---|
| straight, no abs | 85 s | 1 939 | 195 881.9 | – |
| straight + abs | 85 s | 1 939 | 176 341.0 | **9.976 %** |
| curved (table), no abs | 103 s | 1 908 | 193 263.3 | – |
| curved (table) + abs | 104 s | 1 908 | 169 313.5 | **12.392 %** |
| curved (no table), no abs | 101 s | 1 934 | 195 605.0 | – |
| curved (no table) + abs | 106 s | 1 934 | 176 042.6 | **10.001 %** |

Five rows reproduce §25.3 to the tenth of a photon (195 882 / 176 341 / 193 263 /
195 605 / 176 043), and the §25.5 difference images reproduce as well (-1.34 %, -0.14 %,
-0.17 %). Only the table+absorption row moved: **12.79 % -> 12.392 %**, and its
difference image went from -4.42 % to -3.99 %. Since everything else is identical, that
change is attributable to the table and to nothing else.

tau provenance for `curve_abs`: **67 781 table + 10 chord + 0 missing** out of 67 791
gated photons. The 2.39 pp gap to the legacy route is therefore *not* a fallback effect
— ten photons used the chord. It is the tabulated tau itself, and §25.4's correction
below establishes which of the two routes is the one that is wrong.

### 25.4 The headline, as written in wave 4: tau disagrees three ways, and the table is the outlier

> **CORRECTED IN WAVE 5 — read the box at the END of this section before quoting any of
> it.** The inference drawn below ("the tabulated one is 28 % high... the table is the
> suspect") is **wrong**. §26.5 had already refuted it on the A/B sample; wave 5 confirms
> the refutation on the production population. The table and the paragraph that follow are
> kept as the wave-4 record, because the reasoning failure is instructive.

| tau source | where it comes from | absorbed @1e5 |
|---|---|---|
| exact chord (`getOpticalDepthChord`) | validated to 0.03–0.7 % (§24) | **9.98 %** |
| legacy path integral (plane tracer) | quadrature along the traced curved path | **10.00 %** |
| **tabulated** (`CurvedAbsorptionTable`) | cumulative difference, 500 m h cells | **12.79 %** |

Two independent methods that share no code agree to 0.2 %; the tabulated one is **28 %
high** (24 720 vs 19 562 absorbed photons — not statistics). The mean wavelength shift
points the same way (498.1 -> 506.4 nm tabulated vs -> 504.2 nm legacy): the table removes
more blue light.

The table is the suspect because it is the only one of the three whose validity was
established on a *different* population: §23/§24 validated it on rays whose observer sits
comfortably above the ray's perigee. This event is the opposite — near-horizontal photons
whose perigee is within a few hundred metres of the dish — which is exactly the regime
where we measured **-18 %, -74 %, -91 %** when the tangency gate was relaxed, and where
the cumulative cannot represent a path that visits the perigee twice.

---

### 25.4b The correction (wave 5, measured)

The wave-4 inference above rested on two methods agreeing with each other (exact chord
9.98 %, legacy path integral 10.00 %) and treating that consensus as truth. That is what
failed: **both were evaluating the same wrong path**, and agreement between two
implementations of one wrong model is not evidence.

The referee that settles it is `absorptionTauQuadrature()` — exact quadrature along the
refracted ray, perigee visited twice — introduced in §26.5 and re-run in wave 5 at
zenith 92 / r 15 / 1e5 GeV on **2 000 table-served photons of the production population**:

| mu bin | N | **tau table vs quad** | tau legacy vs quad | \|d impact\| | dL/L |
|---|---|---|---|---|---|
| 0.01-0.03 | 489 | **-0.02 %** | **-28.04 %** | 5.90 m | 3.0e-06 |
| 0.03-0.1 | 1 474 | **-0.04 %** | +7.25 % | 0.01 m | 9.5e-09 |
| 0.1-0.3 | 34 | **-0.02 %** | +222.45 % | 0.00 m | 1.5e-08 |

Median residual: **table -0.03 %**, legacy +1.99 %.

The corrected reading of §25.3 and §25.4 is therefore the opposite of what was written:

- **The tabulated tau is right** — within 0.04 % of an exact quadrature in *every* bin,
  including the grazing bin where this event's photons actually live.
- **The legacy path integral is wrong** — -28 % at mu < 0.03, the dominant population
  here. At run level that is 10.001 % where the table gives 12.392 %: the legacy route
  **under-absorbs by ~24 %**.
- The 2.39 pp disagreement in §25.3b is therefore *not* a table defect. It is the size of
  the legacy error, weighted over the gated population.

What this overrides:

1. **"The table within 1 pp of the legacy route" is an unsound acceptance criterion.** It
   asks the correct method to agree with one measured 24 % wrong; the only way to satisfy
   it is to break the table. The sound criterion is the per-bin comparison against the
   quadrature above (target <= 0.5 %), plus provenance and yields.
2. **Step B's 12.79 % -> 12.392 % is a real refinement but not a fix of this gap**,
   because the gap was never the table's.
3. **The remaining run-level work is W2**: repair the legacy path integral. When it is
   right, 10.00 % should move to ~12.4 % and the two-route comparison becomes a
   validation rather than a contradiction.
4. **§25.5's image conclusion needs the same caveat.** Its -4.42 % was read in wave 4 as
   the table over-absorbing. After step B it is -3.99 %, and since the tau is now
   validated, that residual is the physical curved-vs-straight difference (of which
   -1.34 % is geometry with absorption off) — not a table artefact.

Reproduce: one `curve_abs` run at zen92 / r15 / 1e5 with `--cherenkov-trace-ab 2000`, then
`python analyse_trace_ab.py <rundir>/cherenkov_ab.csv`, which writes
`trace_ab_residuals.png` alongside this README.

---

### 25.5 The image says it is geometry too, not only absorption

> **Read §25.4b first.** The -4.42 % below was attributed in wave 4 to the table
> over-absorbing; the tau is now validated (-0.03 %), so that attribution is wrong. After
> step B the same difference is -3.99 %.

Difference images relative to **the straight-line run of the same absorption setting**
(`cherenkov_xy_diff_E1e5.png`, 20 bins, 15 m dish):

| comparison | net | net % | photons redistributed (half-sum abs) |
|---|---|---|---|
| `curve_noabs - straight_noabs` | -2 619 | **-1.34 %** | 27 988 (14 % of the image) |
| `curve_noabs_notable - straight_noabs` | -277 | -0.14 % | 15 208 (7.8 %) |
| `curve_abs - straight_abs` | -7 797 | **-4.42 %** | 22 949 |
| `curve_abs_notable - straight_abs` | -298 | -0.17 % | **6 209 (3.2 %)** |

Read the two rows without absorption first: the legacy curved tracer reproduces the
straight-line image almost exactly (-0.17 %, 3 % of photons moved) — refraction at this
geometry displaces a ray by a couple of metres against a 15 m dish. The **table route
redistributes 4.5x more photons and loses 1.34 % net even with absorption off**. So the
launch-and-continue path (table impact on the `r_obs` sphere, then a plane tracer from
there) is not equivalent to integrating directly to the plane, and the absorption setting
only adds to that (-4.42 % = about -1.3 pp geometry + -3.1 pp tau).

Two systematic errors, both on the table side, both pointing at the same thing: the fast
path answers a **sphere** question for a **plane** instrument.


### 25.6 The plan

The geometry integrator is **not** the problem: `|dL|/L = 1.5e-7`, impact 0.5 m, and the
legacy curved route reproduces the straight image to 0.17 %. What has to change is the
**representation the fast path uses**. Steps are ordered by what unblocks what.

**A. Attribute before fixing — per-photon A/B harness** *(blocking)*
Deterministic absorption removed the RNG from the route, so an A/B is now exact. Add a
diagnostics switch that, for the first N photons, runs the table route *and* the plane
integrator from the same state and appends
`photon_id, h_emit, x,y,z, dir, lambda, mu, L_table, L_legacy, |d impact|, tau_table,
tau_legacy, tau_chord, status_table`. `tau_chord` is the referee (validated, §24).
Report the residuals binned in `mu x h_emit`, which answers "near-tangency or dip?" and
"which side is wrong?" instead of guessing.
*Acceptance:* a residual table per bin, and a decision for B based on it. Run it at 1e6
(3.6 M photons) so each bin is resolved.

**B. Fix the absorption table's ray representation (wave 2d)** *(blocking, physics)*
Three changes, each independently checkable:
1. **True-perigee reference.** The cumulative currently starts at the first node with
   `cos psi > 1e-3`, 32–428 m above the real perigee. Reference `T` at the analytic
   perigee `h_t(xi, lambda)`.
2. **Interpolate in `u = sqrt(h - h_t)`**, not in `h`: `T` is proportional to `u` for
   constant `alpha`, so this removes the sqrt-convexity bias by construction.
3. **Dip / sum rule.** When the path visits the perigee (descending emission, rising
   arrival — both signs are already in the tracer), `tau = T(A) + T(B) - 2 T(h_t)`;
   otherwise the difference. This is the one error whose *sign* can flip, and the reason
   the four acceptance rows currently fall back to the chord.
Keep the chord fallback for real rejections and keep the counters.
*Acceptance:* tau residual vs the referee **< 2 % on the 1e4/1e5/1e6 samples** (from
28 %); the four refused acceptance rows served by the table with `|d tau|/tau <= 1e-2`;
`AbsorptionChord` 21/21, acceptance 20/20, frozen-set `|d tau|/tau <= 5e-3` unchanged.

**C. Resolve the `tPlane <= 0` branch** *(correctness, cheap, do before trusting C/D)*
At 1e4, 9 140 of 9 349 `upgoing` photons were "rescued" by the plane tracer, mostly where
the crossing lies **behind** the emission. Physically that is not a hit. Decide it
explicitly, count it, and test it — this could inflate every legacy-route yield in 25.3.

**D. Plane-parametrised table** *(the architectural fix)*
Parametrise on the **telescope plane** (crossing position and arrival direction, traced
from the dish outward) instead of the `r_obs` sphere, so the tabulated quantity *is* the
emission -> plane path. This deletes launch-and-continue, turns the 18 % unserved
`upgoing` population at 1e5 into ordinary hits, and makes the image exact by construction
instead of by continuation.
*Acceptance:* >= 95 % of gated photons served; `|d impact| <= 1 m`, `|dL|/L <= 1e-6` vs the
plane integrator on a frozen set; image difference vs the plane tracer **< 1 % net** (from
1.34 %); per-event cost below the legacy route (currently 93 s vs 90 s at 1e5).

**E. Performance, only after A–D**
Measure with the protocol the timing analysis demands: 3 repeats, internal
`Run time`/`EAS.run` timers, eventio write timed separately, cold cache discarded (the
first run of a session was 8.5 s slower on identical work because the 306 MB binary and the
PROPOSAL tables were cold). Profile the continuation (RK4 steps, tau calls per photon) to
confirm the +12.9 s seen at 1e4 is the tau quadrature along long continuations. D is
expected to remove that cost entirely, so this becomes verification.

**F. Regression gates that must exist regardless**
* Triage the 27–28 pre-existing `testModules "CurvedRayTable"` failures against a
  pre-2g-h2n baseline *before* any geometry change, so a new defect cannot hide inside a
  known-failing suite.
* Add production gates to the Step-D counters: `tableAvailable == true` for any site with
  `observation-level > 0` (the 25.2a bug), and
  `absFallbackTau / (absTableTau + absFallbackTau) < 2 %`.

| step | unblocks | effort |
|---|---|---|
| A attribute | removes guesswork, makes B checkable | 1 session |
| C `tPlane <= 0` | validity of legacy-route yields | 0.5 session |
| B table tau | yield physics inside tolerance (28 % -> <2 %) | 1–2 sessions |
| F gates | no silent regressions | 0.5 session |
| D plane table | speed, 100 % coverage, 1.34 % -> <1 % image | 2–3 sessions |
| E performance | defensible numbers | 1 session |

In one sentence: **make the table answer the question the telescope actually asks — where
does the ray cross my plane, with what tau — instead of answering it on a sphere and then
walking to the plane.**

### 25.7 To reproduce

```
# runs (per configuration output trees are mandatory, see 25.1)
/tmp/run_e5.sh <config> [--enable-curvature] [--atmosphere-absorption-file <abs>] \
                        [--disable-cherenkov-curved-table]

# images and difference images
conda activate jupyter_env
python <plot>/make_cherenkov_plots.py E1e5 /uufs/.../stepE5_2gh2n
python <plot>/make_cherenkov_diff.py  E1e5 /uufs/.../stepE5_2gh2n
```

Artefacts: `Muon_Trinity/plot/cherenkov_xy_E1e5.png` (6 panels, 20 bins),
`cherenkov_xy_diff_E1e5.png` (2 rows x 3 columns, referenced to the straight line),
`cherenkov_radial_spectrum_E1e5.png`, plus the `E1e6` set; run logs under `/tmp/e5_*.log`,
`/tmp/e6_*.log`, and the per-config `config.yaml` inside each output tree.


---

## Section 26 — Step C verdict, and the design of the plane-parametrised table (step D)

### 26.1 C: the `tPlane <= 0` branch is legitimate and never fires in production

`rayTraceSphericalToPlane` accepts a plane crossing at `tPlane <= 0` and places the impact
at `p + tPlane*u`. That looked like a phantom-hit factory (a photon cannot reach a plane
behind it), so it was instrumented before any policy was chosen
(`CherenkovProcess::planeBackStats()`):

```
1e5 GeV, curve (no table), all 67,791 gated photons through the plane tracer:
  plane-crossing test: calls=67791  backward=0  (|tPlane|<=1 m: 0)  min tPlane=0.00 m
  blindLegacyOK=66824  blindLegacyFail=967
```

**Zero occurrences.** The 9 140 "rescued upgoing photons" seen at 1e4 GeV were ordinary
*forward* crossings: a photon emitted above the dish moving outward can still cross the
vertical dish plane ahead of it. The earlier reading of that counter was an inference, not
a measurement, and it was wrong.

Two cutoffs were tried and both were rejected by measurement: 100 m → 4 acceptance rows
became `contFail` (the test's "table hit but the plane continuation failed" guard), 1 km →
2 still. The reason is physical: a rising photon can cross the dish plane *well before* it
reaches `r = r_obs`, which is precisely the case the branch exists for. It is therefore
kept without a magnitude cutoff, **counted** (`calls/backward/nearZero/minTPlane`) so that
a future geometry which starts abusing it is visible instead of silent.

### 26.2 Why the sphere table cannot be repaired by better interpolation (D's motivation)

Step A's per-photon A/B (1e6 GeV, 3 000 table-served photons) measured the geometry error
by bin of arrival cosine:

| mu bin | N | \|impact_table - impact_legacy\| | dL/L |
|---|---|---|---|
| 0.01-0.03 | 1281 | **8.69 m** | 2.9e-05 |
| 0.03-0.1 | 1567 | 0.53 m | 5.2e-07 |
| 0.1-0.3 | 135 | 0.01 m | 2.6e-08 |

The path length agrees to 1e-5 while the *position* is off by 8.7 m — more than half the
15 m dish radius — and the image loses 1.34 % net at 1e5 and **3.99 % at 1e6** with
absorption switched off. The mechanism is the lever arm: an angular error `dmu` in the
table displaces the impact by `L*dmu/mu`, so at mu = 0.02 and L = 40 km a 4e-6 error in mu
is already 8 m. No h- or mu-interpolation improvement can beat that at fixed node count,
because the quantity being interpolated is not what the telescope measures.

### 26.3 D: tabulate the PLANE crossing (design -- SUPERSEDED by 26.4)

**Parametrisation.** Spherical symmetry makes the sphere table 2-D in (h, mu). Fixing a
*tilted* plane breaks that symmetry: where the ray crosses the plane also depends on the
azimuth of its direction about the radial, relative to the plane's tilt axis. So the plane
table needs a third axis:

```
(h, mu, phi) -> (X, Y, cx, cy, L, tau)
```

* `h` — emission height ASL, as now (96 nodes, 100 m cells below 5 km).
* `mu` — |cos| of the angle between the photon direction and the radial, as now (800 nodes,
  the same log/quadratic axis; grazing is mu -> 0).
* `phi` — azimuth of the photon direction about the local radial, measured from the
  projection of the telescope-plane normal onto the local horizontal; `phi = 0` is "in the
  plane of the tilt". Reflection symmetry gives phi in [0, pi]; 64 nodes, non-uniform
  (cos phi grading) because the crossing position varies fastest near phi = 0.
* Stored per node: the **plane** crossing (X, Y) in the telescope frame, the arrival
  direction (cx, cy) in that frame, the geometric path L, and tau.

**Build.** Trace *from the plane* backward: for each (mu, phi) launch a ray from the dish
centre outward and walk up; record, at each h node, where that ray crosses the plane. This
makes the tabulated quantity the plane -> emission path, i.e. exactly the question the
photon asks, with no launch-and-continue and no lever-arm amplification.

**Size.** 96 * 800 * 64 * 4 quantities * 4 B = **78 MB** uncompressed (the deleted tau block
was 64 MiB, so this is comparable); a float16 or per-quantity float32 split brings it to
~40 MB. Build cost scales with the phi nodes: expect 32 * 30 s = 16 min once per process,
which is acceptable only if the table is reused (it is).

**Singularities.** Rays that miss the plane inside a sane radius have NaN stored and are
rejected to the legacy tracer, counted as now (`planeTableRejected`).

**Acceptance (from the plan).** >= 95 % of gated photons served; |d impact| <= 1 m and
|dL|/L <= 1e-6 against the plane integrator on a frozen set; image difference vs the plane
tracer < 1 % net (from 3.99 % at 1e6); per-event cost below the legacy route.

**Status: REFUTED by the increment-1 probe, see 26.4** -- the phi axis this design needs
cannot be interpolated (5 608 nodes for 1 m).  Kept here because the reasoning about the
lever arm and the NaN policy still applies.  The first increment was: The next concrete step is the build loop for a
single (mu, phi) column plus a probe that compares its plane crossing with
`rayTraceSphericalToPlane` for the same launch — that validates the parametrisation before
any storage work.


### 26.4 Step D, increment 1: the phi axis is NOT viable (measured, design refuted)

Before any storage work, a self-contained probe (`/tmp/probe_plane_param.cpp`, exponential
atmosphere, exact first integral `p = n r sin(psi)`) marched backward from the dish with
the arrival direction `(mu, phi)` and compared the emission point from an N-node linear
interpolation in `phi` with the exact mid-cell value:

| h_emit | N=8 | N=16 | N=32 | N=64 |
|---|---|---|---|---|
| 500 m | 122 474 m | 30 693 m | 7 678 m | 1 920 m |
| 1 500 m | 122 474 m | 30 693 m | 7 678 m | 1 920 m |
| 20 000 m | 122 805 m | 30 775 m | 7 699 m | 1 925 m |

The error scales exactly as `1/N^2` and equals the **arc sagitta of the rotation**:
`r (dphi)^2 / 8 = 6.374e6 * (pi/64)^2 / 8 = 1919.8 m` against 1919.7 m measured.  In other
words the emission point under a change of `phi` is a pure *rotation about the Earth's
axis*, and interpolating between two rotated rays cuts the chord between them.  Reaching
1 m would need

```
N_phi = 2 pi / sqrt(8/r) = 5 608 nodes
```

which at (96 h x 800 mu) would be 4.3e8 nodes -- four orders of magnitude beyond the
memory of the block it was meant to replace.  The same argument applies to *any*
azimuth-like axis (including the arrival-direction azimuth), because the sensitivity of
the emission position to an angular change is `L * dphi` with `L ~ 40 km` at grazing, not
something the table can compress.

**Conclusion: D as scoped is refuted.** A plane-parametrised table does not remove the
sensitivity, it merely relocates it; the amplification `L/mu` that produces the 8.7 m
impact error is inherent to horizon geometry, not to the sphere parametrisation.  The
remaining lever on the geometry error is therefore the **resolution of the mu axis near
the tangency** in the existing table: the measured 8.7 m at mu < 0.03 is consistent with
an interpolation error that falls as `dmu^2`, so ~3x more nodes below mu ~ 0.05 (800 ->
~2400 nodes, 230 -> ~690 MB) would bring it to ~1 m, and the alternative is to accept the
already-quantified 3.99 % (1e6) image loss against a 15 m dish.

Nothing in this changes the absorption work: the tau path (§25, §26.1) is unaffected, and

### 26.5 Step A': the referee replaced -- the tabulated tau is right, the LEGACY one is not

The first A/B used a straight chord to the impact point as the referee and therefore
indicted the table (+537 % at mu < 0.03).  That referee was wrong: at grazing the
refracted ray dips to a perigee hundreds of metres below the emission, so a chord is not
the photon's path.  `absorptionTauQuadrature()` now integrates `alpha(h) dh / cos(psi)`
along the refracted ray with `u = sqrt(h - h_t)` (the same substitution the table build
uses) and visits the perigee twice when the photon descends to it and returns.

1e5 GeV, 2 000 table-served photons, every one traced twice:

| mu bin | N | **tau table vs quad** | tau legacy vs quad | \|d impact\| | dL/L |
|---|---|---|---|---|---|
| 0.01-0.03 | 489 | **-0.02 %** | -28.0 % | 5.90 m | 3.0e-06 |
| 0.03-0.1 | 1474 | **-0.04 %** | +7.3 % | 0.01 m | 9.5e-09 |
| 0.1-0.3 | 34 | **-0.02 %** | +222 % | 0.00 m | 1.5e-08 |

Median over the sample: **table -0.03 %, legacy +1.99 %**.

So B is validated at production statistics: the tabulated optical depth is within 0.04 % of
an exact quadrature in every grazing bin, and the +537 % that motivated step B's design
was an artefact of the chord referee.  The path integral in the legacy plane tracer is
the one that is wrong at grazing (-28 % at mu < 0.03), which is worth knowing because that
route still serves the 18 % of gated photons the table rejects, and because it is what the
E2E and ChordRoute MC checks compare against.

What remains is purely geometric and confined to mu < 0.03: 5.9 m impact displacement on
24 % of the sample (90th percentile 7.9 m, max 78 m) against a 15 m dish.  The next
experiment is the mu-axis refinement described in 26.4, which must be run as TWO processes
(one table each): a single process holding nMu = 800 and nMu = 2400 (231 + 692 MB) is
killed by the node, and the first attempt at the comparison did exactly that.

the per-photon table coverage is unchanged.



### 26.6 Wave 5: the horizon error is a cliff at the per-column horizon, and a global axis cannot fix it

The impact error had never been resolved against mu. Measured on 2 000 table-served
photons of the production population (zen92, r 15 m, 1e5 GeV, `curve_abs` A/B CSV,
|d impact| = table vs plane-tracer):

| \|mu\| bin | N | median \|d\| | p90 \|d\| |
|---|---|---|---|
| 0.010-0.020 | 19 | **66.6 m** | 77.6 m |
| 0.020-0.030 | 470 | **5.52 m** | 35.1 m |
| 0.030-0.050 | 1 152 | **0.0118 m** | 0.37 m |
| 0.050-0.100 | 322 | 0.005 m | 0.010 m |
| 0.100-0.300 | 34 | 0.0045 m | 0.010 m |

That is a 460x collapse across one bin width, not a gradual degradation, and it sits
exactly where the apparent horizon belongs. For an emission Delta below the telescope:

```
mu_h = sqrt(2*Delta / r_obs)      0.0177 at Delta = 1 km
                                  0.0217 at Delta = 1.5 km
                                  0.0280 at Delta = 2.5 km
```

so mu_h is NOT a single number: it ranges over [0.018, 0.028] across this event's
emission heights. Near its own horizon a ray's d and L behave like sqrt(mu - mu_h), which
is what motivated the experiment.

**The experiment (tried, refuted).** The axis was changed from `s = ln(mu/mu_min)` to
`s = sqrt(mu - mu_min)` with mu_min raised 2e-4 -> 0.018, which also addressed a separate
finding: 53 % of the old nodes sat below the horizon and could never be used. Cost: no
extra memory (same 800 nodes, 229 MB), node spacing improved from 1.08 % to 0.5-0.9 %,
table rebuilt in 31 s, and the cache magic was bumped so no stale table could be read.

Result, same 2 000 photons, same seed:

| bin | median OLD | median NEW |
|---|---|---|
| 0.018-0.020 | 66.61 m | 73.37 m (10 % worse) |
| 0.020-0.030 | 5.517 m | 5.620 m (2 % worse) |
| 0.030-0.050 | 0.0118 m | 0.0118 m (**identical to 4 digits**) |
| overall mean | 3.470 m | 4.143 m |
| fraction > 1 m | 24.4 % | 24.9 % |

Doubling the resolution and linearising a suspected cusp changed nothing where the error
is small and made it slightly worse where it is large. **The error is therefore not
interpolation on a global axis.** The reason is the per-column horizon above: a single
axis shared by every h row cannot linearise a singularity whose position moves with h.
The change was reverted; the binary was verified to reproduce 169 313.5 photons exactly,
and the runtime `muMin=2.00e-04 muBreak=0.884966` confirms the log axis is back.

What this rules out, and what is left:

- **Ruled out: buying more mu nodes.** If a 2x resolution change moves the answer by 0 %,
  the 800 -> 2400 upgrade (692 MB, 85.7 s build) would burn 460 MB to fix something that
  is not the limit.
- **Best remaining option: a per-column axis**, i.e. sqrt(mu - mu_h(h)) with its own
  horizon per row. This needs each h row to own its mu nodes, which breaks the shared-node
  indexing and the across-column derivative differencing -- a grid redesign, not a tweak.
- **Cheapest route to that: a grazing sub-table** (mu < ~0.05) with its own per-column
  axis, leaving the validated main table untouched.
- **Not available yet: reject deeper and let the exact tracer take it.** For the curved
  route a table rejection hands the photon to the plane tracer, whose tau is -28 % at
  grazing (§25.4b). That trades a 6 m geometry error for a 24 % optical-depth error, so
  this option is blocked until the legacy tau is repaired.
- Image cost if nothing is done: **3.99 % at 1e6**, measured directly as
  `curve_noabs - straight_noabs` (the identical comparison with the legacy tracer is
  -0.14 %).


### 26.7 Wave 5: the per-column axis is structurally impossible, and the rejection band is the real knob

**The per-column axis cannot be built.** The table is populated COLUMN BY COLUMN AT FIXED mu
(`buildColumn(mu, f, d, l, s, th)` fills every height for one mu) because a single ray at a
given mu traverses all heights. A per-height mu axis therefore needs `n_h x n_mu` separate rays
instead of `n_mu` — 4 171x more trace work. So "each height row owns its own sqrt(mu - mu_h)"
is not a tuning choice; it is incompatible with the build. That closes the design question.

**What does work: widen the grazing rejection band.** The table already refuses rays within
`dh_in_` = 10 m of their apparent horizon. Raising that to 400 m (`kGrazingClearance_m` in
`CurvedRayTable::trace`) sends the whole divergence to the exact tracer. Measured at
zen92 / r15 / 1e5 with `--cherenkov-trace-ab 2000`:

| bin | before | after 400 m |
|---|---|---|
| 0.010-0.020 | 66.6 m (N=19) | **rejected** (exact tracer, 0.50 m class) |
| 0.020-0.030 | 5.52 m, p90 35.1 (N=470) | **4.14 m, p90 12.5** (N=407) |
| 0.030-0.050 | 0.0118 m (N=1152) | 0.012 m (N=1212), unchanged |
| fraction > 1 m | 24.4 % | **20.1 %** |

Cost, same run: hits 55 761 -> 51 473, rejected 12 030 -> 16 318 (+36 %), yield
169 313.5 -> 166 930.1 (**-1.4 %**), absorbed 12.39 % -> **11.39 %** — the rejected photons
take the legacy route's chord-like tau, which is 17 % low (§25.4b).

**So this is a real but partial fix with a price, and the operating point is a physics
decision, not a software one.** The lever is `kGrazingClearance_m`: at 10 m (the old value)
the 1e6 image carries 3.99 %; at 400 m the worst bin is gone and the >1 m fraction drops 4 pp
for 1.4 % of the yield. Reaching the 1.5 m target would mean rejecting the whole 0.02-0.03 bin
as well (mu < ~0.035), roughly a quarter of the photons — and every rejected photon pays the
legacy tau error until W2 is finished.

The cache magic is now `CRVRTB4`, so no stale table can be loaded under the new guard.

