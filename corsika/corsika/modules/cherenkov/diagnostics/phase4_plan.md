# Phase 4 plan -- telescopic 1 mm at tangency: dual-mode propagation + impact-error control

Status: **plan**. Written against the working tree with Phases 1-3 shipped.

## 1. Goal and meaning of "1 mm"

Certify the **telescope impact coordinate** to 1 mm over the 85-95 deg field, including
tangency. An adaptive estimate targets 1 mm; a *bound* needs validated bounds. Deliver
the target plus a per-ray `certified / not_certified` status; never report the target as
a bound.

In scope: the plane tracer (`rayTracePlaneDopri`, Cherenkov.inl:896; `rayTraceSphericalToPlane`,
:441). Out of scope: `rayTraceSpherical` (Ground/RK4), `CurvedRayTable`, the absorption
tables (they are the REFERENCES here).

## 2. Baseline (measured)

| fact | value | source |
|---|---|---|
| frozen-set impact error, DOPRI5 @1 mm / 20 km cap | 0.118 mm max | phase1_geometry.md 3.3 |
| RK4 production caps, same set | 0.00034 mm | 3.3 |
| production photons above 1 mm (of 2075) | 2 | 3.5 |
| outlier 1 | RK4 wrong: DOPRI@1e-3 == DOPRI@1e-5, RK4 2.07 mm away | 3.5 |
| outlier 2 | impact moves 10 mm between atol 1e-3 and 1e-5 | 3.5 |
| table impact cliff near the horizon | 1.05 m @ mu .029-.030, .244 m @ .031-.033, .011 m above .033 | CurvedRayTable.hpp:1137 |
| rejections (10 m model nodes) | 39 % of attempts | 3.6 |

Reuse: `CurvedRayTable::buildColumn` (same ODE in the tangential variable, mm-validated,
= independent reference); the perigee bisection on A(h)=p and the u=sqrt(h-h_t) grading
(`absorptionTauQuadrature`); `RefractiveIndexTable`'s N, logSlope, cached cell slopes.

## 3. Architecture

State y (plane frame) + impact budget + turning data -> mode selector (hysteresis) ->
{regular: (x,z,u_x,u_z,Sigma) DOPRI5 | tangent: (Phi, s, Sigma, v), r = r_t + v^2} ->
shared dense-output plane crossing -> impact-error accounting and per-ray status.

Both modes call the same `kappaAt`/`N`; switching changes the independent variable only.

## 4. Regular mode: keep the ODE, replace the acceptance rule

Acceptance becomes impact-weighted: e_k = |J(s_{k+1}) . dy_k| <= b_k.

KEY POINT the sketch understates: b is a functional of an EVENT (the plane crossing), so

    J = [db/dy]_{s fixed} + [db/ds] * (ds*/dy),   ds*/dy = -(grad d . dy)/(grad d . f)

dropping the second term underestimates by 1/cos(theta_arr) -- the very factor that makes
the cliff. Implementation: two-pass per-ray adjoint (store per-step stages, cap 20k
steps), fallback to finite differences, plus a cheaper forward-sensitivity variant to be
measured against it.

## 5. Tangent mode

p = n r sin(psi), q(r) = n(r) r, r_t from q(r_t) = p, r = r_t + v^2,
B(v) = (q(r_t+v^2)^2 - p^2)/v^2, B(0) = 2 p q'(r_t):

    dPhi/dv = 2p/((r_t+v^2) sqrt(B)),  ds/dv = 2q/sqrt(B),  dSigma/dv = N ds/dv,  dv/dv = 1

v is carried as a state component so the RHS stays autonomous; the dense output gives
Phi, s, Sigma and v at the same t.

Stability: evaluate B cell-exactly (inside a cell n = 1 + N_i exp(b (r-r_i)), so
q' = 1 + n' r, q'' = b N (2 + r b)) and use the series 2p q' + (q'^2 + p q'') v^2 for
v^2 < 1 m. Branch: sign(v) = sign(dr/ds); cos(psi) = sign(v) sqrt(1 - sin^2(psi)).
Refuse (status) if q' <= 0 (ducting/degenerate) or the r_t bracket fails.
If the detector is reached before the turn, stop there.

## 6. Switching

REGULAR->TANGENT when p/p_crit > 1 - eps_in AND (r - r_t) < 1.5 * distance-to-plane;
TANGENT->REGULAR when v > 0 AND |mu_r| > eps_out, eps_out > eps_in (hysteresis).
Start eps_in from the existing grazingClearanceM_ calibration (164 m of clearance per
0.001 of mu); tune by validation, not assumption. Count switches; flag oscillation.

## 7. Budget, certification, counters

Allocation (mm): integration 0.50, turning-radius + coordinate conversion 0.15, crossing
solve 0.10, atmosphere evaluation 0.10, reserve 0.15. Distribute the 0.50 mm as
b_k = 0.5 mm * h_k / L_remaining; h_new = h clip[f_safety (b_k/e_k)^(1/5), f_min, f_max]
with f_safety re-tuned (0.9 already gives 39 % rejections here).
certified iff sum e_k <= 0.5 mm AND r_t converged AND no oscillation AND crossing
converged; otherwise `not_certified`, still recorded, counted and printed.
New counters: switches, r_t iterations, B series/branch use, adjoint passes,
sum e_k (mean/p99/max), certified/not_certified, aperture-edge verdict mismatches.

## 8. Validation

Extend `testCherenkovTracerBenchmark` with an 85-95 deg sweep at 3 km (turning points
below/at/above the telescope, 1-100 km paths, threshold-straddling rays, aperture edge
at R_tel +- 5 mm, layer crossings, misses). References: CurvedRayTable columns (an
independent implementation of the same invariant) and DOPRI5 at 1e-5/1e-6 m; report the
references' own disagreement so the gate is not reference-limited.
Gates: impact <= 1 mm for every CERTIFIED ray; no not_certified ray among today's hits;
direction <= 1e-7; time <= 0.05 ns; hit set identical; no aperture verdict flips.
Adopt if cost <= 2x today's DOPRI5 (>= 2x RK4 in config C) and no frozen-benchmark
regression.

## 9. Order, size, files

1. Tangent propagator behind `TracerIntegrator::TANGENT` (not default): cell-exact q/B,
r_t solve, regularised loop, reconstruction, dense output. ~450 lines
(`TangentPropagator.hpp` + Cherenkov.inl).
2. Switching + shared crossing logic. ~200 lines.
3. Impact-weighted control + status. ~350 lines (touches Dopri5.hpp, Cherenkov.hpp).
4. Validation/certification study, tuning. ~300 lines of test + 2 runs per config.

## 10. Risks, non-goals, stop condition

Cost creep (two passes) -> forward-sensitivity variant and the >=2x RK4 gate.
Degenerate/ducting turning points are DETECTED AND REFUSED, not supported.
Not a bound unless 8 passes across the whole field.
Non-goals: Ground-mode RK4 path, absorption tau, table contents.
Stop: if the tangent propagator cannot cut the certified step count by >=2x on the
near-tangent set, drop it and keep only the impact controller (which is what fixes the
10 mm ray).

## 11. Deliberate deviations from the sketch

(i) two-pass stored-step adjoint instead of retrace-to-convergence (the trajectory is
deterministic per photon, so the stored adjoint is cheaper and exact to the same order);
(ii) no explicit layer-boundary splitting -- 10 m model nodes make it a 10 m cap (3.6
measured the 39 % cost); the cell-exact B(v) replaces the need in tangent mode.

---

# Step 1 -- status: SHIPPED (diagnostic only)

`TangentPropagator.hpp` + `CherenkovProcess::rayTracePlaneTangent` +
`TracerIntegrator::TANGENT` + `testCherenkovTangentTracer.cpp`.

## What it is

State `(Phi, s, Sigma, v)` with `r = r_t + v^2`; `r_t` from a bracketed bisection on the
monotone `q(r) = n(r) r`; `B(v)` evaluated cell-exactly (series below `v^2 = 1 m`, direct
form above); the same DOPRI5 core, the same 5-point dense-output crossing scan, the same
tau rule and the same plane frame as the regular tracer.

## Numerics validated standalone

```
q(r)  exact to 0.00e+00 rel,  q'(r) to 8.0e-14          (vs the analytic exponential model)
B(0)  = 10607442.215253  vs  2 p q'(r_t) = 10607442.215253
series vs direct B:  1.3e-7 at v=0.1 m, 9.2e-10 at v=0.999 m, 0 beyond
                     (the DIRECT form carries 9.2e-4 of cancellation noise at v=0.001 m --
                      the reason the series exists)
```

## Domain, measured

A ray of invariant `p` turns at `q(r_t) = p`.  Sweeping emission height, distance and aim
shows the mode's population: **perigee between the surface and the dish, reached before
the telescope plane**.  Examples (3.1 km emission, 100 km, aim 0 / -0.1 / -0.2 / -1 deg):
perigee 2796 / 2904 / 2989 / 2812 m, turn at 62 / 50 / 38 / 61 km, all crossing the plane
at 100 km on the way back up.  Rays whose perigee is below the surface (a 45 deg downward
ray, an upward ray) are refused -- that is the contract, not a failure.

## Test result (`testCherenkovTangentTracer`)

```
reference: regular DOPRI5 at atol_pos = 1e-6 m, no step cap
6/6 turning rays served, 3/3 control rays correctly refused (0 silences)
max |impact| 0.009979 mm   (limit 1 mm)
max |dU|     1.77e-10      (limit 1e-7)
max dL       0.000790 mm
max dt       2.1e-5 ns     (limit 0.05 ns)
steps        23 accepted / 0 rejected for the 6 rays  (~4 per ray), h in [1, 21] m of v
```

`h` is in `v`, and at the perigee `ds/dv = 2q/sqrt(B) ~ 3.9 km per unit v`, so a 10 m step
in `v` spans ~39 km of path.  That is the whole point of the coordinate change: the rays
that forced the regular tracer onto small steps are now four steps each, at 10 um.

## A bug the test caught (recorded deliberately)

The first version used `groundZ_m_` as the floor.  That member is the **sphere target of
Ground mode** -- 3 km in the telescope runs -- so the check "hit the ground" fired for
every ray that dipped below the dish, i.e. the entire population the mode exists for, and
8 of 9 rays were refused.  The floor is now `earthRadiusM_ - 100 m`, matching the regular
tracers' `h_cur < -100 m`.  This is why the acceptance test runs against a *converged
regular* tracer rather than against the tangent mode's own estimates.

## Usage rule until step 2

`--cherenkov-tracer tangent` is **diagnostic**.  With no mode selector, every ray without
a turning point is refused and therefore *dropped*; the application logs a warning when it
is selected.  Production stays on `dopri5`.

## What step 2 must add

1. The selector (hysteresis, `eps_in` from the `grazingClearanceM_` calibration) with a
   **fallback to the regular tracer on refusal** -- this is what makes the mode usable.
2. `CurvedRayTable`-column cross-checks for the tangent population (the independent
   invariant-based reference), and the aperture-edge and layer-boundary rays.
3. The impact-weighted controller (step 3) and the per-ray certification status.

---

# Step 2 -- status: SHIPPED (`auto` is the application default)

## What was added

* `TracerIntegrator::AUTO`: try the regularised mode for every ray, and on ANY refusal
  fall back to DOPRI5 **for the same ray**, counted (`tangentFallbacks`).  A ray is never
  lost, which is what makes the mode usable without a mid-ray switch.
* `--cherenkov-tracer auto` (new default) | `dopri5` | `rk4` | `tangent` (diagnostic).
* `tangentRays` now means SERVED, `tangentUncertified` means refused; a ray that enters
  the tangent loop and fails there is counted once, not twice.
* A closed-form necessary condition before the bracket solve:
  `p > q(floor)` (one `q` evaluation).  Without it AUTO pays a 90-iteration bisection for
  every non-turning ray -- measured at ~9 us per refusal, i.e. ~1.8 s per event for the
  43 km muon.

**Deviation from the plan, stated:** no mid-ray switching and therefore no hysteresis.
The regularised coordinate covers the WHOLE path (v runs from negative, through the
perigee, to positive), so a ray never needs two integrations; the safeguard hysteresis
would provide is the counted fallback instead.  Velocity: every switch the plan worried
about is a per-ray decision here, and it is recorded.

## Two bugs the new test caught

1. **AUTO fell back to RK4, not DOPRI5.** The DOPRI dispatch was keyed on
   `tracerIntegrator_ == DOPRI5`, so an AUTO ray that refused tangency fell through to the
   RK4 branch below.  Caught by asserting `tangentRays + dopriRays == gated rays` and by
   the 0.118 mm AUTO-vs-DOPRI5 disagreement that RK4 produces.  Fixed by accepting AUTO in
   the DOPRI dispatch.
2. **The Phi tolerance mapping was optimistic.** `atol_Phi = atol_pos/r0` gave 1.19 mm of
   arc for a 1 mm request (the tolerance is a per-step budget and the accumulated error
   sits at it).  Halved; the observed error follows.

## Test result (`testCherenkovTangentTracer`, 3 cases, 26 assertions)

```
AUTO vs DOPRI5 on a mixed 10-ray set:  verdict flips 0, max |dP| 0.0144 mm, |dU| 2.7e-10
tangent served 4, refused 6, re-traced 6   (fallback exercised)
aperture-edge sweep (aim across r = 12 m, 61 rows): 0 verdict flips
```

The third case compares the regularised mode against an INDEPENDENT invariant quadrature
(fixed-step midpoint in u = sqrt(r-r_t), two resolutions to prove the reference converged).
It finds a **2.1 mm residual in path length and 3.3e-10 rad in swept angle**, while the two
tracers agree with each other to 14 um of impact and 1.3 um of path.  So a converged
reference and two integrators disagree at the 2 mm level and the test does not yet say
which is wrong.  Recorded, gated at 5 mm so a regression is still caught, and carried into
step 3 as its first task -- attributing exactly this kind of residual is what the
impact-error controller is for.

## Production A/B, same binary and seed (43 km muon, telescope, moving dish, config C)

```
                hits   bunch loop   gated rays   tangent served   tangent refused
AUTO            2075     16.688 s       197 652          56 759           140 893
DOPRI5          2075     18.631 s       197 652               -                 -
                10.4 % faster, identical yields

steps 10.13 M accepted / 6.85 M rejected (DOPRI5: 11.92 M / 7.50 M)
f-evals 102.1 M (DOPRI5: 116.7 M), 0 min-step failures, h[min/max] 1e-6 / 2e4 m
```

28.7 % of the gated rays are served by the regularised mode in a configuration whose
photons are emitted from 43 km down to the dish; the rest are refused by the closed-form
condition and re-traced.  The attribution closes exactly (56 759 + 140 893 = 197 652),
which it did not before the launch-and-continue branch was made to increment the served
counter -- that branch returns early and was invisible, a 5 615-ray gap.

`--cherenkov-tracer auto` is now the application default; `dopri5`, `rk4` and the
diagnostic `tangent` remain available.

# Step 3, task 1b -- the DOPRI5 plane tracer's 0.12 mm floor: VERIFIED

Probe `/tmp/dopri_crossing_probe.cpp`, full output `/tmp/p4_dopri_cross.log`. It reproduces
`rayDopri`'s plane mode (same `kappa`, ODE, cap rule, 5-point scan, bisection) rather than
calling it, so the re-landed variant needs no change to production code. Validation of the
reproduction: it returns the benchmark's own numbers **exactly** --

| cap | reproduction, shipped | `testCherenkovTracerBenchmark` |
|---|---|---|
| 5 km | 0.002477 mm | 0.0025 mm |
| 20 km | 0.118073 mm | 0.118 mm |
| 50 km | 0.123544 mm | 0.124 mm |

The reference is re-landed at `atol 1e-6` with no cap, so it cannot hide the defect it is
meant to detect; its own dense-vs-relanded spread is <= 2.9e-4 mm (last column of part 1),
i.e. the reference is clean.

## Part 1 -- step caps, telescope-normal plane, max over the 12-ray set

| cap | shipped (dense crossing) | re-landed | gain |
|---|---|---|---|
| 5 km | 0.002477 mm | 0.000251 mm | 9.9x |
| 20 km | **0.118073 mm** | **0.006426 mm** | 18.4x |
| 50 km | 0.123544 mm | 0.015906 mm | 7.8x |

Per ray at the 20 km cap, the two rays that dominate are the two that cross in the INTERIOR
of a long step: `down_10km` 0.118073 -> 0.000758 mm (h_cross 11.2 km, t 0.74) and
`layer_30km` 0.110154 -> 0.006426 mm (h 13.8 km, t 0.76). Rays that cross near a step
boundary show no gap at all: `tang_100km` 0.000006 -> 0.000011 mm (t 0.026) and
`horiz_700km` at the 5 km cap 0.000004 -> 0.000004 mm (t 0.88). The 0.12 mm floor is those
two rays, and it disappears when the crossing is landed instead of interpolated.

## Part 2 -- incidence angle at the 20 km cap (plane normal tilted off the ray)

`theta` is the tilt of the plane normal from the aim direction, `sin(theta_arr)` the
resulting ray-to-normal sin; 80 deg means the ray runs nearly IN the plane (grazing).

| ray | 0 deg | 15 deg | 45 deg | 80 deg |
|---|---|---|---|---|
| `horiz_200km` | 0.0718 / 0.0033 | 0.0782 / 0.0036 | 0.1247 / 0.0058 | **0.7026 / 0.4010** |
| `layer_30km` | 0.1102 / 0.0064 | 0.1294 / 0.0076 | 0.2792 / 0.0164 | 0.4209 / 0.0244 |
| `down_20km` | 0.0033 / 0.0029 | 0.0038 / 0.0034 | 0.0072 / 0.0064 | 0.0198 / 0.0174 |
| `tang_100km` | 0.000006 / 0.000011 | 0.000006 / 0.000013 | 0.000004 / 0.000017 | 0.000085 / 0.000075 |

(`shipped / re-landed`, mm.)

## Part 3 -- the crossing swept across one step (the decisive evidence)

The plane is translated along its normal, which moves `t_cross` through the step **without
changing the trajectory**: the step sequence depends only on the ray, the cap and the
tolerances, never on the plane. The integration is therefore identical in every row below, so
any variation with `t_cross` can only come from the crossing treatment.

| t_cross | 0.028 | 0.143 | 0.257 | 0.371 | 0.486 | 0.600 | 0.714 | 0.799 | 0.914 |
|---|---|---|---|---|---|---|---|---|---|
| shipped | 0.0063 | 0.0460 | 0.0718 | 0.0518 | 0.0039 | 0.0620 | 0.0866 | 0.0474 | 0.0132 |
| re-landed | 0.0034 | 0.0034 | 0.0033 | 0.0033 | 0.0033 | 0.0033 | 0.0035 | 0.0025 | 0.0029 |

`horiz_200km`, 20 km cap, mm. The shipped error swings **22x** with the crossing position and
collapses toward the step boundaries (`t -> 0, 1`); the re-landed error is flat to **1.4x**.
An integration that is identical in every row cannot produce the first line of that table; an
interpolation whose error vanishes at `t = 0` and `t = 1` produces exactly that shape. (Rows
that jump, e.g. 0.914 -> 0.028, fall in a different step, so the swing is step-to-step as well
as within-step -- each step has its own `h` and local curvature.)

This is the same signature measured directly in the tangent mode (task 1, interpolant profile
inside one step), reached from a completely different direction, which is why the attribution
is not a coincidence of one geometry.


Two findings, and the second is new:

1. **The near-normal conclusion does not transfer.** The shipped impact error grows with
   incidence -- 0.072 -> 0.703 mm for `horiz_200km` -- so at grazing incidence it reaches
   **0.70 mm, within a factor 1.4 of the 1 mm goal**, not the 0.12 mm the benchmark suggests.
2. **Re-landing stops being the lever there.** It buys 21x at 0 deg and 17x at 45 deg, but
   only 1.75x on `horiz_200km` at 80 deg, leaving **0.40 mm**. What remains is the
   integration error pushed through an ill-conditioned impact map (the same `1/cos` cliff the
   `CurvedRayTable` documents at 1.05 mm for mu 0.029-0.030), so the grazing case needs the
   impact-weighted controller, not a crossing fix.

Scope caveat: the sweep rotates the plane about the dish centre, which is the cleanest way to
vary the incidence at a fixed path; it measures conditioning, not a real off-axis telescope.
Production Telescope mode arrives within the Cherenkov apex angle (~1.3 deg) of the normal,
i.e. `sin(theta_arr) <= 0.023`, so the 45-80 deg rows are a stress case. Grazing arises for
the horizontal Ground-mode plane (photons hit it steeply) and for extreme-FoV geometries.

# Step 3, task 1 -- status: CLOSED. The 2.1 mm residual is the dense output.
(Closed before task 1b above; the two share one cause.)



The open item carried out of step 2 ("a converged reference and two integrators disagree by
2.1 mm and the test does not say which is wrong") is **attributed**. Nothing needed to be
changed in either integrator; the defect is a third component that neither of the old
assertions measured. Probe: `/tmp/tangent_residual_probe.cpp` (stages A--D below).

## It was not the leg bookkeeping

`exactAt(v0)` -- the two-leg quadrature `r_t -> r_0` plus `r_t -> r_0`, anchored at the
tracer's `(Phi, s) = (0, 0)` -- returns **0.00e+00 rad / 0.000e+00 mm**. The composition the
test uses (sum of the two legs when the perigee was visited, single leg otherwise) is exact.

## It was not the integration in v

An independent `Dopri5<4>` driving the *same* right-hand side to a **fixed** arrival radius
(no plane, no crossing solve, no dense output), against the quadrature at 4e6 points:

| atol (m) | Phi - ref | s - ref | rel. s | steps |
|---|---|---|---|---|
| 1e-3 | -2.77e-12 rad | -0.0177 mm | 1.4e-10 | 4 |
| 1e-5 | 8.15e-14 | +0.0005 mm | 4.2e-12 | 7 |
| 1e-7 | 2.85e-13 | +0.0018 mm | 1.5e-11 | 15 |
| 1e-9 | 8.79e-15 | +0.00006 mm | 4.6e-13 | 268 |

The quadrature's own self-convergence (2e5 vs 4e6) is 5e-16 rad / 3e-6 mm. So the reference
is converged *and* the v-integration is good to 0.02 mm **at the production tolerance**. The
residual does not move with tolerance, i.e. it is a *bias* rather than truncation, exactly as
the handoff's second probe predicted -- but the bias is in the interpolant, not the integration.

## It is the quartic dense output, read at an interior t

The crossing is located on the interpolant, and `geometric_path`/`Phi` are then read from
`dp.dense(tc)` -- not from the step state. Interpolant minus exact, at the same v, inside the
crossing step (`h = 15.09`, crossing at `t = 0.257`):

| t | 0.000 | 0.100 | **0.250** | 0.500 | 0.750 | 1.000 |
|---|---|---|---|---|---|---|
| dPhi [rad] | -2.8e-12 | +1.1e-10 | **+3.3e-10** | +7.2e-11 | -2.4e-10 | -2.3e-12 |
| ds [mm] | -0.018 | +0.722 | **+2.079** | +0.460 | -1.540 | -0.015 |

Both endpoints are inside 0.02 mm (matching stage A exactly) and the interior bulges to
+-2 mm, changing sign across the step: the O(h^5) signature of Hairer's quartic P-matrix
interpolant. This is the same mechanism the plan already suspected in item 3, now measured.

## Why two independent tracers appeared to agree

The residual is **one displacement along the ray**: `ds = r dPhi` holds to 1.3e-6 m. Measured
in the test: `max |ds - r dPhi| = 0.001343 mm` against `max |ds| = 2.138 mm` (a factor 1600).
An along-ray displacement
* is removed from the impact coordinate to first order -- only `ds sin(theta_arr)` survives,
  `theta_arr` being the angle between the ray and the plane normal (these rays arrive within
  3e-3 to 7.5e-3 rad of it, so 2.14 mm -> **0.006 to 0.014 mm** of impact), and
* cancels in `L`: the tracer then crosses the plane a fraction of 2 mm further along the
  same ray, which is why tangent-vs-DOPRI5 `dL` is 1.4 um and not 2.1 mm.

Case 1 of the test compares exactly those two observables, so a 2 mm defect was invisible to
it. The new gate is the along-ray identity, which a broken leg or a changed interpolant does
break.

## The remedy, measured (not adopted in production)

Re-landing the crossing step ON the plane -- take a single step of size `h*tc` from the step
start and bisect on that size until the endpoint lies on the plane, so the **5th-order state**
is the crossing instead of the interpolant's interior value:

| ray | as shipped: dPhi / ds | re-landed: dPhi / ds |
|---|---|---|
| turn_3100_d100 | 3.303e-10 rad / 2.104 mm | -2.78e-12 rad / -0.0177 mm |
| turn_3100_d100_a10 | 3.160e-10 / 2.013 | -1.07e-12 / -0.0068 |
| turn_3100_d100_a20 | 3.002e-10 / 1.912 | -1.11e-12 / -0.0071 |
| turn_3200_d100 | 3.356e-10 / 2.138 | -3.35e-12 / -0.0213 |

120x better, and exactly at the level of the step endpoints. Cost: the split step re-does one
step (`+7` RHS evaluations, more with the bisection; the cheap version is one Newton correction
on `h` using `d(plane distance)/dh = (ds/dv) cos(theta_arr)`, i.e. ~2 extra steps).

**Decision: not adopted, and the reason is measured, not assumed.** For the plane tracers the
ray arrives within the field angle of the plane normal, so the defect costs
`ds sin(theta_arr)`: 0.0144 mm here. Budget for integration is 0.50 mm (plan section 7), so it
uses 3 % of it. Adopting the split step would spend ~+14 % on the 28.7 % of rays served by the
tangent mode to remove an error that is already 70x inside the goal. Revisit if a target needs
the *path length* to better than 0.02 mm at a grazing crossing, or if the incidence angle grows
(the Ground-mode plane is horizontal and photons hit it steeply, `sin(theta_arr) ~ 1`, which is
where this defect maps 1:1 into impact).

## Consequence for the rest of step 3

* `Dopri5.hpp`'s header claim that the dense output "is what makes a LARGE step safe at a
  surface crossing" is **falsified at the 2 mm level** for an interior crossing. Corrected in
  place. The interpolant is still the right tool for *locating* the crossing (it is smooth and
  monotone there); it is not accurate enough to *be* the crossing state.
* `rayTracePlaneDopri` has the identical pattern (`dp.dense(tc)` -> `x_impact`, `geometric_path`),
  so it is exposed to the same defect, and the 0.12 mm benchmark floor at a 20 km cap is
  consistent with it in magnitude (`~2 mm x sin(theta_arr)`); that link is **not yet measured**
  and is the first thing to check when the benchmark set is extended (validation item 4).
  **UPDATE -- measured, see "Step 3, task 1b" above: the link is CONFIRMED (0.118073 ->
  0.006426 mm at a 20 km cap), and at grazing incidence the shipped error reaches 0.70 mm,
  which is the new risk item.**
* The acceptance rule of plan section 4 (`e_k = |J dy_k| <= b_k`) must therefore be evaluated on
  the **state**, not on the interpolant, or it will certify an error it cannot see.


