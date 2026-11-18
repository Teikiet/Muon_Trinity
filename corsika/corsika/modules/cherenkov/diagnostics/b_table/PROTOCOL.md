# Frozen validation protocol — cumulative-b absorption table

Status: **FROZEN** (amendments 1-5 recorded before any validation run). Production
integration is not authorized; this protocol governs diagnostics only.

---

## 1. Held-out accuracy is separated from boundary/guard development

The near-critical guard is derived **exclusively on development probes**. It is frozen
before any held-out ray is evaluated. Held-out failures are reported, never used to adjust
the guard.

- **Development probes** (guard derivation only): b/b_crit = 0.990, 0.995, 0.998, 0.999,
  0.9995, 0.9998, 0.9999, 0.99995, 1.0, 1.0000005, 1.000002, 1.00001, 1.00005;
  h = 1000, 5000, 20000, 40000 m; lambda = 300, 400, 500, 690 nm; f = 0.9713, 1.0, 1.0347.
  These overlap the previously inspected regression material and are therefore not fresh.
- **Held-out set** (frozen, evaluated only after the guard is fixed): see section 5.

## 2. b_crit is defined per configuration

`b_crit(f, lambda) = min over the ground-bound route of [ n_phase(r; f, lambda) * r ]`,
taken over r from r_obs to the emission radius, using the **same phase-index convention as
the production callback** (n = 1 + N(h)*f, so the wavelength enters only through f).

Two distinct probe families, not one:

- **equality probe**: b/b_crit = 1 exactly (the endpoint root is at the observation sphere);
- **beyond-critical probes**: b/b_crit = 1.0000005, 1.000002, 1.00001, 1.00005, 1.0005,
  1.05, 1.5 — each checked for **launchability independently**
  (b <= n_phase(r_e; f, lambda) * r_e) and, if launchable, integrated to confirm the ray
  turns above the observation sphere.

The previously listed held-out ratios are all strictly below 1 and therefore cannot
exercise either family; that is corrected here.

## 3. Boundary probes are explicit and reproducible

All heights in **metres**, all wavelengths in **nm**, all b ratios dimensionless.

- **Observation-height sequence** (ground boundary): h = 0, 1e-6, 1e-3, 1e-2, 0.1, 1.0,
  10, 100, 1000 m. Negative controls: h = -1e-3, -1.0, -100 m (must be rejected).
- **Table-ceiling sequence**: h = h_max - 1e-3, h_max - 1e-2, h_max, h_max + 1e-3,
  h_max + 1.0, h_max + 100 m, with h_max = 50000 m (the profile's last knot).
- **Wavelength endpoints and outside**: lambda = 295, 299.999, 300, 300.001, 699.999,
  700, 700.001, 705 nm.
- **Wavelength midpoints** means the arithmetic mean of two **adjacent listed native
  knots** (native spacing 5 nm in [300, 700]), e.g. (300+305)/2 = 302.5, (695+700)/2 = 697.5.
- **Near-critical sequence**: b/b_crit = 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 1.0.
- **Just-outside eligibility**: f = 0.9712, 0.9713, 1.0347, 1.0348; lambda as above;
  h above the ceiling; b > b_crit.

Recorded table bounds: h in [0, 50000] m; lambda in [300, 700] nm (native knots, 5 nm
spacing, 81 knots); b in [0, b_crit]. Dataset counts and source/log hashes are recorded in
the README of this directory; every run records the profile md5
(dfb2b7fd0d0c57ee0a704035a6a81e0b) alongside its output.

## 4. Guard procedure (reference-convergence degradation alone is insufficient)

The guard is admissible only if **all four** checks pass on development probes:

1. **Fixed reference uncertainty budget**: the reference must be converged so that its own
transmission uncertainty is <= 1e-5, i.e. 100x below the 1e-3 acceptance budget,
demonstrated by successive refinement on each development probe.
2. **Numerical reachability and root conditioning**: confirm the endpoint root exists and is
   well conditioned; report the conditioning measure and where it degrades.
3. **Floating-point sensitivity**: perturb b, h and lambda by relative 1e-15, 1e-13, 1e-11
   and require the table and reference transmission to move by less than 1e-5.
4. **Combined table error near the proposed cutoff**: |dT| <= 1e-5 on development rays within
   the region being considered for acceptance (i.e. 100x inside the acceptance budget).

If the four checks cannot establish a defensible guard, that region is **rejected to legacy
fallback and reported as unvalidated**. No clamping, no extrapolation.

## 5. Acceptance and timing details

**Accuracy acceptance.** |dT| <= 1e-3 for **every** eligible held-out ray, where dT is the
combined error against the independently converged reference with all approximation
contributions evaluated together (integration, interpolation, storage). Rejected rays are
reported as **counts by reason, separately** — they are not part of the accuracy statement.

**Held-out set (frozen).** h in {1500, 6000, 18000, 32000, 46000} m;
b/b_crit in {0.15, 0.35, 0.6, 0.8, 0.92, 0.97, 0.995, 0.9995}; lambda in {310, 355, 400,
445, 490, 535, 580, 625, 670, 695} nm plus adjacent-knot midpoints {312.5, 357.5, 447.5,
537.5, 627.5, 697.5} nm; f in {0.9713, 0.99, 1.0, 1.01, 1.0347} — **including interior f**,
which the earlier measurements never covered in combination.

**Timing protocol.** >= 7 runs, >= 2 discarded warm-up runs; report median, min, max and
IQR (variability, not the median alone). Machine and compiler recorded with the results.

- **Workload**: no measured production photon distribution is available to this workstream.
  The declared substitute is a synthetic distribution uniform in h in [1000, 31000] m,
  mu in [0.02, 0.92], lambda in [300, 700] nm, with a fixed seed (20261110) and a fixed
  sample size (200000 photons). This is a **declared synthetic workload, not a production
  one**, and any timing conclusion is limited accordingly.
- **Timed scope**: eligibility evaluation + geometry lookup + absorption lookup + dispatch.
- **Survival sampling**: **excluded** from the timed loop (documented explicitly).
- **Initialization and build costs are reported separately** from per-photon latency:
  profile load, absorption-table build, geometry-table build.
