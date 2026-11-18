# Curved-ray absorption studies — archived, closed

Status: **closed**. No further implementation is authorized. Production routing is
unchanged and the **legacy curved-ray tracer remains the absorption route**. These
diagnostics are kept as the record of what was measured and rejected.

Each file is a standalone `main()`; they were built against the framework include paths
(`corsika/` plus the Eigen include from the Conan cache) and run with the shipped profile
as `argv[1]`:

```
modules/data/CHERENKOV/atmosphere/abstable_MODTRAN_bottomup.dat
```

## What production retains

Absorption is applied through the legacy tracer's **curved-ray optical depth**.
`ec458103e` fixed the extinction sign convention (the shipped profiles store the column
increasing with height, so `alpha = +d tau/dh`; the previous ordering returned exactly 0);
the quadrature was sub-stepped to 50 m in `46f422d89`. The CurvedRayTable route is
deliberately excluded whenever absorption is enabled, because integrating tau along the
returned chord is not equivalent to integrating it along the curved ray.

**Observed absorption accuracy of the retained route, on tested rays:** `|d tau| ~ 1e-2`
before the sub-stepping fix (a left-endpoint rule at ~900 m geometric steps), and
`~1e-4` to `~1e-5` on the rays tested afterwards. This is an observation on those rays,
**not a domain-wide guarantee**.

## Designs evaluated and rejected

| diagnostic | representation | outcome |
|---|---|---|
| `diag_layers.cpp` | 25 uniform components, thickness-averaged `alpha_j` | **rejected** — combined `|dT| = 0.1619` (162x over the 1e-3 budget) on the 720-ray frozen set |
| `diag_direct.cpp` | per-photon direct native-profile integration | accurate (`7.9e-05`, 13x inside) at 1.9 MB and 6.1 ms init, but **16–183 ms/photon** — not viable |
| `diag_basis.cpp` | native trapezoidal altitude basis, boundary-corrected | the basis itself is **sound**: identity gate 2.6e-13 same-nodes, pointwise 4.1e-11, and `9.5e-05` against an independently converged reference on probes |
| `diag_coarse.cpp` | that basis with merged components + trilinear `(h, mu, f)` storage | **rejected** — combined `|dT| = 2.1e-2 … 2.7e-1` on the held-out set |

The boundary correction in `diag_basis.cpp` is derived, not tuned: the callback clamps
`h-25` to the profile floor while shrinking its denominator `h_hi - h_lo`, and for a
linear column the two cancel, so `alpha_cb` is **constant** across `[0, 25]` and across
`[49975, 50000]`. The correct boundary basis therefore **saturates at 1.0** (one-sided
skirt). Forcing `sum W_k = L` was never the criterion — the pointwise check is.

## Corrected reading of the coarsening study

- **M=25 fails held-out accuracy** (merge-only `2.685e-03`); it passed on the dev set
  (`8.648e-04`). **Whether a finer grid would succeed is unmeasured.**
- **M=50's combined error (`2.280e-02`) is not attributable to interpolation alone** — it
  also contains quadrature and float-storage contributions.
- **Production memory and build times are estimates**, extrapolated from the 693-ray
  diagnostic grid (1.87 GB / 234 MB; ~6.4–100 h). **No production build was made and no
  total peak memory was measured.**
- **The `0.142 us` lookup timing applies specifically to M=10** on this grid, not to the
  representation in general.
- The two gates were found to be mutually exclusive *on the tested configurations*: those
  that could pass accuracy (M >= 25) exceed the 256 MB cap once the 214 MB geometry table
  is included, and those that fit (M <= 10) fail accuracy. This rejects the tested
  uniform merged-component grids, **not every W-grid representation**.

## Frozen datasets — retained as regression sets

- **Regression set** (720 rays): h in {2500, 7500, 15000, 25000, 45000} m; mu in {0.999,
  0.7, 0.3, 0.12, 0.05, 0.01}; f in {0.9713, 1.0, 1.0347}; lambda in {300, 302.5, 400,
  402.5, 500, 502.5, 690, 700} nm.
- **Held-out set** used for the final coarsening evaluation: h in {3500, 12000, 35000} m;
  mu in {0.85, 0.45, 0.15, 0.03}; f in {0.99, 1.02}; lambda in {305, 420, 610, 512.5} nm.

Any future adaptive-grid, compression or restricted-domain proposal needs its **own
separately frozen study**; these sets must not be reused for tuning.

## Method notes worth carrying forward

- The independent reference must be **converged**, and convergence must be demonstrated
  (step refinement) — successive-step differences are uncertainty estimates, not bounds.
- Different quadrature nodes do **not** by themselves establish an independent
  *trajectory* reference; the trajectory model is shared.
- Isolated passes do not imply a combined pass; report both, and always evaluate on a
  held-out set frozen before tuning (here the dev-set merge result did not transfer).
