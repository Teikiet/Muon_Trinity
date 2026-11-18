# Ray-tracing settings comparison — six configurations, one muon point

Diagnostics only. **No production or geometry code was changed**; the only build action was a
rebuild of the application binary, which was stale. Outputs live under
`/scratch/general/vast/u1520754/chercmp/` (nothing was written into `corsika_results/`).

## Input point

```
--pdg 13            muon
-E 1e4              1e4 GeV
--injection-height 43000        --observation-level 2944
--zenith 90 --azimuth 270
--emcut 0.02 --mucut 0.43 --hadcut 0.5
--cherenkov-projection telescope --telescope-azimuth 270 --telescope-zenith 90
--telescope-x 0 --telescope-y 0 --telescope-z 0 --telescope-radius 15.0
--telescope-position-correction 1
--hadronModel SIBYLL-2.3d --cherenkov-output-format eventio --cherenkov-profile -v info
seeds 600 601 602 603 604
```

**`--pdg 13` is the muon.** An earlier attempt used `-p 5`, which is the **b quark**, and that
produced zero Cherenkov bunches in every geometry, projection, radius, correction and cut
combination tried. The application's own help lists the mapping
(`p=2212, gamma=22, e-=11, nu_e=12, mu-=13, nu_mu=14, tau=15, nu_tau=16`). The zero-hit result was
a wrong primary, not a geometry or projection problem — the proton control (`-p 2212`) produced
92 013 902 bunches and 240 hits, and `--pdg 13` produced 13 179 483 bunches and 2053 hits.

## Binary

`corsika-build/applications/muon_cherenkov_3d`, rebuilt 2026-11-11 14:37 from HEAD `fe4edbea8`.
The previous binary (2026-11-09) lacked `--disable-cherenkov-curved-table`, which configuration C
requires, and both it and the installed copy predated the current sources.

## Configurations

| cfg | geometry | absorption | extra flags |
|---|---|---|---|
| A | straight line | off | — |
| B | straight line | **on** | `--atmosphere-absorption-file` |
| C | curved, **legacy tracer** | off | `--enable-curvature --disable-cherenkov-curved-table` |
| D | curved, **legacy tracer** | **on** | both of the above |
| E | curved, **CurvedRayTable** | off | `--enable-curvature` |
| F | straight line | off | `--disable-cherenkov-fast-timing` |

There is no "curved + table + absorption" cell: absorption-active rays are excluded from the
table route (`Cherenkov.inl:763-766`), so D necessarily exercises the legacy tracer.

## Results — recorded Cherenkov hits per seed

```
seed       A       B       C       D       E       F
600     2079    2053    2074    2044    2074    2079
601      883     878     883     870     883     883
602      700     692     700     692     700     700
603      651     647     650     643     650     651
604      686     680     684     674     684     686
```

The runs are **deterministic** for a given (configuration, seed): F reproduces A **exactly** in
all five seeds. The scatter across seeds (~3x, 651 -> 2079) is genuine physics variance of this
point, not run-to-run noise.

## Paired differences (mean +- SE over the 5 seeds)

The ensemble standard error (+-273 hits) is about 20x larger than any configuration effect, so
the comparison must be **paired per seed**:

```
A-B  absorption, straight line   +9.8 +- 4.1  (+0.98 %)   diffs [26, 5, 8, 4, 6]   all positive
C-D  absorption, curved          +13.6 +- 4.2 (+1.36 %)   diffs [30, 13, 8, 7, 10] all positive
A-C  curvature, no absorption     +1.6 +- 0.9  (+0.16 %)   diffs [5, 0, 0, 1, 2]
A-E  curvature + table, no abs    +1.6 +- 0.9  (+0.16 %)   diffs [5, 0, 0, 1, 2]
C-E  table vs legacy tracer        0.0 +- 0.0  ( 0.00 %)   diffs [0, 0, 0, 0, 0]   EXACTLY zero
A-F  fast vs legacy timing         0.0 +- 0.0  ( 0.00 %)   diffs [0, 0, 0, 0, 0]   EXACTLY zero
```

## Findings

1. **The CurvedRayTable reproduces the legacy curved tracer exactly.** C-E is zero in every
   seed, not merely small: for this point the table route changes nothing about which photons
   are recorded.
2. **Curvature changes the yield by ~0.16 %** (5 hits out of 2079 in the worst seed, zero in
   two seeds). At zenith 90 deg the curvature branch is active (threshold is >60 deg), so this
   is a genuine measurement of its effect at this point: negligible.
3. **Absorption removes ~1.0 % (straight) to ~1.4 % (curved)**, consistently positive in all
   five seeds. Small because the path sits at 43 km altitude, where the air is thin.
4. **The fast/legacy timing choice does not touch the photon sample at all** (exactly zero
   difference in all seeds), which is what the design intends: it changes only the delay
   integral.
5. **The ray-tracing setting barely affects cost.** Per recorded photon, production is
   ~1.44-1.52 ms and the bunch loop ~0.98-1.19 ms in every configuration. The timing bucket is
   the only one that separates: ~0.2 us/photon for the closed-form path versus ~1.2 us/photon
   for the legacy segment loop - about 0.05 % of the per-photon budget. Total wall time is
   ~23.5-26 s in all six configurations, dominated by start-up.

## Caveats

- **Yield conclusions are only valid paired.** Unpaired, the +-273 hit ensemble error hides
  every effect measured here.
- **One input point.** Nothing here generalises to other energies, zenith angles or telescope
  radii; the ~1 % absorption figure is specific to a 43 km-altitude path.
- **The routing counters are not printed by the application.** `CurvedTableStats`
  (bypassed / hits / rejected / legacy-fallback) and `absorptionDrawCount()` exist in the API
  but no code path reports them, so for C/E the yield difference cannot be decomposed into
  "geometry" versus "table rejection". Surfacing them would be a production change.
- **Cost figures are per recorded photon** and the recorded count varies 3x across seeds, so
  the per-photon numbers are more meaningful than the bucket totals.

## Incidentally found (not the cause of anything here)

`Cherenkov.inl:824` comments that `refractiveIndexFunc_` "takes the HEIGHT above ground" and
passes `hEmit`, but the application's function actually takes the **absolute z**
(`muon_cherenkov_3d.cpp:1745`: `alt_km = (abs_z_m - seaLevelZ_m)/1000.0`). That call site
therefore passes a height where a radius is expected. It affects the `f_phase` fed to the
curved-table route, not the yields measured above. Reported, not fixed.

## Reproducing

```
bash matrix.sh          # 30 runs -> /scratch/general/vast/u1520754/chercmp/
python3 analyze.py      # table, paired differences, per-photon costs
```
