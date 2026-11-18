# README for the Cherenkov applications: `muon_cherenkov_3d.cpp`, `muon_cherenkov_test.cpp`, `muon_cherenkov.cpp`

This document is an end-to-end explanation of the CORSIKA 8 Cherenkov
applications in this directory. The **primary subject is
`muon_cherenkov_3d.cpp`**, which is the most complete of the three and the
reference implementation for everything described here. The other two
applications are only mentioned where they *differ* (see §20).

All three simulate air showers (primarily **muon** showers for
atmospheric-Cherenkov applications), record the Cherenkov photons that strike a
detector/telescope plane, and write the photon **bunches** out in one of two
formats:

* **CSV** (default) - plain text, columns documented below.
* **eventio** - a CHASM-compatible binary format readable by optics tools such as GrOptics.

`muon_cherenkov_3d.cpp` adds, on top of the shared Cherenkov machinery:

* a **projection-mode driven stopping plane** (ground / shower-axis / telescope),
* **global (telescope-centred) particle output** via `ObservationPlane::setGlobalOutput`,
* **primary injection overrides** (position and/or direction of the primary,
  decoupled from the nominal detector geometry) plus an
  `injection_geometry.yaml` provenance file,
* an optional **3-D observation-plane particle recorder** that writes
  `particles_3d.parquet` + `particles_3d_slices.yaml` (a time-resolved 3-D
  record of the shower as it develops toward the detector).

A large part of the file is a from-scratch **eventio writer** (the
`namespace eventio_writer { ... }` block) plus a **custom refractive-index /
Cherenkov-propagation handler**. The physics core (particle tracking,
hadronic/EM interaction models, continuous energy loss, thinning, radio
observers, etc.) is assembled from standard CORSIKA 8 framework classes.

---

## 0. Scope and feature matrix

| Application | Size | Role |
|---|---|---|
| `muon_cherenkov_3d.cpp` | 115 042 bytes | **Documented here.** Cherenkov photons (CSV/eventio), projection-mode driven stopping plane, global particle output, primary injection overrides, optional 3-D observation-plane recorder. |
| `muon_cherenkov_test.cpp` | 75 057 bytes | Same Cherenkov path and eventio writer, but **no 3-D recorder and no injection overrides**. The stopping plane is chosen by a `zenith > 80 deg` threshold (vertical plane through the telescope above it, horizontal plane at the core below it), `--zenith` is limited to 92 deg, and the ground particles are written in plane-local coordinates. |
| `muon_cherenkov.cpp` | 73 546 bytes | Earliest variant: the same `zenith > 85 deg` vertical-plane switch, no 3-D recorder, no injection overrides. |

Feature matrix (`X` = present, `-` = absent):

| Feature | `muon_cherenkov_3d.cpp` | `muon_cherenkov_test.cpp` | `muon_cherenkov.cpp` |
|---|---|---|---|
| Cherenkov photons, CSV output | X | X | X |
| Cherenkov photons, eventio output | X | X | X |
| `--cherenkov-projection` (ground/shower/telescope) | X | X | X |
| Dispersion / curvature / absorption options | X | X | X |
| Curved-ray table, default for eligible curved-ground rays (`--disable-cherenkov-curved-table`) | X | X | X |
| Telescope position correction | X | X | X |
| Stopping plane follows the projection mode | X | - | - |
| Global particle output (`setGlobalOutput`) | X | - | - |
| Primary injection overrides (`--injection-*`) | X | - | - |
| `injection_geometry.yaml` | X | - | - |
| 3-D observation-plane recorder (`--observation-steps`) | X | - | - |
| `particles_3d.parquet` / `particles_3d_slices.yaml` | X | - | - |
| `--zenith` up to 100 deg | X | - (max 92 deg) | - |

---

## 1. Purpose and high-level flow

The main program performs, in order:

1. Parse all command-line options (CLI11) and apply the verbosity level.
2. Register the random-number streams and set the seed.
3. Build the *environment* (curved Earth + US Standard atmosphere + geomagnetic field + refractive-index model).
4. Build the *primary particle* (single fixed energy, or energy sampled from a power-law range; direction from zenith/azimuth).
5. Compute the *nominal injection point* on a sphere around the observation point and construct a `ShowerAxis` through it.
6. Apply the optional *primary injection overrides* (`--injection-x/-y/-zenith/-azimuth`), validate the resulting
   starting point against the atmosphere, and (if any override was given) log the ray/plane crossing and write
   `injection_geometry.yaml`.
7. Set up the *Cherenkov subsystem*: a `HitRecorder` lambda, the output file (CSV or eventio), and a `CherenkovProcess`.
8. Load and spline-fit an *atmosphere refractive-index table*; build an `n(z, lambda)` function (phase + group index with dispersion).
9. Choose the *stopping (absorbing) observation plane* from the projection mode, and enable global (telescope-centred)
   particle output on it.
10. Assemble the *output writers* (energy loss, longitudinal profile, production profile, interaction writer, primary
    writer, ground plane, and optionally the 3-D recorder).
11. Assemble the *interaction models* (SIBYLL/QGSJet/EPOS/Pythia, UrQMD/FLUKA low-E, PROPOSAL, Sophia, TAUOLA, BetheBloch, thinning).
12. Optionally build the *3-D slice stack* (parallel to the stopping plane) and the `ThreeDPlaneRecorder`.
13. Run each shower through a full `Cascade` (`run_shower` lambda), then close the output library.
14. Write the eventio file (if selected), close the CSV, and log the total run time.


---

## 2. Headers / includes

The file pulls in essentially the entire CORSIKA 8 framework:

| Group | Headers | Purpose |
|---|---|---|
| Framework core | `Cascade.hpp`, `EnergyMomentumOperations.hpp`, `Logging.hpp`, `PhysicalUnits.hpp`, `InteractionCounter.hpp` | cascade orchestration, units, logging |
| Geometry | `CoordinateSystemRegistry.hpp`, `PhysicalGeometry.hpp`, `Plane.hpp`, `Sphere.hpp` | points, directions, coordinate systems, observation planes |
| Processes | `DynamicInteractionProcess.hpp`, `ProcessSequence.hpp`, `SwitchProcessSequence.hpp` | the process-chain machinery |
| Random | `RNGManager.hpp`, `PowerLawDistribution.hpp` | seeded RNG streams |
| Writers | `EnergyLossWriter`, `InteractionWriter`, `LongitudinalWriter`, `ProductionWriter`, `PrimaryWriter`, `SubWriter` | physics output |
| Output | `OutputManager.hpp`, `SaveBoostHistogram.hpp` | the on-disk output library, interaction histograms |
| Media | `corsika/media/CORSIKA7Atmospheres.hpp`, `corsika/media/Environment.hpp`, `corsika/media/magnetic/GeomagneticModel.hpp`, `corsika/media/refractivity/GladstoneDaleRefractiveIndex.hpp`, `corsika/media/density_and_composition/HomogeneousMedium.hpp`, `corsika/media/LayeredSphericalAtmosphereBuilder.hpp`, `corsika/media/medium/MediumPropertyModel.hpp`, `corsika/media/composition/NuclearComposition.hpp`, `corsika/media/ShowerAxis.hpp`, `corsika/media/magnetic/UniformMagneticField.hpp` | atmosphere, media properties and fields |
| Interactions | `BetheBlochPDG`, `Epos`, `EposLhcr`, `PROPOSAL`, `ParticleCut`, `Pythia8`, `QGSJetII`, `QGSJetIII`, `Sibyll`, `Sophia`, `StackInspector`, `EMThinning`, `LongitudinalProfile`, `ProductionProfile`, `FLUKA`/`UrQMD`, `TAUOLA` | hadronic/EM physics |
| Radio (mostly commented out) | `CoREAS`, `RadioProcess`, `ZHS`, `Observer`, `TimeDomainObserver`, `ObserverCollection`, `TabulatedFlatAtmospherePropagator` | radio emission (inactive here) |
| Cherenkov | `corsika/modules/cherenkov/Cherenkov.hpp` | the Cherenkov process |
| CLI | `CLI/App.hpp`, `CLI/Config.hpp`, `CLI/Formatter.hpp` | CLI11 parser |
| Std/OS | `boost/filesystem`, `<algorithm>`, `<cmath>`, `<fstream>`, `<sstream>`, `<chrono>`, `<cstdlib>`, `<cstring>`, `<iomanip>`, `<limits>`, `<string>`, `<optional>`, `<set>`, `<array>`, `<vector>` | OS/file/format utilities |

> **Note (CORSIKA 8 update):** the media headers now live in sub-directories of
> `corsika/media/` (`magnetic/`, `refractivity/`, `medium/`, `composition/`,
> `density_and_composition/`, `interfaces/`) and all media types are in the
> `corsika::media` namespace - hence the `media::` qualifiers used throughout this
> document. `CoordinateSystemRegistry.hpp` is included for the `CoordinateSystem`
> type.

The 3-D recorder uses `<optional>` (for the `std::optional<DirectionVector>` slice
normal and for the `std::optional<ThreeDPlaneRecorder>` in `main`) and `<set>`
(for the PDG filter).

> **Note:** `InteractionCounter.hpp` is included *first*, with a special
> comment - it uses `boost/histogram`, which can fail if other
> `boost/type_traits` headers were already included. Hence the
> `/* clang-format off */` guard.

The `#ifdef WITH_FLUKA` block selects FLUKA as the low-energy hadronic model
when built with `WITH_FLUKA`, otherwise it falls back to UrQMD.

---

## 3. Type aliases and setup

```cpp
using EnvironmentInterface = media::IRefractiveIndexModel<
    media::IMediumPropertyModel<media::IMagneticFieldModel<media::IMediumModel>>>;
using EnvType = media::Environment<EnvironmentInterface>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;
```

* `EnvType` - the concrete environment carrying all medium properties (refractive index -> medium properties -> magnetic field -> base medium model).
* `StackType` - the particle stack living in that environment.
* `TrackingType` - the trajectory tracker.
* `Particle` - a single particle in the stack.

---

## 4. Random number setup - `registerRandomStreams`

```cpp
long registerRandomStreams(long seed) { ... }
```

Each named component gets its *own* independent random stream. The function
registers streams named:

`cascade`, `qgsjet`, `qgsjetIII`, `sibyll`, `sophia`, `epos`, `epos-lhcr`,
`pythia`, `urqmd`, `fluka`, `proposal`, `thinning`, `tauola`,
`primary_particle`.

* If `seed == 0` the seed is drawn from `std::random_device` (nondeterministic);
  otherwise the user-supplied seed is used. The seed value is logged and then
  applied via `RNGManager<>::getInstance().setSeed(seed)`.

---

## 5. `MyExtraEnv` template

```cpp
template <typename T>
using MyExtraEnv = media::GladstoneDaleRefractiveIndex<
    media::MediumPropertyModel<media::UniformMagneticField<T>>>;
```

This stacks additional medium behavior *on top of* the base environment `T`:

* `UniformMagneticField<T>` - a constant magnetic field,
* `MediumPropertyModel<...>` - generic medium-property handling,
* `GladstoneDaleRefractiveIndex<...>` - a refractive-index model (used with the
  custom `n(z)` handler later).

This alias is passed as the second template argument to the atmosphere builder.


---

## 6. `main()` - command-line interface

The parser is built with **CLI11**. Options are grouped into the visual menus
used by the code: `Primary`, `Config`, `Misc.`, `Library/Output`, `Radio`,
`Thinning`, `Cherenkov`. The tables below list the group each option is
registered in (the CLI11 *default* group is marked `(default)`).

### 6.1 Primary options (group `Primary`)

| Option | Type / default | Range | Meaning |
|---|---|---|---|
| `-Z` | int | 0-26 | atomic number (requires `-A`) |
| `-A` | int | 1-58 | atomic mass number (requires `-Z`) |
| `-p, --pdg` | int | - | PDG code (`p=2212, gamma=22, e-=11, nu_e=12, mu-=13, nu_mu=14, tau=15, nu_tau=16`). Mutually exclusive with `-A`/`-Z`. |
| `-E, --energy` | GeV, default `0` | - | fixed primary energy |
| `--energy_range` | 2 doubles | `> 0` | low/high energy range; energies sampled from a power law |
| `--eslope` | -1.0 | - | spectral index `dN/dE = E^eslope` for sampling |
| `-z, --zenith` | 0.0 deg | 0-100 | primary zenith angle (deg). Upgoing showers (`> 90`) are supported. |
| `-a, --azimuth` | 0.0 deg | 0-360 | primary azimuth angle (deg) |
| `--injection-x` | 0.0 m | - | horizontal offset of the **primary injection point** in global x; `0` = on the nominal shower axis. A non-zero x/y offset switches `--injection-height` to the radial (production-height) interpretation, see §10.2 |
| `--injection-y` | 0.0 m | - | horizontal offset of the primary injection point in global y (see §10.2) |
| `--injection-zenith` | -1.0 deg | -1 to 100 | zenith angle of the **primary's injection direction**; `-1` = use the nominal shower zenith |
| `--injection-azimuth` | -1.0 deg | -1 to 360 | azimuth angle of the primary's injection direction; `-1` = use the nominal shower azimuth |
| `--injection-root` | `auto` | `auto`/`near`/`far` | which intersection of the shower axis with the injection sphere to use: `auto` (near whenever it is upstream, else far), `near`, or `far`. See §9 |

> `-E, --energy` is registered without an explicit `->group(...)` call, so CLI11
> places it in the *default* (ungrouped) section of `--help`; physically it
> belongs to the primary settings.

The `--injection-*` overrides are **primary-only**: they move/repoint the
injected particle but leave the detector geometry (shower axis, stopping plane,
telescope) anchored on the nominal direction. See §10. (`--injection-root` only
selects which root of the injection-sphere quadratic is used; it does not move
the detector.)

### 6.2 Configuration options (group `Config`)

| Option | Default | Range | Meaning |
|---|---|---|---|
| `--emcut` | 0.5e-3 GeV | 1e-6 - 1e13 | min kinetic energy for photons/e+/- in tracking |
| `--hadcut` | 0.3 GeV | 0.02 - 1e13 | min kinetic energy for hadrons |
| `--mucut` | 0.3 GeV | 1e-6 - 1e13 | min kinetic energy for muons |
| `--taucut` | 0.3 GeV | 1e-6 - 1e13 | min kinetic energy for tau leptons |
| `--max-deflection-angle` | 0.2 rad | 1e-8 - 1 | max deflection per tracking step |
| `--track-neutrinos` | off | - | enable neutrino tracking |
| `--observation-level` | 0 m | -1e3 - 1e5 | observation altitude above Earth radius |
| `--injection-height` | 112.75e3 m | -1e3 - 1e6 | radial injection altitude above sea level (polar coordinate). With `--injection-x/-y` offsets it is the true production height and the injection `z` is derived from it (see §10.2) |
| `--observation-steps` | 0 | `>= 0` | number of 3-D observation slices; `0` = no 3-D recorder (see §18) |
| `--obs3d-emin` | 0.0 GeV | `>= 0` | min kinetic energy recorded by the 3-D recorder (`0` = no filter) |
| `--obs3d-particles` | `""` | - | comma-separated PDG codes, or the preset `muons`; empty = record all |

The three `--obs3d-*`/`--observation-steps` options live in the **`Config`**
group (not in a dedicated group), and they are inert unless
`--observation-steps > 0`.


### 6.3 Miscellaneous options (group `Misc.`)

| Option | Default | Meaning |
|---|---|---|
| `--neutrino-interaction-type` | `both` | `neutral`/`NC`, `charged`/`CC`, or `both` |
| `-s, --seed` | 0 | RNG seed (0 = random) |
| `--force-interaction` | off | fix the location of the first interaction |
| `--force-decay` | off | force the primary to decay immediately |
| `--disable-interaction-histograms` | off | skip the interaction-histogram NPZ output |
| `-v, --verbosity` | `info` | `warn`/`info`/`debug`/`trace` (`trace` needs a Debug build) |
| `-M, --hadronModel` | `SIBYLL-2.3d` | HE hadronic model: `SIBYLL-2.3d`, `QGSJet-II.04`, `QGSJet-III`, `EPOS-LHC-R`, `Pythia8` |
| `-T, --hadronModelTransitionEnergy` | 10^1.9 ~ 79.4 GeV | LE/HE hadronic model switch energy |
| `--cherenkov-projection` | `ground` | `ground` (fixed z-plane), `shower` (perpendicular to the shower axis), `telescope` (circular detector) |
| `--cherenkov-output-format` | `csv` | `csv` or `eventio` (CHASM-compatible for GrOptics) |
| `--telescope-x` | 0.0 m | telescope centre x, relative to the observation height at the origin |
| `--telescope-y` | 0.0 m | telescope centre y, relative to the observation height at the origin |
| `--telescope-z` | 0.0 m | telescope centre z, relative to the observation height at the origin |
| `--telescope-radius` | 100.0 m | telescope aperture radius (used by `telescope` projection) |
| `--telescope-azimuth` | -1.0 deg | telescope pointing azimuth (0=+x, 90=+y, 180=-x, 270=-y); `-1` = use the shower azimuth |
| `--telescope-zenith` | -1.0 deg | telescope pointing zenith (0=down, 90=horizontal, 180=up); `-1` = use the shower zenith |
| `--atmosphere-file` | absolute path, see below | atmospheric model file (`Alt[km] rho thick n-1`) |
| `--disable-dispersion` | off | disable wavelength-dependent group-velocity dispersion (benchmarking) |
| `--enable-curvature` | off | curved photon-path correction for large zenith angles (> 60 deg) |
| `--disable-cherenkov-curved-table` | off | use the legacy curved-ray tracer instead of the `CurvedRayTable` lookup. Only affects **eligible** curved-ground rays (curvature enabled, ground projection, no absorption, refractive-index table present and covering the ray); the table is the default for those. See §13.2. |

### 6.4 Library / output options (group `Library/Output`)

| Option | Default | Meaning |
|---|---|---|
| `-N, --nevent` | 1 | number of showers to simulate |
| `-f, --filename` | `corsika_library` (required) | output library path; must **not** already exist (`CLI::NonexistentPath`) |
| `--compress` | off | compress the output directory to a tarball |

### 6.5 Thinning options (group `Thinning`)

| Option | Default | Range | Meaning |
|---|---|---|---|
| `--emthin` | 1e-6 | 0-1 | EM thinning onset (fraction of the primary energy) |
| `--max-weight` | 0 | `>= 0` | max weight (0 => Kobal's optimum x 0.5) |
| `--multithin` | off | - | keep thinned particles (weight = 0) |

### 6.6 Radio options (group `Radio`)

| Option | Default | Meaning |
|---|---|---|
| `--ring` | 0 | concentric star pattern of CoREAS/ZHS observers (radius = ring x 25 m). The observers are constructed but the `RadioProcess` registration is commented out, so no radio output is produced. |

### 6.7 Cherenkov options (group `Cherenkov`)

| Option | Default | Meaning |
|---|---|---|
| `--telescope-position-correction` | 0.1 | fraction of Cherenkov steps after which the telescope is re-centred on the most energetic charged particle (0 disables). Silently forced to 0 unless the projection mode is `telescope`. |
| `--atmosphere-absorption-file` | `""` | optional atmospheric absorption table for Monte-Carlo photon absorption. Text file with header `n_wavelengths n_heights`, then the wavelength array [nm], the height array [m] and the `ecoeff` matrix; use `tools/convert_npz_to_dat.py` to convert CHASM `abstable.npz` files. |
| `--disable-cherenkov-fast-timing` | off | disable the closed-form straight-line timing path and use the legacy 10–100-segment loop (validation reference). |
| `--cherenkov-profile` | off | enable the per-step Cherenkov profiler; prints `Cherenkov profile: total=… production=… bunchloop=… bunches=… recorded=…` at the end of the run. |
| `--cherenkov-timing-segments` | 0 | fixed integration-segment count for the non-fast timing path (0 = heuristic). Use e.g. `10000` as a high-resolution validation reference. |

### 6.8 Default `--atmosphere-file` and shipped data

The default is an **absolute path specific to this installation**:

```
/uufs/chpc.utah.edu/common/home/u1520754/corsika/modules/data/CHERENKOV/atmosphere/atmprof9.dat
```

Shipped tables in `modules/data/CHERENKOV/atmosphere/`:

* atmosphere profiles `atmprof1.dat`, `atmprof2.dat` ... `atmprof36.dat` (CORSIKA-style tables with the
  columns `Alt[km] rho thick n-1`), including the default `atmprof9.dat`;
* absorption tables `abstable_MODTRAN_new.dat` and `abstable_MODTRAN_bottomup.dat` for
  `--atmosphere-absorption-file`.

Other data sub-directories: `modules/data/CHERENKOV/telescope/` (telescope position CSVs),
`modules/data/CHERENKOV/corsika7_data/` (gamma/proton reference data) and
`modules/data/CHERENKOV/examples_input/` (example input cards).

### 6.9 Post-parse validation

After `CLI11_PARSE`:

* the verbosity is applied (with `trace` requiring a Debug build, otherwise the program aborts with an error),
* it is validated that either `--pdg` **or** both `-A` and `-Z` are given,
* `registerRandomStreams` is called with the seed.


---

## 7. Environment and root coordinate system

```cpp
EnvType env;
const CoordinateSystem& rootCS = env.getCoordinateSystem();
Point const center{rootCS, 0_m, 0_m, 0_m};
Point const surface_{rootCS, 0_m, 0_m, constants::EarthRadius::Mean};
media::GeomagneticModel wmm(rootCS, center, corsika_data("GeoMag/WMM.COF"));
```

* `center` - Earth center (origin of `rootCS`).
* `surface_` - a point on the mean-radius sphere (sea level).
* `wmm` - world magnetic model loaded from `GeoMag/WMM.COF`.

The atmosphere is built with Keilhauer's parametrization of the US Standard
Atmosphere:

```cpp
media::create_5layer_atmosphere<EnvironmentInterface, MyExtraEnv>(
    env, media::AtmosphereId::USStdBK, center, 1.000327, surface_,
    media::Medium::AirDry1Atm,
    // Bx=28.6 uT (horizontal), Bz=-46.4 uT (downward) in +z-up convention
    MagneticFieldVector{rootCS, 28.6_uT, 0_uT, -46.4_uT});
```

* `1.000327` is the sea-level phase refractive index (`n-1 = 0.000327`).
* The magnetic field is set to a realistic mid-latitude-like value.

---

## 8. Primary particle construction

The beam code is chosen from PDG or Z/A:

```cpp
// --pdg
beamCode = convert_from_PDG(PDGCode(app["--pdg"]->as<int>()));
// or Z/A
if ((A==1)&&(Z==1)) beamCode = Code::Proton;
else if ((A==1)&&(Z==0)) beamCode = Code::Neutron;
else beamCode = get_nucleus_code(A, Z);
```

`neutrinoPrimary` records whether the primary is a neutrino
(`bool const neutrinoPrimary = is_neutrino(beamCode);`); it is used later to
decide whether the neutrino interaction process is added to the shower sequence
(see §16.2 and §16.5).

Energy is resolved to `eMin/eMax` either from `--energy` (equal min/max) or
from `--energy_range` (swapped with a warning if low > high). If neither is set,
the program logs a critical message and returns 0.

Direction:

```cpp
thetaRad = zenith/180*M_PI;  phiRad = azimuth/180*M_PI;
(nx,ny,nz) = (sin(theta)*cos(phi), sin(theta)*sin(phi), -cos(theta))
propDir = DirectionVector(rootCS, {nx, ny, nz});
```

The `-cos(theta)` z-component points *downward* for a downgoing shower
(zenith < 90 deg). For `zenith > 90 deg` the shower propagates upward; the
injection geometry below handles both cases with the same formula.

`primary_pdg` (used for the eventio header) is the `--pdg` value if given, or
otherwise the CORSIKA-style nucleus code `1000*Z + A`.

---

## 9. Geometry: observation point, injection sphere, shower axis

```cpp
observationHeight = obsLevel*1m + EarthRadius::Mean;
injectionHeight   = injLevel*1m  + EarthRadius::Mean;
showerCore = Point(rootCS, 0, 0, observationHeight);
```

The core is placed at the observation altitude. The **nominal** injection point
lies on a sphere of radius `injectionHeight` around Earth's center, back along
the propagation direction. The code solves the quadratic
`s^2 + 2*obsH*cos(theta)*s + obsH^2 - injH^2 = 0` and selects one of its two
roots:

```cpp
a_coeff  = -observationHeight * cos(thetaRad);
disc     = injH^2 - (sin(thetaRad)*obsH)^2;
tNear    = a_coeff - sqrt(disc);   // intersection closer to the detector
tFar     = a_coeff + sqrt(disc);   // intersection farther from the detector
backwardDir = (-sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta));

bool const useNear = (injectionRoot == "near") ||
                     (injectionRoot == "auto" && tNear > 0_m);
t = useNear ? tNear : tFar;
injectionPos = showerCore + backwardDir * t;
```

`--injection-root` selects the root (`auto`/`near`/`far`, default `auto`):

* downgoing showers (theta <= 90 deg): `tNear <= 0`, so only `tFar` is upstream
  and is used by both `auto` and `far`;
* upgoing Earth-skimming showers (theta > 90 deg, injection sphere below the
  observation level): `tNear > 0` is upstream as well, and `auto` picks the
  **near** root so the primary climbs to the detector without crossing the Earth
  (the far root would place the injection point on the opposite side of the
  planet);
* if the selected root is not upstream of the detector, the program logs an error
  and exits with `EXIT_FAILURE`.

```cpp
media::ShowerAxis const showerAxis{injectionPos, (showerCore - injectionPos) * 1.2, env};
auto const dX = 10_g / square(1_cm);   // longitudinal binning (g/cm2)
```

`media::ShowerAxis` is stretched ~20% beyond the core-to-injection distance because
the longitudinal profiles can extend past the core at large zenith angles. All
longitudinal writers bin along this axis with the density-scaled spacing
`dX = 10 g/cm2`. The shower axis is **anchored on the nominal injection point**
`injectionPos`, even when the primary injection is overridden.


---

## 10. Primary injection overrides (`--injection-*`)

The detector geometry (nominal shower axis, stopping plane, telescope) is always
built from `--zenith`/`--azimuth`/`--injection-height` and is **never** modified
by the injection overrides. The overrides only change where the primary particle
starts and which way it points, which is what makes it possible to study, e.g.,
a primary that misses the nominal axis, or an off-axis telescope event.

### 10.1 Direction override

```cpp
DirectionVector injDir = propDir;                      // nominal direction
bool const directionOverridden =
    (injection_zenith_deg >= 0.0) || (injection_azimuth_deg >= 0.0);
if (directionOverridden) {
  injTheta = (injection_zenith_deg >= 0 ? injection_zenith_deg : --zenith)/180*pi;
  injPhi   = (injection_azimuth_deg >= 0 ? injection_azimuth_deg : --azimuth)/180*pi;
  injDir   = DirectionVector(rootCS, {sin(injTheta)*cos(injPhi),
                                      sin(injTheta)*sin(injPhi),
                                     -cos(injTheta)}).normalized();
}
```

Each angle falls back independently, so `--injection-zenith` alone keeps the
nominal azimuth and vice versa.

### 10.2 Position override

The offsets are applied along the **global horizontal axes** of the root
coordinate system. When an x/y offset is given, `--injection-height` is
interpreted as the **radial altitude above sea level** (a polar coordinate,
i.e. the true production height `h_prod`), and the injection `z` is **derived**
from it so that the point stays on the sphere of radius `R + injection_height`:

```cpp
injE1 = (1, 0, 0);  injE2 = (0, 1, 0);                 // global x, y
injectionPosPrimary = injectionPos;
bool const positionOverridden = (injection_x_m != 0.0) || (injection_y_m != 0.0);
if (positionOverridden) {
  double const xNew = injectionPos.getCoordinates().getX()/1_m + injection_x_m;
  double const yNew = injectionPos.getCoordinates().getY()/1_m + injection_y_m;
  double const rNew = injectionHeight/1_m;             // R + injection_height
  double const rho2 = xNew*xNew + yNew*yNew;
  if (rho2 >= rNew*rNew) { /* error: offset beyond the injection sphere */ }
  double const zNew = std::sqrt(rNew*rNew - rho2);     // keep |P| = R + h
  injectionPosPrimary = Point(rootCS, xNew*1_m, yNew*1_m, zNew*1_m);
}
```

With `--zenith 0` the nominal injection point is `(0, 0, injectionHeight)`, so
the x/y offsets are the **absolute** horizontal coordinates of the production
point and its radial distance is exactly `R + injection_height`. The injection
point therefore always lies on the production sphere `|P| = R + h_prod`, which
is what lets a backtracked muon (specified by `h_prod` and its horizontal
position) land exactly on the intended point. The `rho2 >= rNew*rNew` guard
rejects offsets that would fall beyond the injection sphere.

> **Changed behaviour:** previously the offset was a *flat* Cartesian shift that
> preserved the nominal `z`, so an off-axis point ended up at radius
> `sqrt(x² + y² + (R + h_prod)²) > R + h_prod` (i.e. a slightly too-high
> altitude). The `z` is now recomputed from the radial `--injection-height`.

### 10.3 Validation against the atmosphere (always performed)

The radial distance of `injectionPosPrimary` from the Earth centre is converted
to a height above sea level and checked against the validity of the CORSIKA 7
5-layer atmosphere:

```cpp
constexpr double kAtmosphereTopAboveSeaLevel_m = 112.8e3;
injRadius_m = (injectionPosPrimary - center).getNorm()/1_m;
hInj_m      = injRadius_m - EarthRadius::Mean/1_m;
```

* `hInj_m < 0` -> the start point is **below sea level**: log error, `return EXIT_FAILURE`.
* `hInj_m > 112.8e3` -> the start point is **above the atmosphere** (above the top of the
  constant-density outer layer): log error, `return EXIT_FAILURE`.

The check is unconditional, so it also validates the nominal
`--injection-height` even when no override is given. When an x/y offset is
given, `injectionPosPrimary` lies exactly on the sphere `R + injection_height`
(§10.2), so `hInj_m` equals `--injection-height` and the check validates that
value directly.


### 10.4 Ray/plane crossing diagnostics

When at least one override is given, the primary's straight-line ray is
intersected with the stopping plane and logged:

```cpp
denom = obsPlaneNormal . injDir;
// |denom| < 1e-12  -> warning: ray is parallel to the plane
num   = obsPlaneNormal . (obsPlaneCenter - injectionPosPrimary);
s_m   = (num/denom)/1_m;
// s_m <= 0         -> warning: primary is already past the plane / moving away
hit   = injectionPosPrimary + injDir*(s_m*1_m);
d     = hit - obsPlaneCenter;
hit_x_m = injE1 . d;  hit_y_m = injE2 . d;  hit_r_m = hypot(hit_x_m, hit_y_m);
```

Notes on the output:

* `denom` is used through `magnitude()`, so the **sign** is dropped on purpose -
  a downgoing shower has `denom = -1` and would otherwise be misclassified as
  parallel to the plane.
* `hit_x_m`/`hit_y_m` are measured along the **global x/y axes**, not along the
  plane's own x-axis, so they are directly comparable to `--telescope-x/-y`.
* In `telescope` projection mode the code additionally warns if
  `hit_r_m > --telescope-radius` ("MISSES the telescope") or confirms that the
  crossing is inside the aperture.
* In `ground`/`shower` mode it notes that the recording plane is infinite, so a
  non-parallel ray is always recorded, just at non-zero horizontal (x, y).

### 10.5 `injection_geometry.yaml`

Written into the Cherenkov directory (`dirname(-f)`) whenever position and/or
direction are overridden, with `std::setprecision(9)` and bare doubles (no
`phys::units` suffix):

```yaml
earth_radius_m: <m>
nominal:
  direction: [nx, ny, nz]
  injection_point_m: [x, y, z]
primary:
  direction: [nx, ny, nz]
  injection_point_m: [x, y, z]
  radius_m: <m>            # |injectionPosPrimary - EarthCentre|
  altitude_m: <m>          # radius_m - EarthRadius::Mean
offset_global_m: {x: <m>, y: <m>}
angle_to_axis_deg: <deg>   # angle between nominal and primary direction
observation_plane:
  center_m: [x, y, z]
  normal: [nx, ny, nz]
  x_axis: [1, 0, 0]        # = injE1 (global x)
  y_axis: [0, 1, 0]        # = injE2 (global y)
  crosses_plane: true|false
  crossing_local_m: {x: <m>, y: <m>, r: <m>}   # or: null
  telescope_radius_m: <m>
```

`crossing_local_m` is written as `null` (rather than zeros) when the ray never
crosses the plane, so a missing crossing cannot be mistaken for a crossing
exactly at the plane centre.

> **Caveat:** the file is written **once per run**, before the shower loop, and
> is truncated (not appended). With `-N > 1` it therefore documents the run's
> nominal/primary geometry, not the per-event sampled energies.

The values actually injected are logged per shower:

```
Primary name / PDG / Z-A, Primary Total Energy, Primary Momentum,
Primary Direction:    <injDir.getNorm()>      <- norm only, see caveat
Point of Injection:   <injectionPosPrimary coordinates>
Shower Axis Length:   |showerCore - injectionPos| * 1.2
```

> **Caveat (cosmetic):** `Primary Direction:` prints `injDir.getNorm()`, i.e. the
> norm of the (unit) direction vector - roughly `1` - instead of the direction
> components. The real injected direction is in `injection_geometry.yaml`
> (`primary.direction`) and in `primary.parquet`.

The primary is added to the stack with the overridden values:

```cpp
auto const primaryProperties =
    std::make_tuple(beamCode, eKin, injDir.normalized(), injectionPosPrimary, 0_ns);
stack.addParticle(primaryProperties);
...
primaryWriter.recordPrimary(primaryProperties);
```


---

## 11. Cherenkov subsystem

### 11.1 Ground Z and output path

```cpp
double const groundZ_m_abs = observationHeight / 1_m;   // absolute coordinate
boost::filesystem::path cherDir = basePath.parent_path();  // dir part of -f
if (cherDir.empty()) cherDir = ".";
boost::filesystem::create_directories(cherDir);
```

* The "ground" Z for projection is the absolute observation altitude.
* The Cherenkov files are placed in the directory part of `--filename`
  (e.g. `-f muon_test` -> current dir; `-f /home/user/test_muon` -> `/home/user`).

### 11.2 Output file selection

* **eventio** format -> `cherenkov_hits.dat` in `cherDir`.
* **csv** format -> `cherenkov_hits.csv` in `cherDir`.

For CSV, the header written is:

```
x_m,y_m,z_m,wavelength_nm,z_production_m,dir_x,dir_y,dir_z,
dir_ground_x,dir_ground_y,dir_ground_z,time_ns,weight,n_minus_1,path_length_m
```

The `n_minus_1` column is stored as `n_prod - 1.0` (the refractive-index surplus
at the production point).

The eventio header is pre-filled with the **nominal** shower and telescope
parameters:

```cpp
primary_energy = (--energy > 0) ? --energy
               : (energy_range.size() ? energy_range[0] : 1.0);
tel_z          = observationHeight/1_m + --telescope-z;      // absolute
eventioWriterPtr->set_event_params(thetaRad, phiRad, primary_energy, primary_pdg,
                                   tel_x, tel_y, tel_z, tel_radius);
```

> **Caveat:** `primary_energy` is *not* the per-event sampled energy - with
> `--energy_range` the eventio `EventHeader` always carries the lower bound of
> the range, and the zenith/azimuth are the nominal ones.

### 11.3 The `HitRecorder` lambda

The recorder is a callable receiving one *photon bunch*. Two implementations
are selected at run time:

* **eventio recorder** - converts meters to cm for `(x, y, zem)`, keeps
  `cx=dir_x`, `cy=dir_y`, `time_ns`, `weight` (= no. of photons in the bunch)
  and `wavelength_nm`, and pushes them into the writer via `add_photon(...)`.
* **CSV recorder** - writes the 15 columns above, one line per bunch.

Both increment `cherenkovHitCount`. Unused variables (`z_m`, `dir_z` and the
ground-direction components, `n_prod`, `path_length_m`) are cast to `void` in
the eventio path to suppress warnings.

> **Caveats:**
> * `cherenkovHitCount` is **cumulative over all showers**; the per-shower log
>   line `Total Cherenkov hits recorded: N` therefore shows a running total, not
>   that shower's count.
> * Both output paths accumulate **all** showers into a single file: the CSV has
>   no event column, and the eventio writer merges every shower's bunches into
>   one `EventHeader`/`TelescopeData`/`EventEnd` record (see §19.8).

### 11.4 Projection mode

```cpp
projModeStr in {ground, shower, telescope}
-> ProjectionMode::{Ground, ShowerAxis, Telescope}
```

`--telescope-position-correction` is silently zeroed (with an informational log
message) unless the projection mode is `telescope`.

### 11.5 Telescope parameters and basis

`tel_x`, `tel_y`, `tel_z_relative`, `tel_radius` come from the CLI; the absolute
`Z` is `observationHeight + tel_z_relative`. The pointing is:

```cpp
tel_azimuth_rad = (tel_azimuth_deg >= 0) ? tel_azimuth_deg*pi/180 : phiRad;
tel_zenith_rad  = (tel_zenith_deg  >= 0) ? tel_zenith_deg *pi/180 : thetaRad;
```

i.e. if not specified, the telescope points along the shower direction.
`shower_core_{x,y,z}` are read off the core point for the process constructor.

The observing-plane basis is deliberately identical to the one used inside
`corsika::cherenkov::CherenkovProcess` (see
`corsika/modules/cherenkov/Cherenkov.inl`, `computeTelescopeBasis()`):

```cpp
telBasisZ = ( sin(zen)*cos(az),  sin(zen)*sin(az), -cos(zen) );  // pointing = plane normal
telBasisX = ( -sin(az),           cos(az),          0.       );  // in-plane x-axis (horizontal)
// in-plane y-axis, where needed: telBasisZ x telBasisX
```

`telCenter = Point(rootCS, tel_x*1_m, tel_y*1_m, tel_z*1_m)` is used both for the
`telescope` stopping plane and as the **origin of all particle output
coordinates** (§15.2).


---

## 12. Atmosphere table loading and cubic spline

The file `--atmosphere-file` supplies altitude-resolved `n-1` values (rows
`# Alt[km] rho thick n-1`; `#`-prefixed or empty lines skipped). The data is:

1. **Sorted** by altitude (explicit `std::sort` by index).
2. `ln(n-1)` is precomputed (values `<= 1e-10` clamped to `log(1e-10)` to avoid `log(0)`).
3. A **natural cubic spline** (second-derivative coefficients, tridiagonal solve)
   is built over `ln(n-1)` vs. altitude. This guarantees C^2-continuous,
   oscillation-free interpolation of the refractive index.
4. A **fine lookup table** is pre-interpolated with a step of **0.01 km = 10 m**
   from table min to max altitude, giving refractive-index precision better than
   1e-6 with a fast linear interpolation per photon.

If the coarse table has fewer than 3 points or the spline cannot be built, the
fine table is skipped and the code falls back to the coarse data (or to a
constant `n` when the file cannot be opened at all, with a warning).

The load is reported as:

```
Loaded <n> atmosphere data points from '<path>'
  Altitude range: <min> - <max> km
  n-1 range: <min> - <max>
  Using cubic spline interpolation for smooth refractive index
Created high-resolution lookup table with <n> points (step=0.01 km)
  Refractive index precision: better than 1e-6
```

---

## 13. The refractive-index function `refrIndexFunc`

This lambda returns `std::pair<double,double> = (n_phase, n_group)` given an
absolute Z (m) and a wavelength (nm):

1. **Base `n-1` at altitude**:
   * below the lowest tabulated altitude -> clamp to the first value,
   * above the highest tabulated altitude -> *exponential extrapolation* using the
     scale height of the last two points (with a hard floor at `n-1 = 1e-10`),
   * otherwise -> binary-search bracketing and **linear interpolation of `ln(n-1)`**
     on the (fine) grid.
2. **Dispersion**, unless `--disable-dispersion`:
   ```
   ratio    = (400 / lambda)
   f_phase  = 0.967 + 0.033  * ratio^2.5
   f_group  = 0.967 + 0.1155 * ratio^2.5
   n_phase  = 1 + (n-1)*f_phase
   n_group  = 1 + (n-1)*f_group
   ```
   The phase index drives the Cherenkov angle; the group index drives the arrival
   time.

This gives wavelength-dependent photon speeds and angles (realistic Cherenkov
ring broadening in time and angle).

### 13.1 Closed-form straight-line timing (P1)

Because the altitude and wavelength dependences factorise
(`n(z,λ) − 1 = N(z)·f(λ)`), the group-delay integral along a **straight** photon
path collapses to a closed form:

```
∫ (n_group − 1) ds = f_group(λ) · ∫ N(z) ds ,   ∫ N ds = (A(z_end) − A(z_start)) / u_z
```

where `A(z) = ∫ N dz'` is a vertical column integral precomputed on the same
uniform 10 m grid the runtime lookup uses (`RefractiveIndexTable`, built from
`atmAlt_final`/`atmLnNminus1_final`). This replaces the legacy 10–100-segment
loop exactly. The switch to a midpoint form is on `|Δz| < 10 m` (where the closed
form's `ΔA` would cancel); the midpoint relative error is
`≈ (Δz/H)²/24 ≈ 6.5e-8` at `Δz = 10 m`, `H ≈ 8 km`.

The fast path is enabled automatically when a valid table is supplied **and** the
refractive-index function is separable (checked at construction). It changes
**only** the arrival time: photon positions, directions, wavelengths, weights and
`n_minus_1` are bit-identical to the legacy path, and the times agree with a
10 000-segment numerical integral to the CSV's 6-significant-figure precision.

> **Flat-Earth note:** the altitude used by the refractive-index and absorption
> lookups is the Cartesian `z − R_earth`, **not** the spherical `|r| − R_earth`.
> This is a pre-existing approximation (unchanged by this work) and applies even
> when `--enable-curvature` is recommended above 60°; the closed-form integral is
> exact *with respect to this flat-Earth model*.

**Measured performance** (single 100 GeV muon and 1 TeV tau, **eventio** output,
median of 3 runs, via `--cherenkov-profile`):

| | 100 GeV muon (185 k bunches) | 1 TeV tau (4.1 M bunches) |
|---|---|---|
| `doContinuous` legacy → closed form | 1.36 s → 0.13 s (**10.4×**) | 39.4 s → 3.09 s (**12.8×**) |
| timing sub-part | 1.26 s → 0.045 s (**28×**) | 37.1 s → 1.00 s (**37×**) |
| eventio recorder | ~0.017 s | ~0.34 s |
| cascade `EAS.run` | 1.54 s → 0.26 s | 45.6 s → 9.1 s |
| whole run | 15.0 s → 14.1 s (1.06×) | 58.5 s → 22.4 s (**2.6×**) |

The single-muon run is dominated by ~13 s of one-time model initialisation, so the
end-to-end gain is small there; for a real cascade (tau) it is 2.6×. Note the
recorder cost is strongly output-format dependent: **eventio ≈ 0.02–0.34 s**,
whereas **CSV ≈ 1.4 s** (text formatting), which limits the bunch-loop speedup to
~1.8× with CSV.

### 13.2 Curved-ray table (P7)

`--enable-curvature` replaces the straight-line projection with an exact spherical
(Bouguer) ray trace. That trace integrates an ODE, so it costs ~17 us to ~2 ms per
photon depending on how grazing the ray is. The `CurvedRayTable`
(`corsika/modules/cherenkov/CurvedRayTable.hpp`, header-only, no framework deps)
precomputes the ray geometry on a two-level altitude x mu grid and reduces each
lookup to a bicubic Hermite interpolation, paying a one-off build instead.

**Eligibility.** One predicate, evaluated in `CherenkovProcess::traceCurvedGround`:

```
enableCurvature_ && enableCurvedTable_ && projection mode == Ground
  && !(absorption enabled)          // curved-path absorption is deferred
  && refractiveIndexTable_ valid    // the table is built from it
  && coverage: [r_obs, r_obs + 120 km] within [zMin(), zMax()]
```

Every ineligible ray - and every ray the table *rejects* - goes through the
**unchanged** legacy tracer. Plane projection (`rayTraceSphericalToPlane`) is not
affected. The table is built lazily, once per process, under a `once_flag`; copies
of `CherenkovProcess` share one table and one counter set. If the refractive-index
table does not cover the sampled interval the table is not built at all and the
reason is logged: the RI table clamps below `zMin()` and extrapolates above
`zMax()`, so a mismatch must not silently build.

**Accounting.** `CherenkovProcess::curvedTableStats()` separates five classes,
which satisfy by construction:

```
total = bypassed + attempted
attempted = hits + sum(rejected)
sum(rejected) = legacyFallbackSucceeded + legacyFallbackFailed
finalDrop = bypassedLegacyFailed + legacyFallbackFailed
```

`rejected[]` is indexed by `CurvedRayTable::Status` (`grazing`,
`dispersion-out-of-range`, `upgoing`, `beyond-critical`, `above-table`,
`emission-below-obs`, `invalid-input`, `invalid-table`). A table rejection rescued
by the legacy tracer reproduces the baseline acceptance for that ray; whether a
table **hit** would have failed in the legacy tracer is **not** observable from
these counters - that needs a separate both-routes comparison run.

**Accuracy.** Worst case per observable over the strict suite, against an
independent RK4 reference integrated with ds = 2 m:

| observable | dataset | worst |
|---|---|---|
| impact distance `|dd|` | strict (h, mu) grid, mu in [2e-4, 1] | 4.2 mm |
| impact distance `|dd|` | break region (h ~ 119.9 km, mu ~ 0.888) | 3.2 mm |
| time `|dt|` | strict grid | 1.5e-2 ns |
| arrival angle `|dtheta|` | strict grid | 6.2e-10 rad |

Below `muMin() = 2e-4` the near-horizontal integrand is not resolved and the table
**rejects** rather than extrapolating (the physical 1 PeV Cherenkov range is
mu ~ 0.002-0.06, i.e. 10-300x above the cutoff). `f_phase` outside `[0.96, 1.04]`
is rejected too; the production range is 0.9713-1.0347.

**Cost** (one process, h = 20 km, single run; table build 29.5 s and 214.3 MB
against the project 256 MB cap, once per process):

| route | per photon |
|---|---|
| table lookup, integrated (`traceCurvedGround`) | **0.449 us** |
| legacy curved trace, mu = 0.5 | 16.74 us (37x slower) |
| legacy curved trace, mu = 0.05 (thousands of steps) | 671.6 us (1496x slower) |

The integrated route costs ~0.13 us more than the bare table lookup (0.32 us)
because it also computes the dispersion factors. A full-event benchmark
(`--cherenkov-profile`) has **not** been run for this change.

**Known limitations** (reported in the log, not asserted as accuracy):

* The propagation-plane convention repair made the legacy tracer's *signed* impact
  position correct (it previously landed on the wrong side of the emission, 18.4 km
  at h = 20 km, mu = 0.9). It did **not** change the legacy's accuracy. Against the
  independent reference the legacy's **relative path error** is 3.2e-7 ... 2.213e-1
  (dimensionless; h-independent, mu-only) and its **position residual** is
  0.011-0.081 m at mu = 0.9, ~19 m at mu = 0.1, ~148 m at mu = 0.05 and ~15.6 km at
  mu = 0.01. Both are dominated by the deliberate 1 km straight-line switch;
  **full attribution of the 22 % path error is an open workstream (P3)** and is not
  claimed here.
* The legacy tracer cannot trace mu = 0.001 (its straight-line chord misses the
  ground sphere), nor h = 100 km at mu = 0.01 (`max_steps`), and an exactly vertical
  direction returns `false`. These are counted and logged, not silently dropped.
* Absorption-enabled runs stay on the legacy **tracer**, and that is the **retained
  production route** — the tau-table route is closed, not deferred. Four table
  representations were studied and rejected; the record, the measured numbers and the
  frozen datasets are archived in `corsika/modules/cherenkov/diagnostics/README.md`.
  The optical depth actually used is the legacy tracer's
  **curved-ray** integral. It was previously computed and then discarded in favour
  of a chord-based approximation (`survivesSpherical`); the ownership of that
  decision is now explicit (`AbsorptionState`, one sample per evaluated photon),
  and the quadrature is sub-stepped to 50 m so `|d tau| <= 1e-3`. This is a physics
  change: absorption-enabled curved outputs differ from earlier revisions. The
  table route is excluded when absorption is on.
* **Absorption accuracy of the retained route is an observation on tested rays, not a
  domain-wide guarantee:** `|d tau| ~ 1e-2` before the sub-stepping fix, `~1e-4` to
  `~1e-5` on the rays tested afterwards.
* **Revision dependency:** the ownership change is only valid together with the
  extinction sign correction that followed it (commit `ec458103e`). The shipped
  profiles store the cumulative column from the ground up, i.e. tau **increasing**
  with height, while the callback had assumed the opposite ordering and returned
  `alpha = 0` at every altitude. Used on its own, the ownership change would
  therefore have disabled absorption for ground+curvature runs. Do not cherry-pick
  the ownership change without the sign fix.
* The flat-Earth `z - R_earth` approximation of 13.1 applies unchanged.

---

## 14. `CherenkovProcess` constructor

The process is instantiated with (see the code comments for the parameter order):

```cpp
corsika::cherenkov::CherenkovProcess cherenkovProcess(
    n                 = 1.0003,               // default (sea-level) phase n
    lambdaMin         = 300.0, lambdaMax = 900.0, // wavelength window [nm]
    bunchSize         = 5.0,                  // photons per bunch
    groundZ_m         = groundZ_m_abs,        // absolute ground z
    recorder          = cherenkovRecorder,    // callback per bunch
    mode              = projMode,             // Ground/ShowerAxis/Telescope
    showerAxis        = std::make_shared<corsika::media::ShowerAxis>(showerAxis),
    showerCoreX/Y/Z_m, telescopeX/Y/Z_m, telescopeRadius_m,
    refractiveIndexFn = refrIndexFunc,
    earthRadiusM      = EarthRadius::Mean/1m,
    telescopeAzimuth, telescopeZenith,
    showerAzimuth(phiRad), showerZenith(thetaRad),
    enableDispersion  = !disable_dispersion,
    enableCurvature   = enable_curvature,
    enableCurvedTable = !disable_curved_table,   // default; see §13.2
    absorptionFile    = --atmosphere-absorption-file,
    telescopePosCorr  = --telescope-position-correction);
```

Responsibilities: produce Cherenkov photons from charged particles above
threshold, propagate them (refraction, group delay/dispersion, optionally
curved paths and absorption), and project the arrival positions onto the chosen
detector geometry.

Note that the constructor receives the **nominal** shower parameters
(`thetaRad`/`phiRad`, `shower_core_*`), i.e. the photon projection is unaffected
by `--injection-*`.


---

## 15. Output manager, stopping plane and physics writers

### 15.1 Output manager

```cpp
std::stringstream args;
for (int i = 0; i < argc; ++i) { args << argv[i] << " "; }
OutputManager output(app["--filename"], seed, args.str(), compressOutput);
```

* `--filename` -> output library directory,
* `seed` -> reproducibility metadata,
* `args.str()` -> a full copy of the command line stored as provenance,
* `compressOutput` -> optional tarball.

### 15.2 The stopping (absorbing) observation plane

In `muon_cherenkov_3d.cpp` the absorbing plane **follows the Cherenkov
projection mode**, so the particle record and the photon record always refer to
the same detector surface. There is no zenith-angle threshold.

| `--cherenkov-projection` | Plane center | Plane normal | In-plane x-axis |
|---|---|---|---|
| `telescope` | `telCenter` | `-telBasisZ` (opposite to the pointing direction) | `telBasisX` |
| `shower` | `showerCore` | `(-nx, -ny, -nz)` (back along the shower) | see below |
| `ground` | `showerCore` | `(0, 0, 1)` | `(1, 0, 0)` |

For the `shower` mode the in-plane x-axis is built exactly like in
`Cherenkov.inl` (`ShowerAxis` mode):

```cpp
if (std::abs(nx) < 0.9) b1 = ( 0.0,  nz, -ny);
else                    b1 = (-nz,  0.0,  nx);
obsPlaneXAxis = b1 / |b1|;
```

> **Note on the normal sign:** `ObservationPlane` treats the *positive*-normal
> side as the "catching" side and, in its leap-frog fallback, writes+removes any
> particle found on the negative side. The normal therefore points **back
> towards the incoming shower** (i.e. opposite to the shower propagation
> direction), which is the convention all three branches above follow.

Each branch logs its choice, e.g.:

```
Stopping plane = TELESCOPE observing plane (zenith ... deg, azimuth ... deg): center=(...), normal=(...)
Stopping plane = SHOWER-axis surface (zenith ... deg)
Stopping plane = GROUND surface (zenith ... deg)
```

```cpp
Plane const obsPlane(obsPlaneCenter, obsPlaneNormal);
ObservationPlane<TrackingType, ParticleWriterParquet> observationLevel{
    obsPlane, obsPlaneXAxis,
    true,        // plane should "absorb" particles
    1e-6 * 1_m,  // padding, unused for absorbing planes
    true};       // forwarded to ParticleWriterParquet(printZ = true)
observationLevel.setGlobalOutput(telCenter);
output.add("particles", observationLevel);
```

**Global particle output.** `setGlobalOutput(telCenter)` switches the recorded
coordinates of `particles.parquet` from plane-local to **global, telescope-
centred**: the `x, y, z` columns then hold the crossing position in the root
coordinate system relative to `telCenter`, and `nx, ny, nz` hold the direction
cosines in that same global system. This makes the ground-plane particles
directly comparable with `particles_3d.parquet` (§18), which uses the same
origin.

> **Caveat:** the trailing `true` is the `printZ` argument of
> `ParticleWriterParquet(bool printZ = true)`, so `particles.parquet` **does**
> contain a `z` column. The inline source comment "do not print z-coordinate" is
> misleading.

### 15.3 Registered writers

| Writer | Output name / directory | Files written | Purpose |
|---|---|---|---|
| `EnergyLossWriter dEdX{showerAxis, dX}` | `"energyloss"` -> `<lib>/energyloss/` | `dEdX.parquet` | continuous energy losses along the axis |
| `LongitudinalWriter profile{showerAxis, dX}` + `LongitudinalProfile longprof` | `"profile"` -> `<lib>/profile/` | `profile.parquet` | longitudinal profile vs. depth |
| `ProductionWriter prod_profile{showerAxis, dX}` + `ProductionProfile prodprof` | `"production_profile"` -> `<lib>/production_profile/` | `profile.parquet` | production-height spectra |
| `ObservationPlane observationLevel` | `"particles"` -> `<lib>/particles/` | `particles.parquet` | particles reaching the ground/shower/telescope plane (global, telescope-centred coordinates) |
| `PrimaryWriter primaryWriter(observationLevel)` | `"primary"` -> `<lib>/primary/` | `summary.yaml` only | primary particle record (YAML summary, no parquet) |
| `InteractionWriter inter_writer(showerAxis, observationLevel)` | `"interactions"` -> `<lib>/interactions/` | `interactions.parquet` | hadronic interaction records |
| `ThreeDPlaneRecorder obs3d` (only if `--observation-steps > 0`) | `"particles_3d"` -> `<lib>/particles_3d/` | `particles_3d.parquet`, `particles_3d_slices.yaml` | 3-D slice crossings (§18) |

See §21 for the full on-disk layout.


---

## 16. Interaction models and process chain

### 16.1 High-energy hadronic model (`--hadronModel`)

Built into `DynamicInteractionProcess<StackType> heModel` (wrapped in
`InteractionCounter heCounted`):

| CLI value | Model |
|---|---|
| `SIBYLL-2.3d` | `corsika::sibyll::Interaction` (default) |
| `QGSJet-II.04` | `qgsjetII::Interaction` |
| `QGSJet-III` | `qgsjetIII::Interaction` |
| `EPOS-LHC-R` | `corsika::EPOS_LHCR::Interaction` |
| `Pythia8` | `corsika::pythia8::Interaction` |

SIBYLL is always allocated separately and reused for PROPOSAL photo-hadronic
interactions.

### 16.2 Decays and neutrino interactions

* `decayTauola` (TAUOLA, left-handed helicity) is selected for `Tau+/-` by
  `make_select(IsTauSwitch(), decayTauola, decayPythia)`; everything else uses
  `decayPythia`.
* The neutrino interaction process is built **per shower**, and only when
  `neutrinoPrimary` (i.e. the primary is a neutrino), as
  `corsika::pythia8::SingleUseNeutrinoInteraction neutrinoPrimaryPythia(trackedParticles, NC, CC)`
  followed by `neutrinoPrimaryPythia.reset()`. The NC/CC toggles come from
  `--neutrino-interaction-type`. `SingleUseNeutrinoInteraction` allows exactly one
  interaction per shower, so tracked **secondary** neutrinos do not trigger a
  second interaction, and it is never reused across events (the underlying
  `corsika::pythia8::NeutrinoInteraction` re-initialises Pythia on every call and
  fails on reuse).

### 16.3 Low-energy / EM

* `sophia` - hadronic photon interactions in the resonance region.
* `emCascade` - PROPOSAL EM cascade; the photon, e+/-, mu+/- and tau+/- production
  thresholds are set to `min(emcut, hadcut, mucut, taucut)` via
  `set_energy_production_threshold`.
* Continuous loss: BetheBlochPDG for hadrons, PROPOSAL `ContinuousProcess`
  otherwise - selected with `make_select(EMHadronSwitch(), ...)`.
* Low-energy hadron model: UrQMD (or FLUKA with `WITH_FLUKA`).

### 16.4 Energy cuts / thinning

```cpp
cut = ParticleCut(emcut, emcut, hadcut, mucut, taucut, !track_neutrinos, dEdX);
EMThinning thinning{emthinfrac * primaryTotalEnergy, maxWeight, !multithin};
StackInspector<StackType> stackInspect(10000, false, primaryTotalEnergy);
```

### 16.5 Ordered process sequence

The sequence is assembled inside the per-shower loop. The neutrino process
(built in §16.2) is prepended **only for neutrino primaries**:

```cpp
auto sequence = make_sequence(
    stackInspect,          // StackInspector
    neutrinoPrimaryPythia, // primary neutrino interactions (only if neutrinoPrimary)
    hadronSequence,        // LE(UrQMD/FLUKA) / HE(heCounted) by energy
    decaySequence,         // TAUOLA/Pythia decays
    emCascade,             // PROPOSAL EM
    prodprof,              // production profile
    emContinuous,          // continuous losses (BetheBloch/PROPOSAL)
    cherenkovProcess,      // <-- Cherenkov photon production
    longprof,              // longitudinal profile
    observationLevel,      // absorbing record plane
    inter_writer,          // record interactions
    thinning,              // EM thinning
    cut);                  // particle cuts
```

With `--observation-steps > 0` the 3-D recorder is inserted **immediately before
the absorbing plane**, so it sees the particles just before they are absorbed:

```cpp
auto sequence = make_sequence(stackInspect, neutrinoPrimaryPythia, hadronSequence,
                              decaySequence, emCascade, prodprof, emContinuous,
                              cherenkovProcess, longprof, *obs3d, observationLevel,
                              inter_writer, thinning, cut);
```

Because the presence of the neutrino process changes the *type* of the sequence,
the code selects between the variants with plain `if (neutrinoPrimary)` /
`if (obs3d)` statements, and the (identical) shower body is run through a
`run_shower(sequence)` lambda (see §17). For non-neutrino primaries the sequence
is built **without** `neutrinoPrimaryPythia`, so tracked secondary neutrinos
never interact.

---

## 17. Per-shower loop and finalization

The per-shower body is factored into a `run_shower` lambda that takes the
process sequence, so the (up to four) sequence variants above share one
implementation:

```cpp
auto run_shower = [&](auto& sequence) {
  TrackingType tracking(app["--max-deflection-angle"]->as<double>());
  StackType stack;
  Cascade EAS(env, tracking, sequence, output, stack);
  ...
};
```

For `i_shower = 1..nevent`:

1. `cherenkovProcess.resetTelescopePositionCorrectionForEvent()` - reset the
   re-centering state for the new event.
2. **Sample the energy** with `PowerLawDistribution(eSlope, eMin, eMax)` from the
   `primary_particle` stream. If `eMin == eMax` (fixed energy) it returns `eMin`
   directly.
   ```cpp
   eKin = primaryTotalEnergy - get_mass(beamCode);
   ```
3. **Thinning weights** from the primary energy; `maxWeight` is the CLI value or
   `0.5 * emthinfrac * E` when `--max-weight` is 0.
4. Fresh `StackInspector{10000, false, primaryTotalEnergy}`.
5. Build the neutrino process (only if `neutrinoPrimary`, §16.2), assemble
   `sequence` (§16.5), and run the body through `run_shower(sequence)`, which
   creates `TrackingType tracking(...)`, `StackType stack`, and
   `Cascade EAS(env, tracking, sequence, output, stack)`.
6. `stack.clear()`.
7. Log the primary parameters (name, PDG/Z-A, energy, momentum, direction,
   injection point, shower-axis length) - see §10.5 for the caveats.
8. Add the primary: `(beamCode, eKin, injDir.normalized(), injectionPosPrimary, 0_ns)`.
9. Optionally `EAS.forceInteraction()` / `EAS.forceDecay()`.
10. `primaryWriter.recordPrimary(primaryProperties)`.
11. `EAS.run()` - run the full shower.
12. **Energy budget** check:
    ```cpp
    Efinal = dEdX.getEnergyLost() + observationLevel.getEnergyGround();
    ```
    log the relative difference vs. `primaryTotalEnergy` (energy-conservation
    diagnostic).
13. **Interaction histograms**, unless `--disable-interaction-histograms`:
    `inthist_lab_<i>.npz` and `inthist_cms_<i>.npz` under
    `<--filename>/interaction_hist/`.
14. Log `Total Cherenkov hits recorded` (cumulative, see §11.3).
15. If `--telescope-position-correction > 0` and the mode is `telescope`, append
    the per-event re-centering result to
    `telescope_position_correction.yaml` (in `cherDir`, append mode, one YAML
    document per event):
    ```yaml
    ---
    event: <i>
    offset_local_m: {x: <m>, y: <m>, z: <m>}
    has_update: true|false
    last_update_pid: <pdg>
    ```
16. Flush the CSV (if applicable).

### Finalization

After the loop:

* `output.endOfLibrary()` closes the output library;
* **eventio:** `eventioWriterPtr->write()` writes `cherenkov_hits.dat`
  (exceptions are caught and logged as warnings);
* **CSV:** `cherenkovCsvPtr->close()`;
* the total wall time is logged as `Run time: {:.3f} s`, and `main` returns
  `EXIT_SUCCESS`.


---

## 18. The 3-D observation-plane recorder (`muon_cherenkov_3d.cpp`)

`muon_cherenkov_3d.cpp` extends the Cherenkov program with an optional
**three-dimensional, time-resolved record of the particle shower**. For every
particle crossing each of a stack of evenly spaced observation planes it stores
the position (plane-local *and* global), the direction, the kinetic energy, the
weight, the arrival time and the slice index - enabling a 3-D "movie" of the
shower as it develops toward the detector.

### 18.1 Enabling and insertion into the sequence

Set `--observation-steps N` with `N > 0`. With the default `N = 0` the recorder
is **not** constructed, no extra file is produced, and the run reduces to the
plain Cherenkov simulation (only the projection-driven stopping plane, the global
particle output and the injection overrides remain as features specific to this
application).

When enabled, the recorder is registered as the output `"particles_3d"` and
inserted into the process sequence immediately **before** the absorbing
detector plane `observationLevel` (§16.5), so it records the particles just
before they are absorbed. It is itself **never** absorbing.

### 18.2 The two classes

Both are defined inside `muon_cherenkov_3d.cpp`, next to the `eventio_writer`
namespace:

| Class | Role |
|---|---|
| `ThreeDParticleWriter : BaseOutput` | Owns the `ParquetStreamer` for `particles_3d.parquet`, the shower id and the **output origin** (`setOutputOrigin(Point)`); writes one row per crossing (`write(...)`) and emits the side table `particles_3d_slices.yaml` (via `startOfLibrary`). |
| `ThreeDPlaneRecorder : ContinuousProcess<ThreeDPlaneRecorder>, ThreeDParticleWriter` | Owns the `N` planes (centers, normals, in-plane axes, slice altitudes); `getMaxStepLength` limits the tracking step to the nearest sampled plane; `doContinuous` records every plane crossed within the step (sign-change detection). |

Constructor:

```cpp
ThreeDPlaneRecorder(std::vector<Point> const& centers,
                    std::vector<DirectionVector> const& normals,
                    std::vector<DirectionVector> const& xaxes,
                    std::vector<DirectionVector> const& yaxes,
                    std::vector<double> const& hSlice,        // slice altitudes [m]
                    Point const& earthCenter, LengthType const rEarth,
                    double const eminGeV, std::set<Code> const& filter,
                    Point const& outputOrigin);
```

`outputOrigin` is passed as `telCenter`, so the global columns of
`particles_3d.parquet` use the **same telescope-centred origin** as
`particles.parquet` (§15.2).

### 18.3 Slice geometry: a uniform stack parallel to the stopping plane

The slice stack is built in `main` from the *stopping plane* that was selected
from `--cherenkov-projection` (so the 3-D planes are parallel to the same
detector surface the photons and ground particles refer to):

```cpp
sInj_signed_m = (injectionPosPrimary - obsPlaneCenter) . obsPlaneNormal / 1_m;
sInj_m        = |sInj_signed_m|;
Delta_m       = sInj_m / N;                       // spacing along the normal
stackSign     = (sInj_signed_m >= 0.0) ? +1.0 : -1.0;

sliceNormal = obsPlaneNormal.normalized();
sliceX      = obsPlaneXAxis.normalized();
sliceY      = sliceNormal.cross(sliceX);

for (k = 0 .. N-1) {
  s_m  = stackSign * k * Delta_m;
  ck   = obsPlaneCenter + sliceNormal * (s_m * 1_m);   // slice centre
  hSlice[k] = ((ck - EarthCentre).getNorm() - EarthRadius::Mean) / 1_m;
}
```

Properties of this construction:

* **Slice 0 lies exactly on the stopping surface** (`s_m = 0`), i.e. on the
  absorbing `observationLevel` plane, and slice `N-1` is one spacing short of the
  injection point, so **no slice coincides with the primary's starting point**.
* The stack marches **towards the injection point** (`stackSign` from the sign of
  the signed injection distance), which is the upstream side for both the
  `telescope`/`shower` surfaces (normal points back toward the shower) and the
  `ground` surface.
* All planes share one normal (`sliceNormal`) and one in-plane x-axis
  (`sliceX`, which is the stopping plane's own x-axis); the y-axis is
  `sliceNormal x sliceX`, i.e. the same convention as `ObservationPlane`.
* `hSlice[k]` is the slice-centre altitude above sea level (`R_E` is
  `constants::EarthRadius::Mean`, `EarthCentre` is the root origin).

> **Caveat (coincident planes):** because slice 0 sits exactly on the absorbing
> plane, the 3-D recorder and `ObservationPlane` are stacked with zero gap. The
> framework documents (issue #397) that two `ObservationPlane`s should not be
> placed exactly on top of each other; the recorder is a different class and does
> not absorb, but if you see duplicated/odd records at slice 0 this is the reason
> to look at.

> **Abort condition:** if the primary injection point lies *in* the stopping
> plane (`sInj_m == 0`, e.g. `--telescope-zenith 90` with `--zenith` unset), the
> program logs an error and exits with `EXIT_FAILURE`, because no stack can be
> built. Give a non-zero `--injection-x/-y`, or set `--zenith` to the shower
> direction.

On success the configuration is logged:

```
3-D observation recorder active: <N> slices parallel to the stopping plane
(spacing=<Delta> m, altitude span=<h_N-1 - h_0> m)
```


### 18.4 Crossing detection and step limiting (O(1) fast path)

Because the slices form a **uniform stack of parallel planes**, the whole stack
is described by one normal and one signed spacing, which the constructor
derives once and verifies:

```cpp
sliceNormal_ = normals_[0].normalized();
dStep_m_     = sliceNormal_ . (centers_[1] - centers_[0]) / 1_m;   // signed
analyticOk_  = (dStep_m_ != 0.0);
for (k = 1 .. N-1) {                       // verify uniformity
  dk = sliceNormal_ . (centers_[k] - centers_[0]) / 1_m;
  nk = sliceNormal_ . normals_[k];         // = 1 for parallel planes
  if (|dk - k*dStep_m_| > 1e-6 || |nk - 1.0| > 1e-9) analyticOk_ = false;
}
```

With the signed distances `d_k(p) = n . (p - c_k)` this gives the exact relation
`d_k(p) = d_0(p) - k*dStep`, which is what makes the fast path possible.

**`getMaxStepLength`** (called by the tracker) returns the length of the step to
the *nearest* sampled plane, or `inf` when there is nothing to hit:

* `inf` if the trajectory duration is not positive;
* analytic path:
  * `q = d_0 / dStep_m_` is the **fractional slice index** of the particle, and
    `qdot = (n . v)/dStep_m_` gives the direction of travel;
  * `inf` if the particle moves parallel to the planes (`n . v == 0`) or if the
    estimated next slice index lies far outside the stack
    (`kCenter < -3` or `> N+3`);
  * otherwise a small window of candidate planes
    (`int(kCenter) + dk`, `dk = -2..2`) is intersected with the trajectory and the
    smallest admissible length is returned. The window (rather than a single
    plane) makes this robust against the straight-line index estimate being off
    by one; and since `doContinuous` records on the *sign change* of the signed
    distance, an imperfect estimate here can never lose or duplicate a record.
  * crossings with `tEntry <= epsS_` (`epsS_ = 1e-12 s`) or `tEntry > duration`
    are ignored;
* **fallback**: if `analyticOk_` is false (single-slice or non-uniform stacks), the
  original per-plane loop over all `N` planes is used.

**`doContinuous`** records every plane crossed during the step:

```cpp
if (ekinGeV < eminGeV_) return ProcessReturn::Ok;            // energy filter
if (!filter_.empty() && filter_.count(pid) == 0) return Ok;  // PDG filter

// analytic path: restrict the loop to the contiguous slice range that can
// possibly be crossed, in O(1) instead of scanning all N planes
lo = min(dPre0, dPost0);  hi = max(dPre0, dPost0);           // d_0 at pre/post
ka = (dStep > 0) ? lo/dStep : hi/dStep;
kb = (dStep > 0) ? hi/dStep : lo/dStep;
kLo = max(ceil(ka), 0);  kHi = min(floor(kb), N-1);
if (kHi < 0 || kLo > N-1) return Ok;                         // crosses no slice

for (k = kLo .. kHi) {
  dPre  = normals_[k] . (pPre  - centers_[k]);
  dPost = normals_[k] . (pPost - centers_[k]);
  crossed = (dPre > 0 && dPost <= 0) || (dPre < 0 && dPost >= 0);
  if (!crossed) continue;
  frac = dPre / (dPre - dPost);                 // linear interpolation in the step
  pos  = pPre + (pPost - pPre) * frac;
  t    = tPre + (tPost - tPre) * frac;
  xl   = xAxes_[k] . (pos - centers_[k]);       // plane-local coordinates
  yl   = yAxes_[k] . (pos - centers_[k]);
  h    = (pos - earthCenter_).getNorm() - rEarth_;
  write(pid, ekinPost, xl, yl, pos.x, pos.y, pos.z,
        nx, ny, nz, t, weight, h, k);
}
return ProcessReturn::Ok;                       // never absorbs
```

Summary of the semantics:

* a particle is recorded **at most once per plane per step**, at the
  interpolated crossing point (linear interpolation of position and time within
  the step);
* the recorded direction is the **post-step** direction and the recorded kinetic
  energy the **post-step** kinetic energy, while the weight is taken from the
  pre-step particle;
* the recorder always returns `ProcessReturn::Ok`, i.e. it never absorbs, kills
  or modifies a particle;
* `epsS_ = 1e-12 s` protects against zero-length/duplicate crossings at the very
  beginning of a step.

### 18.5 Output: `<lib>/particles_3d/particles_3d.parquet`

Column order (the leading `shower` column is added automatically by
`ParquetStreamer::initStreamer`, all others by `ThreeDParticleWriter`):

| Column | Parquet type | Meaning |
|---|---|---|
| `shower` | int32 (UINT_32) | shower id |
| `slice` | int32 | slice index `k` (`0 .. N-1`) |
| `pdg` | int32 | particle PDG code |
| `kinetic_energy` | float | kinetic energy (GeV) |
| `x`, `y` | float | plane-local lateral offsets from the slice centre (m) |
| `X`, `Y`, `Z` | double | position in global coordinates **relative to the output origin** (= `telCenter`) (m) |
| `nx`, `ny`, `nz` | float | global direction cosines |
| `time` | double | arrival time (s), measured from primary injection at `t = 0` |
| `weight` | float | particle weight |
| `h` | float | altitude above sea level at the crossing, `|pos - EarthCentre| - R_E` (m) |


### 18.6 Output: `<lib>/particles_3d/particles_3d_slices.yaml`

A companion side table, written by `ThreeDPlaneRecorder::startOfLibrary` into the
same directory as the parquet file, fully describes the slice stack (all values
in metres, printed with `setprecision(9)` as plain doubles):

```yaml
earth_radius_m: 6371000
output_origin_m: [x, y, z]          # = telCenter, origin of the X/Y/Z columns
slices:
  - {slice: 0, h_m: <m>, center_m: [x, y, z], normal: [nx, ny, nz],
     x_axis: [xx, xy, xz], y_axis: [yx, yy, yz]}
  - {slice: 1, ...}
  ...
```

It is intended for validation (check the `h` column of the parquet against
`h_m`) and for referencing each slice's local origin and axes.

### 18.7 Filters

| Option | Effect |
|---|---|
| `--obs3d-emin GEV` | skip particles whose **post-step** kinetic energy is below `GEV` (`0` = no filter) |
| `--obs3d-particles LIST` | record only particles whose PDG code is in the comma-separated `LIST`; whitespace around the tokens is trimmed; the preset `muons` expands to `{mu+, mu-}` (`Code::MuMinus`, `Code::MuPlus`); empty string = no particle filter |

`--obs3d-particles` is parsed only when `--observation-steps > 0`; tokens are
converted with `convert_from_PDG(PDGCode(value))`, and a token that cannot be
parsed as an integer makes the program log an error and exit with
`EXIT_FAILURE`.

### 18.8 Validation recipe and caveats

* **Slice altitudes:** compare the `h` column of `particles_3d.parquet` with
  `h_m` of the corresponding slice in `particles_3d_slices.yaml`.
* **Slice 0 vs. the ground plane:** slice 0 coincides with the absorbing plane,
  so for `ground` projection the slice-0 rows of `particles_3d.parquet` should
  reproduce the corresponding rows of `particles.parquet` (same telescope-centred
  origin, same direction cosines), apart from the energy/PDG filters and the fact
  that the recorder sees the particle *before* it is absorbed/padded past the
  plane.
* **Geometry:** `x_axis`, `y_axis`, `normal` and `center_m` of each slice let you
  reconstruct the local frames, and `output_origin_m` plus `earth_radius_m`
  convert the `X/Y/Z` and `h` columns back to absolute positions.
* **Caveat (altitude vs. slice centre):** `h_m` is the altitude of the slice
  *centre*. For the `shower` and `telescope` projections the planes are not
  geopotential surfaces, so the crossing altitude `h` varies across the plane and
  will differ from `h_m` away from the centre. Only for `ground` projection
  (horizontal planes) is `h` essentially constant across a slice.
* **Caveat (single file, all showers):** one `particles_3d.parquet` and one
  `particles_3d_slices.yaml` are written per output library; use the `shower`
  column to separate events (`startOfShower` stores the id, `startOfLibrary`
  resets it to 0).
* **Caveat (filters are applied before the recorder's geometry):** the energy and
  PDG filters are evaluated per step, using the post-step kinetic energy.

### 18.9 Example

```bash
./muon_cherenkov_3d -p 13 -E 100 -z 0 -a 0 -f my3d \
    --observation-steps 10 --obs3d-particles muons
```

records muons crossing 10 slices and produces
`my3d/particles_3d/particles_3d.parquet` plus
`my3d/particles_3d/particles_3d_slices.yaml` (and, as always, the Cherenkov file
`cherenkov_hits.csv` next to the library directory, i.e. in `.`). With
`--observation-steps 0` (or omitted) no 3-D output is produced.


---

## 19. The `eventio_writer` namespace (CHASM-compatible binary writer)

This namespace implements a from-scratch writer for the **eventio** container
format used by IACT tools. Eventio is a block-based format with specific byte
packing - all integers are little-endian.

### 19.1 Constants and type codes

```cpp
static const std::array<uint8_t, 4> SYNC_MARKER = {0x37, 0x8a, 0x1f, 0xd4};
enum class IACTType : uint32_t {
  RunHeader=1200, InputCard=1212, AtmosphericProfile=1216,
  TelescopeDefinition=1201, EventHeader=1202, ArrayOffsets=1203,
  Longitudinal=1211, TelescopeData=1204, Photons=1205, CameraLayout=1206,
  TriggerTime=1207, PhotoElectrons=1208, EventEnd=1209, RunEnd=1210 };
```

The `0x37 0x8a 0x1f 0xd4` pattern is the eventio top-level **sync marker**
preceding every top-level object.

### 19.2 Low-level byte writers

* `write_int32` / `write_int16` / `write_float32` / `write_bytes` - raw
  little-endian streaming.
* `str_to_float(const char* s)` - copies 4 bytes of a literal (e.g. `"RUNH"`)
  into a `float` via `memcpy`, so 4-character "words" written in C7-style float
  arrays match the expected eventio interpretation.

### 19.3 Top-level object header layout

```
4 bytes : SYNC_MARKER
4 bytes : type_word  = type | (user<<16) | (extended<<17) | (version<<20)
4 bytes : id          (int32)
4 bytes : length_word = content_size_lo(30 bits) | (only_subobjects<<30)
[4 bytes]: extension_word = content_size >> 30      (only if extended)
```

`write_toplevel_header` decides `extended` when the content size exceeds
`0x3FFFFFFF` (~ 1 GiB), setting bit 17 and appending the upper-size word.
Header = 16 bytes normally, 20 if extended.

### 19.4 Sub-object header (inside a top-level object)

Same but **without** the sync marker:

```
4 bytes : type_word
4 bytes : id
4 bytes : length_word
[4 bytes]: extension_word (only if extended)
```

Total = 12 bytes normally, 16 if extended (`write_subobject_header`).

### 19.5 `PhotonBunch` and writer

Eventio "bunches" carry 8 floats (CHASM layout) written by
`write_photon_bunch`:

| Field | Meaning | Stored as |
|---|---|---|
| `x`, `y` | impact position on the detector | cm (already x100) |
| `cx`, `cy` | direction cosines (x, y) | unit-less |
| `time` | arrival time | ns |
| `zem` | emission height | cm |
| `n_photons` | number of photons in this bunch (= `weight`) | count |
| `wavelength` | photon wavelength | nm |


### 19.6 Block classes

**`RunHeaderBlock`** (type 1200, id 102)
* Content size **1096** bytes (matches C7/reference).
* Payload: `count=273`, then the `"RUNH"` word, the run number, a date
  placeholder, and 270 zero floats (words 4..273).

**`InputCardBlock`** (type 1212, id 0)
* Content = just `n_strings=0` (4 bytes). Readers only need `n_strings` to
  parse; the full 2748-byte reference card content is not reproduced.

**`EventHeaderBlock`** (type 1202, id = event number)
* Content size **1096**, `count=273`.
* Words: `"EVTH"`, event number, particle ID (`6.0` = muon), total energy in GeV,
  two placeholder floats, first-interaction height, momentum x/y/z, zenith (rad),
  azimuth (rad), then zeros (words 13..273).

**`TelescopeDefBlock`** (type 1201, id 0)
* `content_size = 4 + ntel*16`.
* Payload: `ntel`, then per telescope `(x, y, z, radius)` as 4 floats.
* Written with the single configured telescope position/radius.

**`CherenkovEvent`**
* Per-event container: the bunch vector plus all header metadata (zenith,
  azimuth, primary energy/pdg, telescope geometry, event/run numbers).

**`CherenkovEventioWriter`**
* `add_photon(...)` appends one `PhotonBunch`.
* `set_event_params(...)` records shower/telescope parameters.
* `write()` streams the whole run in this order:

### 19.7 `write()` - block order

1. **RunHeader** (C7 format)
2. **InputCard** (`n_strings = 0`)
3. **TelescopeDefinition** (CHASM-compatible)
4. **EventHeader** (shower parameters)
5. **ArrayOffsets** (type 1203): `n_arrays=1`, `time_offset=0`, per-array
   `(x=0, y=0)` offsets - content 16 bytes.
6. **TelescopeData** (type 1204):
   * If bunches exist: `only_subobjects=true` header whose length is the Photons
     sub-object (header 12/16 bytes depending on extension) plus the Photons
     content. Inside, a **Photons** (type 1205) sub-object with
     `arr(int16)=0`, `tel(int16)=0`, `n_photons(float)=sum of bunch.n_photons`,
     `n_bunches(int32)`, then every `PhotonBunch` (8 floats each).
   * If no bunches: an empty TelescopeData (`content_size=0`,
     `only_subobjects=true`).
7. **EventEnd** (type 1209, id = event number): `count=273`, `"EVTE"`, event
   number, zeros.
8. **RunEnd** (type 1210, id 102): `content_size=16`, `count=3`, `"RUNE"`,
   run number, `n_events=1`.

The result is a single-run, single-event CHASM-style file designed for ingestion
by tools such as GrOptics.

### 19.8 Caveats of the eventio path

* **All showers are merged into one event record.** `write()` is called once
  after the shower loop and streams whatever `CherenkovEvent` accumulated, so
  with `-N > 1` every shower's bunches end up in a single `EventHeader` /
  `TelescopeData` / `EventEnd` triple (`event_no` stays 1, `RunEnd` reports
  `n_events=1`).
* **The particle-ID word is hard-coded** to `6.0` (muon) in `EventHeaderBlock`,
  and the `pdg` argument of `set_event_params` is stored in `CherenkovEvent` but
  never written.
* **The header energy is nominal** (see §11.2): with `--energy_range` the written
  energy is the lower bound of the range, not the sampled event energy.
* The file is written **after** `output.endOfLibrary()`; a failure to open/write
  it is caught, logged as a warning, and does not abort the run.


---

## 20. Differences of the sibling applications

All three applications share:

* the environment / atmosphere setup and the magnetic field,
* the Cherenkov subsystem: CSV/eventio output, the `HitRecorder` lambdas, the
  atmosphere-table spline and `refrIndexFunc`, dispersion / curvature /
  absorption, and the telescope position correction,
* the full interaction, decay and thinning process chain,
* the physics writers (`energyloss`, `profile`, `production_profile`,
  `particles`, `primary`, `interactions`).

They differ in the following points:

| Aspect | `muon_cherenkov_3d.cpp` | `muon_cherenkov_test.cpp` | `muon_cherenkov.cpp` |
|---|---|---|---|
| Stopping-plane selection | follows `--cherenkov-projection` (no zenith threshold) | `zenith > 80 deg`: vertical plane through the telescope; otherwise horizontal plane at the core | `zenith > 85 deg`: vertical plane; otherwise horizontal |
| Plane center (inclined case) | `telCenter` in `telescope` mode, `showerCore` in `shower`/`ground` mode | telescope position | telescope position |
| In-plane x-axis (inclined case) | `telBasisX` (`telescope`) or the shower-axis construction (`shower`) | `(-ny, nx, 0)/sqrt(nx^2+ny^2)`, horizontal and perpendicular to the shower projection | analogous |
| Ground particles | global, telescope-centred (`setGlobalOutput`) | plane-local | plane-local |
| `--injection-x/-y/-zenith/-azimuth` | yes | no | no |
| `injection_geometry.yaml` | yes | no | no |
| 3-D recorder (`--observation-steps`, `--obs3d-*`) | yes | no | no |
| `particles_3d.parquet` / `particles_3d_slices.yaml` | yes | no | no |
| `--zenith` accepted range | 0-100 deg | 0-92 deg | 0-92 deg |
| Per-shower body | factored into a `run_shower` lambda; sequence chosen with `if (neutrinoPrimary)` / `if (obs3d)` | factored into a `run_shower` lambda; sequence chosen with `if (neutrinoPrimary)` | factored into a `run_shower` lambda; sequence chosen with `if (neutrinoPrimary)` |
| Neutrino process | per-shower `SingleUseNeutrinoInteraction`, only for neutrino primaries | same | same |
| Source size | 115 042 bytes | 75 057 bytes | 73 546 bytes |

> The old geometry description ("altitude shells below 80 deg, axis-perpendicular
> planes above, 5 mm slice-0 offset") belongs to a **previous version** of the
> 3-D recorder. The current implementation always builds a uniform stack
> parallel to the stopping plane (§18.3) and has no 5 mm offset.

---

## 21. Output file layout

`OutputManager::add(name, output)` creates a sub-directory `<lib>/<name>/` for
every registered output, and `startOfLibrary` is called with that directory, so
the parquet files live one level below `--filename`:

```
<lib>/                              # = -f / --filename
├── config.yaml                     # top-level config (command line, version)
├── summary.yaml                    # written in the OutputManager destructor
├── energyloss/dEdX.parquet
├── profile/profile.parquet
├── production_profile/profile.parquet
├── particles/particles.parquet     # global, telescope-centred coordinates
├── primary/summary.yaml            # PrimaryWriter stores a YAML summary, no parquet
├── interactions/interactions.parquet
├── particles_3d/                   # only with --observation-steps > 0
│   ├── particles_3d.parquet
│   └── particles_3d_slices.yaml
├── interaction_hist/               # unless --disable-interaction-histograms
│   ├── inthist_lab_<i>.npz
│   └── inthist_cms_<i>.npz
└── <name>/config.yaml, <name>/summary.yaml   # per-output metadata
```

Files written **next to** the library directory (into `dirname(-f)`, not into
`<lib>`), because they are not registered outputs:

| File | When |
|---|---|
| `cherenkov_hits.csv` (default) or `cherenkov_hits.dat` (eventio) | always |
| `injection_geometry.yaml` | only when `--injection-x/-y/-zenith/-azimuth` is used |
| `telescope_position_correction.yaml` | only when `--telescope-position-correction > 0` and the mode is `telescope` |


---

## 22. Practical example invocations

Simulate one 100 GeV muon, vertical (`-z 0`), record Cherenkov photons in CSV
with a 1 m telescope aperture:

```bash
./muon_cherenkov_3d -p 13 -E 100 -z 0 -a 0 -f mylib \
    --cherenkov-projection telescope --telescope-radius 1 \
    --cherenkov-output-format csv
```

Horizontal muons with the telescope-plane observation surface and eventio
output (CHASM-compatible), fixed seed, five events:

```bash
./muon_cherenkov_3d -p 13 -E 112 -z 89.9 -a 0 -f horiz \
    --cherenkov-output-format eventio -s 42 -N 5 \
    --cherenkov-projection telescope --telescope-radius 1
```

3-D record of the muon component through 10 slices:

```bash
./muon_cherenkov_3d -p 13 -E 100 -z 0 -a 0 -f my3d \
    --observation-steps 10 --obs3d-particles muons
```

3-D record with an energy threshold, plus an off-axis primary injection (the
detector geometry stays on the nominal axis; the primary is displaced by
(+200 m, +100 m) in global x/y and injected 2 deg off-axis):

```bash
./muon_cherenkov_3d -p 13 -E 100 -z 20 -a 0 -f offaxis \
    --observation-steps 12 --obs3d-emin 0.1 \
    --injection-x 200 --injection-y 100 --injection-zenith 22 \
    --cherenkov-projection telescope --telescope-radius 50 \
    --telescope-azimuth 0 --telescope-zenith 20
```

Notes:

* `-f` sets the base output path; the library is created there (it must **not**
  exist yet), and the Cherenkov files / YAML provenance files go to its
  directory part.
* `-s 0` (default) yields a random, nondeterministic seed.
* Without `--injection-*` the file `injection_geometry.yaml` is not written and
  the primary is injected on the nominal shower axis.

---

## 23. Key physics parameters to tune

| Parameter | Where to set | Effect |
|---|---|---|
| Cherenkov wavelength window | hard-coded 300-900 nm in the process constructor | which photons are produced |
| Photons per bunch | `bunchSize = 5.0` (hard-coded) | output granularity / weight distribution |
| Dispersion on/off | `--disable-dispersion` | Cherenkov ring timing & angle vs. lambda |
| Curvature correction | `--enable-curvature` | needed at zenith > 60 deg |
| Curved-ray table | `--disable-cherenkov-curved-table` | on by default for eligible curved-ground rays; 0.449 us/photon vs 16.7-672 us for the legacy trace (13.2) |
| Refractive-index model | `--atmosphere-file` | drives the Cherenkov angle / timing |
| Absorption | `--atmosphere-absorption-file` | Monte-Carlo attenuation of photons |
| Projection | `--cherenkov-projection` | ground / shower / telescope surfaces |
| Telescope | `--telescope-x/-y/-z/-radius/-azimuth/-zenith` | detector aperture & pointing |
| Telescope re-centering | `--telescope-position-correction` | fraction of steps used to follow the brightest charged particle (`telescope` mode only) |
| Energy cuts | `--emcut`/`--hadcut`/`--mucut`/`--taucut` | tracking threshold, CPU cost |
| Hadronic model | `-M` + `-T` | high-energy physics + LE/HE switch energy |
| Observation geometry | `--observation-level`, `--injection-height` | detector altitude and primary start altitude |
| 3-D slices | `--observation-steps` | number of parallel observation planes (CPU and output size) |
| 3-D energy filter | `--obs3d-emin` | min kinetic energy of recorded particles |
| 3-D particle filter | `--obs3d-particles` | restrict to PDG codes (or the `muons` preset) |
| Primary injection | `--injection-x/-y/-zenith/-azimuth` | displace/repoint the primary without moving the detector |

---

## 24. Consolidated caveats

* **Eventio:** all showers are merged into a single event record; the particle-ID
  word is hard-coded to muon (6); the header energy/zenith/azimuth are the
  *nominal* ones, not the per-event sampled values (§11.2, §19.8).
* **CSV:** one file for all showers, no event column (§11.3).
* **`cherenkovHitCount`** is cumulative over showers, so the per-shower log line
  is a running total (§11.3).
* **Slice 0 of the 3-D recorder coincides with the absorbing plane** (§18.3);
  the geometry is a uniform parallel stack, not altitude shells, and there is no
  5 mm offset.
* **`h_m` in `particles_3d_slices.yaml` is the slice-centre altitude**; for
  non-horizontal slices the crossing altitude `h` varies across the plane (§18.8).
* **`Primary Direction:` log line prints the norm** of the direction vector, not
  its components (§10.5).
* **The `--injection-*` overrides do not move the detector geometry**: the
  shower axis, stopping plane and telescope remain nominal (§10).
* **`--injection-x/-y` use the radial (polar) `--injection-height`**: with a
  non-zero offset the injection `z` is derived so the point sits on the sphere
  `R + injection_height` (see §10.2). Before this change the offset preserved the
  nominal `z`, giving a slightly wrong altitude for off-axis injections.
* **`injection_geometry.yaml` is written once per run**, before the shower loop
  (§10.5).
* **`particles.parquet` contains a `z` column** even though the source comment
  says otherwise (§15.2).
* **`--atmosphere-file` defaults to an absolute, machine-specific path** (§6.8).
* **`-f/--filename` must not already exist** (`CLI::NonexistentPath`), and
  `OutputManager` throws if the output directory is present.
* **Radio observers** are constructed but not registered, so no radio output is
  produced (§6.6).
* **The neutrino interaction process is added to the sequence only for neutrino
  primaries** (`is_neutrino(beamCode)`); it is a per-shower
  `SingleUseNeutrinoInteraction`, so tracked secondary neutrinos do not interact
  and the process is never reused across events (§16.2, §16.5).
