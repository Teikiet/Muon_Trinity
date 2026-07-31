# Muon Flux Analysis for the Muon Trinity Demonstrator

This repository contains a full end-to-end analysis pipeline for estimating the atmospheric muon flux detected by the **Muon Trinity Demonstrator** telescope. The analysis proceeds from raw Geant4/CORSIKA simulation output through trigger-efficiency modelling, effective-area estimation, and final muon rate prediction.

---

## Table of Contents

1. [Overview](#overview)
2. [Dependencies](#dependencies)
3. [Configuration](#configuration)
4. [Pipeline Steps](#pipeline-steps)
   - [1. Data Loading & Preprocessing](#1-data-loading--preprocessing)
   - [2. Summary Statistics & Coverage Bands](#2-summary-statistics--coverage-bands)
   - [3. Empirical 2D Trigger-Efficiency Map](#3-empirical-2d-trigger-efficiency-map)
   - [4. Per-Incidence-Bin Logistic Fits](#4-per-incidence-bin-logistic-fits)
   - [5. Global Models of Logistic Parameters](#5-global-models-of-logistic-parameters)
   - [6. Cross-Group Parameter Correlations](#6-cross-group-parameter-correlations)
   - [7. Combined Statistical-Systematic Uncertainty](#7-combined-statistical-systematic-uncertainty)
   - [8. Efficiency Grid on (Zenith, Azimuth, Height, Energy)](#8-efficiency-grid-on-zenith-azimuth-height-energy)
   - [9. Effective Area from Digitized Trigger-Rate Data](#9-effective-area-from-digitized-trigger-rate-data)
   - [10. Muon Flux Calculation](#10-muon-flux-calculation)
5. [Key Results](#key-results)
6. [Outputs](#outputs)

---

## Overview

The goal is to compute the expected muon trigger rate at the Muon Trinity Demonstrator site (2944 m altitude, Dugway Proving Ground, Utah). The pipeline:

1. Loads thousands of simulated muon showers (CORSIKA8 + Trinity Simulation Chain) organised by energy, seed, and telescope position.
2. Builds an empirical trigger-efficiency model $\epsilon(N_\gamma, \theta_{\text{inc}})$ — the probability that a shower with `photon_count` $N_\gamma$ and incidence angle $\theta_{\text{inc}}$ deposits $\geq 20$ photo-electrons (PE) in any pixel.
3. Fits a logistic (sigmoid) function per incidence bin and then models the logistic parameters ($b$, $x_0$, $w$) as smooth functions of incidence angle.
4. Combines the model with a pre-computed **MCEq** atmospheric muon flux grid (zenith 85°–90°) to obtain the detected muon rate by integrating $\Phi(E, \theta, \phi) \times \epsilon(E, H, \theta, \phi)$ over energy, emission height, zenith, and azimuth.
5. Propagates uncertainties through bootstrapping, binning-systematic variations, and analytic Jacobian error propagation.

---

## Dependencies

- Python ≥ 3.8
- `numpy`, `scipy`, `pandas`, `matplotlib`
- `sympy` (symbolic Jacobian)
- `joblib` (parallel computation)
- `pickle` (load pre-computed MCEq grid)

---

## Configuration

Key configurable parameters at the top of the script:

| Parameter        | Value    | Description                                              |
|------------------|----------|----------------------------------------------------------|
| `PE_THRESHOLD`   | 20       | Minimum photo-electrons in a pixel to count as triggered |
| `E_MIN_GEV`      | 1 × 10¹  | Minimum muon energy (GeV)                                |
| `E_MAX_GEV`      | 1 × 10⁶  | Maximum muon energy (GeV)                                |
| `SEEDS`          | [1,2,3,4]| Random seeds used per energy point                       |
| `PID`            | 13       | Particle ID (muon)                                       |
| `RADIUS`         | 5        | Telescope radial offset (m)                              |
| `BASE_DIR`       | (scratch)| Path to simulation output directory                      |

---

## Pipeline Steps

### 1. Data Loading & Preprocessing

**Function**: `load_scan()`, `get_available_energies()`

- Scans `BASE_DIR` for directories matching the pattern `Muon_pid13_E{energy}_R5`.
- For each energy and each of the 4 random seeds, reads the CSV output of the simulation scan.
- Drops rows where `file_found == 0` (missing simulation files).
- Computes derived quantities:
  - **Radial distance**: $r = \sqrt{tel\_x^2 + tel\_z^2}$
  - **Detection flag**: `detected = (max_pe >= PE_THRESHOLD)`
  - **Slant depth** through the atmosphere (km) — geometric line-of-sight from emission point to observer.
  - **Incidence angle** relative to the telescope pointing direction $\hat{n}_{\text{tel}}$:

  \[
  $\cos\theta_{\text{inc}} = \sqrt{1 - \cos^2\theta_z - \cos^2(\phi - 180^\circ)}$
  \]

  where $\theta_z$ is the shower zenith and $\phi$ is the shower azimuth.

### 2. Summary Statistics & Coverage Bands

For each simulated energy bin, the code computes the mean and percentile coverage bands (50 %, 80 %, 95 %, 100 %) of:

- **Photon count** (`cph_photon_count`)  
- **Max PE per pixel**  
- **File size** (MB)  
- **Runtime** (hours converted from seconds)

Plots are generated with `plot_with_bands()`, showing log-log curves of each quantity vs. muon energy with nested shaded percentile bands.

Total storage consumption and the photon-per-MB ratio are also reported.

### 3. Empirical 2D Trigger-Efficiency Map

**Location**: Section "Empirical trigger efficiency: P(Max_PE >= threshold | photon_count, incidence)"

**Scope**: Only showers with `R == 0` (telescope at origin) are used.

- **Inputs**: `photon_count` $(N_\gamma)$, `Incidence` $(\theta_\text{inc})$, `Max_PE`
- **Bins**: 50 × 50 histogram in ($\log_{10} N_\gamma$,$\theta_\text{inc}$).
- **Efficiency**: $\epsilon(N_\gamma,\theta_\text{inc}) = N_\text{hit} / N_\text{total}$ for bins with ≥5 events.
- A light Gaussian smooth (σ = 0.8) is applied to fill sparse bins.
- Saved as an interpolator (`trigger_eff()`) using `RegularGridInterpolator` with bilinear interpolation; values are clamped to [0,1].
- **Validation**: The efficiency is validated by slicing the incidence range into 7 sub-bands and overlaying empirical points with the interpolator prediction.

### 4. Per-Incidence-Bin Logistic Fits

**Model**: For each incidence bin, trigger efficiency vs. $\log_{10} N_\gamma$ is modelled as a logistic (sigmoid):

\[
$\epsilon(\log_{10} N_\gamma) = \frac{b}{1 + \exp\left(-\frac{\log_{10} N_\gamma - x_0}{w}\right)}$
\]

where:
- $b$ = saturation efficiency (upper asymptote)
- $x_0$ = $\log_{10} N_\gamma$ at which $\epsilon = b/2$
- $w$ = width (steepness) of the transition

**Fitting method**: Maximum-likelihood with binomial likelihood (NLL minimisation via L-BFGS-B). Clopper-Pearson 68 % confidence intervals are used for data uncertainties.

**Binning systematics**: The analysis is repeated `N_BIN_VARIATIONS = 20` times with randomly chosen numbers of photon-count bins (20–200) and incidence bins (5–100). This captures the sensitivity of fitted parameters to the binning choice.

Quality cuts for a fitted bin:
- Chi²/dof < 10
- Relative error on $x_0$ < 1.0 and on $w$ < 1.0
- At least one binning variation yields a valid fit per incidence bin.

### 5. Global Models of Logistic Parameters

The incidence-dependent logistic parameters $(b(\theta), x_0(\theta), w(\theta))$ are fit globally across all incidence bins with quality-passed fits:

| Parameter | Model Selected | Equation |
|-----------|---------------|----------|
| $b(\theta)$ | `sigmoid_fixed` | $b(\theta) = B_\infty + \displaystyle\frac{B_0 - B_\infty}{1 + e^{k(\theta - \theta_m)}}$<br>with $B_0=1.0$, $B_\infty=0.0$ |
| $x_0(\theta)$ | `log_sigmoid` | $x_0(\theta) = x_\min + \displaystyle\frac{x_\max - x_\min}{1 + e^{-k(\theta - \theta_m)}}$ |
| $w(\theta)$ | `sigmoid` | $w(\theta) = w_\infty + \displaystyle\frac{w_0 - w_\infty}{1 + e^{k(\theta - \theta_m)}}$ |

The choice of `sigmoid_fixed` for $b$ forces the saturation efficiency to drop from 1.0 to 0.0 across a sharp transition at some midpoint angle.

### 6. Cross-Group Parameter Correlations

For each quality incidence bin, the ensemble of $(b, x_0, w)$ triplets from the binning variations provides an empirical 3×3 covariance. These are pooled across bins to obtain a representative correlation matrix $R_{b,x_0,w}$, used to populate off-diagonal blocks of the joint parameter covariance matrix.

### 7. Combined Statistical-Systematic Uncertainty

The final trigger-efficiency uncertainty surface combines two complementary approaches:

| Region | Method |
|--------|--------|
| **Empirical region** (within the ($\log_{10}N_\gamma,\theta_\text{inc}$) range of triggered data) | Bootstrap resampling of the raw (photon_count, incidence) pairs → 1000 replicas → 68 % percentile band ($N_\text{boot}=1000$, random binning & smoothing per replica) |
| **Extrapolation region** (outside data coverage) | Analytic Jacobian error propagation: $\sigma_{\text{jac}}^2 = J(\theta^*) \cdot \Sigma_{\text{joint}} \cdot J(\theta^*)^T$, where the Jacobian $J$ is computed symbolically with SymPy and $\Sigma_{\text{joint}}$ is the full (b,x0,w) joint covariance matrix |

**Coverage blending**: The empirical and Jacobian uncertainties are blended using a smooth coverage mask $C(\log_{10}N_\gamma, \theta_\text{inc})$:

\[
\begin{align*}
\text{lo} &= \mu - \sqrt{C^2 \cdot (\mu - \text{lo}_\text{emp})^2 + (1-C)^2 \cdot (z\sigma_\text{jac})^2} \\
\text{hi} &= \mu + \sqrt{C^2 \cdot (\text{hi}_\text{emp} - \mu)^2 + (1-C)^2 \cdot (z\sigma_\text{jac})^2}
\end{align*}
\]

where $z$ is the Gaussian quantile for the chosen CL ($68\%$ → $z\approx1.0$) and $\mu$ is the blended median surface (empirical median in data region, parametric model elsewhere).

The final interpolators `triger_eff_ci()` return (nominal, low, high) for any query ($N_\gamma$,$\theta_\text{inc}$) pair via fast `map_coordinates` lookup.

### 8. Efficiency Grid on (Zenith, Azimuth, Height, Energy)

The trigger-efficiency model is promoted to a 4-dimensional grid:

**Axes**:
- Zenith angle ($\theta_z$: 10 bins, linear from data min to max)
- Azimuth offset from telescope pointing ($\phi$: 5 bins, $0$–$10^\circ$)
- Emission height ($H$: 10 bins, log-spaced from data min to max)
- Muon energy ($E$: 20 bins, log-spaced from 10 GeV to 1 PeV)

**Procedure**:

1. For each (E, H, Zenith, Azimuth) cell, collect the photon counts of all showers falling into that bin.
2. For each cell, evaluate `t_rigger_eff_ci()` at each azimuth offset (5 phi values) and average over azimuth using trapezoidal integration.
3. Bootstrap uncertainty from shower statistics is propagated in quadrature with the model-side uncertainty.
4. **H-fill**: Forward-fill NaN values along the height axis (emissions at greater heights tend to be undetectable; later heights may be populated).
5. **Slant-depth extrapolation**: For cells that remain empty after forward-fill, a pre-computed (Energy, Slant Depth) → photon-count lookup table is used to find photon counts from showers at similar slant depths and energies. The efficiency is then evaluated at those counts.
6. **Distance-transform fill**: Any remaining NaNs are filled by nearest-neighbour interpolation (EDT).
7. **Edge extrapolation**: Height range is extended from (data min/max) to (0 m, 120 km) and zenith range to the telescope field of view.

The result is a `RegularGridInterpolator` over (zenith, azimuth, height, log₁₀E).

### 9. Effective Area from Digitized Trigger-Rate Data

The central effective area is derived from digitised data of the trigger rate as a function of distance from the detector centre.

- **Input**: 13 data points $(R_i, \text{trigger\_rate}_i)$ from a published plot.
- **Processing**: Subtracted baseline, normalised to max = 1.
- **Fit**: Weighted linear regression $\text{eff}(R) = mR + b$ (weights $=1/\sigma^2$).
- **Effective area**:

  \[
  A_{\text{eff}} = 2\pi \int_0^\infty \epsilon(R) \cdot R \, dR
  \]

- **Uncertainty**: Monte Carlo propagation (1000 samples from the ($m,b$) covariance matrix).

Results:
- Central estimate: ~**29.5 cm²**
- MC mean ± std: values reported in the output
- Raw geometric area: $\pi \times (500\ \text{m})^2 \approx 7.85 \times 10^5\ \text{m}^2$

### 10. Muon Flux Calculation

**MCEq grid**: A pre-computed pickle file contains the differential muon flux $\frac{d^2N}{dE\,dX\,d\Omega}$ for 100 zenith bins in 85°–90°, including emission height $X$ (slant depth) and muon energy $E$.

**Flux integration**:

For each zenith slice in 85°–91.8°:

\[
\frac{dN}{dt} = A_{\text{eff}} \int \!\! \int \!\! \int \!\! \int
\Phi(E, X, \theta, \phi) \; \epsilon(E, H, \theta, \phi) \;
dE \; dX \; d\Omega
\]

where:
- $d\Omega = \cos\theta \, d\theta \, d\phi$ (solid angle element)
- Height is converted to slant depth using an analytic 5-layer exponential atmosphere model (`convert_altitude_to_depth()`)
- The integration order is: energy (trapezoidal), height/slant (trapezoidal), azimuth (trapezoidal with $\cos\phi$ weighting), and finally zenith (trapezoidal with $\cos^2\theta$ weighting accounting for the telescope's pointing response).

**Zenith extrapolation to 91.8°**: Since the MCEq grid only covers up to 90°, the 90° flux slice is copied to 10 new zenith angles linearly spaced in $\cos\theta$ between 90° and 91.8°. Efficiency is recomputed at these new angles using the 4D interpolator. This extends coverage to the full field of view of the telescope (pointing at 90° zenith, FoV ≈ ±1.8°).

**Rate per energy bin**: The function `Muon_Rate_per_E_Calcution()` returns the rate spectrum $d^2N / (dt\, dE)$ by skipping the energy integration step.

### Outputs & Final Visualisation

- **Trigger-efficiency surfaces** saved to:
  ```
  /scratch/general/vast/u1520754/data_Muon_Trinity/triger_eff_surfces.npz
  ```
  Contains: BOOT_MED, BOOT_LO, BOOT_HI, BOOT_GLOBAL, SIGMA_JAC, and grid definitions.

- **Plots saved**:
  - `../plot/eff_r.png` — effective area fit (trigger rate vs. radius)
  - `../plot/Muon_rate.png` — final muon rate spectrum with twin axes (rate and efficiency)

- **Final rate output**: ~0.3 muons per day (≈ 0.1 muons per 8-hour observing night), with uncertainty bounds from the lo/hi surfaces.

---

## Key Results

| Quantity | Value | Notes |
|----------|-------|-------|
| Energy range | 10 GeV – 1 PeV | Simulated muon energies |
| PE threshold | 20 PE | Detection criterion |
| Effective area (central) | ≈ 29.5 cm² | Derived from digitised data |
| Raw geometric area | ≈ 7.85 × 10⁵ m² | π × (500 m)² |
| Predicted muon rate | ≈ 0.3 day⁻¹ | Integrated over full spectrum |
| Per-night rate | ≈ 0.1 night⁻¹ | For an 8-hour observation |
| Efficiency vs. energy | Peaks ∼10⁻³–10⁻² | At low energy; drops steeply above 10⁵ GeV |
| Total simulation storage | Reported in TB | Sum of all file sizes across energies |
| Photon-per-MB ratio | Reported | Measure of simulation efficiency |

---

## Reference

