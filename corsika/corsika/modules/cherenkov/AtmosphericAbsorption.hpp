#pragma once

/*
 * AtmosphericAbsorption.hpp
 * 
 * This class handles atmospheric absorption of Cherenkov photons using
 * pre-computed optical depth tables (similar to CHASM's abstable format).
 * 
 * The absorption table contains cumulative optical depth τ(λ, h) from altitude h
 * to the top of the atmosphere. For a photon traveling from h_emit to h_det:
 *   - Vertical optical depth: Δτ = τ(λ, h_emit) - τ(λ, h_det)
 *   - Slant path correction: τ_slant = Δτ / cos(θ)  where θ is zenith angle
 *   - Transmission probability: T = exp(-τ_slant)
 * 
 * File format (text):
 *   Line 1: n_wavelengths n_heights
 *   Line 2: wavelength_0 wavelength_1 ... wavelength_(n-1)  [nm]
 *   Line 3: height_0 height_1 ... height_(m-1)  [m]
 *   Lines 4+: ecoeff values, one row per wavelength, columns are heights
 * 
 * Usage:
 *   AtmosphericAbsorption abs("abstable.dat");
 *   double transmission = abs.getTransmission(wavelength_nm, h_emit_m, h_det_m, cos_zenith);
 */

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <random>

namespace corsika::cherenkov {

class AtmosphericAbsorption {
private:
    std::vector<double> wavelengths_;  // wavelengths in nm
    std::vector<double> heights_;      // heights in m
    std::vector<std::vector<double>> ecoeff_;  // optical depth τ(λ, h), indexed as [wavelength_idx][height_idx]
    /// Per-slab extinction alpha_k(λ) = dτ/dh of the model, [wavelength_idx][slab].
    /// Built once at load; the chord integral (getOpticalDepthChord) walks it.
    std::vector<std::vector<double>> alphaSlab_;
    bool loaded_{false};
    
public:
    AtmosphericAbsorption() = default;
    
    explicit AtmosphericAbsorption(const std::string& filename) {
        load(filename);
    }
    
    /**
     * Load absorption table from file.
     * Supports text format (.dat, .txt) with the format described above.
     */
    void load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("AtmosphericAbsorption: Cannot open file: " + filename);
        }
        
        // Read header: n_wavelengths n_heights
        size_t n_wavelengths, n_heights;
        file >> n_wavelengths >> n_heights;
        
        if (n_wavelengths == 0 || n_heights == 0) {
            throw std::runtime_error("AtmosphericAbsorption: Invalid header in file: " + filename);
        }
        
        // Read wavelengths
        wavelengths_.resize(n_wavelengths);
        for (size_t i = 0; i < n_wavelengths; ++i) {
            file >> wavelengths_[i];
        }
        
        // Read heights
        heights_.resize(n_heights);
        for (size_t i = 0; i < n_heights; ++i) {
            file >> heights_[i];
        }
        
        // Read ecoeff values (one row per wavelength)
        ecoeff_.resize(n_wavelengths);
        for (size_t iw = 0; iw < n_wavelengths; ++iw) {
            ecoeff_[iw].resize(n_heights);
            for (size_t ih = 0; ih < n_heights; ++ih) {
                file >> ecoeff_[iw][ih];
            }
        }
        
        if (file.fail()) {
            throw std::runtime_error("AtmosphericAbsorption: Error reading file: " + filename);
        }
        
        loaded_ = true;
        buildSlabExtinction();
    }
    
    bool isLoaded() const { return loaded_; }

  private:
    /// Index of the cell of a monotone (usually uniform) grid containing x.
    static size_t gridCell(std::vector<double> const& g, double x) {
        size_t const n = g.size();
        if (n < 2) return 0;
        double const d0 = g[1] - g[0];
        bool const uniform = (n > 2) && (d0 > 0.0) &&
                             (std::abs((g[n - 1] - g[0]) / d0 - double(n - 1)) < 1.0e-6);
        if (uniform) {
            double const t = std::ceil((x - g[0]) / d0) - 1.0;
            if (t <= 0.0) return 0;
            if (t >= double(n - 2)) return n - 2;
            return static_cast<size_t>(t);
        }
        size_t best = 0;
        for (size_t i = 0; i + 1 < n; ++i) {
            best = i;
            if (g[i + 1] >= x) break;
        }
        return best;
    }

    /// Per-slab extinction alpha_k = (tau_{k+1} - tau_k)/(h_{k+1} - h_k) at every
    /// wavelength node.  The table's tau is piecewise LINEAR in altitude (that is
    /// exactly what the bilinear lookup means), so this is the extinction of the model
    /// rather than a finite-difference approximation -- and the chord integral below is
    /// then exact for the model, not merely accurate.
    void buildSlabExtinction() {
        alphaSlab_.assign(wavelengths_.size(), {});
        for (size_t iw = 0; iw < wavelengths_.size(); ++iw) {
            auto& a = alphaSlab_[iw];
            a.assign(heights_.size() > 1 ? heights_.size() - 1 : 0, 0.0);
            for (size_t k = 0; k + 1 < heights_.size(); ++k) {
                double const dh = heights_[k + 1] - heights_[k];
                a[k] = (dh > 0.0)
                           ? std::max(0.0, (ecoeff_[iw][k + 1] - ecoeff_[iw][k]) / dh)
                           : 0.0;
            }
        }
    }

  public:
    
    /**
     * Interpolate optical depth at given wavelength and height.
     * Uses bilinear interpolation.
     * 
     * @param lambda_nm Wavelength in nm
     * @param h_m Height in meters
     * @return Optical depth τ(λ, h) from height h to top of atmosphere
     */
    double getOpticalDepth(double lambda_nm, double h_m) const {
        if (!loaded_ || wavelengths_.empty() || heights_.empty()) {
            return 0.0;  // No absorption
        }
        
        // Clamp to table bounds
        lambda_nm = std::max(wavelengths_.front(), std::min(wavelengths_.back(), lambda_nm));
        h_m = std::max(heights_.front(), std::min(heights_.back(), h_m));
        
        // Cell indices.  This is the hot path of every tau integration: one alpha()
        // evaluation is two of these calls, a 100 km grazing ray makes ~20k alpha
        // evaluations, and the previous linear scans cost ~70 comparisons each.
        // Both shipped grids are uniform (5 nm in wavelength, 1 km in height), so the
        // cell is arithmetic; the scan is kept only as a fallback for a non-uniform
        // grid.  (Deliberately stateless: the table is shared across threads, so no
        // mutable cache index.)
        auto cellOf = [](std::vector<double> const& grid, double const x) -> size_t {
            size_t const n = grid.size();
            if (n < 2) return 0;
            double const d0 = grid[1] - grid[0];
            bool const uniform =
                (n > 2) && (d0 > 0.0) &&
                (std::abs((grid[n - 1] - grid[0]) / d0 - double(n - 1)) < 1.0e-6);
            if (uniform) {
                double const t = std::ceil((x - grid[0]) / d0) - 1.0;
                if (t <= 0.0) return 0;
                if (t >= double(n - 2)) return n - 2;
                return static_cast<size_t>(t);
            }
            size_t best = 0;
            for (size_t i = 0; i + 1 < n; ++i) {
                best = i;
                if (grid[i + 1] >= x) break;
            }
            return best;
        };
        size_t const iw = cellOf(wavelengths_, lambda_nm);
        size_t const ih = cellOf(heights_, h_m);
        
        // Bilinear interpolation
        size_t iw1 = std::min(iw + 1, wavelengths_.size() - 1);
        size_t ih1 = std::min(ih + 1, heights_.size() - 1);
        
        double fw = 0.0;
        if (iw1 > iw && wavelengths_[iw1] != wavelengths_[iw]) {
            fw = (lambda_nm - wavelengths_[iw]) / (wavelengths_[iw1] - wavelengths_[iw]);
        }
        
        double fh = 0.0;
        if (ih1 > ih && heights_[ih1] != heights_[ih]) {
            fh = (h_m - heights_[ih]) / (heights_[ih1] - heights_[ih]);
        }
        
        // Interpolate along height first, then wavelength
        double tau_w0 = ecoeff_[iw][ih] * (1.0 - fh) + ecoeff_[iw][ih1] * fh;
        double tau_w1 = ecoeff_[iw1][ih] * (1.0 - fh) + ecoeff_[iw1][ih1] * fh;
        
        return tau_w0 * (1.0 - fw) + tau_w1 * fw;
    }
    
    /**
     * Calculate transmission probability for a photon traveling from h_emit to h_det.
     * 
     * The optical depth along the slant path is:
     *   τ_slant = (τ(λ, h_emit) - τ(λ, h_det)) / cos(θ)
     * 
     * where θ is the angle from vertical (zenith angle of the photon path).
     * 
     * @param lambda_nm Wavelength in nm
     * @param h_emit_m Emission altitude in meters
     * @param h_det_m Detection altitude in meters
     * @param cos_zenith Cosine of the zenith angle of the photon path (|cos_zenith| used)
     * @return Transmission probability in range [0, 1]
     */
    double getTransmission(double lambda_nm, double h_emit_m, double h_det_m, double cos_zenith) const {
        if (!loaded_) {
            return 1.0;  // No absorption file loaded - full transmission
        }
        
        // Get optical depths at emission and detection heights
        double tau_emit = getOpticalDepth(lambda_nm, h_emit_m);
        double tau_det = getOpticalDepth(lambda_nm, h_det_m);
        
        // Vertical optical depth between emission and detection
        double delta_tau_vertical = std::abs(tau_emit - tau_det);
        
        // Apply slant path correction
        // Use |cos_zenith| to handle both upward and downward going photons
        double abs_cos_zenith = std::abs(cos_zenith);
        if (abs_cos_zenith < 0.01) {
            // For nearly horizontal paths, use a minimum to avoid division by zero
            // This corresponds to air mass ~100
            abs_cos_zenith = 0.01;
        }
        
        double tau_slant = delta_tau_vertical / abs_cos_zenith;
        
        // Transmission probability
        return std::exp(-tau_slant);
    }
    
    /**
     * Transmission for a straight chord in spherical geometry -- EXACT for the
     * table's slab model.
     *
     * The chord is fixed by its three side lengths (r_emit, r_det, path_length), from
     * which getOpticalDepthChord() derives both zenith cosines and the dip test, and
     * the optical depth is the exact slab integral; no air-mass or Chapman
     * approximation is involved.  The previous implementation is kept, unchanged, as
     * getTransmissionSphericalChapman() so the two can be compared.
     *
     * @param lambda_nm Wavelength in nm
     * @param h_emit_m Emission altitude in meters (above sea level)
     * @param h_det_m Detection altitude in meters (above sea level)
     * @param path_length_m Geometric chord length in meters
     * @param earth_radius_m Earth's radius in meters (default 6.371e6)
     * @return Transmission probability in range [0, 1]
     */
    double getTransmissionSpherical(double lambda_nm, double h_emit_m, double h_det_m,
                                    double path_length_m,
                                    double earth_radius_m = 6.371e6) const {
        // The chord length is passed through unchanged: the three side lengths
        // (r_emit, r_det, L) fix the triangle, so both zenith cosines and the dip test
        // are unique.  Passing the emission ANGLE instead is ambiguous -- two chords
        // solve it -- which is the bug getOpticalDepthChord() now cannot have.
        return getTransmissionChord(lambda_nm, h_emit_m, h_det_m, path_length_m, earth_radius_m);
    }

    /**
     * Previous implementation, kept for A/B comparison: air-mass multiplication of the
     * vertical optical depth with a Chapman correction for grazing angles.  Retained
     * only so the tests can quote the bias the exact integral removed; production calls
     * getTransmissionSpherical().
     */
    double getTransmissionSphericalChapman(double lambda_nm, double h_emit_m, double h_det_m,
                                           double path_length_m,
                                           double earth_radius_m = 6.371e6) const {
        if (!loaded_) {
            return 1.0;  // No absorption file loaded - full transmission
        }
        
        // Get optical depths at emission and detection heights
        double tau_emit = getOpticalDepth(lambda_nm, h_emit_m);
        double tau_det = getOpticalDepth(lambda_nm, h_det_m);
        
        // Vertical optical depth between emission and detection
        double delta_tau_vertical = std::abs(tau_emit - tau_det);
        
        // Calculate vertical distance
        double const dh = std::abs(h_emit_m - h_det_m);
        
        if (dh < 1.0 || path_length_m < 1.0) {
            // Very short path - use simple approximation
            return std::exp(-delta_tau_vertical);
        }
        
        // Compute air mass factor with spherical geometry correction
        // For flat Earth: air_mass = path_length / dh = 1/cos(zenith)
        // For curved Earth: need to account for ray passing through varying density shells
        
        double const cos_zenith_geom = dh / path_length_m;
        double const sin_zenith_geom = std::sqrt(std::max(0.0, 1.0 - cos_zenith_geom * cos_zenith_geom));
        
        // Atmospheric scale height (typical value for lower atmosphere)
        double const H_scale = 8400.0;  // meters
        
        // Average altitude and radius
        double const h_avg = 0.5 * (h_emit_m + h_det_m);
        double const r_avg = earth_radius_m + h_avg;
        
        // Chapman function parameter: X = r/H
        double const X = r_avg / H_scale;
        
        // Compute air mass using modified Chapman approximation
        // For moderate zenith angles (< 75°), use geometric air mass
        // For grazing angles, apply Chapman correction
        double air_mass;
        
        if (cos_zenith_geom > 0.26) {  // zenith < 75°
            // Standard geometric air mass with small curvature correction
            // Rozenberg formula for moderate angles
            air_mass = 1.0 / (cos_zenith_geom + 0.025 * std::exp(-11.0 * cos_zenith_geom));
        } else {
            // Grazing incidence - use full Chapman function approximation
            // Chapman function: Ch(X, θ) for exponential atmosphere
            // Simplified Kasten & Young (1989) formula adapted for finite paths
            double const y = std::sqrt(0.5 * X) * std::abs(cos_zenith_geom);
            
            if (y < 8.0) {
                // Full Chapman function approximation
                // Ch(X, θ) ≈ sqrt(π*X/2) * exp(X*(1-sin(θ))) * erfc(y)
                // Using Abramowitz & Stegun approximation for erfc
                double const t = 1.0 / (1.0 + 0.3275911 * y);
                double const erfc_y = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 
                                     + t * (-1.453152027 + t * 1.061405429)))) * std::exp(-y * y);
                
                double const chapman = std::sqrt(M_PI * X / 2.0) * 
                                       std::exp(X * (1.0 - sin_zenith_geom)) * erfc_y;
                
                // Clamp to reasonable range
                air_mass = std::min(chapman, 40.0);
            } else {
                // Asymptotic limit for very grazing angles
                air_mass = 1.0 / std::max(0.025, cos_zenith_geom);
            }
        }
        
        // Apply air mass to get slant optical depth
        double const tau_slant = delta_tau_vertical * air_mass;
        
        // Transmission probability
        return std::exp(-tau_slant);
    }
    
    /**
     * Monte Carlo decision: should this photon survive atmospheric absorption?
     * 
     * @param lambda_nm Wavelength in nm
     * @param h_emit_m Emission altitude in meters
     * @param h_det_m Detection altitude in meters
     * @param cos_zenith Cosine of the zenith angle of the photon path
     * @param rng Random number generator
     * @return true if photon survives, false if absorbed
     */
    template<typename RNG>
    bool survives(double lambda_nm, double h_emit_m, double h_det_m, double cos_zenith, RNG& rng) const {
        if (!loaded_) {
            return true;  // No absorption - all photons survive
        }
        
        double transmission = getTransmission(lambda_nm, h_emit_m, h_det_m, cos_zenith);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return uniform(rng) < transmission;
    }
    
    /**
     * Monte Carlo decision for spherical Earth geometry.
     * Uses Chapman function for accurate air mass calculation at large zenith angles.
     * 
     * @param lambda_nm Wavelength in nm
     * @param h_emit_m Emission altitude in meters
     * @param h_det_m Detection altitude in meters
     * @param path_length_m Geometric path length in meters
     * @param earth_radius_m Earth's radius in meters
     * @param rng Random number generator
     * @return true if photon survives, false if absorbed
     */
    template<typename RNG>
    bool survivesSpherical(double lambda_nm, double h_emit_m, double h_det_m, 
                           double path_length_m, double earth_radius_m, RNG& rng) const {
        if (!loaded_) {
            return true;  // No absorption - all photons survive
        }
        
        double transmission = getTransmissionSpherical(lambda_nm, h_emit_m, h_det_m, 
                                                        path_length_m, earth_radius_m);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        return uniform(rng) < transmission;
    }
    
    /**
     * EXACT optical depth along the straight chord from the emission height to the
     * detection height, for this table's model.
     *
     * The stored quantity is the cumulative column tau(h) = int_0^h alpha dh',
     * piecewise linear in altitude, so alpha is CONSTANT inside each 1 km slab.  With
     * the impact parameter p = r_emit sin(theta) (conserved along a chord) the path
     * length inside a slab is a difference of
     *     x(h) = sqrt(r(h)^2 - p^2),
     * hence
     *     tau = int_{r_min}^{r_emit} alpha dx + sum over legs,
     * which is exact to within the table's own interpolation.  r_min is p when the
     * chord dips (the foot of the perpendicular falls strictly inside the segment) and
     * the lower endpoint otherwise.
     *
     * WHY THE ARGUMENT IS THE CHORD LENGTH AND NOT THE EMISSION ANGLE.
     * The plane is fixed by (r_emit, r_det) and the angle at the emission -- but that
     * pair is ambiguous: BOTH a near and a far chord solve it, and they differ by the
     * whole width of the "flat" middle.  Passing the emission angle and assuming the
     * far intersection made a nearly horizontal chord of 10 km be integrated as 126 km
     * (measured: tau = 5.08 instead of 0.367, 12.8x, on the acceptance set's
     * 3100 m -> 3000 m / 10 km row), i.e. a photon with 69 % survival was absorbed
     * 99.4 % of the time.  The three side lengths (r_emit, r_det, L) determine the
     * triangle uniquely, so taking L removes the ambiguity by construction: the zenith
     * cosine at each end follows from the law of cosines and the dip test becomes a
     * test of both base angles.
     *
     * This replaces getTransmissionSphericalChapman()'s air-mass approximation on the
     * straight-line route; the two are compared by testCherenkovAbsorption* and the
     * measured bias is recorded there.
     *
     * @param path_length_m Geometric chord length between the two ends [m].  Must be
     *        the REAL distance between them; <= 0 means "unknown" and falls back to
     *        the vertical column difference.
     */
    double getOpticalDepthChord(double lambda_nm, double h_emit_m, double h_det_m,
                                double path_length_m, double earth_radius_m) const {
        if (!loaded_ || alphaSlab_.empty() || heights_.size() < 2) return 0.0;

        lambda_nm = std::max(wavelengths_.front(), std::min(wavelengths_.back(), lambda_nm));
        double const hLo = heights_.front(), hHi = heights_.back();
        double const hE = std::max(hLo, std::min(hHi, h_emit_m));
        double const hD = std::max(hLo, std::min(hHi, h_det_m));
        double const rE = earth_radius_m + hE;
        double const rD = earth_radius_m + hD;
        double const L = path_length_m;
        if (!(L > 1.0)) {
            // Degenerate/absent chord (a caller that only knows the heights): the
            // vertical column difference is the only defined answer.
            return std::abs(getOpticalDepth(lambda_nm, hE) - getOpticalDepth(lambda_nm, hD));
        }
        // Zenith cosine at each end.  The law of cosines gives the angle at that vertex
        // from the INWARD radial, so the zenith cosine is its NEGATIVE (the same sign
        // trap that doubled the column at the zenith before it was caught).
        double czE = -(rE * rE + L * L - rD * rD) / (2.0 * rE * L);
        double czD = -(rD * rD + L * L - rE * rE) / (2.0 * rD * L);
        czE = std::max(-1.0, std::min(1.0, czE));
        czD = std::max(-1.0, std::min(1.0, czD));
        double const p = rE * std::sqrt(std::max(0.0, 1.0 - czE * czE));
        double const discD = rD * rD - p * p;
        if (!(discD >= 0.0)) {
            // Cannot happen for a real triangle; keep it finite rather than NaN.
            return std::abs(getOpticalDepth(lambda_nm, hE) - getOpticalDepth(lambda_nm, hD));
        }

        // The foot of the perpendicular from the geocentre lies strictly inside the
        // segment iff the chord is DESCENDING at both ends: descending at E means it
        // heads towards the foot, and descending still at D means it has not passed it.
        // Then the closest approach p is part of the path and the two legs are added;
        // otherwise one endpoint is the closest approach and the sum degenerates to
        // their difference, because cum(rMin) = 0.
        double rMin = std::min(rE, rD);
        if (czE < 0.0 && czD < 0.0) rMin = p;

        auto cumFromMin = [&](size_t iw, double rTop) -> double {
            double acc = 0.0;
            auto const& a = alphaSlab_[iw];
            for (size_t k = 0; k < a.size(); ++k) {
                double const r0 = earth_radius_m + heights_[k];
                double const r1 = earth_radius_m + heights_[k + 1];
                if (!(r1 > rMin)) continue; // entirely below the closest approach
                if (!(r0 < rTop)) break;    // entirely above the top
                double const lo = std::max(r0, rMin);
                double const hi = std::min(r1, rTop);
                if (!(hi > lo)) continue;
                double const xl = std::sqrt(std::max(0.0, lo * lo - p * p));
                double const xh = std::sqrt(std::max(0.0, hi * hi - p * p));
                acc += a[k] * (xh - xl);
            }
            return acc;
        };

        size_t const iw = gridCell(wavelengths_, lambda_nm);
        size_t const iw1 = std::min(iw + 1, wavelengths_.size() - 1);
        double fw = 0.0;
        if (iw1 > iw && wavelengths_[iw1] != wavelengths_[iw]) {
            fw = (lambda_nm - wavelengths_[iw]) / (wavelengths_[iw1] - wavelengths_[iw]);
        }
        double const t0 = cumFromMin(iw, rE) + cumFromMin(iw, rD);
        if (!(fw > 0.0) || iw1 == iw) return t0;
        double const t1 = cumFromMin(iw1, rE) + cumFromMin(iw1, rD);
        return t0 * (1.0 - fw) + t1 * fw;
    }

    /// Transmission from the exact chord optical depth (see above).
    double getTransmissionChord(double lambda_nm, double h_emit_m, double h_det_m,
                               double path_length_m, double earth_radius_m) const {
        if (!loaded_) return 1.0;
        return std::exp(-getOpticalDepthChord(lambda_nm, h_emit_m, h_det_m, path_length_m,
                                             earth_radius_m));
    }

    /**
     * Get range of wavelengths in the table
     */
    std::pair<double, double> getWavelengthRange() const {
        if (wavelengths_.empty()) return {0.0, 0.0};
        return {wavelengths_.front(), wavelengths_.back()};
    }
    
    /**
     * Get range of heights in the table
     */
    std::pair<double, double> getHeightRange() const {
        if (heights_.empty()) return {0.0, 0.0};
        return {heights_.front(), heights_.back()};
    }
};

} // namespace corsika::cherenkov
