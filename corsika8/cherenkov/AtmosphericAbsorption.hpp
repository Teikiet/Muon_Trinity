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
    }
    
    bool isLoaded() const { return loaded_; }
    
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
        
        // Find wavelength indices
        size_t iw = 0;
        for (size_t i = 0; i < wavelengths_.size() - 1; ++i) {
            if (wavelengths_[i + 1] >= lambda_nm) {
                iw = i;
                break;
            }
            iw = i;
        }
        
        // Find height indices
        size_t ih = 0;
        for (size_t i = 0; i < heights_.size() - 1; ++i) {
            if (heights_[i + 1] >= h_m) {
                ih = i;
                break;
            }
            ih = i;
        }
        
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
     * Calculate transmission probability for curved (spherical) Earth geometry.
     * 
     * Uses the Chapman function approximation for air mass in a spherical atmosphere.
     * More accurate than flat-Earth for zenith angles > 60°.
     * 
     * The air mass factor accounts for the increased path length through denser
     * lower atmosphere when the ray curves around the Earth.
     * 
     * @param lambda_nm Wavelength in nm
     * @param h_emit_m Emission altitude in meters (above sea level)
     * @param h_det_m Detection altitude in meters (above sea level)
     * @param path_length_m Geometric path length in meters
     * @param earth_radius_m Earth's radius in meters (default 6.371e6)
     * @return Transmission probability in range [0, 1]
     */
    double getTransmissionSpherical(double lambda_nm, double h_emit_m, double h_det_m, 
                                     double path_length_m, double earth_radius_m = 6.371e6) const {
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
