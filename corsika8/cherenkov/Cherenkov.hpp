#pragma once

#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>

#include <functional>
#include <memory>
#include <utility>

namespace corsika::cherenkov {

  // Projection mode for Cherenkov photons
  enum class ProjectionMode {
    Ground,      // Project to fixed ground z-plane (fails for horizontal showers)
    ShowerAxis,  // Project to plane perpendicular to shower axis
    Telescope    // Project to telescope surface perpendicular to shower axis with finite radius
  };

  class CherenkovProcess : public corsika::ContinuousProcess<CherenkovProcess> {
    /// Uniform distribution for absorption survival draws
    std::uniform_real_distribution<double> absorptionUniform_{0.0, 1.0};
  public:
    // HitRecorder stores projected hit position (ground coords for Ground/Telescope modes,
    // shower-plane coords for ShowerAxis mode)
    using HitRecorder = std::function<void(
      double x_m,              // projected x coordinate
      double y_m,              // projected y coordinate
      double z_m,              // projected z coordinate
      double lambda_nm,        // wavelength
      double z_prod_m,         // production height
      double dir_x,            // photon direction x (projection coordinates)
      double dir_y,            // photon direction y (projection coordinates)
      double dir_z,            // photon direction z (projection coordinates)
      double dir_ground_x,     // photon direction x in ground coordinates
      double dir_ground_y,     // photon direction y in ground coordinates
      double dir_ground_z,     // photon direction z in ground coordinates
      double time_ns,          // arrival time
      double weight,           // bunch weight (number of photons this bunch represents)
      double n_prod,           // refractive index at production height
      double path_length_m)>;  // optical path length in meters

  private:
    double refractiveIndex_;
    double fLambdaMin;
    double fLambdaMax;
    double fBunchSize;
    double groundZ_m_;
    double earthRadiusM_;  // Earth's radius in meters (for height above sea level calculation)
    ProjectionMode projectionMode_;
    // Fixed shower axis (primary particle direction) used for all projections
    std::shared_ptr<corsika::ShowerAxis const> primaryShowerAxis_;
    // Shower core coordinates (for ShowerAxis projection mode)
    double showerCoreX_m_;
    double showerCoreY_m_;
    double showerCoreZ_m_;
    // Telescope parameters (for Telescope projection mode)
    double telescopeX_m_;
    double telescopeY_m_;
    double telescopeZ_m_;
    double telescopeRadius_m_;
    // Telescope pointing direction (azimuth and zenith angles in radians)
    double telescopeAzimuth_rad_;
    double telescopeZenith_rad_;
    // Shower pointing direction (as backup for telescope if not specified)
    double showerAzimuth_rad_;
    double showerZenith_rad_;
    // Optional altitude and wavelength-dependent refractive index: takes absolute z [m] and wavelength [nm]
    // Returns pair<n_phase, n_group> for refraction and timing calculations respectively
    // When enableDispersion is false, both values are equal (no dispersion correction)
    std::function<std::pair<double, double>(double, double)> refractiveIndexFunc_;
    bool enableDispersion_{};  // apply wavelength-dependent group delay when true
    bool enableCurvature_{};   // account for curved photon paths in stratified atmosphere when true
    // Atmospheric absorption table for Monte Carlo absorption of Cherenkov photons
    AtmosphericAbsorption absorptionTable_;
    bool enableAbsorption_{};  // apply atmospheric absorption when true
    HitRecorder recorder_;
    mutable std::shared_ptr<spdlog::logger> logger_;
    /// Dedicated RNG for atmospheric absorption decisions.
    /// Isolated from the main Cherenkov RNG to avoid disturbing
    /// photon generation and particle propagation sequences.
    mutable std::mt19937_64 absorptionRNG_;
    // Extinction coefficient function: returns extinction per meter at altitude and wavelength
    std::function<double(double, double)> extinctionCoeffFunc_;

    // Helper function to compute local atmospheric scale height from refractive index gradient
    double computeLocalScaleHeight(double z) const;

    // Spherical ray tracing to tilted plane using Bouguer's formula: n(r)*r*sin(theta) = const
    // Returns true if ray reaches the plane; outputs impact point and optical/geometric path lengths
    bool rayTraceSphericalToPlane(double x0, double y0, double z0_altitude,
           double dx, double dy, double dz,
           double plane_x, double plane_y, double plane_z,
           double norm_x, double norm_y, double norm_z,
           double& x_impact, double& y_impact, double& z_impact,
           double& optical_path, double& geometric_path,
           double lambdaNm,
           double& accumulated_tau) const;

    // Spherical ray tracing using Bouguer's formula: n(r)*r*sin(theta) = const
    // Returns true if ray reaches ground; outputs ground x/y (Cartesian), optical and geometric path lengths
    bool rayTraceSpherical(double x0, double y0, double z0_altitude,
                 double dx, double dy, double dz,
                 double& x_ground, double& y_ground, double& z_ground,
                 double& optical_path, double& geometric_path,
                 double lambdaNm,
                 double& accumulated_tau) const;

  public:
    CherenkovProcess(double refractiveIndex,
                     double lambdaMin   = 300.0,
                     double lambdaMax   = 900.0,
                     double bunchSize   = 5.0,
                     double groundZ_m   = 0.0,
                     HitRecorder recorder = HitRecorder{},
                     ProjectionMode mode = ProjectionMode::Ground,
                     std::shared_ptr<corsika::ShowerAxis const> primaryShowerAxis = nullptr,
                     double showerCoreX_m = 0.0,
                     double showerCoreY_m = 0.0,
                     double showerCoreZ_m = 0.0,
                     double telescopeX_m = 0.0,
                     double telescopeY_m = 0.0,
                     double telescopeZ_m = 0.0,
                     double telescopeRadius_m = 0.0,
                     std::function<std::pair<double, double>(double, double)> refractiveIndexFunc = {},
                     double earthRadiusM = 6.371e6,
                     double telescopeAzimuth_rad = -1.0,
                     double telescopeZenith_rad = -1.0,
                     double showerAzimuth_rad = 0.0,
                     double showerZenith_rad = 0.0,
                     bool enableDispersion = true,
                     bool enableCurvature = false,
                     const std::string& absorptionFile = "");

    template <typename TParticle, typename TTrack>
    corsika::units::si::LengthType
    getMaxStepLength(TParticle const& particle, TTrack const& track) const;

    // REQUIRED by ContinuousProcess
    template <typename TParticle>
    corsika::ProcessReturn doContinuous(corsika::Step<TParticle>& step,
                                        bool const stepLimit) const;
  };

} // namespace corsika::cherenkov

#include <corsika/modules/cherenkov/Cherenkov.inl>
