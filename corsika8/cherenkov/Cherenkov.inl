#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/random/RNGManager.hpp>

namespace corsika::cherenkov {

  using namespace corsika::units::si;

  inline CherenkovProcess::CherenkovProcess(
      double refractiveIndex,
      double lambdaMin,
      double lambdaMax,
      double bunchSize,
      double groundZ_m,
      HitRecorder recorder,
      ProjectionMode mode,
      std::shared_ptr<corsika::ShowerAxis const> primaryShowerAxis,
      double showerCoreX_m,
      double showerCoreY_m,
      double showerCoreZ_m,
      double telescopeX_m,
      double telescopeY_m,
      double telescopeZ_m,
      double telescopeRadius_m,
      std::function<std::pair<double, double>(double, double)> refractiveIndexFunc,
      double earthRadiusM,
      double telescopeAzimuth_rad,
      double telescopeZenith_rad,
      double showerAzimuth_rad,
      double showerZenith_rad,
      bool enableDispersion,
      bool enableCurvature,
        const std::string& absorptionFile)
      : refractiveIndex_(refractiveIndex)
      , fLambdaMin(lambdaMin)
      , fLambdaMax(lambdaMax)
      , fBunchSize(bunchSize)
      , groundZ_m_(groundZ_m)
      , earthRadiusM_(earthRadiusM)
      , projectionMode_(mode)
      , primaryShowerAxis_(primaryShowerAxis)
      , showerCoreX_m_(showerCoreX_m)
      , showerCoreY_m_(showerCoreY_m)
      , showerCoreZ_m_(showerCoreZ_m)
      , telescopeX_m_(telescopeX_m)
      , telescopeY_m_(telescopeY_m)
      , telescopeZ_m_(telescopeZ_m)
      , telescopeRadius_m_(telescopeRadius_m)
      , telescopeAzimuth_rad_(telescopeAzimuth_rad)
      , telescopeZenith_rad_(telescopeZenith_rad)
      , showerAzimuth_rad_(showerAzimuth_rad)
      , showerZenith_rad_(showerZenith_rad)
      , refractiveIndexFunc_(std::move(refractiveIndexFunc))
      , enableDispersion_(enableDispersion)
      , enableCurvature_(enableCurvature)
      , enableAbsorption_(false)
      , recorder_(std::move(recorder))
      , logger_(corsika::get_logger("CherenkovProcess"))
      , absorptionRNG_(std::mt19937_64(98765))
  {

    // Load atmospheric absorption table if file is provided
    if (!absorptionFile.empty()) {
      try {
        absorptionTable_.load(absorptionFile);
        enableAbsorption_ = true;
        auto [wl_min, wl_max] = absorptionTable_.getWavelengthRange();
        auto [h_min, h_max] = absorptionTable_.getHeightRange();
        CORSIKA_LOGGER_INFO(
            logger_,
            "Atmospheric absorption enabled: file={}, wavelength range=[{}, {}] nm, height range=[{}, {}] m",
            absorptionFile, wl_min, wl_max, h_min, h_max);
      } catch (const std::exception& e) {
        CORSIKA_LOGGER_ERROR(logger_, "Failed to load absorption file '{}': {}", absorptionFile, e.what());
        enableAbsorption_ = false;
      }
    }
    
    // Build extinction coefficient function if absorption table is loaded
    if (absorptionTable_.isLoaded()) {
      extinctionCoeffFunc_ = [this](double h_m, double lambda_nm) -> double {
        double const dh = 50.0;
        // Clamp to table bounds
        auto height_range = absorptionTable_.getHeightRange();
        double h_lo = std::max(height_range.first, h_m - dh * 0.5);
        double h_hi = std::min(height_range.second, h_m + dh * 0.5);
        double actual_dh = h_hi - h_lo;
        if (actual_dh <= 0.0) return 0.0;
        double tau_lo = absorptionTable_.getOpticalDepth(lambda_nm, h_lo);
        double tau_hi = absorptionTable_.getOpticalDepth(lambda_nm, h_hi);
        return std::max(0.0, (tau_lo - tau_hi) / actual_dh);
      };
    }

    std::string modeStr = (mode == ProjectionMode::Ground ? "Ground" : 
                          (mode == ProjectionMode::ShowerAxis ? "ShowerAxis" : "Telescope"));
    CORSIKA_LOGGER_INFO(
        logger_,
        "CherenkovProcess initialized: n={}, lambda=[{},{}] nm, groundZ={} m, earthRadius={} m, mode={}, absorption={}",
        refractiveIndex_, fLambdaMin, fLambdaMax, groundZ_m_, earthRadiusM_, modeStr, enableAbsorption_ ? "on" : "off");
    
    if (mode == ProjectionMode::Telescope) {
      CORSIKA_LOGGER_INFO(
          logger_,
          "Telescope surface: position=({}, {}, {}) m, radius={} m, azimuth={:.4f} rad, zenith={:.4f} rad",
          telescopeX_m_, telescopeY_m_, telescopeZ_m_, telescopeRadius_m_, 
          telescopeAzimuth_rad_, telescopeZenith_rad_);
    }
  }

  // Helper function to compute local atmospheric scale height from refractive index gradient
  inline double CherenkovProcess::computeLocalScaleHeight(double z) const {
    if (!refractiveIndexFunc_) return 8000.0;  // default
    
    double const dz = 100.0;  // small altitude step in meters
    // Use reference wavelength 400nm for scale height calculation (phase index)
    double const n1 = refractiveIndexFunc_(z, 400.0).first;
    double const n2 = refractiveIndexFunc_(z + dz, 400.0).first;
    double const N1 = n1 - 1.0;
    double const N2 = n2 - 1.0;
    
    // Guard against invalid refractivity values
    if (N1 <= 1e-12 || N2 <= 1e-12 || N2 >= N1) {
      return 8000.0;  // fallback for extreme altitudes
    }
    
    // Guard against N2 being too close to N1 (where log(N2/N1) → 0)
    // This would cause unreasonably large scale heights or division by near-zero
    double const ratio = N2 / N1;
    double const log_ratio = std::log(ratio);
    
    // If magnitude of logarithm is too small, the gradient is too flat
    if (std::abs(log_ratio) < 1e-6) {
      return 8000.0;  // fallback when gradient is effectively zero
    }
    
    // H = -dz / ln(N2/N1)
    return -dz / log_ratio;
  }

  // Ray trace photon in a spherically stratified atmosphere to a tilted plane (Bouguer's formula)
  // n(r) * r * sin(theta) = constant along the ray
  // Uses FULL SPHERICAL GEOMETRY (Earth-centered coordinates) for accurate Bouguer invariant
  // Properly handles super-inclined/grazing rays
  inline bool CherenkovProcess::rayTraceSphericalToPlane(
      double x0, double y0, double z0_altitude,
      double dx, double dy, double dz,
      double plane_x, double plane_y, double plane_z,
      double norm_x, double norm_y, double norm_z,
      double& x_impact, double& y_impact, double& z_impact,
      double& optical_path, double& geometric_path,
      double lambdaNm,
      double& accumulated_tau) const {

    // Require an altitude-dependent refractive index
      if (!refractiveIndexFunc_) {
        CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: no refractiveIndexFunc");
      return false;
    }

    double const R_E = earthRadiusM_;
    // Thread-local RNG for absorption culling (independent of physics RNG)
    static thread_local std::mt19937 localRNG(
        std::hash<std::thread::id>{}(std::this_thread::get_id()) ^ 
        static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    static thread_local std::uniform_real_distribution<double> unitDist(0.0, 1.0);

    // Draw survival threshold ONCE before tracing
    double const survivalThreshold = (enableAbsorption_ && extinctionCoeffFunc_) 
        ? unitDist(localRNG) : 0.0;

    // Position relative to Earth center (x0,y0,z0 are absolute coordinates in rootCS)
    double const px = x0;
    double const py = y0;
    double const pz = z0_altitude;
    double const r0 = std::sqrt(px * px + py * py + pz * pz);

    // Radial unit vector at emission point (points outward/up from Earth center)
    double const r_hat_x = px / r0;
    double const r_hat_y = py / r0;
    double const r_hat_z = pz / r0;

    // CRITICAL: cos(theta) = dot(direction, radial) for SPHERICAL geometry
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

    // Refractive index at emission altitude (above sea level)
    double const h0 = r0 - R_E;
    double const n0 = refractiveIndexFunc_(h0, lambdaNm).first;  // phase index for Bouguer

    // Ray invariant (Bouguer's formula): n * r * sin(theta) = constant
    // THIS IS NOW CORRECT FOR SPHERICAL GEOMETRY
    double const ray_inv = n0 * r0 * sin_theta0;


    // Helper: compute signed distance from point to plane
    auto distanceToPlane = [&](double x, double y, double z) -> double {
      return (x - plane_x) * norm_x + (y - plane_y) * norm_y + (z - plane_z) * norm_z;
    };

    // Tangential unit vector in propagation plane (perpendicular to radial)
    double tang_x = dx - cos_theta0 * r_hat_x;
    double tang_y = dy - cos_theta0 * r_hat_y;
    double tang_z = dz - cos_theta0 * r_hat_z;
    double const tang_mag = std::sqrt(tang_x * tang_x + tang_y * tang_y + tang_z * tang_z);

    double t_hat_x, t_hat_y, t_hat_z;
    if (tang_mag > 1e-12) {
      t_hat_x = tang_x / tang_mag;
      t_hat_y = tang_y / tang_mag;
      t_hat_z = tang_z / tang_mag;
    } else {
      // Nearly vertical ray - build an arbitrary tangential direction via cross product
      // Choose a reference vector not parallel to r_hat, then t_hat = r_hat x ref
      double ref_x = 0.0;
      double ref_y = 0.0;
      double ref_z = 1.0;
      if (std::abs(r_hat_z) > 0.9) {
        ref_x = 1.0;
        ref_y = 0.0;
        ref_z = 0.0;
      }

      t_hat_x = r_hat_y * ref_z - r_hat_z * ref_y;
      t_hat_y = r_hat_z * ref_x - r_hat_x * ref_z;
      t_hat_z = r_hat_x * ref_y - r_hat_y * ref_x;

      double const t_norm = std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
      t_hat_x /= t_norm;
      t_hat_y /= t_norm;
      t_hat_z /= t_norm;
    }

    // Compute azimuthal sign: determines which way around the propagation plane x_plane evolves
    // Use t_hat_y as reference: if tangential component has positive y, then x_plane evolves positive
    double const az_sign = (t_hat_y >= 0.0) ? 1.0 : -1.0;

    // 2D coordinates in the propagation plane (Earth-centered Cartesian)
    double x_plane = 0.0;
    double z_plane = r0;
    double r_cur = r0;
    
    // Current 3D position in world coordinates
    double x_cur = x0;
    double y_cur = y0;
    double z_cur = z0_altitude;

    // Initial direction: +1 for upward, -1 for downward
    int direction = (cos_theta0 > 0.0) ? 1 : -1;

    optical_path = 0.0;
    geometric_path = 0.0;

    int const max_steps = 10000;
    double const max_altitude = 100000.0;  // 100 km maximum

    double r_prev = r_cur;
    double x_prev = x_cur;
    double y_prev = y_cur;
    double z_prev = z_cur;
    double optical_prev = 0.0;
    double geometric_prev = 0.0;
    double dist_prev = distanceToPlane(x_cur, y_cur, z_cur);

    for (int step = 0; step < max_steps; ++step) {
      double const h_cur = r_cur - R_E;
      
      // Check if we've escaped atmosphere
      if (h_cur > max_altitude) {
        double const photon_theta_rad = std::acos(std::max(-1.0, std::min(1.0, -dz)));
        double const photon_theta_deg = photon_theta_rad * 180.0 / M_PI;
        CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: escaped atmosphere at step {}, h_cur={:.1f}m > max_altitude={:.1f}m, photon zenith=({:.2f} deg)", 
                           step, h_cur, max_altitude, photon_theta_deg);
        return false;
      }

      // Get both phase and group refractive indices
      auto const [n_phase, n_group] = refractiveIndexFunc_(h_cur, lambdaNm);

      // Current angle from Bouguer's formula (SPHERICAL) using phase index
      double sin_theta_cur = ray_inv / (n_phase * r_cur); //angle with respect to local radial direction at current position
      
      // Check for turning point
      if (sin_theta_cur > 1.0) {
        // Revert to previous state before turning point
        r_cur = r_prev;
        x_plane = x_prev;
        z_plane = z_prev;
        x_cur = x_prev;
        y_cur = y_prev;
        z_cur = z_prev;
        optical_path = optical_prev;
        geometric_path = geometric_prev;
        // Flip direction at turning point
        direction = -direction;
        if (step % 100 == 0) {
          CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: turning point at step {}, reverted and flipped direction to {}", 
                             step, direction);
        }
        continue;  // Skip rest of iteration with reverted state
      }
      
      double const cos_theta_cur = std::sqrt(std::max(0.0, 1.0 - sin_theta_cur * sin_theta_cur));
      double const tan_theta_cur = (cos_theta_cur > 1e-12) ? sin_theta_cur / cos_theta_cur : 1e12;


      // Compute signed distance from current 3D position to telescope plane
      double const cos_alpha_d = z_plane / r_cur;
      double const sin_alpha_d = x_plane / r_cur;
      double const cur3d_x = r_cur * (cos_alpha_d * r_hat_x + sin_alpha_d * t_hat_x);
      double const cur3d_y = r_cur * (cos_alpha_d * r_hat_y + sin_alpha_d * t_hat_y);
      double const cur3d_z = r_cur * (cos_alpha_d * r_hat_z + sin_alpha_d * t_hat_z);
      double const dist_to_plane = (cur3d_x - plane_x) * norm_x + (cur3d_y - plane_y) * norm_y + (cur3d_z - plane_z) * norm_z;
      double const abs_dist_to_plane = std::abs(dist_to_plane);

      // --- Switch to straight-line propagation if within 1000m of the plane ---
      if (abs_dist_to_plane <= 1000.0 && direction == -1) {
        // Reconstruct the current 3D direction from the 2D propagation plane state
        double const cos_alpha_c = z_plane / r_cur;
        double const sin_alpha_c = x_plane / r_cur;

        // Local radial direction in 3D (outward from Earth center)
        double const r3d_x = cos_alpha_c * r_hat_x + sin_alpha_c * t_hat_x;
        double const r3d_y = cos_alpha_c * r_hat_y + sin_alpha_c * t_hat_y;
        double const r3d_z = cos_alpha_c * r_hat_z + sin_alpha_c * t_hat_z;

        // Local tangential direction in 3D (perpendicular to radial, in propagation plane)
        double const t3d_x = -sin_alpha_c * r_hat_x + cos_alpha_c * t_hat_x;
        double const t3d_y = -sin_alpha_c * r_hat_y + cos_alpha_c * t_hat_y;
        double const t3d_z = -sin_alpha_c * r_hat_z + cos_alpha_c * t_hat_z;

        // Current photon direction in 3D: radial component (downward) + tangential component
        double const dir3d_x = -cos_theta_cur * r3d_x + sin_theta_cur * t3d_x;
        double const dir3d_y = -cos_theta_cur * r3d_y + sin_theta_cur * t3d_y;
        double const dir3d_z = -cos_theta_cur * r3d_z + sin_theta_cur * t3d_z;

        // Current 3D position (Earth-centered)
        double const pos3d_x = r_cur * (cos_alpha_c * r_hat_x + sin_alpha_c * t_hat_x);
        double const pos3d_y = r_cur * (cos_alpha_c * r_hat_y + sin_alpha_c * t_hat_y);
        double const pos3d_z = r_cur * (cos_alpha_c * r_hat_z + sin_alpha_c * t_hat_z);

        // Analytical straight-line intersection with plane: (p + t*d - plane)·n = 0
        double const num = (plane_x - pos3d_x) * norm_x + (plane_y - pos3d_y) * norm_y + (plane_z - pos3d_z) * norm_z;
        double const den = dir3d_x * norm_x + dir3d_y * norm_y + dir3d_z * norm_z;
        if (std::abs(den) < 1e-10) {
          // Ray is parallel to plane
          return false;
        }
        double const t_hit = num / den;
        if (t_hit < 0.0) {
          // Intersection is behind current position
          return false;
        }

        // Impact point
        x_impact = pos3d_x + t_hit * dir3d_x;
        y_impact = pos3d_y + t_hit * dir3d_y;
        z_impact = pos3d_z + t_hit * dir3d_z;

        // Accumulate optical and geometric path for the straight-line segment
        double const h_cur_plane = r_cur - R_E;
        double const h_impact = std::sqrt(x_impact * x_impact + y_impact * y_impact + z_impact * z_impact) - R_E;
        double const h_mid = 0.5 * (h_cur_plane + h_impact);
        auto const [n_phase_mid, n_group_mid] = refractiveIndexFunc_(h_mid, lambdaNm);
        double const seg_len = t_hit;
        geometric_path += seg_len;
        optical_path += n_group_mid * seg_len;

        // Accumulate extinction along the straight-line segment (sub-steps)
        int const n_tau_steps = std::max(1, static_cast<int>(seg_len / 50.0));
        double const ds_tau = seg_len / static_cast<double>(n_tau_steps);
        for (int ts = 0; ts < n_tau_steps; ++ts) {
          double const s_mid = (static_cast<double>(ts) + 0.5) * ds_tau;
          double const tau_x = pos3d_x + s_mid * dir3d_x;
          double const tau_y = pos3d_y + s_mid * dir3d_y;
          double const tau_z = pos3d_z + s_mid * dir3d_z;
          double const tau_r = std::sqrt(tau_x * tau_x + tau_y * tau_y + tau_z * tau_z);
          double const tau_h = tau_r - R_E;
          if (enableAbsorption_ && extinctionCoeffFunc_) {
            accumulated_tau += extinctionCoeffFunc_(tau_h, lambdaNm) * ds_tau;
          }
        }

        // Early termination: if transmission < random threshold, photon is already dead
        if (enableAbsorption_ && extinctionCoeffFunc_) {
          if (std::exp(-accumulated_tau) < survivalThreshold) {
            return false;
          }
        }

        return true;
      }

      // --- Normal curved ray tracing step (when > 500m from plane) ---
      double const H_scale = computeLocalScaleHeight(h_cur);
      double const dx_step_base = H_scale / 20.0;
      double dx;
      if (sin_theta_cur > 0.99) {
        dx = H_scale / 100.0;
      } else {
        dx = std::min(dx_step_base, 1000.0);
      }

      // Store previous state for interpolation
      r_prev = r_cur;
      x_prev = x_cur;
      y_prev = y_cur;
      z_prev = z_cur;
      optical_prev = optical_path;
      geometric_prev = geometric_path;
      double const dist_at_step_start = dist_prev;

      // Local angular position in Earth-centered plane
      double const cos_alpha = z_plane / r_cur;
      double const sin_alpha = x_plane / r_cur;

      // Coordinate increments for this radial step
      double dx_sign = -1.0;
      double dz_sign = -1.0;
      if (static_cast<double>(direction) > 0) {
        //dz is positive for upward rays
        dz_sign = 1.0;
      }
      double const dr = std::abs(dx * cos_theta_cur/(sin_theta_cur * cos_alpha + cos_theta_cur * sin_alpha));
      double const dx_step = dx_sign * dx; //abs(dr *  (tan_theta_cur * cos_alpha + sin_alpha));
      double const dz_step = dz_sign * std::abs(dr * (cos_alpha - tan_theta_cur * sin_alpha));

      // Path length element
      double const ds = std::sqrt(dx_step * dx_step + dz_step * dz_step);
      geometric_path += ds;

      // Accumulate optical path using group index (for timing)
      optical_path += n_group * ds;

      // --- FIX: Accumulate extinction tau (main loop) ---
      if (enableAbsorption_ && extinctionCoeffFunc_) {
        double const h_mid = 0.5 * (r_prev - R_E + r_cur - R_E);
        double const kext = extinctionCoeffFunc_(h_mid, lambdaNm);
        accumulated_tau += kext * ds;
        // Early termination: if transmission < random threshold, photon is already dead
        if (std::exp(-accumulated_tau) < survivalThreshold) {
          return false;
        }
      }

      // Update position in the propagation plane
      // Apply azimuthal sign to allow x_plane to be both positive and negative
      x_plane += az_sign * dx_step;
      z_plane += dz_step;
      r_cur = std::sqrt(x_plane * x_plane + z_plane * z_plane);

      // Convert 2D plane coordinates back to 3D world coordinates
      double const cos_alpha_step = z_plane / r_cur;
      double const sin_alpha_step = x_plane / r_cur;
      double const gx = r_cur * (cos_alpha_step * r_hat_x + sin_alpha_step * t_hat_x);
      double const gy = r_cur * (cos_alpha_step * r_hat_y + sin_alpha_step * t_hat_y);
      double const gz = r_cur * (cos_alpha_step * r_hat_z + sin_alpha_step * t_hat_z);

      x_cur = gx;
      y_cur = gy;
      z_cur = gz;  // absolute z in Earth-centered coordinates

      // Check for plane crossing: sign change in distance
      double const dist_cur = distanceToPlane(x_cur, y_cur, z_cur);

      if (dist_at_step_start * dist_cur <= 0.0) {
        // Plane crossing detected within this step!
        // Interpolate to find exact crossing point and accumulated paths
        double alpha;
        if (std::abs(dist_at_step_start - dist_cur) > 1e-12) {
          alpha = dist_at_step_start / (dist_at_step_start - dist_cur);
        } else {
          alpha = 0.5;
        }
        
        // Clamp alpha to [0, 1]
        alpha = std::max(0.0, std::min(1.0, alpha));

        // Interpolate position at plane crossing
        x_impact = x_prev + alpha * (x_cur - x_prev);
        y_impact = y_prev + alpha * (y_cur - y_prev);
        z_impact = z_prev + alpha * (z_cur - z_prev);

        // CRITICAL: Interpolate accumulated paths at crossing
        // This ensures correct arrival time accounting for curved ray path
        optical_path = optical_prev + alpha * (optical_path - optical_prev);
        geometric_path = geometric_prev + alpha * (geometric_path - geometric_prev);

        return true;
      }
      dist_prev = dist_cur;

      // Plane-aware escape check: only abort if the ray moves away from Earth
      // AND the distance to the plane is growing without a sign change.
      if (direction < 0 && r_cur > r_prev + 1e-6) {
        if (dist_prev * dist_cur > 0.0 && std::abs(dist_cur) > std::abs(dist_prev)) {
          CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: ray diverging from plane at step {}, r_prev={:.1f}m, r_cur={:.1f}m",
                             step, r_prev, r_cur);
          return false;
        }
      }
    }

    CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: FAILED - max_steps={} reached without hitting plane", max_steps);
    return false;  // Max steps reached or ray never crossed plane
  }

  // Ray trace photon in a spherically stratified atmosphere with full 3D support
  // Handles both upward and downward-going rays using Bouguer's invariant
  // Detects turning points for upward rays and integrates until hitting ground
  inline bool CherenkovProcess::rayTraceSpherical(
    double x0, double y0, double z0_altitude,
    double dx, double dy, double dz,
    double& x_ground, double& y_ground, double& z_ground,
    double& optical_path, double& geometric_path,
    double lambdaNm,
    double& accumulated_tau) const {

    // Require an altitude-dependent refractive index; otherwise fall back to planar logic
    if (!refractiveIndexFunc_) {
      return false;
    }

    double const R_E = earthRadiusM_;
    // groundZ_m_ is absolute z (Earth-centered coordinate)
    double const r_ground = groundZ_m_;
    double const h_ground = r_ground - R_E;

    // Position relative to Earth center (x0,y0,z0 are absolute coordinates in rootCS)
    double const px = x0;
    double const py = y0;
    double const pz = z0_altitude;
    double const r0 = std::sqrt(px * px + py * py + pz * pz);

    // Already at or below ground
    if (r0 <= r_ground) {
      x_ground = x0;
      y_ground = y0;
      z_ground = z0_altitude;  // absolute z in Earth-centered coordinates
      optical_path = 0.0;
      geometric_path = 0.0;
      return true;
    }

    // Radial unit vector at emission point (points outward/up from Earth center)
    double const r_hat_x = px / r0;
    double const r_hat_y = py / r0;
    double const r_hat_z = pz / r0;

    // cos(theta) = dot(direction, radial)
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

    // Refractive index at emission altitude
    double const h0 = r0 - R_E;
    double const n0 = refractiveIndexFunc_(h0, lambdaNm).first;  // phase index for Bouguer

    // Ray invariant (Bouguer's formula): n * r * sin(theta) = constant
    double const ray_inv = n0 * r0 * sin_theta0;

    // Check if ray can physically reach ground (no total internal reflection)
    double const n_ground_val = refractiveIndexFunc_(h_ground, lambdaNm).first;  // phase index
    double const sin_theta_ground = ray_inv / (n_ground_val * r_ground);
    if (sin_theta_ground > 1.0 && cos_theta0 < 0.0) {
      // Only reject if ray is initially going downward and hits total internal reflection
      return false;  // Ray cannot reach ground
    }

    // Tangential unit vector in propagation plane (perpendicular to radial)
    double tang_x = dx - cos_theta0 * r_hat_x;
    double tang_y = dy - cos_theta0 * r_hat_y;
    double tang_z = dz - cos_theta0 * r_hat_z;
    double const tang_mag = std::sqrt(tang_x * tang_x + tang_y * tang_y + tang_z * tang_z);

    double t_hat_x, t_hat_y, t_hat_z;
    if (tang_mag > 1e-12) {
      t_hat_x = tang_x / tang_mag;
      t_hat_y = tang_y / tang_mag;
      t_hat_z = tang_z / tang_mag;
    } else {
      // Nearly vertical ray - build an arbitrary tangential direction via cross product
      // Choose a reference vector not parallel to r_hat, then t_hat = r_hat x ref
      double ref_x = 0.0;
      double ref_y = 0.0;
      double ref_z = 1.0;
      if (std::abs(r_hat_z) > 0.9) {
        ref_x = 1.0;
        ref_y = 0.0;
        ref_z = 0.0;
      }

      t_hat_x = r_hat_y * ref_z - r_hat_z * ref_y;
      t_hat_y = r_hat_z * ref_x - r_hat_x * ref_z;
      t_hat_z = r_hat_x * ref_y - r_hat_y * ref_x;

      double const t_norm = std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
      t_hat_x /= t_norm;
      t_hat_y /= t_norm;
      t_hat_z /= t_norm;
    }

    // Compute azimuthal sign: determines which way around the propagation plane x_plane evolves
    // Use t_hat_y as reference: if tangential component has positive y, then x_plane evolves positive
    double const az_sign = (t_hat_y >= 0.0) ? 1.0 : -1.0;

    // 2D coordinates in the propagation plane (Earth-centered Cartesian)
    double x_plane = 0.0;
    double z_plane = r0;
    double r_cur = r0;
    
    // Initial direction: +1 for upward, -1 for downward
    int direction = (cos_theta0 > 0.0) ? 1 : -1;

    optical_path = 0.0;
    geometric_path = 0.0;
    accumulated_tau = 0.0;

    int const max_steps = 10000; 
    double const max_altitude = 100000.0;  // 100 km maximum

    double x_prev = x_plane;
    double z_prev = z_plane;
    double r_prev = r_cur;
    double optical_prev = 0.0;
    double geometric_prev = 0.0;

    for (int step = 0; step < max_steps; ++step) {
      double const h_cur = r_cur - R_E;
      // Check if we've escaped atmosphere
      if (h_cur > max_altitude) {
        double const photon_theta_rad = std::acos(std::max(-1.0, std::min(1.0, -dz)));
        double const photon_theta_deg = photon_theta_rad * 180.0 / M_PI;
        CORSIKA_LOGGER_WARN(logger_, "rayTraceSpherical: escaped atmosphere at step {}, h_cur={:.1f}m > max_altitude={:.1f}m, photon zenith=({:.2f} deg)", 
                           step, h_cur, max_altitude, photon_theta_deg);
        return false;
      }

      auto const [n_phase, n_group] = refractiveIndexFunc_(h_cur, lambdaNm);
      double sin_theta_cur = ray_inv / (n_phase * r_cur);
      // Check for turning point
      if (sin_theta_cur > 1.0) {
        r_cur = r_prev;
        x_plane = x_prev;
        z_plane = z_prev;
        optical_path = optical_prev;
        geometric_path = geometric_prev;
        direction = -direction;
        continue;
      }
      double const cos_theta_cur = std::sqrt(std::max(0.0, 1.0 - sin_theta_cur * sin_theta_cur));
      double const tan_theta_cur = (cos_theta_cur > 1e-12) ? sin_theta_cur / cos_theta_cur : 1e12;


      // --- Switch to straight-line propagation if within 1000m of ground ---
      double const dist_to_ground = r_cur - r_ground;
      if (dist_to_ground <= 1000.0 && direction == -1) {
        // Reconstruct the current 3D direction from the 2D propagation plane state
        double const cos_alpha_c = z_plane / r_cur;
        double const sin_alpha_c = x_plane / r_cur;

        // Local radial direction in 3D (outward from Earth center)
        double const r3d_x = cos_alpha_c * r_hat_x + sin_alpha_c * az_sign * t_hat_x;
        double const r3d_y = cos_alpha_c * r_hat_y + sin_alpha_c * az_sign * t_hat_y;
        double const r3d_z = cos_alpha_c * r_hat_z + sin_alpha_c * az_sign * t_hat_z;

        // Local tangential direction in 3D (perpendicular to radial, in propagation plane)
        double const t3d_x = -sin_alpha_c * r_hat_x + cos_alpha_c * az_sign * t_hat_x;
        double const t3d_y = -sin_alpha_c * r_hat_y + cos_alpha_c * az_sign * t_hat_y;
        double const t3d_z = -sin_alpha_c * r_hat_z + cos_alpha_c * az_sign * t_hat_z;

        // Current photon direction in 3D: radial component (downward) + tangential component
        double const dir3d_x = -cos_theta_cur * r3d_x + sin_theta_cur * t3d_x;
        double const dir3d_y = -cos_theta_cur * r3d_y + sin_theta_cur * t3d_y;
        double const dir3d_z = -cos_theta_cur * r3d_z + sin_theta_cur * t3d_z;

        // Current 3D position (Earth-centered)
        double const pos3d_x = r_cur * (cos_alpha_c * r_hat_x + sin_alpha_c * az_sign * t_hat_x);
        double const pos3d_y = r_cur * (cos_alpha_c * r_hat_y + sin_alpha_c * az_sign * t_hat_y);
        double const pos3d_z = r_cur * (cos_alpha_c * r_hat_z + sin_alpha_c * az_sign * t_hat_z);

        // Analytical straight-line intersection with ground sphere: |pos + t*dir|^2 = r_ground^2
        double const pos_dot_dir = pos3d_x * dir3d_x + pos3d_y * dir3d_y + pos3d_z * dir3d_z;
        double const pos_dot_pos = pos3d_x * pos3d_x + pos3d_y * pos3d_y + pos3d_z * pos3d_z;
        double const discriminant = pos_dot_dir * pos_dot_dir - (pos_dot_pos - r_ground * r_ground);
        if (discriminant < 0.0) {
          // Ray does not intersect the ground sphere — miss
          return false;
        }
        double const sqrt_disc = std::sqrt(discriminant);
        double t_hit = -pos_dot_dir - sqrt_disc;
        if (t_hit < 0.0) {
          t_hit = -pos_dot_dir + sqrt_disc;
        }
        if (t_hit < 0.0) {
          // Both intersections are behind us
          return false;
        }

        // Ground intersection point in 3D
        double const gx = pos3d_x + t_hit * dir3d_x;
        double const gy = pos3d_y + t_hit * dir3d_y;
        double const gz = pos3d_z + t_hit * dir3d_z;

        // Accumulate optical path and geometric path for the straight-line segment
        double const h_mid = 0.5 * (h_cur + h_ground);
        auto const [n_phase_mid, n_group_mid] = refractiveIndexFunc_(h_mid, lambdaNm);
        geometric_path += t_hit;
        optical_path += n_group_mid * t_hit;

        // Accumulate extinction along the straight-line segment (sub-steps)
        int const n_tau_steps = std::max(1, static_cast<int>(t_hit / 50.0));
        double const ds_tau = t_hit / static_cast<double>(n_tau_steps);
        for (int ts = 0; ts < n_tau_steps; ++ts) {
          double const s_mid = (static_cast<double>(ts) + 0.5) * ds_tau;
          double const tau_x = pos3d_x + s_mid * dir3d_x;
          double const tau_y = pos3d_y + s_mid * dir3d_y;
          double const tau_z = pos3d_z + s_mid * dir3d_z;
          double const tau_r = std::sqrt(tau_x * tau_x + tau_y * tau_y + tau_z * tau_z);
          double const tau_h = tau_r - R_E;
          if (enableAbsorption_ && extinctionCoeffFunc_) {
            accumulated_tau += extinctionCoeffFunc_(tau_h, lambdaNm) * ds_tau;
          }
        }

        x_ground = gx;
        y_ground = gy;
        z_ground = gz;
        return true;
      }

      // --- Normal curved ray tracing step (when > 1000m from ground) ---
      double const H_scale = computeLocalScaleHeight(h_cur);
      double const dx_step_base = H_scale / 20.0;
      double dx;
      if (sin_theta_cur > 0.99) {
        dx = H_scale / 100.0;
      } else {
        dx = std::min(dx_step_base, 1000.0);
      }

      // Store previous state for interpolation
      x_prev = x_plane;
      z_prev = z_plane;
      r_prev = r_cur;
      optical_prev = optical_path;
      geometric_prev = geometric_path;

      double const cos_alpha = z_plane / r_cur;
      double const sin_alpha = x_plane / r_cur;

      double dx_sign = -1.0;
      double dz_sign = -1.0;
      if (static_cast<double>(direction) > 0) {
        dz_sign = 1.0;
      }
      double const dr = std::abs(dx * cos_theta_cur/(sin_theta_cur * cos_alpha + cos_theta_cur * sin_alpha));
      double const dx_step = dx_sign * dx;
      double const dz_step = dz_sign * std::abs(dr * (cos_alpha - tan_theta_cur * sin_alpha));

      double const ds = std::sqrt(dx_step * dx_step + dz_step * dz_step);
      geometric_path += ds;
      optical_path += n_group * ds;

      // --- Accumulate extinction optical depth along the curved ray ---
      double const h_mid = 0.5 * (r_prev - R_E + r_cur - R_E);
      double const kext = extinctionCoeffFunc_ ? extinctionCoeffFunc_(h_mid, lambdaNm) : 0.0;
      accumulated_tau += kext * ds;

      // Update position in the propagation plane
      x_plane += az_sign * dx_step;
      z_plane += dz_step;
      r_cur = std::sqrt(x_plane * x_plane + z_plane * z_plane);

      // Detect crossing of ground radius and interpolate
      if (r_cur <= r_ground && direction < 0) {
        double alpha = 0.0;
        if (std::abs(r_prev - r_cur) > 1e-12) {
          alpha = (r_prev - r_ground) / (r_prev - r_cur);
        }
        alpha = std::max(0.0, std::min(1.0, alpha));

        // Partial step to ground
        double ds_partial = alpha * ds;
        double h_mid_partial = 0.5 * (r_prev - R_E + (r_ground - R_E));
        double kext_partial = extinctionCoeffFunc_ ? extinctionCoeffFunc_(h_mid_partial, lambdaNm) : 0.0;
        accumulated_tau += kext_partial * ds_partial;

        x_plane = x_prev + alpha * (x_plane - x_prev);
        z_plane = z_prev + alpha * (z_plane - z_prev);
        r_cur = r_ground;
        optical_path = optical_prev + alpha * (optical_path - optical_prev);
        geometric_path = geometric_prev + alpha * (geometric_path - geometric_prev);

        // Convert plane coordinates back to 3D world coordinates (Step 9)
        double const gx = x_plane * t_hat_x + z_plane * r_hat_x;
        double const gy = x_plane * t_hat_y + z_plane * r_hat_y;
        double const gz = x_plane * t_hat_z + z_plane * r_hat_z;

        x_ground = gx;
        y_ground = gy;
        z_ground = gz;  // absolute z in Earth-centered coordinates
        return true;
      }

      // Safety: ray curving away from Earth when moving downward
      if (direction < 0 && r_cur > r_prev + 1e-6) {
        return false;
      }
    }

    return false;  // Max steps reached or ray escaped
  }

  // --- maximum step length for continuous effects ---
  template <typename TParticle, typename TTrack>
  inline corsika::units::si::LengthType
  CherenkovProcess::getMaxStepLength(TParticle const& particle,
                                      TTrack const& /*track*/) const {
    
    auto const E    = particle.getEnergy();
    auto const mass = particle.getMass();

    if (particle.getChargeNumber() == 0) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }

    double const m_over_E = mass / E;      // dimensionless
    double const beta2    = 1. - m_over_E * m_over_E;

    if (beta2 <= 0.) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }

    double const beta = std::sqrt(beta2);
    
    // Use refractive index for threshold check. For altitude-dependent case,
    // we use the constant value as a conservative estimate.
    // Actual per-point threshold check happens in doContinuous().
    double const n_threshold = refractiveIndex_;
    double const betaThreshold = 1.0 / n_threshold;

    if (beta < betaThreshold) {
      return std::numeric_limits<double>::infinity() * 1_m;
    }

    // For now: allow arbitrary step lengths when above threshold
    // Step will be limited by other processes (decay, interaction, etc.)
    return std::numeric_limits<double>::infinity() * 1_m;
  }

  // --- apply continuous Cherenkov photon production ---
  template <typename TParticle>
  inline corsika::ProcessReturn
  CherenkovProcess::doContinuous(corsika::Step<TParticle>& step,
                                  bool const /*stepLimit*/) const {

    auto const& projectile = step.getParticlePre();

    auto const E    = projectile.getEnergy();
    auto const mass = projectile.getMass();
    auto const q    = projectile.getChargeNumber();

    double const m_over_E = mass / E;
    double const beta2    = 1. - m_over_E * m_over_E;

    if (beta2 <= 0. || q == 0) {
      return corsika::ProcessReturn::Ok;
    }

    double const beta          = std::sqrt(beta2);
    // Local refractive index at production start height (absolute z in meters)
    // Use reference wavelength 400nm for threshold check (phase index)
    double const n_pre = refractiveIndexFunc_ ?
      refractiveIndexFunc_(step.getPositionPre().getZ(step.getPositionPre().getCoordinateSystem()) / 1_m, 400.0).first
      : refractiveIndex_;
    double const betaThreshold = 1.0 / n_pre;

    if (beta < betaThreshold) {
      return corsika::ProcessReturn::Ok;
    }

    // Calculate step length
    auto const posPre = step.getPositionPre();
    auto const posPost = step.getPositionPost();
    auto const displacement = posPost - posPre;
    double const stepLength_cm = displacement.getNorm() / 1_cm;

    if (stepLength_cm <= 0.) {
      return corsika::ProcessReturn::Ok;
    }

    double const alpha = 1.0 / 137.036;
    
    double const lambda_min_cm = fLambdaMin * 1e-7;
    double const lambda_max_cm = fLambdaMax * 1e-7;

    // Integrate photon production: dN/dxdλ = 2πα q² sin²θ(z,λ) / λ²
    // With dispersion, sin²θ = 1 - 1/(β² n(z,λ)²) depends on both altitude and wavelength
    double nPhotonsMean = 0.0;
    
    if (refractiveIndexFunc_ && enableDispersion_) {
      // Full 2D integration over position and wavelength when dispersion is enabled
      int const nPosSamples = 10;      // spatial samples along step
      int const nLambdaSamples = 10;   // wavelength samples
      double const dLambda_cm = (lambda_max_cm - lambda_min_cm) / nLambdaSamples;
      
      for (int i = 0; i < nPosSamples; ++i) {
        double const frac = (i + 0.5) / nPosSamples;
        auto const posIntegrate = posPre + displacement * frac;
        double const z_integrate = posIntegrate.getZ(posIntegrate.getCoordinateSystem()) / 1_m;
        
        // Wavelength integral: ∫ sin²θ(λ) / λ² dλ
        double wavelengthIntegral = 0.0;
        for (int j = 0; j < nLambdaSamples; ++j) {
          double const lambda_cm = lambda_min_cm + (j + 0.5) * dLambda_cm;
          double const lambda_nm = lambda_cm * 1e7;  // convert cm to nm
          
          // Use phase index for Cherenkov angle calculation (sin²θ)
          double const n_integrate = refractiveIndexFunc_(z_integrate, lambda_nm).first;
          double const sinTheta2 = 1.0 - 1.0 / (beta2 * n_integrate * n_integrate);
          
          if (sinTheta2 > 0.) {
            wavelengthIntegral += sinTheta2 / (lambda_cm * lambda_cm) * dLambda_cm;
          }
        }
        
        if (wavelengthIntegral > 0.) {
          double const photonsPerCm = 2.0 * M_PI * alpha * q * q * wavelengthIntegral;
          nPhotonsMean += photonsPerCm * stepLength_cm / nPosSamples;
        }
      }
    } else if (refractiveIndexFunc_) {
      // Altitude-dependent n but no dispersion: integrate over position only
      int const nSamples = 20;
      for (int i = 0; i < nSamples; ++i) {
        double const frac = (i + 0.5) / nSamples;
        auto const posIntegrate = posPre + displacement * frac;
        double const z_integrate = posIntegrate.getZ(posIntegrate.getCoordinateSystem()) / 1_m;
        // When dispersion disabled, wavelength parameter is ignored; use phase index
        double const n_integrate = refractiveIndexFunc_(z_integrate, 400.0).first;
        double const sinTheta2_integrate = 1.0 - 1.0 / (beta2 * n_integrate * n_integrate);
        
        if (sinTheta2_integrate > 0.) {
          double const photonsPerCm_integrate = 2.0 * M_PI * alpha * q * q * sinTheta2_integrate *
                                                 (1.0 / lambda_min_cm - 1.0 / lambda_max_cm);
          nPhotonsMean += photonsPerCm_integrate * stepLength_cm / nSamples;
        }
      }
    } else {
      // Constant refractive index case (original code)
      double const sinTheta2 = 1.0 - 1.0 / (beta2 * refractiveIndex_ * refractiveIndex_);
      if (sinTheta2 > 0.) {
        double const photonsPerCm = 2.0 * M_PI * alpha * q * q * sinTheta2 *
                                     (1.0 / lambda_min_cm - 1.0 / lambda_max_cm);
        nPhotonsMean = photonsPerCm * stepLength_cm;
      }
    }

    if (nPhotonsMean <= 0.) {
      return corsika::ProcessReturn::Ok;
    }

    auto& rng = RNGManager<>::getInstance().getRandomStream("cascade");
    std::poisson_distribution<int> poissonDist(nPhotonsMean);
    int nPhotons = poissonDist(rng);

    if (nPhotons <= 0) {
      return corsika::ProcessReturn::Ok;
    }

    int nBunches = std::max(1, static_cast<int>(nPhotons / fBunchSize));
    
    // Calculate weight per bunch (number of photons each bunch represents)
    double const weightPerBunch = static_cast<double>(nPhotons) / static_cast<double>(nBunches);

    // Counters for debugging telescope mode
    int hitCounter = 0;
    int missCounter = 0;
    int parallelCounter = 0;  // photons traveling parallel to telescope plane

    if (!recorder_) {
      return corsika::ProcessReturn::Ok;
    }

    auto const dir = projectile.getDirection();
    auto const& cs = posPre.getCoordinateSystem();

    double const px = dir.getX(cs);
    double const py = dir.getY(cs);
    double const pz = dir.getZ(cs);
    
    // Pre-compute perpendicular vectors to particle direction (independent of emission point)
    double perp1_x, perp1_y, perp1_z;
    if (std::abs(px) < 0.9) {
      perp1_x = 0.0;
      perp1_y = pz;
      perp1_z = -py;
    } else {
      perp1_x = py;
      perp1_y = -px;
      perp1_z = 0.0;
    }
    double const pnorm1 = std::sqrt(perp1_x*perp1_x + perp1_y*perp1_y + perp1_z*perp1_z);
    perp1_x /= pnorm1;
    perp1_y /= pnorm1;
    perp1_z /= pnorm1;

    double const perp2_x = py*perp1_z - pz*perp1_y;
    double const perp2_y = pz*perp1_x - px*perp1_z;
    double const perp2_z = px*perp1_y - py*perp1_x;

    // Pre-compute shower coordinate basis vectors for ShowerAxis/Telescope modes
    double sa_x_n = 0, sa_y_n = 0, sa_z_n = 0;
    double b1_x = 0, b1_y = 0, b1_z = 0;
    double b2_x = 0, b2_y = 0, b2_z = 0;
    bool useShowerCoords = false;

    // Pre-compute telescope coordinate basis vectors for Telescope mode
    double tel_z_x = 0, tel_z_y = 0, tel_z_z = 0;  // Z-axis (pointing direction)
    double tel_x_x = 0, tel_x_y = 0, tel_x_z = 0;  // X-axis (parallel to ground)
    double tel_y_x = 0, tel_y_y = 0, tel_y_z = 0;  // Y-axis (perpendicular to both)
    bool useTelescopeCoords = false;

    if ((projectionMode_ == ProjectionMode::ShowerAxis || 
         projectionMode_ == ProjectionMode::Telescope) && primaryShowerAxis_) {
      
      auto const& showerDir = primaryShowerAxis_->getDirection();
      double const sa_x = showerDir.getX(cs);
      double const sa_y = showerDir.getY(cs);
      double const sa_z = showerDir.getZ(cs);
      double const sa_norm = std::sqrt(sa_x*sa_x + sa_y*sa_y + sa_z*sa_z);
      
      if (sa_norm >= 1e-10) {
        useShowerCoords = true;
        sa_x_n = sa_x / sa_norm;
        sa_y_n = sa_y / sa_norm;
        sa_z_n = sa_z / sa_norm;

        // Create orthonormal basis on the shower plane
        if (std::abs(sa_x_n) < 0.9) {
          b1_x = 0.0;
          b1_y = sa_z_n;
          b1_z = -sa_y_n;
        } else {
          b1_x = -sa_z_n;
          b1_y = 0.0;
          b1_z = sa_x_n;
        }
        double const norm1 = std::sqrt(b1_x*b1_x + b1_y*b1_y + b1_z*b1_z);
        b1_x /= norm1;
        b1_y /= norm1;
        b1_z /= norm1;

        b2_x = sa_y_n * b1_z - sa_z_n * b1_y;
        b2_y = sa_z_n * b1_x - sa_x_n * b1_z;
        b2_z = sa_x_n * b1_y - sa_y_n * b1_x;
      }
    }

    // For Telescope mode: compute telescope coordinate system
    if (projectionMode_ == ProjectionMode::Telescope) {
      useTelescopeCoords = true;
      
      // Z-axis points in the telescope direction (azimuth, zenith)
      // azimuth: 0 = +x, pi/2 = +y, pi = -x, 3pi/2 = -y
      // zenith: 0 = -z (downward), pi/2 = horizontal, pi = +z (upward)
      double const cos_zen = std::cos(telescopeZenith_rad_);
      double const sin_zen = std::sin(telescopeZenith_rad_);
      double const cos_az = std::cos(telescopeAzimuth_rad_);
      double const sin_az = std::sin(telescopeAzimuth_rad_);
      
      tel_z_x = sin_zen * cos_az;
      tel_z_y = sin_zen * sin_az;
      tel_z_z = -cos_zen;  // Negative because zenith 0 points down
      
      // X-axis is parallel to ground and perpendicular to Z-axis
      // Ground is the xy-plane, so we need a horizontal vector perpendicular to Z
      // Take horizontal component of Z: (sin_zen*cos_az, sin_zen*sin_az, 0)
      // and rotate 90 degrees in xy-plane: (-sin_az, cos_az, 0)
      double const horiz_norm = std::sqrt(tel_z_x*tel_z_x + tel_z_y*tel_z_y);
      if (horiz_norm > 1e-10) {
        // Perpendicular to Z in xy-plane, pointing "right" when looking in telescope direction
        tel_x_x = -sin_az;
        tel_x_y = cos_az;
        tel_x_z = 0.0;
      } else {
        // Telescope pointing straight up or down: use x-axis direction
        tel_x_x = 1.0;
        tel_x_y = 0.0;
        tel_x_z = 0.0;
      }
      
      // Y-axis completes the right-handed system: Y = Z × X
      tel_y_x = tel_z_y * tel_x_z - tel_z_z * tel_x_y;
      tel_y_y = tel_z_z * tel_x_x - tel_z_x * tel_x_z;
      tel_y_z = tel_z_x * tel_x_y - tel_z_y * tel_x_x;
    }

    std::uniform_real_distribution<double> posDist(0., 1.);
    std::uniform_real_distribution<double> phiDist(0., 2.0 * M_PI);
    
    // Pre-compute constants for 1/λ² wavelength sampling (inverse transform method)
    double const lambda_min_inv = 1.0 / fLambdaMin;
    double const lambda_max_inv = 1.0 / fLambdaMax;
    double const lambda_inv_range = lambda_min_inv - lambda_max_inv;

    // Note: cosTheta_C and sinTheta_C will be recalculated per emission point
    // inside the bunch loop below to account for altitude variation

    for (int i = 0; i < nBunches; ++i) {
      double const frac = posDist(rng);
      auto const posEmit = posPre + displacement * frac;

      double const x0 = posEmit.getX(cs) / 1_m;
      double const y0 = posEmit.getY(cs) / 1_m;
      double const z0 = posEmit.getZ(cs) / 1_m;

      // Sample wavelength from 1/λ² distribution early (needed for ray tracing)
      // p(λ) ∝ 1/λ², so λ = 1 / (1/λ_min - u*(1/λ_min - 1/λ_max))
      double const u_wavelength = posDist(rng);
      double const lambdaNm = 1.0 / (lambda_min_inv - u_wavelength * lambda_inv_range);

      // Get altitude and wavelength-dependent refractive index at emission point
      // Use phase index for Cherenkov angle calculation
      double const n_local = refractiveIndexFunc_ ?
        refractiveIndexFunc_(z0, lambdaNm).first : refractiveIndex_;
      
      // Recalculate Cherenkov angle for this emission point
      double const cosTheta_C = 1.0 / (beta * n_local);
      
      // Check if particle is above Cherenkov threshold at this altitude
      if (cosTheta_C >= 1.0) {
        continue;  // Below threshold at this altitude
      }
      double const sinTheta_C = std::sqrt(1.0 - cosTheta_C * cosTheta_C);

      double const phi = phiDist(rng);

      // Photon direction in world coordinates
      double const dx = px * cosTheta_C + sinTheta_C * (perp1_x * std::cos(phi) + perp2_x * std::sin(phi));
      double const dy = py * cosTheta_C + sinTheta_C * (perp1_y * std::cos(phi) + perp2_y * std::sin(phi));
      double const dz = pz * cosTheta_C + sinTheta_C * (perp1_z * std::cos(phi) + perp2_z * std::sin(phi));

      // Direction to store (will be transformed for shower/telescope modes)
      double dx_out = dx;
      double dy_out = dy;
      double dz_out = dz;

      // Transform direction to shower coordinates if applicable
      if (useShowerCoords) {
        // Project direction onto shower coordinate basis
        // dx_shower = d · b1 (transverse component 1)
        // dy_shower = d · b2 (transverse component 2)  
        // dz_shower = d · sa (along shower axis, positive = towards ground)
        dx_out = dx * b1_x + dy * b1_y + dz * b1_z;
        dy_out = dx * b2_x + dy * b2_y + dz * b2_z;
        dz_out = dx * sa_x_n + dy * sa_y_n + dz * sa_z_n;
      }

      // Project to detection plane/surface
      double xg, yg, zg;
      double det_x_world = 0.0, det_y_world = 0.0, det_z_world = 0.0;  // world coords for travel time integration
      bool used_ray_trace = false;
      double ray_optical_path = 0.0;   // only valid when used_ray_trace
      double ray_geom_path = 0.0;      // only valid when used_ray_trace

      double ray_accumulated_tau = 0.0;
      if (projectionMode_ == ProjectionMode::Ground) {
        double ray_accumulated_tau = 0.0;
        if (enableCurvature_ && refractiveIndexFunc_) {
          bool hit = rayTraceSpherical(x0, y0, z0, dx, dy, dz,
                                       xg, yg, zg,
                                       ray_optical_path, ray_geom_path,
                                       lambdaNm,
                                       ray_accumulated_tau);
          if (hit) {
            used_ray_trace = true;
            det_x_world = xg;
            det_y_world = yg;
            det_z_world = zg;
          } else {
            // Curved ray tracing failed - skip this photon entirely
            continue;
          }
        } else if (enableCurvature_ && !refractiveIndexFunc_) {
          CORSIKA_LOGGER_WARN(logger_, "enableCurvature=true but refractiveIndexFunc is not set - using straight-line");
        }

        // Fallback: straight-line intersection with flat ground plane (only if curvature disabled)
        if (!used_ray_trace) {
          if (std::abs(dz) < 1e-10) {
            continue;
          }
          double const t_int = (groundZ_m_ - z0) / dz;
          if (t_int <= 0.) {
            continue;
          }
          xg = x0 + t_int * dx;
          yg = y0 + t_int * dy;
          zg = groundZ_m_;
          det_x_world = xg;
          det_y_world = yg;
          det_z_world = zg;
          ray_geom_path = 0.0;  // not used
          ray_optical_path = 0.0;
          ray_accumulated_tau = 0.0;
        }
      } else if (projectionMode_ == ProjectionMode::Telescope && useTelescopeCoords) {
        // For Telescope mode, project onto the telescope surface
        // The telescope surface is perpendicular to tel_z at distance 0 from telescopeX/Y/Z

        // Fast pre-check: does straight-line path come within 10x telescope radius?
        // This avoids expensive curved ray tracing for photons that will clearly miss
        bool passesPreCheck = false;
        if (enableCurvature_ && refractiveIndexFunc_) {
          // Straight-line intersection with telescope plane
          double const num_check =
              (telescopeX_m_ - x0) * tel_z_x +
              (telescopeY_m_ - y0) * tel_z_y +
              (telescopeZ_m_ - z0) * tel_z_z;

          double const den_check =
              dx * tel_z_x +
              dy * tel_z_y +
              dz * tel_z_z;

          if (std::abs(den_check) >= 1e-10) {
            double const t_check = num_check / den_check;
            if (t_check > 0.) {
              // Straight-line intersection point
              double const x_check = x0 + t_check * dx;
              double const y_check = y0 + t_check * dy;
              double const z_check = z0 + t_check * dz;

              // Distance from telescope center
              double const rx_check = x_check - telescopeX_m_;
              double const ry_check = y_check - telescopeY_m_;
              double const rz_check = z_check - telescopeZ_m_;

              // Project onto telescope plane
              double const xg_check = rx_check * tel_x_x + ry_check * tel_x_y + rz_check * tel_x_z;
              double const yg_check = rx_check * tel_y_x + ry_check * tel_y_y + rz_check * tel_y_z;

              // Check if within 10x radius (safety factor for curved paths)
              double const dist_check = std::sqrt(xg_check*xg_check + yg_check*yg_check);
              if (dist_check <= 10.0 * telescopeRadius_m_) {
                passesPreCheck = true;
              }
            }
          }
        } else {
          // If curvature disabled, always do the check (will fall through to straight-line code)
          passesPreCheck = true;
        }

        // Skip expensive curved ray tracing if straight-line check failed
        if (!passesPreCheck) {
          missCounter++;
          continue;
        }

        // Try spherical ray tracing when curvature is enabled and pre-check passed
        if (enableCurvature_ && refractiveIndexFunc_) {
          bool hit = rayTraceSphericalToPlane(
              x0, y0, z0,
              dx, dy, dz,
              telescopeX_m_, telescopeY_m_, telescopeZ_m_,
              tel_z_x, tel_z_y, tel_z_z,
              det_x_world, det_y_world, det_z_world,
              ray_optical_path, ray_geom_path,
              lambdaNm,
              ray_accumulated_tau);
          
          if (hit) {
            used_ray_trace = true;
          } else {
            // Curved ray tracing failed - skip this photon entirely
            continue;
          }
        } else if (enableCurvature_ && !refractiveIndexFunc_) {
          CORSIKA_LOGGER_WARN(logger_, "enableCurvature=true but refractiveIndexFunc is not set - using straight-line");
        }
        
        // Fallback: straight-line ray-plane intersection (only if curvature disabled)
        if (!used_ray_trace) {
          double const num =
              (telescopeX_m_ - x0) * tel_z_x +
              (telescopeY_m_ - y0) * tel_z_y +
              (telescopeZ_m_ - z0) * tel_z_z;

          double const den =
              dx * tel_z_x +
              dy * tel_z_y +
              dz * tel_z_z;

          if (std::abs(den) < 1e-10) {
            parallelCounter++;
            continue;
          }

          double const t_int = num / den;
          if (t_int <= 0.) {
            missCounter++;
            continue;
          }

          // 3D intersection point in world coordinates
          det_x_world = x0 + t_int * dx;
          det_y_world = y0 + t_int * dy;
          det_z_world = z0 + t_int * dz;
          ray_geom_path = 0.0;
          ray_optical_path = 0.0;
        }

        // Vector from telescope center to intersection point
        double const rx = det_x_world - telescopeX_m_;
        double const ry = det_y_world - telescopeY_m_;
        double const rz = det_z_world - telescopeZ_m_;

        // Project onto telescope coordinate system
        xg = rx * tel_x_x + ry * tel_x_y + rz * tel_x_z;
        yg = rx * tel_y_x + ry * tel_y_y + rz * tel_y_z;
        zg = 0.0;  // On the telescope surface

        // Apply radius cut
        double const dist_from_center = std::sqrt(xg*xg + yg*yg);
        if (dist_from_center > telescopeRadius_m_) {
          missCounter++;
          continue;
        }
        
        hitCounter++;
        
        // Transform photon direction from world to telescope coordinates
        dx_out = dx * tel_x_x + dy * tel_x_y + dz * tel_x_z;
        dy_out = dx * tel_y_x + dy * tel_y_y + dz * tel_y_z;
        dz_out = dx * tel_z_x + dy * tel_z_y + dz * tel_z_z;

      } else {
        // ShowerAxis mode

        if (!useShowerCoords) {
          // Fallback to ground projection
          if (std::abs(dz) < 1e-10) continue;
          double const t_int = (groundZ_m_ - z0) / dz;
          if (t_int <= 0.) continue;
          xg = x0 + t_int * dx;
          yg = y0 + t_int * dy;
          zg = groundZ_m_;
          det_x_world = xg;
          det_y_world = yg;
          det_z_world = zg;
        } else {
          // Reference point on the plane depends on mode
          double plane_x, plane_y, plane_z;
          plane_x = showerCoreX_m_;
          plane_y = showerCoreY_m_;
          plane_z = showerCoreZ_m_;

          // Ray-plane intersection
          double const num =
              (plane_x - x0) * sa_x_n +
              (plane_y - y0) * sa_y_n +
              (plane_z - z0) * sa_z_n;

          double const den =
              dx * sa_x_n +
              dy * sa_y_n +
              dz * sa_z_n;

          if (std::abs(den) < 1e-10) continue;

          double const t_int = num / den;
          if (t_int <= 0.) continue;

          // 3D intersection point
          double const xi = x0 + t_int * dx;
          double const yi = y0 + t_int * dy;
          double const zi = z0 + t_int * dz;
          det_x_world = xi;
          det_y_world = yi;
          det_z_world = zi;

          // Vector from plane reference point to intersection
          double const rx = xi - plane_x;
          double const ry = yi - plane_y;
          double const rz = zi - plane_z;

          // Project onto plane coordinates (centered at plane reference point)
          xg = rx * b1_x + ry * b1_y + rz * b1_z;
          yg = rx * b2_x + ry * b2_y + rz * b2_z;
          zg = 0.0;
        }
      }

      // Calculate photon time
      auto const timePre = step.getTimePre();
      auto const timePost = step.getTimePost();
      auto const timeEmit = timePre + (timePost - timePre) * frac;

      // Compute arrival time
      constexpr double c_light = 299792458.0; // m/s
      double totalTime_ns;

      // Variables needed for timing calculation (declared here for scope)
      double path_len = 0.0;    // geometric path length
      double delay_ns = 0.0;    // atmospheric delay

      if (used_ray_trace) {
        // Ray tracer computed actual curved paths through atmosphere
        // ray_optical_path = ∫ n(s) ds  (accounts for refractive index along curved ray)
        // ray_geom_path = ∫ ds          (physical path length along curved ray)
        //
        // Travel time = optical_path / c
        // This correctly accounts for both geometry (curved path) and refractive slowdown
        double const t_propagation_ns = (ray_optical_path / c_light) * 1e9;
        
        totalTime_ns = (timeEmit / 1_ns) + t_propagation_ns;
        
      } else {
        // Straight-line approximation with refractive index integration
        double const dx_path = det_x_world - x0;
        double const dy_path = det_y_world - y0;
        double const dz_path = det_z_world - z0;
        path_len = std::sqrt(dx_path*dx_path + dy_path*dy_path + dz_path*dz_path);

        double path_correction = 1.0;
        
        // Check if curvature correction is needed (large zenith angle)
        double const h_horiz = std::sqrt((det_x_world - x0)*(det_x_world - x0) + 
                                          (det_y_world - y0)*(det_y_world - y0));
        double const z_drop = std::abs(z0 - det_z_world);
        
        if (enableCurvature_ && z_drop > 100.0 && h_horiz > 0.1) {
          double const theta_rad = std::atan2(h_horiz, z_drop);
          double const theta_deg = theta_rad * 180.0 / M_PI;
          static constexpr double kCurvatureThresholdDeg = 60.0;
          
          if (theta_deg > kCurvatureThresholdDeg) {
            double const z_avg = 0.5 * (z0 + det_z_world);
            double const H = computeLocalScaleHeight(z_avg);
            double const n_avg = refractiveIndexFunc_ ? refractiveIndexFunc_(z_avg, lambdaNm).first : refractiveIndex_;
            double const N_avg = n_avg - 1.0;
            double const abs_dN_dz = N_avg / H;
            double const tan_theta = std::tan(theta_rad);
            double const delta_theta = abs_dN_dz * z_drop * tan_theta;
            path_correction = 1.0 + (delta_theta * delta_theta) / 24.0;
          }
        }
        
        // Geometric time: straight-line distance at vacuum speed of light
        double const t_geom_ns = (path_len / c_light) * 1e9;

        delay_ns = 0.0;  // atmospheric delay only

        if (path_len < 1e-2) {
          // Very short path: single evaluation at emission altitude using group index
          auto const [n_phase_emit, n_group_emit] = refractiveIndexFunc_ ? 
              refractiveIndexFunc_(z0, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
          double const N_group = n_group_emit - 1.0;
          delay_ns = (path_len * path_correction * N_group / c_light) * 1e9;
        } else {
          double const uz_path = dz_path / path_len;
          double const delta_z = std::abs(dz_path);

          int segments = 10;
          if (delta_z < 100.0) {
            segments = 10;
          } else if (delta_z < 1000.0) {
            segments = 20;
          } else if (delta_z < 10000.0) {
            segments = 50;
          } else {
            segments = 100;
          }

          static constexpr double kNearHorizontalUz = 1e-2;  // |uz| threshold for near-horizontal treatment

          // Near-horizontal paths: use average altitude with single segment
          if (std::abs(uz_path) < kNearHorizontalUz) {
            double const z_avg = 0.5 * (z0 + det_z_world);
            auto const [n_phase_avg, n_group_avg] = refractiveIndexFunc_ ? 
                refractiveIndexFunc_(z_avg, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
            double const N_group = n_group_avg - 1.0;
            delay_ns = (path_len * path_correction * N_group / c_light) * 1e9;
          } else {
            double const ds = path_len / static_cast<double>(segments);
            double delay_m = 0.0;  // integral of (n-1) ds in meters
            for (int si = 0; si < segments; ++si) {
              double const s_mid = (static_cast<double>(si) + 0.5) * ds;
              double const z_mid = z0 + s_mid * uz_path;
              // Use group index directly for timing
              auto const [n_phase_mid, n_group_mid] = refractiveIndexFunc_ ? 
                  refractiveIndexFunc_(z_mid, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
              double const N_group = n_group_mid - 1.0;
              delay_m += N_group * ds;
            }
            // Apply curvature correction to delay path as well
            delay_ns = (delay_m * path_correction / c_light) * 1e9;
          }
        }

        totalTime_ns = (timeEmit / 1_ns) + t_geom_ns + delay_ns;
      }

      // Convert production height from Earth-centered z to altitude above sea level
      // z0 is in Earth-centered coordinates (Earth radius + altitude)
      double const z0_above_sea_level = z0 - earthRadiusM_;

      // Apply atmospheric absorption (Monte Carlo)
      if (enableAbsorption_) {
        double const h_det_m = det_z_world - earthRadiusM_;
        bool photon_survives = true;
        static thread_local std::uniform_real_distribution<double> absorptionUniform(0.0, 1.0);
        if (used_ray_trace && ray_accumulated_tau > 0.0) {
          // Photon already survived during ray tracing - no additional check needed
          photon_survives = true;
        } else {
          // Fallback: post-hoc method for straight-line or when ray trace didn't compute tau
          double const dx_path = det_x_world - x0;
          double const dy_path = det_y_world - y0;
          double const dz_path = det_z_world - z0;
          double const path_length = std::sqrt(dx_path*dx_path + dy_path*dy_path + dz_path*dz_path);
          if (enableCurvature_) {
            photon_survives = absorptionTable_.survivesSpherical(
                lambdaNm, z0_above_sea_level, h_det_m, path_length, earthRadiusM_, absorptionRNG_);
          } else {
            double const cos_zenith = (path_length > 1e-6) ? std::abs(dz_path) / path_length : 1.0;
            photon_survives = absorptionTable_.survives(
                lambdaNm, z0_above_sea_level, h_det_m, cos_zenith, absorptionRNG_);
          }
        }
        if (!photon_survives) {
          continue;  // Photon absorbed
        }
      }

      // Record the bunch once with weight (number of photons it represents)
      // For Ground mode: directions are in world coordinates
      // For ShowerAxis/Telescope modes: directions are in shower coordinates
      // Pass optical_path for timing calculations (accounts for refractive index)
      double const pathLengthForRecording = used_ray_trace ? ray_optical_path : 
                                             (path_len + (delay_ns * c_light / 1e9));
      recorder_(xg, yg, zg, lambdaNm, z0_above_sea_level, dx_out, dy_out, dz_out,
            dx, dy, dz,
            totalTime_ns, weightPerBunch, n_local, pathLengthForRecording);
    }

    return corsika::ProcessReturn::Ok;
  }

} // namespace corsika::cherenkov
