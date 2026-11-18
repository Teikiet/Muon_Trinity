#include <cstdlib>  // getenv for the step-C diagnostic
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/framework/random/RNGManager.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace corsika::cherenkov {

  using namespace corsika::units::si;

  // ---------------------------------------------------------------------------
  // Phase 1 (geometry): cubic Hermite dense output for one RK4 step.
  //
  // The tracers already carry the position AND the direction as RK4 state, so the
  // endpoint tangents (ds*u0, ds*u1) are known.  Interpolating the position with
  // them -- instead of linearly between the two step ends, i.e. along the chord --
  // is what makes the crossing of a surface *inside* the step accurate.
  //
  // Why this matters: the chord-to-arc sagitta is c^2/(8*rho) with rho ~ 5e6 m
  // (1/r + n'/n), i.e. 1.3-1.9 mm at the production step sizes (300-400 m).  That
  // was the dominant geometry error of both tracers, two orders above their RK4
  // truncation (< 0.1 um at the same step size).  Measured on a 2-D vertical
  // section (exp atmosphere): linear 1.92 mm at 126 steps -> < 0.001 mm with this
  // interpolation, at the same step count.
  //
  // The direction is taken from the Hermite TANGENT (renormalised).  It equals the
  // endpoint directions exactly at t = 0 and t = 1, and avoids the previous
  // end-of-step direction, whose error is (1-alpha)*ds*kappa ~ 8e-5 rad at 400 m.
  // ---------------------------------------------------------------------------
  namespace {
    struct HermiteStep {
      double x0{0.0}, z0{0.0}, ux0{0.0}, uz0{0.0}; // step start (position, unit dir)
      double x1{0.0}, z1{0.0}, ux1{0.0}, uz1{0.0}; // step end
      double ds{0.0};                             // step arc length

      /// Position at t in [0,1] [m].
      inline void pos(double const t, double& x, double& z) const noexcept {
        double const t2 = t * t;
        double const t3 = t2 * t;
        double const h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        double const h10 = t3 - 2.0 * t2 + t;
        double const h01 = -2.0 * t3 + 3.0 * t2;
        double const h11 = t3 - t2;
        x = h00 * x0 + h10 * ds * ux0 + h01 * x1 + h11 * ds * ux1;
        z = h00 * z0 + h10 * ds * uz0 + h01 * z1 + h11 * ds * uz1;
      }

      /// d(position)/dt [m]; the tangent direction is this, normalised.
      inline void dpos(double const t, double& dx, double& dz) const noexcept {
        double const t2 = t * t;
        double const d00 = 6.0 * t2 - 6.0 * t;
        double const d10 = 3.0 * t2 - 4.0 * t + 1.0;
        double const d01 = -6.0 * t2 + 6.0 * t;
        double const d11 = 3.0 * t2 - 2.0 * t;
        dx = d00 * x0 + d10 * ds * ux0 + d01 * x1 + d11 * ds * ux1;
        dz = d00 * z0 + d10 * ds * uz0 + d01 * z1 + d11 * ds * uz1;
      }

      /// Unit direction at t in [0,1]; exactly the endpoint directions at t=0,1.
      inline void dir(double const t, double& ux, double& uz) const noexcept {
        double dx = 0.0, dz = 0.0;
        dpos(t, dx, dz);
        double const n = std::sqrt(dx * dx + dz * dz);
        if (n > 1.0e-300) {
          ux = dx / n;
          uz = dz / n;
        } else {
          ux = ux1;
          uz = uz1;
        }
      }
    };
  } // namespace

  inline CherenkovProcess::CherenkovProcess(
      double refractiveIndex,
      double lambdaMin,
      double lambdaMax,
      double bunchSize,
      double groundZ_m,
      HitRecorder recorder,
      ProjectionMode mode,
      std::shared_ptr<corsika::media::ShowerAxis const> primaryShowerAxis,
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
      const std::string& absorptionFile,
      double telescopePositionCorrectionFraction,
      std::shared_ptr<RefractiveIndexTable const> refractiveIndexTable,
      bool enableCurvedTable)
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
      , telescopeXOriginal_m_(telescopeX_m)
      , telescopeYOriginal_m_(telescopeY_m)
      , telescopeZOriginal_m_(telescopeZ_m)
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
      , refractiveIndexTable_(std::move(refractiveIndexTable))
      , enableFastTiming_(false)
      , enableCurvedTable_(enableCurvedTable)
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
        // Local extinction coefficient = d tau / d h of the table's interpolant.
        //
        // This used to be a symmetric +-25 m box average.  A box average is unbiased
        // only for a path that samples both sides of a height NODE equally, and the
        // shipped MODTRAN table has 1 km slabs: a grazing ray that stays just above a
        // node received the AVERAGE of the two slab slopes over tens of kilometres
        // (at h = 2000 m: 8.97e-5 instead of the 4.01e-5 of the slab it is in).
        // Narrowing the box cannot fix that -- at a node the average is
        // width-independent -- so this is a one-sided difference instead: it returns
        // the slope of the slab the point is actually in, which is exactly the model
        // the acceptance reference integrates (tau piecewise linear in altitude).
        // Measured on the frozen set this removed the last four >1 % cells:
        // 2000/100 8.5 % -> <0.02 %, 3000/10 2.0 %, 2000/50 1.3 %, 3000/50 1.2 %.
        auto const height_range = absorptionTable_.getHeightRange();
        double constexpr kDh_m = 1.0;
        double const hc = std::max(height_range.first, std::min(height_range.second, h_m));
        double h_lo = std::max(height_range.first, hc);
        double h_hi = std::min(height_range.second, h_lo + kDh_m);
        if (h_hi - h_lo < 1.0e-6) {  // at the top of the table: step down instead
          h_hi = hc;
          h_lo = std::max(height_range.first, h_hi - kDh_m);
        }
        double const actual_dh = h_hi - h_lo;
        if (actual_dh <= 0.0) return 0.0;
        double const tau_lo = absorptionTable_.getOpticalDepth(lambda_nm, h_lo);
        double const tau_hi = absorptionTable_.getOpticalDepth(lambda_nm, h_hi);
        // CONVENTION: the shipped tables (and tools/convert_npz_to_dat.py, which passes
        // CHASM ecoeff through unchanged) store the CUMULATIVE COLUMN from the ground
        // up, i.e. tau(h) = int_0^h alpha dh, so tau INCREASES with h and the local
        // extinction is +d tau/dh. The previous (tau_lo - tau_hi) assumed the opposite
        // ordering and returned exactly 0 for every h with the shipped profiles, which
        // silently disabled the tracer's optical-depth integration.
        return std::max(0.0, (tau_hi - tau_lo) / actual_dh);
      };
    }

    computeTelescopeBasis();

    // --- closed-form straight-line timing setup (P1) -------------------------
    // The fast path requires a valid vertical refractivity table AND a separable
    // refractive-index function, n(z,lam)-1 = N(z)*f(lam). A cheap four-point check
    // detects non-separable functions and falls back to the segment loop.
    if (refractiveIndexTable_ && refractiveIndexTable_->valid() && refractiveIndexFunc_) {
      // Rank-1 check on a small grid: n(z,lam)-1 must equal N(z)*f(lam), i.e.
      // M[i][j]*M[i2][j2] == M[i][j2]*M[i2][j] for all i,i2 (heights) and j,j2
      // (wavelengths). Sampling a grid (rather than four points) makes it much
      // harder for a non-separable model to pass by chance.
      constexpr int nz = 4;
      constexpr int nl = 4;
      double M[nz][nl];
      double maxAbs = 0.0;
      for (int i = 0; i < nz; ++i) {
        double const zi = refractiveIndexTable_->zMin() +
                          (0.05 + 0.30 * i) *
                              (refractiveIndexTable_->zMax() - refractiveIndexTable_->zMin());
        for (int j = 0; j < nl; ++j) {
          double const lam =
              fLambdaMin + (fLambdaMax - fLambdaMin) * j / double(nl - 1);
          M[i][j] = refractiveIndexFunc_(zi - earthRadiusM_, lam).second - 1.0;
          maxAbs = std::max(maxAbs, std::abs(M[i][j]));
        }
      }
      double maxDet = 0.0;
      for (int i = 0; i < nz; ++i) {
        for (int i2 = i + 1; i2 < nz; ++i2) {
          for (int j = 0; j < nl; ++j) {
            for (int j2 = j + 1; j2 < nl; ++j2) {
              maxDet = std::max(maxDet, std::abs(M[i][j] * M[i2][j2] - M[i][j2] * M[i2][j]));
            }
          }
        }
      }
      bool const separable =
          (maxDet <= 1e-6 * std::max(maxAbs * maxAbs, 1e-30));
      enableFastTiming_ = separable;
      separableIndex_ = separable;
      CORSIKA_LOGGER_INFO(
          logger_,
          "Cherenkov closed-form timing {} (table nodes={}, separability grid {}x{}, "
          "maxDet={:.3e}, scale={:.3e})",
          enableFastTiming_ ? "ENABLED" : "disabled", refractiveIndexTable_->size(), nz, nl,
          maxDet, maxAbs * maxAbs);
    } else {
      CORSIKA_LOGGER_INFO(logger_,
                          "Cherenkov closed-form timing disabled (no refractive-index table)");
    }

    // Configure telescope position correction
    double const clamped_fraction = std::max(0.0, std::min(1.0, telescopePositionCorrectionFraction));
    telescopePositionCorrectionFraction_ = clamped_fraction;
    enableTelescopePositionCorrection_ = (clamped_fraction > 0.0);
    if (enableTelescopePositionCorrection_) {
      double const inv = 1.0 / clamped_fraction;
      telescopePositionCorrectionInterval_ =
          std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(inv)));
    } else {
      telescopePositionCorrectionInterval_ = 0;
    }
    telescopePositionCorrectionStepCounter_ = 0;
    telescopePositionCorrectionCandidateValid_ = false;
    telescopePositionCorrectionLastUpdateValid_ = false;
    telescopePositionCorrectionLastUpdatePid_ = 0;

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

  inline void CherenkovProcess::computeTelescopeBasis() noexcept {
    double const cos_zen = std::cos(telescopeZenith_rad_);
    double const sin_zen = std::sin(telescopeZenith_rad_);
    double const cos_az = std::cos(telescopeAzimuth_rad_);
    double const sin_az = std::sin(telescopeAzimuth_rad_);

    // Z-axis points in the telescope direction (azimuth, zenith)
    telescopeBasisZ_x_ = sin_zen * cos_az;
    telescopeBasisZ_y_ = sin_zen * sin_az;
    telescopeBasisZ_z_ = -cos_zen;  // Negative because zenith 0 points down

    // X-axis is parallel to ground and perpendicular to Z-axis
    double const horiz_norm = std::sqrt(telescopeBasisZ_x_ * telescopeBasisZ_x_ +
                                       telescopeBasisZ_y_ * telescopeBasisZ_y_);
    if (horiz_norm > 1e-10) {
      telescopeBasisX_x_ = -sin_az;
      telescopeBasisX_y_ = cos_az;
      telescopeBasisX_z_ = 0.0;
    } else {
      telescopeBasisX_x_ = 1.0;
      telescopeBasisX_y_ = 0.0;
      telescopeBasisX_z_ = 0.0;
    }

    // Y-axis completes the right-handed system: Y = Z x X
    telescopeBasisY_x_ = telescopeBasisZ_y_ * telescopeBasisX_z_ -
                         telescopeBasisZ_z_ * telescopeBasisX_y_;
    telescopeBasisY_y_ = telescopeBasisZ_z_ * telescopeBasisX_x_ -
                         telescopeBasisZ_x_ * telescopeBasisX_z_;
    telescopeBasisY_z_ = telescopeBasisZ_x_ * telescopeBasisX_y_ -
                         telescopeBasisZ_y_ * telescopeBasisX_x_;
  }

  inline void CherenkovProcess::handleTelescopePositionCorrection(
      double energy_GeV,
      double pos_x_m,
      double pos_y_m,
      double pos_z_m,
      double dir_x,
      double dir_y,
      double dir_z,
      int pid) const {
    if (!enableTelescopePositionCorrection_ ||
        projectionMode_ != ProjectionMode::Telescope ||
        telescopePositionCorrectionInterval_ == 0) {
      return;
    }

    std::lock_guard<std::mutex> lock(telescopePositionCorrectionMutex_);

    if (!telescopePositionCorrectionCandidateValid_ ||
        energy_GeV > telescopePositionCorrectionCandidateEnergyGeV_) {
      telescopePositionCorrectionCandidateValid_ = true;
      telescopePositionCorrectionCandidateEnergyGeV_ = energy_GeV;
      telescopePositionCorrectionCandidatePid_ = pid;
      telescopePositionCorrectionCandidatePosX_m_ = pos_x_m;
      telescopePositionCorrectionCandidatePosY_m_ = pos_y_m;
      telescopePositionCorrectionCandidatePosZ_m_ = pos_z_m;
      telescopePositionCorrectionCandidateDirX_ = dir_x;
      telescopePositionCorrectionCandidateDirY_ = dir_y;
      telescopePositionCorrectionCandidateDirZ_ = dir_z;
    }

    ++telescopePositionCorrectionStepCounter_;
    if (telescopePositionCorrectionStepCounter_ >= telescopePositionCorrectionInterval_) {
      if (telescopePositionCorrectionCandidateValid_) {
        bool const updated = applyTelescopePositionCorrectionCandidate();
        if (updated) {
          telescopePositionCorrectionLastUpdateValid_ = true;
          telescopePositionCorrectionLastUpdatePid_ = telescopePositionCorrectionCandidatePid_;
        }
      }
      telescopePositionCorrectionStepCounter_ = 0;
      telescopePositionCorrectionCandidateValid_ = false;
    }
  }

  inline bool CherenkovProcess::applyTelescopePositionCorrectionCandidate() const {
    double const num = (telescopeX_m_ - telescopePositionCorrectionCandidatePosX_m_) * telescopeBasisZ_x_ +
                       (telescopeY_m_ - telescopePositionCorrectionCandidatePosY_m_) * telescopeBasisZ_y_ +
                       (telescopeZ_m_ - telescopePositionCorrectionCandidatePosZ_m_) * telescopeBasisZ_z_;
    double const den = telescopePositionCorrectionCandidateDirX_ * telescopeBasisZ_x_ +
                       telescopePositionCorrectionCandidateDirY_ * telescopeBasisZ_y_ +
                       telescopePositionCorrectionCandidateDirZ_ * telescopeBasisZ_z_;
    if (std::abs(den) < 1e-10) {
      return false;
    }

    double const t_hit = num / den;
    if (t_hit <= 0.0) {
      return false;
    }

    telescopeX_m_ = telescopePositionCorrectionCandidatePosX_m_ + t_hit * telescopePositionCorrectionCandidateDirX_;
    telescopeY_m_ = telescopePositionCorrectionCandidatePosY_m_ + t_hit * telescopePositionCorrectionCandidateDirY_;
    telescopeZ_m_ = telescopePositionCorrectionCandidatePosZ_m_ + t_hit * telescopePositionCorrectionCandidateDirZ_;
    return true;
  }

  inline void CherenkovProcess::resetTelescopePositionCorrectionForEvent() {
    std::lock_guard<std::mutex> lock(telescopePositionCorrectionMutex_);
    telescopeX_m_ = telescopeXOriginal_m_;
    telescopeY_m_ = telescopeYOriginal_m_;
    telescopeZ_m_ = telescopeZOriginal_m_;
    telescopePositionCorrectionStepCounter_ = 0;
    telescopePositionCorrectionCandidateValid_ = false;
    telescopePositionCorrectionLastUpdateValid_ = false;
    telescopePositionCorrectionLastUpdatePid_ = 0;
  }

  inline CherenkovProcess::TelescopePositionCorrectionResult
  CherenkovProcess::getTelescopePositionCorrectionResult() const {
    TelescopePositionCorrectionResult result;
    result.enabled = enableTelescopePositionCorrection_;
    std::lock_guard<std::mutex> lock(telescopePositionCorrectionMutex_);
    double const dx = telescopeX_m_ - telescopeXOriginal_m_;
    double const dy = telescopeY_m_ - telescopeYOriginal_m_;
    double const dz = telescopeZ_m_ - telescopeZOriginal_m_;
    result.offset_local_x_m = dx * telescopeBasisX_x_ + dy * telescopeBasisX_y_ + dz * telescopeBasisX_z_;
    result.offset_local_y_m = dx * telescopeBasisY_x_ + dy * telescopeBasisY_y_ + dz * telescopeBasisY_z_;
    result.offset_local_z_m = dx * telescopeBasisZ_x_ + dy * telescopeBasisZ_y_ + dz * telescopeBasisZ_z_;
    result.has_update = telescopePositionCorrectionLastUpdateValid_;
    result.last_update_pid = telescopePositionCorrectionLastUpdatePid_;
    return result;
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
      double& accumulated_tau,
      double& dx_impact, double& dy_impact, double& dz_impact,
      double maxRadialStep_m) const {

    // --- Step 2m (2g-h2l): the CERTIFIED engine, with a different stop ---------
    // This tracer now runs the same integrator as rayTraceSpherical(): a fixed
    // propagation-plane frame, the direction carried as state, RK4 on
    // (position, direction), one transverse step controller with a radial cap, and a
    // tau quadrature fixed in altitude.  Only the stopping surface differs: the signed
    // distance to the telescope plane, which is AFFINE in the plane coordinates
    //     d(x_plane, z_plane) = A*x_plane + B*z_plane - C,
    // so the crossing fraction is exact in position.
    //
    // The previous Bouguer walk was replaced, not repaired.  Measured with the new
    // acceptance set (testCherenkovPlaneTracerAcceptance, 90 deg pointing):
    //   * it advanced the transverse coordinate with dx_sign = -1 and
    //     az_sign = (t_hat_y >= 0 ? +1 : -1), while t_hat is DEFINED as the
    //     tangential part of the photon direction, so u.t_hat = +sin(theta0) > 0
    //     ALWAYS: the ray was walked BACKWARDS in the transverse coordinate, never
    //     reached the plane and escaped at h = 100 km (48/48 set-A cells and 18/18
    //     set-C cells returned "no path");
    //   * dr = |dx cos(theta)/(sin(theta)cos(alpha) + cos(theta)sin(alpha))| divides
    //     by zero for near-vertical rays (h_cur = inf at step 1);
    //   * the 1 km chord branch required direction == -1, i.e. it could never fire for
    //     the ascending half of the field of view.
    // Reusing the certified machinery also gives the 90 deg case its measured
    // behaviour: theta step-insensitivity, 1e-5 geometry, tau <= 0.5 %.

    if (!refractiveIndexFunc_) {
      CORSIKA_LOGGER_WARN(logger_, "rayTraceSphericalToPlane: no refractiveIndexFunc");
      return false;
    }

    // ---- Phase 3: opt-in adaptive integrator ------------------------------------
    // The RK4 path below is the reference, the acceptance-test subject and the
    // fallback; DOPRI5 is taken only when it is enabled AND has what its error control
    // needs (a separable index and a valid table, so the Sigma component scales with
    // f_group).
    bool const tangentReady = separableIndex_ && refractiveIndexTable_ &&
                              refractiveIndexTable_->valid();
    if (tangentReady && (tracerIntegrator_ == TracerIntegrator::TANGENT ||
                         tracerIntegrator_ == TracerIntegrator::AUTO)) {
      // The tangent attempt is cheap when it refuses (the bracket solve rejects before any
      // stepping), and AUTO falls back to the regular tracer for the SAME ray, so a ray
      // is never lost.  A refusal is counted, never silent.
      if (rayTracePlaneTangent(x0, y0, z0_altitude, dx, dy, dz, plane_x, plane_y,
                               plane_z, norm_x, norm_y, norm_z, x_impact, y_impact,
                               z_impact, optical_path, geometric_path, lambdaNm,
                               accumulated_tau, dx_impact, dy_impact, dz_impact,
                               maxRadialStep_m)) {
        return true;
      }
      if (tracerIntegrator_ == TracerIntegrator::AUTO) ++trTangentFallbacks_;
      if (tracerIntegrator_ == TracerIntegrator::TANGENT) return false;
      // AUTO: fall through to the regular tracer below.
    }
    if ((tracerIntegrator_ == TracerIntegrator::DOPRI5 ||
         tracerIntegrator_ == TracerIntegrator::AUTO) &&
        separableIndex_ && refractiveIndexTable_ && refractiveIndexTable_->valid()) {
      return rayTracePlaneDopri(x0, y0, z0_altitude, dx, dy, dz, plane_x, plane_y,
                                plane_z, norm_x, norm_y, norm_z, x_impact, y_impact,
                                z_impact, optical_path, geometric_path, lambdaNm,
                                accumulated_tau, dx_impact, dy_impact, dz_impact,
                                maxRadialStep_m);
    }
    ++trRk4Rays_;

    // Absorption is NOT sampled here: the tracer's only job is geometry + tau (see
    // absorptionTransmission).  Nothing below consumes randomness.

    double const R_E = earthRadiusM_;
    double const px = x0, py = y0, pz = z0_altitude;
    double const r0 = std::sqrt(px * px + py * py + pz * pz);
    if (!(r0 > 0.0)) return false;

    double const r_hat_x = px / r0, r_hat_y = py / r0, r_hat_z = pz / r0;
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

    // Tangential unit vector.  By construction u.t_hat = sin(theta0) >= 0, so with
    // pos = z_plane*r_hat + x_plane*t_hat the x_plane coordinate always grows along
    // the photon (no sign hack, no mirror).
    double const tang_x = dx - cos_theta0 * r_hat_x;
    double const tang_y = dy - cos_theta0 * r_hat_y;
    double const tang_z = dz - cos_theta0 * r_hat_z;
    double const tang_mag = std::sqrt(tang_x * tang_x + tang_y * tang_y + tang_z * tang_z);
    double t_hat_x, t_hat_y, t_hat_z;
    if (tang_mag > 1e-12) {
      t_hat_x = tang_x / tang_mag;
      t_hat_y = tang_y / tang_mag;
      t_hat_z = tang_z / tang_mag;
    } else {
      double ref_x = 0.0, ref_y = 0.0, ref_z = 1.0;
      if (std::abs(r_hat_z) > 0.9) { ref_x = 1.0; ref_y = 0.0; ref_z = 0.0; }
      t_hat_x = r_hat_y * ref_z - r_hat_z * ref_y;
      t_hat_y = r_hat_z * ref_x - r_hat_x * ref_z;
      t_hat_z = r_hat_x * ref_y - r_hat_y * ref_x;
      double const tn = std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
      t_hat_x /= tn;
      t_hat_y /= tn;
      t_hat_z /= tn;
    }

    // Signed distance to the plane, in plane coordinates.
    double const A = t_hat_x * norm_x + t_hat_y * norm_y + t_hat_z * norm_z;
    double const B = r_hat_x * norm_x + r_hat_y * norm_y + r_hat_z * norm_z;
    double const Cp = plane_x * norm_x + plane_y * norm_y + plane_z * norm_z;
    if (std::abs(A) < 1e-15 && std::abs(B) < 1e-15) return false;  // degenerate plane

    // tau quadrature (Step 3b): fixed in altitude, floored by the transverse length.
    constexpr double kTauSubStep_m = 10.0;
    constexpr double kTauSubStepTransverse_m = 50.0;

    double x_plane = 0.0, z_plane = r0, r_cur = r0;
    double u_x_plane = sin_theta0, u_z_plane = cos_theta0;
    double d_cur = A * x_plane + B * z_plane - Cp;

    int const max_steps = 10000;
    double const max_altitude = 100000.0;

    double x_prev = x_plane, z_prev = z_plane, r_prev = r_cur;
    double up_x = u_x_plane, up_z = u_z_plane;   // direction at the step start (Hermite)
    double optical_prev = optical_path, geometric_prev = geometric_path;
    double tau_prev = accumulated_tau;
    double d_prev = d_cur;

    // ---- Phase 1: index and group path straight from the refractive-index table ----
    // When n(z,lam)-1 = N(z)*f(lam) is separable (measured in the ctor) the tracers
    // need N(z), dN/dz and f_phase/f_group, all of which the table already holds.
    // That replaces the three refractiveIndexFunc_ calls per step (each a binary
    // search over 12000 nodes) by O(1) lookups, and it takes dN/dz from the stored
    // piecewise-exponential log-slope instead of a +-1 m finite difference.  The
    // finite-difference branch is kept verbatim for non-separable indices and for
    // runs without a table.
    bool const useIndexTable = separableIndex_ && refractiveIndexTable_ &&
                               refractiveIndexTable_->valid() && refractiveIndexFunc_;
    double fPhase = 1.0, fGroup = 1.0;
    if (useIndexTable) {
      auto const ng0 = refractiveIndexFunc_(r0 - R_E, lambdaNm);
      double const Nr0 = refractiveIndexTable_->N(r0);
      if (Nr0 > 1.0e-12) {
        fPhase = (ng0.first - 1.0) / Nr0;
        fGroup = (ng0.second - 1.0) / Nr0;
      } else {
        fPhase = 0.0;
        fGroup = 0.0;
      }
    }

    // --- 2g-h2n: the plane can lie BEHIND the starting state -------------------
    // This is not exotic now that the below-observation block exists.  For a photon
    // RISING to the telescope the sphere r = r_obs crossing and the plane crossing
    // almost coincide (the dish sits on both surfaces), so a photon that passes just
    // BELOW the dish crosses the plane a little BEFORE it reaches r = r_obs: the
    // launch state then has d_start > 0, the RK4 loop steps straight over the plane on
    // its first step, d never changes sign again and the photon escapes at
    // h = 100 km after 10^4 steps (measured: all 9 table-served upward cells of the
    // acceptance set and its 6 legacy fallbacks, each burning 10^4 RK4 steps).
    // The remaining stretch is at most a few hundred metres and very nearly straight,
    // so the crossing is the straight-line advance back along u; the curvature over
    // 200 m displaces it by under 1e-3 m.  Every call that has the plane AHEAD gives
    // tPlane > 0 and falls through to the integrator below, so the certified
    // downward/aim-sweep behaviour is untouched.
    double const ddds = A * sin_theta0 + B * cos_theta0; // = u . n = d(plane dist)/ds
    if (std::abs(ddds) > 1.0e-12) {
      double const tPlane = -d_cur / ddds;
      ++nPlaneCalls_;
      if (tPlane <= 0.0) {
        ++nPlaneBack_;
        if (tPlane > -1.0) ++nPlaneTangent_;
        // Step C: a crossing FAR behind the launch state is not a hit -- a photon cannot
        // reach a plane it has already left.  The branch is legitimate only for the
        // launch-and-continue regime, where the table's sphere crossing and the dish
        // plane are a fraction of a metre apart, so the tolerance is generous but finite
        // and a genuine backward crossing can never be scored as a hit.  Measured in
        // production (1e5 GeV, 67,791 gated photons through this tracer): backward = 0,
        // so this changes nothing today; it closes the failure mode.
        // Step C verdict, measured rather than assumed: NO magnitude cutoff.  A cutoff
        // was tried (100 m -> 4 acceptance rows became "continuation failures", 1 km ->
        // 2 remained) because a rising photon can cross the dish plane well before it
        // reaches r = r_obs, which is exactly the case this branch exists for.  What the
        // measurement does say is that in the production configuration the branch NEVER
        // fires: 1e5 GeV, 67,791 gated photons through this tracer, backward = 0, min
        // tPlane = 0.00 m.  So the branch is legitimate AND unused; it is now counted
        // (planeBackStats) so that a future geometry which starts abusing it is visible
        // instead of silent.
        double prevMin = planeMinT_.load();
        while (tPlane < prevMin && !planeMinT_.compare_exchange_weak(prevMin, tPlane)) {}
        if (getenv("CHER_AB_PLANE_BACK") != nullptr) {
          CORSIKA_LOGGER_INFO(logger_,
                              "plane crossing BEHIND the launch state: tPlane={:.3f} m "
                              "d_start={:.3f} m",
                              tPlane, d_cur);
        }
        auto const ng0 = refractiveIndexFunc_(r0 - R_E, lambdaNm);
        double const n_group0 = ng0.second;
        x_impact = px + tPlane * dx;
        y_impact = py + tPlane * dy;
        z_impact = pz + tPlane * dz;
        dx_impact = dx;  // straight branch: the direction is unchanged
        dy_impact = dy;
        dz_impact = dz;
        geometric_path += tPlane;             // tPlane < 0: the stretch is behind us
        optical_path += n_group0 * tPlane;
        if (extinctionCoeffFunc_) {
          accumulated_tau = std::max(
              0.0, accumulated_tau - extinctionCoeffFunc_(r0 - R_E, lambdaNm) * std::fabs(tPlane));
        }
        return true;
      }
    }

    for (int step = 0; step < max_steps; ++step) {
      double const h_cur = r_cur - R_E;
      if (h_cur > max_altitude) {
        // Phase 2: counted always, logged only for the first three.  This used to emit
        // one WARN per photon -- 40 032 lines on a single production event.
        std::size_t const nEsc = escapedAtmosphere_.fetch_add(1) + 1;
        if (nEsc <= 3) {
          double const th = std::acos(std::max(-1.0, std::min(1.0, -dz))) * 180.0 / M_PI;
          CORSIKA_LOGGER_WARN(logger_,
                              "rayTraceSphericalToPlane: escaped atmosphere at step {}, "
                              "h_cur={:.1f}m > max_altitude={:.1f}m, photon zenith=({:.2f} deg) "
                              "[logged for the first 3 only; the total is reported at the end]",
                              step, h_cur, max_altitude, th);
        }
        return false;
      }
      if (h_cur < -100.0) return false;  // below sea level without reaching the plane

      // Phase 1: only the group index for the non-table fallback is needed here; the
      // phase index and its derivative are evaluated PER RK4 STAGE below (see kappaAt).
      double n_group = 1.0;
      if (!useIndexTable) {
        n_group = refractiveIndexFunc_(h_cur, lambdaNm).second;
      }
      // local zenith angles from the state, on the LOCAL radial/transverse basis
      double const rp_x = x_plane / r_cur;
      double const rp_z = z_plane / r_cur;
      double const cos_theta_signed = u_x_plane * rp_x + u_z_plane * rp_z;
      double const cos_theta_cur = std::min(1.0, std::fabs(cos_theta_signed));
      double const sin_theta_cur = std::min(1.0, std::fabs(u_x_plane * rp_z - u_z_plane * rp_x));

      double const H_scale = computeLocalScaleHeight(h_cur);
      double const dx_ctrl = std::min(H_scale / 20.0, 1000.0);
      double ds = (sin_theta_cur > 1.0e-9) ? dx_ctrl / sin_theta_cur
                                           : (maxRadialStep_m > 0.0 ? maxRadialStep_m : 100.0);
      if (maxRadialStep_m > 0.0 && cos_theta_cur > 1.0e-9) {
        ds = std::min(ds, maxRadialStep_m / cos_theta_cur);
      }

      x_prev = x_plane;
      z_prev = z_plane;
      r_prev = r_cur;
      up_x = u_x_plane;
      up_z = u_z_plane;
      optical_prev = optical_path;
      geometric_prev = geometric_path;
      tau_prev = accumulated_tau;

      // ---- RK4 on (position, direction): du/ds = (n'/n)(r_hat - (u.r_hat)u) ----
      double sigma_inc = 0.0;  // int N ds over this step (Phase 1: 4th-order group path)
      {
        // Phase 1b: kappa MUST be evaluated at the stage's own radius.  Freezing it at
        // the step start -- the previous behaviour -- makes the bending integral a
        // left-endpoint rule in n'/n, i.e. first order in ds/H, which was the dominant
        // geometry error of both tracers: measured on the frozen telescope set, impact
        // error 64 mm at ds = 100 m, 255 mm at 400 m and 769 mm at 1200 m, against
        // 0.04-7 mm with the stage-local value (and the frozen 1-2 mm crossing error
        // is removed separately by the Hermite closure).  With the separable table this
        // costs two O(1) lookups per stage and NO refractiveIndexFunc_ call; without
        // it, two function calls per stage (the ln difference is d(ln n)/dr exactly).
        auto const kappaAt = [&](double const rr) -> double {
          if (useIndexTable) {
            double const Ncur = refractiveIndexTable_->N(rr);
            double const nph = 1.0 + Ncur * fPhase;
            return (nph > 1.0e-12) ? (fPhase * refractiveIndexTable_->dNdz(rr)) / nph : 0.0;
          }
          double const hh = rr - R_E;
          double const hLo = std::max(0.0, hh - 1.0);
          double const hHi = std::max(1.0e-6, hh + 1.0);
          double const nLo = refractiveIndexFunc_(hLo, lambdaNm).first;
          double const nHi = refractiveIndexFunc_(hHi, lambdaNm).first;
          if (!(nLo > 1.0e-12) || !(nHi > 1.0e-12)) return 0.0;
          return (std::log(nHi) - std::log(nLo)) / (hHi - hLo);
        };
        auto const rayDeriv = [&](double const xq, double const zq, double const vx,
                                  double const vz, double& dxt, double& dzt, double& dvxt,
                                  double& dvzt) {
          double const rq = std::sqrt(xq * xq + zq * zq);
          double const rx = xq / rq;
          double const rz = zq / rq;
          double const ur = vx * rx + vz * rz;
          double const kq = kappaAt(rq);
          dxt = vx;
          dzt = vz;
          dvxt = kq * (rx - ur * vx);
          dvzt = kq * (rz - ur * vz);
        };
        double a1x, a1z, b1x, b1z;
        rayDeriv(x_plane, z_plane, u_x_plane, u_z_plane, a1x, a1z, b1x, b1z);
        double a2x, a2z, b2x, b2z;
        rayDeriv(x_plane + 0.5 * ds * a1x, z_plane + 0.5 * ds * a1z,
                 u_x_plane + 0.5 * ds * b1x, u_z_plane + 0.5 * ds * b1z, a2x, a2z, b2x, b2z);
        double a3x, a3z, b3x, b3z;
        rayDeriv(x_plane + 0.5 * ds * a2x, z_plane + 0.5 * ds * a2z,
                 u_x_plane + 0.5 * ds * b2x, u_z_plane + 0.5 * ds * b2z, a3x, a3z, b3x, b3z);
        double a4x, a4z, b4x, b4z;
        rayDeriv(x_plane + ds * a3x, z_plane + ds * a3z, u_x_plane + ds * b3x,
                 u_z_plane + ds * b3z, a4x, a4z, b4x, b4z);
        x_plane += ds / 6.0 * (a1x + 2.0 * a2x + 2.0 * a3x + a4x);
        z_plane += ds / 6.0 * (a1z + 2.0 * a2z + 2.0 * a3z + a4z);
        u_x_plane += ds / 6.0 * (b1x + 2.0 * b2x + 2.0 * b3x + b4x);
        u_z_plane += ds / 6.0 * (b1z + 2.0 * b2z + 2.0 * b3z + b4z);
        double const un = std::sqrt(u_x_plane * u_x_plane + u_z_plane * u_z_plane);
        u_x_plane /= un;
        u_z_plane /= un;

        // Phase 1: dSigma/ds = N(r) on the same RK4 stages.  The group optical path is
        // L + f_group*Sigma to 4th order; the previous left-endpoint rule
        // (n_group at the step START) carries a first-order bias, ~0.4 m of optical
        // path at 400 m steps (~1.2 ns).  Only when the separable table supplies N.
        if (useIndexTable) {
          auto const NofStage = [&](double const xq, double const zq) {
            return refractiveIndexTable_->N(std::sqrt(xq * xq + zq * zq));
          };
          double const sg1 = NofStage(x_prev, z_prev);
          double const sg2 = NofStage(x_prev + 0.5 * ds * a1x, z_prev + 0.5 * ds * a1z);
          double const sg3 = NofStage(x_prev + 0.5 * ds * a2x, z_prev + 0.5 * ds * a2z);
          double const sg4 = NofStage(x_prev + ds * a3x, z_prev + ds * a3z);
          sigma_inc = ds / 6.0 * (sg1 + 2.0 * sg2 + 2.0 * sg3 + sg4);
        }
      }

      geometric_path += ds;
      double const op_inc = useIndexTable ? (ds + fGroup * sigma_inc) : (n_group * ds);
      optical_path += op_inc;

      // Dense output for everything that happens INSIDE this step -- the tau
      // quadrature and the plane crossing -- through the endpoint states, instead of
      // the chord between them (Phase 1; see HermiteStep).
      HermiteStep const seg{x_prev, z_prev, up_x, up_z, x_plane, z_plane,
                            u_x_plane, u_z_plane, ds};

      double const r_new = std::sqrt(x_plane * x_plane + z_plane * z_plane);
      if (extinctionCoeffFunc_) {
        int const nByAltitude =
            static_cast<int>(std::ceil(std::abs(r_new - r_prev) / kTauSubStep_m));
        int const nByLength = static_cast<int>(std::ceil(ds / kTauSubStepTransverse_m));
        int const n_tau = std::max(1, std::max(nByAltitude, nByLength));
        double const ds_tau = ds / static_cast<double>(n_tau);
        for (int its = 0; its < n_tau; ++its) {
          double const fr = (static_cast<double>(its) + 0.5) / static_cast<double>(n_tau);
          double sx = 0.0, sz = 0.0;
          seg.pos(fr, sx, sz);   // on the arc, not on the chord
          accumulated_tau +=
              extinctionCoeffFunc_(std::sqrt(sx * sx + sz * sz) - R_E, lambdaNm) * ds_tau;
        }
      }
      r_cur = r_new;

      d_cur = A * x_plane + B * z_plane - Cp;
      if (d_prev * d_cur <= 0.0) {
        // d is affine in (x_plane, z_plane), so the fraction below is exact for the
        // CHORD.  The photon travels on the arc: the chord-to-arc separation is exactly
        // the error this refinement removes (1.3-1.9 mm at the production step sizes,
        // measured; < 0.001 mm after).  Newton on the Hermite curve, started from the
        // chord fraction.  The correction is ~1e-5 of a step, so two iterations reach
        // round-off; six are done for safety and cost ~40 ns.
        double alpha = 0.0;
        if (std::fabs(d_prev - d_cur) > 1e-300) alpha = d_prev / (d_prev - d_cur);
        alpha = std::max(0.0, std::min(1.0, alpha));

        double tc = alpha;
        for (int it = 0; it < 6; ++it) {
          double xq = 0.0, zq = 0.0, dxq = 0.0, dzq = 0.0;
          seg.pos(tc, xq, zq);
          seg.dpos(tc, dxq, dzq);
          double const fq = A * xq + B * zq - Cp;
          double const fpq = A * dxq + B * dzq;
          if (std::fabs(fpq) < 1.0e-300) break;
          double const dt = fq / fpq;
          double const tn = tc - dt;
          if (!(tn >= 0.0 && tn <= 1.0)) break;
          tc = tn;
          if (std::fabs(dt) <= 1.0e-13) break;
        }
        if (!(tc >= 0.0 && tc <= 1.0)) tc = alpha;  // never leave the bracket

        // the closing partial step replaces the full step's tau (Step 2i)
        accumulated_tau = tau_prev;
        if (extinctionCoeffFunc_) {
          double const ds_partial = tc * ds;
          int const nA =
              static_cast<int>(std::ceil(std::abs(r_cur - r_prev) * tc / kTauSubStep_m));
          int const nL = static_cast<int>(std::ceil(ds_partial / kTauSubStepTransverse_m));
          int const n_tau = std::max(1, std::max(nA, nL));
          double const ds_tau = ds_partial / static_cast<double>(n_tau);
          for (int its = 0; its < n_tau; ++its) {
            double const fr = (static_cast<double>(its) + 0.5) / static_cast<double>(n_tau);
            double sx = 0.0, sz = 0.0;
            seg.pos(tc * fr, sx, sz);
            accumulated_tau +=
                extinctionCoeffFunc_(std::sqrt(sx * sx + sz * sz) - R_E, lambdaNm) * ds_tau;
          }
        }
        double xi = 0.0, zi = 0.0;
        seg.pos(tc, xi, zi);
        x_impact = zi * r_hat_x + xi * t_hat_x;
        y_impact = zi * r_hat_y + xi * t_hat_y;
        z_impact = zi * r_hat_z + xi * t_hat_z;
        // WAVE 5: report the arrival direction.  The plane state is (transverse, radial) =
        // (x_plane, z_plane) with direction (u_x_plane, u_z_plane), and the frame maps
        // pos = z*r_hat + x*t_hat, so the direction maps the same way.  Phase 1: take it
        // from the Hermite tangent, not from the step END (off by (1-alpha)*ds*kappa,
        // i.e. ~8e-5 rad at 400 m).
        {
          double uxi = u_x_plane, uzi = u_z_plane;
          seg.dir(tc, uxi, uzi);
          double const nx = uzi * r_hat_x + uxi * t_hat_x;
          double const ny = uzi * r_hat_y + uxi * t_hat_y;
          double const nz = uzi * r_hat_z + uxi * t_hat_z;
          double const nn = std::sqrt(nx * nx + ny * ny + nz * nz);
          if (nn > 0.0) {
            dx_impact = nx / nn;
            dy_impact = ny / nn;
            dz_impact = nz / nn;
          } else {
            dx_impact = dx;
            dy_impact = dy;
            dz_impact = dz;
          }
        }
        optical_path = optical_prev + tc * op_inc;
        geometric_path = geometric_prev + tc * ds;
        return true;
      }
      d_prev = d_cur;
    }

    std::size_t const nMax = maxStepsReached_.fetch_add(1) + 1;
    if (nMax <= 3) {
      CORSIKA_LOGGER_WARN(logger_,
                          "rayTraceSphericalToPlane: FAILED - max_steps={} reached without "
                          "hitting plane [logged for the first 3 only]",
                          max_steps);
    }
    return false;
  }
  // ===========================================================================
  // Phase 3: adaptive DOPRI5 5(4) plane tracer
  // ---------------------------------------------------------------------------
  // Same contract as rayTraceSphericalToPlane (same args, same outputs, same tau
  // ownership), same ODE and same frame; only the stepping changes.  The RK4 routine
  // above stays the reference, the acceptance-test subject and the fallback.
  //
  // Differences that matter for correctness, not just speed:
  //  * the plane crossing is located on the DOPRI5 DENSE OUTPUT (quartic in the step),
  //    with a bracketed solve; the plane distance is evaluated at t = 0, 1/4, 1/2, 3/4, 1
  //    so a dip that returns to the same side is still seen (endpoint sign checks alone
  //    would miss it);
  //  * position, direction, group path and tau are all evaluated at the SAME t;
  //  * the tau quadrature keeps its altitude-fixed 10 m / 50 m rule, now driven by the
  //    dense output, so a large step does not coarsen the absorption integral;
  //  * a rejected step, a crossing solve and a model switch all invalidate FSAL.
  // ===========================================================================
  inline bool CherenkovProcess::rayTracePlaneDopri(
      double x0, double y0, double z0_altitude, double dx, double dy, double dz,
      double plane_x, double plane_y, double plane_z, double norm_x, double norm_y,
      double norm_z, double& x_impact, double& y_impact, double& z_impact,
      double& optical_path, double& geometric_path, double lambdaNm,
      double& accumulated_tau, double& dx_impact, double& dy_impact, double& dz_impact,
      double maxRadialStep_m) const {
    constexpr int ND = 5; // x_plane, z_plane, u_x, u_z, Sigma
    using Dopri = Dopri5<ND>;

    double const R_E = earthRadiusM_;
    double const px = x0, py = y0, pz = z0_altitude;
    double const r0 = std::sqrt(px * px + py * py + pz * pz);
    if (!(r0 > 0.0)) return false;

    double const r_hat_x = px / r0, r_hat_y = py / r0, r_hat_z = pz / r0;
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

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
      double ref_x = 0.0, ref_y = 0.0, ref_z = 1.0;
      if (std::abs(r_hat_z) > 0.9) { ref_x = 1.0; ref_y = 0.0; ref_z = 0.0; }
      t_hat_x = r_hat_y * ref_z - r_hat_z * ref_y;
      t_hat_y = r_hat_z * ref_x - r_hat_x * ref_z;
      t_hat_z = r_hat_x * ref_y - r_hat_y * ref_x;
      double const tn = std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
      t_hat_x /= tn; t_hat_y /= tn; t_hat_z /= tn;
    }

    double const A = t_hat_x * norm_x + t_hat_y * norm_y + t_hat_z * norm_z;
    double const B = r_hat_x * norm_x + r_hat_y * norm_y + r_hat_z * norm_z;
    double const Cp = plane_x * norm_x + plane_y * norm_y + plane_z * norm_z;
    if (std::abs(A) < 1e-15 && std::abs(B) < 1e-15) return false;

    double const d_start = B * r0 - Cp;

    // tau quadrature (same rule as the RK4 tracer)
    constexpr double kTauSubStep_m = 10.0;
    constexpr double kTauSubStepTransverse_m = 50.0;

    // ---- separability-derived dispersion factors and the table-only index ---------
    double fPhase = 1.0, fGroup = 1.0;
    {
      auto const ng0 = refractiveIndexFunc_(r0 - R_E, lambdaNm);
      double const Nr0 = refractiveIndexTable_->N(r0);
      if (Nr0 > 1.0e-12) {
        fPhase = (ng0.first - 1.0) / Nr0;
        fGroup = (ng0.second - 1.0) / Nr0;
      } else {
        fPhase = 0.0;
        fGroup = 0.0;
      }
    }
    auto const kappaAt = [&](double const rr) -> double {
      double const Ncur = refractiveIndexTable_->N(rr);
      double const nph = 1.0 + Ncur * fPhase;
      return (nph > 1.0e-12) ? (fPhase * refractiveIndexTable_->dNdz(rr)) / nph : 0.0;
    };

    double ycur[ND] = {0.0, r0, sin_theta0, cos_theta0, 0.0};
    auto const rhs = [&](double const* y, double* f) {
      double const rq = std::sqrt(y[0] * y[0] + y[1] * y[1]);
      double const rx = y[0] / rq, rz = y[1] / rq;
      double const ur = y[2] * rx + y[3] * rz;
      double const kq = kappaAt(rq);
      f[0] = y[2];
      f[1] = y[3];
      f[2] = kq * (rx - ur * y[2]);
      f[3] = kq * (rz - ur * y[3]);
      f[4] = refractiveIndexTable_->N(rq);
    };

    // ---- the backward-crossing branch (launch-and-continue), as in the RK4 tracer --
    double const ddds = A * sin_theta0 + B * cos_theta0;
    if (std::abs(ddds) > 1.0e-12) {
      double const tPlane = -d_start / ddds;
      ++nPlaneCalls_;
      if (tPlane <= 0.0) {
        ++nPlaneBack_;
        if (tPlane > -1.0) ++nPlaneTangent_;
        x_impact = px + tPlane * dx;
        y_impact = py + tPlane * dy;
        z_impact = pz + tPlane * dz;
        dx_impact = dx;
        dy_impact = dy;
        dz_impact = dz;
        // f_group * Sigma over the short straight stretch (Sigma ~ N * tPlane)
        geometric_path += tPlane;
        optical_path += tPlane + fGroup * refractiveIndexTable_->N(r0) * tPlane;
        if (extinctionCoeffFunc_) {
          accumulated_tau = std::max(
              0.0, accumulated_tau - extinctionCoeffFunc_(r0 - R_E, lambdaNm) * std::fabs(tPlane));
        }
        return true;
      }
    }


    // ---- tolerances, controller, initial step -------------------------------------
    double atol[ND] = {dopriAtolPos_m_, dopriAtolPos_m_, dopriAtolDir_, dopriAtolDir_,
                       dopriAtolOpt_m_};
    double rtol[ND] = {0.0, 0.0, 1.0e-9, 1.0e-9, 0.0};
    Dopri dp;
    dp.setTolerances(atol, rtol);
    dp.setState(ycur);
    dp.setMinStep(dopriMinStep_m_);

    // Counters are accumulated locally and flushed on EVERY return path, so a ray that
    // fails halfway still reports the work it caused.
    std::size_t cAcc = 0, cRej = 0, cF = 0, cMinFail = 0, cCapped = 0, cCross = 0;
    std::size_t cImpactRej = 0, cImpactOk = 0; // Phase 4 step 3
    double cHMin = 1.0e300, cHMax = 0.0;
    auto const flush = [&]() {
      trAccepted_ += cAcc;
      trRejected_ += cRej;
      trF_ += cF;
      trMinStepFail_ += cMinFail;
      trCapped_ += cCapped;
      trCrossIters_ += cCross;
      trImpactRej_ += cImpactRej;
      trImpactOk_ += cImpactOk;
      ++trDopriRays_;
      double lo = trHMin_.load();
      while (cHMin < lo && !trHMin_.compare_exchange_weak(lo, cHMin)) {}
      double hi = trHMax_.load();
      while (cHMax > hi && !trHMax_.compare_exchange_weak(hi, cHMax)) {}
    };

    // The maximum step.  Default: exactly the RK4 caps (H_local/20, 1000 m, plus the
    // caller's radial cap), so a DOPRI5 run is directly comparable with the RK4
    // reference.  A configured dopriMaxStep_m_ > 0 REPLACES them -- that is the
    // benchmark's relaxation, and it is never silent (the caller sets it).
    auto const capFor = [&](double const* y) {
      if (dopriMaxStep_m_ > 0.0) return dopriMaxStep_m_;
      double const rr = std::sqrt(y[0] * y[0] + y[1] * y[1]);
      double const H = computeLocalScaleHeight(rr - R_E);
      double cap = std::min(H / 20.0, 1000.0);
      if (maxRadialStep_m > 0.0) { // the RK4 radial cap
        double const rx = y[0] / rr, rz = y[1] / rr;
        double const cosT = std::fabs(y[2] * rx + y[3] * rz);
        if (cosT > 1.0e-9) cap = std::min(cap, maxRadialStep_m / cosT);
      }
      return cap;
    };

    double h = std::min(1.0, capFor(ycur)); // conservative start
    double sTotal = 0.0;                  // accepted arc length of this call
    double const sigmaPrev = ycur[4];
    double const geoPrev = geometric_path, opPrev = optical_path;
    int const maxSteps = 200000;

    // ---- Phase 4, step 3: impact-weighted acceptance -------------------------------
    // A position error dx at the crossing moves the hit point by
    //     dX = dx - u (n.dx)/(n.u),
    // so the component ALONG THE PLANE NORMAL is amplified by 1/(n.u).  The per-component
    // tolerance controls |dx| and cannot see that amplification, which is exactly why the
    // shallow rays carried the largest impact error.  Everything is in the (t_hat, r_hat)
    // frame the ODE lives in: n has components (A, B) there, and u = (y[2], y[3]).
    // Returns |n.u|; sets eImpact = |dX|.
    auto const impactError = [&](double const* yq, double const* dErr) {
      double const cn = A * yq[2] + B * yq[3];      // n . u
      double const ndx = A * dErr[0] + B * dErr[1]; // n . dx
      double const udx = yq[2] * dErr[0] + yq[3] * dErr[1];
      // s = the arc-length shift of the crossing caused by dx; n.u -> 0 means the ray runs
      // in the plane and there is no crossing to speak of, so the ratio is clamped.  Such a
      // ray is refused by the ddds test above, never certified here.
      double const s = (std::fabs(cn) > 1.0e-6) ? ndx / cn : 0.0;
      double const m2 = dErr[0] * dErr[0] + dErr[1] * dErr[1] - 2.0 * s * udx + s * s;
      return std::sqrt(std::max(0.0, m2));
    };

    // Extinction over a sub-interval of the last accepted step, on the dense output,
    // with the SAME altitude-fixed rule as the RK4 tracer (10 m of radius, floored by
    // 50 m of transverse length).  A large step therefore does not coarsen tau.
    auto const tauOverStep = [&](double const t0, double const t1, double const hStep) {
      if (!extinctionCoeffFunc_) return 0.0;
      double y0q[ND], y1q[ND];
      dp.dense(t0, y0q);
      dp.dense(t1, y1q);
      double const rA = std::sqrt(y0q[0] * y0q[0] + y0q[1] * y0q[1]);
      double const rB = std::sqrt(y1q[0] * y1q[0] + y1q[1] * y1q[1]);
      double const dsSeg = (t1 - t0) * hStep;
      int const nA = static_cast<int>(std::ceil(std::abs(rB - rA) / kTauSubStep_m));
      int const nL = static_cast<int>(std::ceil(dsSeg / kTauSubStepTransverse_m));
      int const nT = std::max(1, std::max(nA, nL));
      double const dt = (t1 - t0) / double(nT);
      double const dsTau = dsSeg / double(nT);
      double acc = 0.0;
      for (int i = 0; i < nT; ++i) {
        double yq[ND];
        dp.dense(t0 + (double(i) + 0.5) * dt, yq);
        double const rr = std::sqrt(yq[0] * yq[0] + yq[1] * yq[1]);
        acc += extinctionCoeffFunc_(rr - R_E, lambdaNm) * dsTau;
      }
      return acc;
    };


    std::size_t fPrev = dp.nF;
    for (int step = 0; step < maxSteps; ++step) {
      double const* yCur = dp.state();
      double yPrev[ND];
      for (int i = 0; i < ND; ++i) yPrev[i] = yCur[i]; // needed to undo / to re-land
      double const rCur = std::sqrt(yCur[0] * yCur[0] + yCur[1] * yCur[1]);
      double const hCur = rCur - R_E;
      if (hCur > 100000.0 || hCur < -100.0) { // escaped / below sea level
        ++escapedAtmosphere_;
        flush();
        return false;
      }

      double const cap = capFor(yCur);
      double const hTry = std::min(h, cap);
      if (hTry < h - 1.0e-12) ++cCapped;
      double hNext = hTry;
      Dopri5Result const res = dp.tryStep(rhs, hTry, hNext);
      cF += dp.nF - fPrev;
      fPrev = dp.nF;

      if (res == Dopri5Result::MinStepFailed) {
        ++cMinFail;
        flush();
        return false;
      }
      if (res == Dopri5Result::Rejected) {
        ++cRej;
        h = hNext; // FSAL stays valid: the state did not move
        continue;
      }

      // ---- Phase 4, step 3: is this step good enough FOR THE HIT? --------------
      // The ODE-level test (err <= 1) has passed.  An impact budget additionally requires
      // the PROJECTED error to fit this step's share of the budget, so that the accumulated
      // hit error is bounded by impactBudget_m_ instead of by whatever 1/cos the geometry
      // happens to apply.  Share = h/L_remaining (equidistribution), with L_remaining read
      // off the plane distance and the exact rate d(plane distance)/ds = n.u.
      if (impactBudget_m_ > 0.0) {
        double const hAcc2 = dp.lastStepSize();
        double const eHit = impactError(dp.state(), dp.lastErrorVector());
        double const cn2 = A * dp.state()[2] + B * dp.state()[3];
        double const dPl = std::fabs(A * dp.state()[0] + B * dp.state()[1] - Cp);
        double const rem =
            std::max(hAcc2, (std::fabs(cn2) > 1.0e-6) ? dPl / std::fabs(cn2) : 1.0e9);
        double const lim = impactBudget_m_ * (hAcc2 / rem);
        if (eHit > lim) {
          // the embedded estimate scales as h^5, so the admissible size is h (lim/e)^(1/5)
          double const ratio = lim / (eHit > 0.0 ? eHit : 1.0);
          double const hNew =
              hAcc2 * std::max(0.1, std::min(1.0, 0.9 * std::pow(ratio, 0.2)));
          if (hNew >= dopriMinStep_m_) {
            dp.setState(yPrev); // undo; setState invalidates FSAL, so k0 is recomputed
            h = hNew;
            ++cImpactRej;
            continue;
          }
          ++impactUncertified_; // accepted but NOT certified: counted, never silent
        }
        ++cImpactOk;
      }

      ++cAcc;

      // ---- accepted: look for the plane INSIDE this step ------------------------
      // The plane distance is sampled at 0, 1/4, 1/2, 3/4, 1 on the dense output.  A dip
      // that returns to the same side is therefore visible, which an endpoint-only sign
      // test would miss (the "grazing approaches" requirement).
      double yEnd[ND];
      dp.dense(1.0, yEnd);
      double const ts[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
      double ds5[5];
      for (int k = 0; k < 5; ++k) {
        double yq[ND];
        dp.dense(ts[k], yq);
        ds5[k] = A * yq[0] + B * yq[1] - Cp;
      }
      int kBr = -1;
      for (int k = 0; k < 4; ++k)
        if (ds5[k] * ds5[k + 1] <= 0.0) { kBr = k; break; }

      double const hAcc = dp.lastStepSize();
      if (kBr >= 0) {
        // ---- bracketed solve on the dense output; position, direction, group path and
        // ---- tau are all taken at the SAME tc --------------------------------
        double lo = ts[kBr], hi = ts[kBr + 1];
        double fLo = ds5[kBr];
        for (int it = 0; it < 60; ++it) {
          double const mid = 0.5 * (lo + hi);
          double yq[ND];
          dp.dense(mid, yq);
          double const fMid = A * yq[0] + B * yq[1] - Cp;
          if (fLo * fMid <= 0.0) {
            hi = mid;
          } else {
            lo = mid;
            fLo = fMid;
          }
          ++cCross;
          if (!(hi - lo > 1.0e-13)) break;
        }
        double const tc = 0.5 * (lo + hi);
        double yq[ND];
        dp.dense(tc, yq);
        double hUse = tc * hAcc; // arc length consumed by the part of this step that counts

        // ---- Phase 4, step 3: optional crossing correction -- land ON the plane -----
        // The interpolant is the right tool to LOCATE the crossing, but its INTERIOR value
        // is off by O(h^5): measured +2.08 mm at t = 0.25 of a 15-unit step while both
        // endpoints are inside 0.02 mm (phase4_plan.md, step 3 task 1).  Re-take the
        // crossing step so it ENDS on the plane, and read the 5th-order STATE.  Newton on
        // the step length, with d(plane distance)/d(arc length) = n.u exactly.
        if (crossingRelanding_) {
          double const cn = A * dp.state()[2] + B * dp.state()[3];
          double const hStart = tc * hAcc;
          double hSub = hUse;
          bool landed = false;
          for (int it = 0; it < 8; ++it) {
            Dopri5<ND> d2;
            d2.setTolerances(atol, rtol);
            d2.setMinStep(dopriMinStep_m_);
            double ys[ND];
            for (int k = 0; k < ND; ++k) ys[k] = yPrev[k];
            d2.setState(ys);
            double hn2 = hSub;
            if (d2.tryStep(rhs, hSub, hn2) != Dopri5Result::Accepted) {
              hSub *= 0.5;
              continue;
            }
            double const f = A * d2.state()[0] + B * d2.state()[1] - Cp;
            if (std::fabs(f) <= 1.0e-9) {
              for (int k = 0; k < ND; ++k) yq[k] = d2.state()[k];
              hUse = hSub;
              landed = true;
              break;
            }
            if (!(std::fabs(cn) > 1.0e-9)) break;
            double const hNextLand = hSub - f / cn;
            if (!(hNextLand > 0.0) || !(hNextLand < 4.0 * hStart + 1.0)) break;
            hSub = hNextLand;
          }
          if (!landed) ++impactUncertified_; // fell back to the interpolant, and counted
        }

        accumulated_tau += tauOverStep(0.0, tc, hAcc);

        x_impact = yq[1] * r_hat_x + yq[0] * t_hat_x;
        y_impact = yq[1] * r_hat_y + yq[0] * t_hat_y;
        z_impact = yq[1] * r_hat_z + yq[0] * t_hat_z;
        {
          double const uxi = yq[2], uzi = yq[3];
          double const nx = uzi * r_hat_x + uxi * t_hat_x;
          double const ny = uzi * r_hat_y + uxi * t_hat_y;
          double const nz = uzi * r_hat_z + uxi * t_hat_z;
          double const nn = std::sqrt(nx * nx + ny * ny + nz * nz);
          if (nn > 0.0) {
            dx_impact = nx / nn;
            dy_impact = ny / nn;
            dz_impact = nz / nn;
          } else {
            dx_impact = dx;
            dy_impact = dy;
            dz_impact = dz;
          }
        }
        geometric_path = geoPrev + sTotal + hUse;
        optical_path = opPrev + (sTotal + hUse) + fGroup * (yq[4] - sigmaPrev);
        flush();
        return true;
      }

      // ---- no crossing: accumulate the whole step ------------------------------
      accumulated_tau += tauOverStep(0.0, 1.0, hAcc);
      sTotal += hAcc;
      geometric_path = geoPrev + sTotal;
      optical_path = opPrev + sTotal + fGroup * (yEnd[4] - sigmaPrev);
      cHMin = std::min(cHMin, hAcc);
      cHMax = std::max(cHMax, hAcc);
      h = hNext; // FSAL hand-off is armed by the stepper itself
    }

    std::size_t const nMax = maxStepsReached_.fetch_add(1) + 1;
    if (nMax <= 3) {
      CORSIKA_LOGGER_WARN(logger_,
                          "rayTracePlaneDopri: step limit reached without hitting the "
                          "plane [logged for the first 3 only]");
    }
    flush();
    return false;
  }
  // ===========================================================================
  // Phase 4, step 1: regularised (turning-point) plane tracer
  // ---------------------------------------------------------------------------
  // Same contract as rayTraceSphericalToPlane / rayTracePlaneDopri: same args, same
  // outputs, same tau ownership.  Only the independent variable changes: the state is
  //
  //     y = (Phi, s, Sigma, v),    r = r_t + v^2,
  //
  // with r_t the turning radius (q(r_t) = p, q = n r) and p = n r sin(psi) the conserved
  // invariant.  v is carried as a state component so the RHS stays autonomous and the
  // dense output delivers (Phi, s, Sigma, v) at the same t:
  //
  //     dPhi/dv = 2 p / [ (r_t+v^2) sqrt(B) ],  ds/dv = 2 q / sqrt(B),
  //     dSigma/dv = N ds/dv,                    dv/dv = 1,
  //     B(v) = ( q(r_t+v^2)^2 - p^2 ) / v^2,    B(0) = 2 p q'(r_t)  (finite).
  //
  // Refusals (returns false, counted) rather than silent approximations:
  //   * no bracketed turning radius, or q'(r_t) <= 0 (degenerate/ducting) -> TangentGeometry;
  //   * the ray reaches the ground before the plane;
  //   * the ray reaches the top of the atmosphere;
  //   * the step guard / minimum step is hit.
  // ===========================================================================
  inline bool CherenkovProcess::rayTracePlaneTangent(
      double x0, double y0, double z0_altitude, double dx, double dy, double dz,
      double plane_x, double plane_y, double plane_z, double norm_x, double norm_y,
      double norm_z, double& x_impact, double& y_impact, double& z_impact,
      double& optical_path, double& geometric_path, double lambdaNm,
      double& accumulated_tau, double& dx_impact, double& dy_impact, double& dz_impact,
      double maxRadialStep_m) const {
    constexpr int ND = 4; // Phi, s, Sigma, v
    using Dopri = Dopri5<ND>;

    double const R_E = earthRadiusM_;
    double const px = x0, py = y0, pz = z0_altitude;
    double const r0 = std::sqrt(px * px + py * py + pz * pz);
    if (!(r0 > 0.0)) return false;

    double const r_hat_x = px / r0, r_hat_y = py / r0, r_hat_z = pz / r0;
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

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
      double ref_x = 0.0, ref_y = 0.0, ref_z = 1.0;
      if (std::abs(r_hat_z) > 0.9) { ref_x = 1.0; ref_y = 0.0; ref_z = 0.0; }
      t_hat_x = r_hat_y * ref_z - r_hat_z * ref_y;
      t_hat_y = r_hat_z * ref_x - r_hat_x * ref_z;
      t_hat_z = r_hat_x * ref_y - r_hat_y * ref_x;
      double const tn = std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
      t_hat_x /= tn; t_hat_y /= tn; t_hat_z /= tn;
    }

    double const A = t_hat_x * norm_x + t_hat_y * norm_y + t_hat_z * norm_z;
    double const Bn = r_hat_x * norm_x + r_hat_y * norm_y + r_hat_z * norm_z;
    double const Cp = plane_x * norm_x + plane_y * norm_y + plane_z * norm_z;
    if (std::abs(A) < 1e-15 && std::abs(Bn) < 1e-15) return false;

    constexpr double kTauSubStep_m = 10.0;
    constexpr double kTauSubStepTransverse_m = 50.0;

    // ---- dispersion factors and the tangent geometry ---------------------------
    double fPhase = 1.0, fGroup = 1.0;
    {
      auto const ng0 = refractiveIndexFunc_(r0 - R_E, lambdaNm);
      double const Nr0 = refractiveIndexTable_->N(r0);
      if (Nr0 > 1.0e-12) {
        fPhase = (ng0.first - 1.0) / Nr0;
        fGroup = (ng0.second - 1.0) / Nr0;
      } else {
        fPhase = 0.0;
        fGroup = 0.0;
      }
    }
    TangentGeometry const geo(*refractiveIndexTable_, fPhase);

    // ---- the backward-crossing branch (launch-and-continue), as in the others -----
    double const d_start = Bn * r0 - Cp;
    double const ddds = A * sin_theta0 + Bn * cos_theta0;
    if (std::abs(ddds) > 1.0e-12) {
      double const tPlane = -d_start / ddds;
      if (tPlane <= 0.0) {
        ++nPlaneCalls_;
        ++nPlaneBack_;
        if (tPlane > -1.0) ++nPlaneTangent_;
        x_impact = px + tPlane * dx;
        y_impact = py + tPlane * dy;
        z_impact = pz + tPlane * dz;
        dx_impact = dx;
        dy_impact = dy;
        dz_impact = dz;
        geometric_path += tPlane;
        optical_path += tPlane + fGroup * geo.N(r0) * tPlane;
        if (extinctionCoeffFunc_) {
          accumulated_tau = std::max(
              0.0, accumulated_tau - extinctionCoeffFunc_(r0 - R_E, lambdaNm) * std::fabs(tPlane));
        }
        ++trTangentRays_; // served by the launch-and-continue branch
        return true;
      }
    }

    // ---- invariant, turning radius and the starting v ---------------------------
    double const p = TangentGeometry::invariant(geo.q(r0), cos_theta0);
    // The floor of the tangent mode is the EARTH'S SURFACE, not the observation level:
    // groundZ_m_ is the sphere target of Ground mode (3 km in the telescope runs), and
    // using it here refused every ray that dips below the dish -- which is exactly the
    // population the mode exists for.  The floor matches the regular tracers' -100 m.
    double const rGround = R_E - 100.0;
    double rt = 0.0;
    if (!(p > 0.0)) return false; // radial ray: nothing to regularise
    // CHEAP NECESSARY CONDITION, O(1): a turning point above the floor exists only if
    // p = q(r0) sin(psi) > q(floor).  Without this, AUTO pays a 90-iteration bracket
    // solve for every ray of a shower that is not near-horizontal -- which is the common
    // case (a 43 km muon descends at ~3 deg, so every photon fails the test).
    if (!(p > geo.q(rGround))) {
      ++trTangentUncertified_;
      return false;
    }
    if (!geo.turningRadius(p, std::max(1.0, rGround - 900.0), r0, rt)) {
      ++trTangentUncertified_; // degenerate or un-bracketed: report, do not guess
      return false;
    }
    // A turning radius below the ground does NOT invalidate the ray: it may cross the
    // telescope plane first.  The loop stops at r <= r_ground instead ("do not force the
    // ray through the turn" is enforced there, not here).

    // Starting v: sign from the radial direction cosine (dr/ds), magnitude from the
    // radial distance to the turning point.  v > 0 means "already past the turn".
    double const v0 = ((cos_theta0 >= 0.0) ? 1.0 : -1.0) *
                      std::sqrt(std::max(0.0, r0 - rt));
    double const rTop = R_E + 100000.0;
    double const vMax = std::sqrt(std::max(0.0, rTop - rt));

    double ycur[ND] = {0.0, 0.0, 0.0, v0};
    auto const rhs = [&](double const* y, double* f) {
      double const v = y[3];
      double const r = rt + v * v;
      double const qv = geo.q(r);
      double const Bv = geo.B(rt, p, v);
      double const sB = std::sqrt(std::max(Bv, 1.0e-300));
      f[0] = 2.0 * p / (r * sB); // dPhi/dv
      f[1] = 2.0 * qv / sB;      // ds/dv
      f[2] = geo.N(r) * f[1];    // dSigma/dv
      f[3] = 1.0;                // dv/dv
    };

    // Tolerances.  Phi is in RADIANS and the impact lever arm is r, so its absolute
    // tolerance is atol_pos/r0 -- 1e-9 rad at the Earth's surface would be 6 mm of
    // impact and must not be used.  v is exact (dv/dv = 1), so it is not controlled.
    // Factor 1/2: with atol_Phi = atol_pos/r0 exactly, the measured arc error at the
    // crossing was 1.19 mm for a 1 mm request (the tolerance is a per-step budget, and
    // the accumulated error sits at it), so the request is halved to make the stated
    // millimetre conservative.
    double atol[ND] = {0.5 * dopriAtolPos_m_ / std::max(1.0, r0),
                       0.5 * dopriAtolPos_m_, dopriAtolOpt_m_, 1.0e30};
    double rtol[ND] = {0.0, 0.0, 0.0, 0.0};
    Dopri dp;
    dp.setTolerances(atol, rtol);
    dp.setState(ycur);
    dp.setMinStep(dopriMinStep_m_);

    std::size_t cAcc = 0, cRej = 0, cF = 0, cMinFail = 0, cCapped = 0, cCross = 0;
    double cHMin = 1.0e300, cHMax = 0.0;
    auto const flush = [&]() {
      trAccepted_ += cAcc;
      trRejected_ += cRej;
      trF_ += cF;
      trMinStepFail_ += cMinFail;
      trCapped_ += cCapped;
      trCrossIters_ += cCross;
      double lo = trHMin_.load();
      while (cHMin < lo && !trHMin_.compare_exchange_weak(lo, cHMin)) {}
      double hi = trHMax_.load();
      while (cHMax > hi && !trHMax_.compare_exchange_weak(hi, cHMax)) {}
    };

    // Step cap expressed as a RADIAL advance: dr = 2 v h + h^2 = capR  =>  h = -v + sqrt(v^2+capR).
    auto const capFor = [&](double const v) {
      double capR;
      if (dopriMaxStep_m_ > 0.0) {
        capR = dopriMaxStep_m_;
      } else {
        capR = std::min(computeLocalScaleHeight(rt + v * v - R_E) / 20.0, 1000.0);
        if (maxRadialStep_m > 0.0) capR = std::min(capR, maxRadialStep_m);
      }
      return -v + std::sqrt(v * v + capR);
    };

    // Extinction over a sub-interval of the last accepted step, from the dense output,
    // with the same altitude-fixed 10 m / 50 m rule as the other tracers.
    auto const tauOverStep = [&](double const t0, double const t1) {
      if (!extinctionCoeffFunc_) return 0.0;
      double y0q[ND], y1q[ND];
      dp.dense(t0, y0q);
      dp.dense(t1, y1q);
      double const rA = rt + y0q[3] * y0q[3];
      double const rB = rt + y1q[3] * y1q[3];
      double const dsSeg = y1q[1] - y0q[1];
      int const nA = static_cast<int>(std::ceil(std::abs(rB - rA) / kTauSubStep_m));
      int const nL = static_cast<int>(std::ceil(std::abs(dsSeg) / kTauSubStepTransverse_m));
      int const nT = std::max(1, std::max(nA, nL));
      double const dt = (t1 - t0) / double(nT);
      double const dsTau = dsSeg / double(nT);
      double acc = 0.0;
      for (int i = 0; i < nT; ++i) {
        double yq[ND];
        dp.dense(t0 + (double(i) + 0.5) * dt, yq);
        double const rr = rt + yq[3] * yq[3];
        acc += extinctionCoeffFunc_(rr - R_E, lambdaNm) * dsTau;
      }
      return acc;
    };

    // (Phi, v) at parameter t of the current step -> plane coordinates and plane distance.
    auto const planeDist = [&](double const t, double& dOut, double& xpOut,
                               double& zpOut) {
      double yq[ND];
      dp.dense(t, yq);
      double const v = yq[3];
      double const r = rt + v * v;
      double const xp = r * std::sin(yq[0]);
      double const zp = r * std::cos(yq[0]);
      dOut = A * xp + Bn * zp - Cp;
      xpOut = xp;
      zpOut = zp;
    };
    double const geoPrev = geometric_path, opPrev = optical_path;
    double h = std::min(1.0, capFor(v0));
    std::size_t fPrev = dp.nF;
    int const maxSteps = 400000;

    for (int step = 0; step < maxSteps; ++step) {
      double const vCur = dp.state()[3];
      if (rt + vCur * vCur <= rGround) { // hit the ground before the plane
        ++trTangentUncertified_;
        flush();
        return false;
      }
      if (vCur > vMax) { // above the atmosphere without a crossing
        ++escapedAtmosphere_;
        flush();
        return false;
      }
      double const hTry = std::min(std::min(h, capFor(vCur)), std::max(1.0e-6, vMax - vCur));
      if (hTry < h - 1.0e-12) ++cCapped;
      double hNext = hTry;
      Dopri5Result const res = dp.tryStep(rhs, hTry, hNext);
      cF += dp.nF - fPrev;
      fPrev = dp.nF;
      if (res == Dopri5Result::MinStepFailed) {
        ++cMinFail;
        ++trTangentUncertified_;
        flush();
        return false;
      }
      if (res == Dopri5Result::Rejected) {
        ++cRej;
        h = hNext;
        continue;
      }
      ++cAcc;
      cHMin = std::min(cHMin, dp.lastStepSize());
      cHMax = std::max(cHMax, dp.lastStepSize());

      // ---- plane crossing inside this step (5-point scan, as in the DOPRI tracer) ---
      double const ts[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
      double ds5[5];
      for (int k = 0; k < 5; ++k) {
        double xpd = 0.0, zpd = 0.0;
        planeDist(ts[k], ds5[k], xpd, zpd);
      }
      int kBr = -1;
      for (int k = 0; k < 4; ++k)
        if (ds5[k] * ds5[k + 1] <= 0.0) { kBr = k; break; }

      if (kBr >= 0) {
        double lo = ts[kBr], hi = ts[kBr + 1];
        double fLo = ds5[kBr];
        for (int it = 0; it < 60; ++it) {
          double const mid = 0.5 * (lo + hi);
          double dq = 0.0, xpd = 0.0, zpd = 0.0;
          planeDist(mid, dq, xpd, zpd);
          if (fLo * dq <= 0.0) {
            hi = mid;
          } else {
            lo = mid;
            fLo = dq;
          }
          ++cCross;
          if (!(hi - lo > 1.0e-13)) break;
        }
        double const tc = 0.5 * (lo + hi);
        double yq[ND];
        dp.dense(tc, yq);
        double const v = yq[3];
        double const r = rt + v * v;
        double const Phi = yq[0];
        double const xp = r * std::sin(Phi);
        double const zp = r * std::cos(Phi);

        accumulated_tau += tauOverStep(0.0, tc);

        x_impact = zp * r_hat_x + xp * t_hat_x;
        y_impact = zp * r_hat_y + xp * t_hat_y;
        z_impact = zp * r_hat_z + xp * t_hat_z;
        {
          // direction: sin(psi) = p/q(r), cos(psi) = sign(v) sqrt(1 - sin^2)
          double const sPsi = geo.sinPsi(r, p);
          double const cPsi = ((v >= 0.0) ? 1.0 : -1.0) *
                              std::sqrt(std::max(0.0, 1.0 - sPsi * sPsi));
          double const uxp = cPsi * std::sin(Phi) + sPsi * std::cos(Phi);
          double const uzp = cPsi * std::cos(Phi) - sPsi * std::sin(Phi);
          double const nx = uzp * r_hat_x + uxp * t_hat_x;
          double const ny = uzp * r_hat_y + uxp * t_hat_y;
          double const nz = uzp * r_hat_z + uxp * t_hat_z;
          double const nn = std::sqrt(nx * nx + ny * ny + nz * nz);
          if (nn > 0.0) {
            dx_impact = nx / nn;
            dy_impact = ny / nn;
            dz_impact = nz / nn;
          } else {
            dx_impact = dx;
            dy_impact = dy;
            dz_impact = dz;
          }
        }
        geometric_path = geoPrev + yq[1];
        optical_path = opPrev + yq[1] + fGroup * yq[2];
        ++trTangentRays_; // SERVED: a failed attempt is counted only as uncertified
        flush();
        return true;
      }

      // ---- no crossing: accumulate and continue ------------------------------------
      accumulated_tau += tauOverStep(0.0, 1.0);
      geometric_path = geoPrev + dp.state()[1];
      optical_path = opPrev + dp.state()[1] + fGroup * dp.state()[2];
      h = hNext;
    }

    std::size_t const nMax = maxStepsReached_.fetch_add(1) + 1;
    if (nMax <= 3) {
      CORSIKA_LOGGER_WARN(logger_,
                          "rayTracePlaneTangent: step limit reached without hitting the "
                          "plane [logged for the first 3 only]");
    }
    ++trTangentUncertified_;
    flush();
    return false;
  }




  // Ray trace photon in a spherically stratified atmosphere with full 3D support
  // Handles both upward and downward-going rays using Bouguer's invariant
  // Detects turning points for upward rays and integrates until hitting ground
  // --- absorption: deterministic transmission, never a draw --------------------
  // Returns exp(-tau) with tau from the state-appropriate source.  Absorption is applied
  // by multiplying the bunch weight (see processPhotons), so:
  //   * the recorded yield is unbiased -- exp(-tau) is exactly the expected fraction of
  //     a bunch's photons that survive, and a bunch is a set of photons on ONE path
  //     (same lambda, same trajectory), so one exp(-tau) is the right factor for all of
  //     them;
  //   * the binomial variance of the old survival sampling is removed (variance
  //     reduction, not a bias);
  //   * nothing here consumes randomness, so the photon stream stays paired between
  //     configurations -- which is what makes an A/B yield comparison meaningful.
  inline double CherenkovProcess::absorptionTransmission(AbsorptionState state,
                                                         double accumulatedTau, double lambdaNm,
                                                         double hEmitAboveSea,
                                                         double hDetAboveSea, double pathLength,
                                                         double cosZenith) const {
    if (!enableAbsorption_) return 1.0;
    if (state == AbsorptionState::OpticalDepthAvailable) {
      return std::exp(-std::max(0.0, accumulatedTau));
    }
    // NotEvaluated: no tau was integrated along the ray -> exact chord integral.
    // arcLength <= 0 means "only the heights are known"; the helper then falls back to
    // the vertical column difference.
    double const tau = absorptionTable_.getOpticalDepthChord(
        lambdaNm, hEmitAboveSea, hDetAboveSea, pathLength, earthRadiusM_);
    (void)cosZenith;
    return std::exp(-std::max(0.0, tau));
  }

  // --- step A': the quadrature referee ----------------------------------------
  // tau = int alpha(h) ds along the REFRACTED ray, with ds = dh / cos(psi) and
  // cos(psi) = sqrt(1 - (p / (n(h)(R+h)))^2), p = xi * n(r_obs) r_obs.  The integral is
  // taken in u = sqrt(h - h_t), where h_t is the true perigee (A(h_t) = p), so the sqrt
  // singularity is removed by the substitution rather than sampled.  When the photon
  // descends at emission and its perigee lies below the observation level, the path
  // visits the perigee and BOTH legs are accumulated.
  inline double CherenkovProcess::absorptionTauQuadrature(double lambdaNm, double hEmit,
                                                          double hObs, double xi,
                                                          bool descendsAtEmit,
                                                          int nSamples) const {
    if (!extinctionCoeffFunc_ || !refractiveIndexTable_) return 0.0;
    auto const& ri = *refractiveIndexTable_;
    auto Aof = [&](double h) {
      double const z = earthRadiusM_ + h;
      return z + ri.N(z) * z;  // n r = (1 + N) r
    };
    double const p = xi * Aof(hObs);
    double h_t = 0.0;
    if (Aof(0.0) < p) {
      double lo = 0.0, hi = 50000.0;
      for (int it = 0; it < 80; ++it) {
        double const mid = 0.5 * (lo + hi);
        (Aof(mid) < p) ? (lo = mid) : (hi = mid);
      }
      h_t = 0.5 * (lo + hi);
    }
    auto cosPsi = [&](double h) {
      double const s = p / Aof(h);
      return (s < 1.0) ? std::sqrt(std::max(0.0, 1.0 - s * s)) : 0.0;
    };
    auto leg = [&](double a, double b) {
      if (!(b > a)) return 0.0;
      double const ua = std::sqrt(std::max(0.0, a - h_t));
      double const ub = std::sqrt(std::max(1.0e-9, b - h_t));
      int const n = std::max(16, nSamples);  // WAVE 5: was a hard 20000 (see convergence test)
      double acc = 0.0;
      for (int i = 0; i < n; ++i) {
        double const u = ua + (ub - ua) * (double(i) + 0.5) / double(n);
        double const h = h_t + u * u;
        double const c = cosPsi(h);
        if (!(c > 0.0)) continue;
        acc += extinctionCoeffFunc_(h, lambdaNm) * 2.0 * u * ((ub - ua) / double(n)) / c;
      }
      return acc;
    };
    double const hLo = std::min(hEmit, hObs), hHi = std::max(hEmit, hObs);
    if (descendsAtEmit && h_t < hLo) return leg(h_t, hEmit) + leg(h_t, hObs);
    return leg(hLo, hHi);
  }

  // --- step A: per-photon A/B diagnostics (see Cherenkov.hpp) ------------------
  inline void CherenkovProcess::enableTraceAB(std::size_t count, std::string path) const {
    std::lock_guard<std::mutex> lock(abMutex_);
    abMaxRows_ = count;
    abRows_ = 0;
    abPath_ = std::move(path);
    if (count == 0 || abPath_.empty()) return;
    std::ofstream out(abPath_);  // truncate and write the header
    out << "id,h_emit_m,mu,lambda_nm,status_table,L_table_m,L_legacy_m,dL_rel,d_impact_m,"
           "x_table,y_table,x_legacy,y_legacy,tau_table,tau_legacy,"
           "tau_chord_table,tau_chord_legacy,legacy_ok,tau_ref_quad,tau_ref_quad_lo\n";
    abLegacyRows_ = 0;
    std::ofstream leg(abPath_ + ".legacy.csv");
    leg << "id,h_emit_m,mu_derived,lambda_nm,tau_tracer,tau_table,tau_chord,h_impact_m,perigee\n";
  }

  inline void CherenkovProcess::writeTraceAB_(TraceABRow const& r) const {
    std::lock_guard<std::mutex> lock(abMutex_);
    std::ofstream out(abPath_, std::ios::app);
    if (!out) return;
    out.precision(10);
    out << abRows_ << ',' << r.hEmit << ',' << r.mu << ',' << r.lambdaNm << ','
        << r.statusTable << ',' << r.Ltable << ',' << r.Llegacy << ',' << r.dLrel << ','
        << r.dImpact << ',' << r.xTable << ',' << r.yTable << ',' << r.xLegacy << ','
        << r.yLegacy << ',' << r.tauTable << ',' << r.tauLegacy << ',' << r.tauChordTable
        << ',' << r.tauChordLegacy << ',' << r.legacyOk << ',' << r.tauRefQuad << ','
        << r.tauRefQuadLo << '\n';
  }

  // --- curved ground tracing: table by default, legacy on rejection/ineligible ---
  inline bool CherenkovProcess::traceCurvedGround(
      double x0, double y0, double z0, double dx, double dy, double dz, double lambdaNm,
      double& x_ground, double& y_ground, double& z_ground, double& optical_path,
      double& geometric_path, double& accumulated_tau, bool& usedTable) const {
    usedTable = false;
    if (!curvedState_) curvedState_ = std::make_shared<CurvedTableState>();
    CurvedTableState& st = *curvedState_;

    // ---- phase 1: configuration/physics prerequisites (cheap, no table) ------
    // 2g-h2m: tau is now TABULATED (CurvedRayTable stores int alpha ds at a few
    // wavelength nodes), so absorption no longer excludes the fast path.  The table
    // is built with the production extinction function and returns the ray's own
    // optical depth, which the caller turns into a weight factor (no sampling).
    bool const absorptionActive = enableAbsorption_ && extinctionCoeffFunc_;
    // 2g-h2m: Telescope mode is NOT served by this table.  The table returns the
    // crossing of the r = r_obs SPHERE, the telescope needs the crossing of its
    // PLANE, and at grazing incidence the two are unrelated: a photon at
    // cos(theta) = 0.0095 crosses the plane 1.4 m below the sphere, which is 147 m of
    // lateral travel (measured, 3100 m / 10 km: plane crossing 1.39 m from the dish = a
    // hit, sphere crossing 146.8 m away = a miss).  The lever is 1/cos(theta) ~ 100 at
    // horizon pointing, so a plane-parametrised table (traced from the dish) is needed
    // instead of a reuse.  Ground mode is unaffected: there the sphere IS the target.
    bool const prerequisites = enableCurvature_ && enableCurvedTable_ &&
                               (projectionMode_ == ProjectionMode::Ground ||
                                projectionMode_ == ProjectionMode::Telescope) &&
                               refractiveIndexTable_ && refractiveIndexTable_->valid();

    // WAVE 5: the plane tracer reports the arrival direction too, so the legacy route can
    // build its tau with the SAME (xi, perigee) test the table route uses.
    double ground_dir_x = 0.0, ground_dir_y = 0.0, ground_dir_z = 0.0;
    auto legacy = [&]() {
      if (projectionMode_ == ProjectionMode::Telescope) {
        // The Telescope route must fall back to the PLANE tracer, not the sphere one.
        return rayTraceSphericalToPlane(x0, y0, z0, dx, dy, dz, telescopeX_m_, telescopeY_m_,
                                        telescopeZ_m_, telescopeBasisZ_x_, telescopeBasisZ_y_,
                                        telescopeBasisZ_z_, x_ground, y_ground, z_ground,
                                        optical_path, geometric_path, lambdaNm, accumulated_tau,
                                        ground_dir_x, ground_dir_y, ground_dir_z);
      }
      return rayTraceSpherical(x0, y0, z0, dx, dy, dz, x_ground, y_ground, z_ground, optical_path,
                               geometric_path, lambdaNm, accumulated_tau);
    };

    if (!prerequisites) {
      ++st.bypassed;
      bool const ok = legacy();
      // ---- WAVE 5: give the legacy route a REFRACTED-PATH tau ---------------------
      // Measured on the legacy route: its own integral equals the exact CHORD to its own
      // impact to 0.998-1.000 in every mu bin, which leaves the mean tau 17 % low
      // (0.2103 vs the table's 0.2461) and the absorbed fraction 9.81 % instead of 12.14 %.
      // The absorption table IS the refracted-path integral (0.03 % vs an exact quadrature),
      // so when absorption is on, tau comes from it.
      //
      // The endpoints follow the TABLE route's convention exactly -- emission height to the
      // OBSERVER height (not the plane impact, which spans 2.3-5.0 km for a 3 km telescope),
      // xi from the conserved invariant, perigee flag from the emission's radial direction.
      // The first attempt used the plane-impact height and a different perigee test and moved
      // the wrong way; that is what this fixes.
      if (abMaxRows_ > 0 && ok && absorptionActive && extinctionCoeffFunc_ && refractiveIndexTable_ &&
          refractiveIndexTable_->valid()) {
        if (!legacyAbsTable_) {
          std::call_once(legacyAbsOnce_, [&]() {
            RefractiveIndexTable const& ril = *refractiveIndexTable_;
            std::vector<double> absLam;
            for (double l = 300.0; l <= 620.0 + 1.0e-9; l += 5.0) absLam.push_back(l);
            auto dispF = [this](double lam) {
              double const Nr = refractiveIndexTable_->N(groundZ_m_);
              if (!refractiveIndexFunc_ || !(Nr > 1.0e-12)) return 1.0;
              double const hRef = groundZ_m_ - earthRadiusM_;
              return (refractiveIndexFunc_(hRef, lam).first - 1.0) / Nr;
            };
            auto t = std::make_shared<CurvedAbsorptionTable>(ril, earthRadiusM_, groundZ_m_,
                                                            extinctionCoeffFunc_, absLam, dispF);
            if (t->valid()) legacyAbsTable_ = std::move(t);
          });
        }
        if (legacyAbsTable_) {
          double const tauTracer = accumulated_tau;
          double const rEmit = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
          double const rImp =
              std::sqrt(x_ground * x_ground + y_ground * y_ground + z_ground * z_ground);
          double const vRadEmit = (rEmit > 0.0) ? (x0 * dx + y0 * dy + z0 * dz) / rEmit : 0.0;
          // WAVE 5: xi and the perigee flag come from the ARRIVAL state, exactly as the table
          // route does.  sin(theta_arr) is the table's second axis, and a path that descends
          // at emission yet arrives ASCENDING has been through its perigee, so both legs count.
          // The earlier attempt used the emitted direction instead: xi came out wrong and the
          // perigee branch over-triggered (28 % vs the 12 % target).
          double const vRadArr =
              (rImp > 0.0) ? (x_ground * ground_dir_x + y_ground * ground_dir_y +
                              z_ground * ground_dir_z) / rImp
                           : 0.0;
          double const xi = std::sqrt(std::max(0.0, 1.0 - vRadArr * vRadArr));
          double const hA = rEmit - earthRadiusM_;
          // ENDPOINT: the tracer's own impact height, not the observer level.  The observer
          // level is the TABLE route's convention only because that route then adds the
          // sphere->plane continuation leg; the legacy tracer traces straight to the plane, so
          // its tau must be evaluated between the same two heights it actually crossed.
          double const hB = rImp - earthRadiusM_;
          bool const visitsPerigee = (vRadEmit < 0.0 && vRadArr > 0.0);
          auto const ta = legacyAbsTable_->tauBetween(hA, hB, xi, lambdaNm, visitsPerigee);
          // WAVE 5 — STILL NOT THE DEFAULT.  With the arrival direction threaded through and
          // xi / perigee taken from it, the substitution gives absorbed = 26.76 % against the
          // table's 12.14 % (it was 28.16 % with the emitted-direction xi).  So the perigee
          // flag was NOT the missing piece; the residual is in the xi convention (the table
          // route's |sin(theta_obs)| comes from the table's own reconstruction, which is not
          // the same quantity as sin(theta) computed from the plane tracer's final state).
          // Gated behind --cherenkov-trace-ab until that is understood; the validated default
          // (legacy keeps its own tau) is unchanged.
          if (abMaxRows_ > 0 && ta.ok()) {
            accumulated_tau = ta.tau;
            // WAVE 5 probe: the substituted tau against the tracer's, per photon.
            ++absorptionWeights_.substCount;
            absorptionWeights_.substTracerSum += tauTracer;
            absorptionWeights_.substTableSum += ta.tau;
            double const dsub = std::abs(ta.tau - tauTracer);
            if (dsub > absorptionWeights_.substDiffMaxAbs) {
              absorptionWeights_.substDiffMaxAbs = dsub;
            }
          } else if (!ta.ok()) {
            ++st.absFallbacks;
            ++absorptionWeights_.substRejected;
          }
          if (abMaxRows_ > 0) {
            double const tauChord =
                absorptionTable_.getOpticalDepthChord(lambdaNm, hA, rImp - earthRadiusM_,
                                                      std::sqrt((x_ground - x0) * (x_ground - x0) +
                                                                (y_ground - y0) * (y_ground - y0) +
                                                                (z_ground - z0) * (z_ground - z0)),
                                                      earthRadiusM_);
            std::lock_guard<std::mutex> lock(abMutex_);
            if (abLegacyRows_ < abMaxRows_) {
              std::ofstream out(abPath_ + ".legacy.csv", std::ios::app);
              if (out) {
                out.precision(10);
                out << abLegacyRows_ << ',' << hA << ','
                    << std::sqrt(std::max(0.0, 1.0 - xi * xi)) << ',' << lambdaNm << ','
                    << tauTracer << ',' << (ta.ok() ? ta.tau : -1.0) << ',' << tauChord << ','
                    << hB << ',' << int(visitsPerigee) << '\n';
                ++abLegacyRows_;
              }
            }
          }
        }
      }
      if (ok) { ++st.bypassedLegacySucceeded; } else { ++st.bypassedLegacyFailed; ++st.finalDrop; }
      return ok;
    }

    // ---- phase 2: lazy, synchronized initialization -------------------------
    std::call_once(st.once, [&]() {
      st.attempted = true;
      RefractiveIndexTable const& ri = *refractiveIndexTable_;
      double const rObs = groundZ_m_;
      // The table samples N(r_obs + h) for h in [0, hTop], so coverage is checked on
      // that radial interval - BOTH bounds, with r_obs (not the Earth radius) as the
      // origin, because the RI table is indexed by ABSOLUTE z: below zMin() it silently
      // clamps and above zMax() it silently extrapolates.
      //
      // hTop is therefore CLAMPED to the coverage instead of being a fixed 120 km.  A
      // fixed 120 km made the build refuse for EVERY observation level above sea level:
      // a 3 km ASL telescope asks for [6374000, 6494000] while the shipped 0-120 km
      // profile tops out at 6490990, so the table never existed, traceCurvedGround
      // silently bypassed to the legacy tracer for every photon (measured on the 2g-h2n
      // Step D run: tableAvailable=false, 3804/3804 photons blind-legacy, 2071 of them
      // dropped) and the acceptance set could not see it because it builds its own
      // 0-130 km table.  Emissions above the clamped top are rejected as AboveTable and
      // take the documented fallback, which is correct: there is no Cherenkov up there.
      double const hTop = std::min(120000.0, ri.zMax() - rObs);
      if (!(rObs >= ri.zMin() - 1.0) || !(hTop > 1000.0)) {
        st.available = false;
        st.unavailableReason = "refractive-index table does not cover [" + std::to_string(rObs) +
                               ", " + std::to_string(rObs + 120000.0) + "] (has [" +
                               std::to_string(ri.zMin()) + ", " + std::to_string(ri.zMax()) +
                               "], usable top " + std::to_string(rObs + hTop) + ")";
        return;
      }
      if (hTop < 119999.0) {
        CORSIKA_LOGGER_INFO(logger_,
                            "CurvedRayTable top clamped to the refractive-index coverage: "
                            "hTop={:.0f} m (requested 120000, RI table tops out at {:.0f} m)",
                            hTop, ri.zMax());
      }
      // tau no longer lives in this table (wave 4b): CurvedAbsorptionTable owns it, with
      // 65 wavelength nodes instead of 5 and no 64 MiB block here.
      // ---- wave 5: reuse an on-disk build when it matches this configuration ----
      // The build is 31.3 s and 229 MB; the file also carries a key over the build
      // parameters AND 32 sampled N/dNdz pairs, so a different atmosphere file or
      // observation level is a MISS (rebuild), never a wrong answer.
      std::string const cacheDir = tableCacheDir_;
      std::string const geomFile =
          cacheDir.empty() ? std::string{} : cacheDir + "/curved_ray_table.bin";
      auto const geomKey =
          CurvedRayTable::buildKey(ri, earthRadiusM_, rObs, hTop, 800, 10.0, 100.0, 30000.0, 0.02);
      std::shared_ptr<CurvedRayTable> tab;
      if (!geomFile.empty()) {
        auto cand = std::make_shared<CurvedRayTable>();
        std::string why;
        auto const t0 = std::chrono::steady_clock::now();
        if (cand->loadFromFile(geomFile, ri, geomKey, &why)) {
          double const dt =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
          CORSIKA_LOGGER_INFO(logger_, "CurvedRayTable LOADED from cache: {:.2f} s, {:.0f} MB ({})",
                              dt, cand->memoryMB(), geomFile);
          tab = std::move(cand);
        } else {
          CORSIKA_LOGGER_INFO(logger_, "CurvedRayTable cache miss ({}): rebuilding", why);
        }
      }
      if (!tab) {
      auto const tBuild0 = std::chrono::steady_clock::now();
      tab = std::make_shared<CurvedRayTable>(ri, earthRadiusM_, rObs, hTop, 800, 10.0, 100.0,
                                                 30000.0, 0.02);
      double const buildS =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - tBuild0).count();
      CORSIKA_LOGGER_INFO(logger_,
                          "CurvedRayTable build: {:.1f} s, {:.0f} MB (one-time, "
                          "process lifetime)",
                          buildS, tab->memoryMB());
      if (tab->valid() && !geomFile.empty()) {
        if (tab->saveToFile(geomFile, geomKey)) {
          CORSIKA_LOGGER_INFO(logger_, "CurvedRayTable cached to {}", geomFile);
        } else {
          CORSIKA_LOGGER_WARN(logger_, "CurvedRayTable cache write FAILED ({}); continuing",
                              geomFile);
        }
      }
      }  // if (!tab)
      if (!tab->valid()) {
        st.available = false;
        st.unavailableReason = "CurvedRayTable build failed";
        return;
      }
      tab->setGrazingClearance(grazingClearanceM_);  // WAVE 5 lookup guard (not in the cache)
      st.table = std::move(tab);
      st.available = true;
      CORSIKA_LOGGER_INFO(logger_,
                          "CurvedRayTable built: hNodes={} muNodes={} memory={:.1f} MB "
                          "muMin={:.2e} muBreak={:.6f}",
                          st.table->hNodes(), st.table->muNodes(), st.table->memoryMB(),
                          st.table->muMin(), st.table->muBreak());

      // ---- wave 2c: the absorption table, (h, xi, lambda) --------------------
      // 65 x 56 x 256 = 3.6 MB against 67 MB for the geometry table's five tau nodes, and
      // measured 3.8e-4 against an independent quadrature (the in-table version is
      // 3.4e-3).  Built with the same extinction function the geometry table uses, and
      // with the dispersion model of this configuration, because xi = sin(theta_arrival)
      // is only exact when p_crit carries the photon's own f(lambda).
      if (absorptionActive) {
        std::vector<double> absLam;
        for (double l = 300.0; l <= 620.0 + 1.0e-9; l += 5.0) absLam.push_back(l);
        // The dispersion callback is needed twice (key + build), so name it once.
        auto const dispF = [this](double lam) {
          double const Nr = refractiveIndexTable_->N(groundZ_m_);
          if (!refractiveIndexFunc_ || !(Nr > 1.0e-12)) return 1.0;
          double const hRef = groundZ_m_ - earthRadiusM_;
          return (refractiveIndexFunc_(hRef, lam).first - 1.0) / Nr;
        };
        // ---- wave 5: same on-disk cache, keyed on sampled alpha() and f(lambda) ----
        std::string const absFile =
            cacheDir.empty() ? std::string{} : cacheDir + "/curved_abs_table.bin";
        auto const absKey = CurvedAbsorptionTable::buildKey(earthRadiusM_, rObs, absLam,
                                                            extinctionCoeffFunc_, dispF);
        std::shared_ptr<CurvedAbsorptionTable> aTab;
        if (!absFile.empty()) {
          auto cand = std::make_shared<CurvedAbsorptionTable>();
          std::string why;
          auto const t0 = std::chrono::steady_clock::now();
          if (cand->loadFromFile(absFile, ri, absKey, &why)) {
            double const dt =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            CORSIKA_LOGGER_INFO(logger_,
                                "CurvedAbsorptionTable LOADED from cache: {:.2f} s, {:.2f} MB ({})",
                                dt, cand->memoryMB(), absFile);
            aTab = std::move(cand);
          } else {
            CORSIKA_LOGGER_INFO(logger_, "CurvedAbsorptionTable cache miss ({}): rebuilding", why);
          }
        }
        if (!aTab) {
          auto const tA0 = std::chrono::steady_clock::now();
          aTab = std::make_shared<CurvedAbsorptionTable>(ri, earthRadiusM_, rObs,
                                                         extinctionCoeffFunc_, absLam, dispF);
          double const aS =
              std::chrono::duration<double>(std::chrono::steady_clock::now() - tA0).count();
          CORSIKA_LOGGER_INFO(logger_,
                              "CurvedAbsorptionTable built: {:.1f} s, {:.2f} MB, h={} xi={} "
                              "lam={} filled={}/{}",
                              aS, aTab->memoryMB(), aTab->hNodes(), aTab->xiNodes(),
                              aTab->lambdaNodes(), aTab->filledNodes(),
                              aTab->hNodes() * aTab->xiNodes() * aTab->lambdaNodes());
          if (aTab->valid() && !absFile.empty()) {
            if (aTab->saveToFile(absFile, absKey)) {
              CORSIKA_LOGGER_INFO(logger_, "CurvedAbsorptionTable cached to {}", absFile);
            } else {
              CORSIKA_LOGGER_WARN(logger_,
                                  "CurvedAbsorptionTable cache write FAILED ({}); continuing",
                                  absFile);
            }
          }
        }
        if (aTab->valid()) {
          st.absTable = aTab;
        } else {
          CORSIKA_LOGGER_WARN(logger_,
                              "CurvedAbsorptionTable build failed: keeping the geometry "
                              "table's tau");
        }
      }
    });

    if (!st.available) {
      ++st.bypassed;
      bool const ok = legacy();
      if (ok) { ++st.bypassedLegacySucceeded; } else { ++st.bypassedLegacyFailed; ++st.finalDrop; }
      return ok;
    }

    // ---- phase 3: the lookup -----------------------------------------------
    // Dispersion factors follow the production convention n = 1 + N * f, i.e.
    // f = (n - 1) / N, with N from the same refractive-index table (indexed by the
    // ABSOLUTE z = radius). refractiveIndexFunc_ takes the HEIGHT above ground.
    double const rEmit = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
    // refractiveIndexFunc_ takes HEIGHT ABOVE SEA LEVEL (that is the convention used by
    // rayTraceSpherical: h = r - earthRadiusM_).  This used the height above the
    // observation level, so with a 3 km detector the emission was described 3 km too
    // high and the dispersion factor came out ~1.45 instead of 1: every ray was then
    // rejected as DispersionOutOfRange (measured: 22/22 eligible cells in the Telescope
    // acceptance set, and the same silent bypass in Ground mode whenever groundZ != R).
    double const hEmit = rEmit - earthRadiusM_;
    double const Nemit = refractiveIndexTable_->N(rEmit);
    std::pair<double, double> ng = refractiveIndexFunc_
                                       ? refractiveIndexFunc_(hEmit, lambdaNm)
                                       : std::make_pair(refractiveIndex_, refractiveIndex_);
    double const fPhase = (Nemit > 1e-12) ? (ng.first - 1.0) / Nemit : 1.0;
    double const fGroup = (Nemit > 1e-12) ? (ng.second - 1.0) / Nemit : 1.0;

    ++st.attemptedLookups;
    CurvedRayTable::Result const r =
        st.table->trace(x0, y0, z0, dx, dy, dz, fPhase, fGroup);
    // ---- wave 2c/4b: tau from the absorption table, exact chord as last resort -----
    // xi = sin(theta_arrival) is the photon's own Bouguer invariant normalised by
    // p_crit(lambda), which is exactly the table's second axis.  The table cannot serve
    // a ray whose observer sits inside the tangency of that column (T diverges on the
    // far side of it); those rays - measured 4 of 37 on the frozen acceptance set - take
    // the exact straight-chord integral instead, and are COUNTED (absFallbackTau), never
    // silently zero.  With the chord's ambiguity removed (it is now parametrised by the
    // chord length) that path is good to 0.03-0.7 % on the same four rows.
    double tauRoute = 0.0;
    double xiTable = -1.0;
    bool visPerigeeTable = false;
    if (st.absTable) {
      double const xi = std::fabs(std::sin(r.theta_obs));
      // Does the photon's path descend to the perigee and climb back?  Emitting with a
      // negative radial component and arriving with a positive one means the perigee was
      // visited, and then BOTH legs accumulate tau (the table's cumulative is referenced
      // at the perigee).  The table cannot know this from two heights; the tracer can.
      double const rImp = std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
      double const vRadEmit = (rEmit > 0.0) ? (x0 * dx + y0 * dy + z0 * dz) / rEmit : 0.0;
      double const vRadArr = (rImp > 0.0) ? (r.x * r.dx + r.y * r.dy + r.z * r.dz) / rImp : 0.0;
      bool const visitsPerigee = (vRadEmit < 0.0 && vRadArr > 0.0);
      xiTable = xi;
      visPerigeeTable = visitsPerigee;
      auto const ta = st.absTable->tauBetween(rEmit - earthRadiusM_,
                                              groundZ_m_ - earthRadiusM_, xi, lambdaNm,
                                              visitsPerigee);
      if (ta.ok()) {
        tauRoute = ta.tau;
        ++st.absHits;
        if (visitsPerigee) ++st.absSumRule;
      } else {
        ++st.absFallbacks;
        // Endpoints of the chord the table refused: the emission and the sphere impact.
        double const cx = r.x - x0, cy = r.y - y0, czc = r.z - z0;
        double const chordL = std::sqrt(cx * cx + cy * cy + czc * czc);
        tauRoute = absorptionTable_.getOpticalDepthChord(lambdaNm, hEmit,
                                                        std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z) -
                                                            earthRadiusM_,
                                                        chordL, earthRadiusM_);
      }
    } else if (absorptionActive) {
      // Absorption is on but the absorption table was not built (its WARN is emitted at
      // build time).  Nothing is invented: tau stays 0 and is counted SEPARATELY from
      // the chord fallback, so "0" can never be mistaken for "small".
      ++st.absMissing;
    }
    if (r.status == CurvedRayTable::Status::Hit &&
        projectionMode_ == ProjectionMode::Telescope) {
      // 2g-h2n: launch-and-continue.  The table answers for the r = r_obs SPHERE; the
      // telescope needs the crossing of its PLANE.  Projecting the sphere impact onto the
      // plane is wrong (a different point), but the state the table returns -- impact
      // POSITION and arrival DIRECTION -- is a launch point for the certified plane
      // tracer.  For photons that hit the dish the remaining path is short (147 m at
      // 3100 m / 10 km, ~R otherwise), so this costs one to three RK4 steps instead of
      // the whole 82 us trace, and nothing is dropped: there is no pre-check gate on
      // this route, so rays that curve into the dish cannot be lost.
      ++st.hits;
      // split by branch: emission above r_obs is answered by the outward block, below
      // it by the inward one (the population the below-observation work exists for)
      if (std::sqrt(x0 * x0 + y0 * y0 + z0 * z0) >= groundZ_m_) {
        ++st.tableHitAbove;
      } else {
        ++st.tableHitBelow;
      }
      ++st.attemptedLookups;
      double op = r.optical_path_m;
      double gp = r.geometric_path_m;
      double tau = tauRoute; // absorption table, or the counted chord last resort
      double dcx = 0.0, dcy = 0.0, dcz = 0.0;
      bool const cont = rayTraceSphericalToPlane(
          r.x, r.y, r.z, r.dx, r.dy, r.dz, telescopeX_m_, telescopeY_m_, telescopeZ_m_,
          telescopeBasisZ_x_, telescopeBasisZ_y_, telescopeBasisZ_z_, x_ground, y_ground,
          z_ground, op, gp, lambdaNm, tau, dcx, dcy, dcz, 100.0);
      optical_path = op;
      geometric_path = gp;
      accumulated_tau = tau;
      usedTable = true;

      // ---- step A: per-photon A/B against the plane integrator ---------------
      // Diagnostic only: the production result above is returned unchanged, the second
      // trace is discarded.  tau_table is the TABULATED value before this continuation
      // added its own segment; the two chord columns are the exact integrals (§24) to
      // each route's own impact point, which separates an optical-depth error from a
      // geometry error.
      if (abRows_ < abMaxRows_) {
        double xl = 0.0, yl = 0.0, zl = 0.0, opl = 0.0, gpl = 0.0, taul = 0.0;
        double dlx = 0.0, dly = 0.0, dlz = 0.0;
        bool const okL = rayTraceSphericalToPlane(
            x0, y0, z0, dx, dy, dz, telescopeX_m_, telescopeY_m_, telescopeZ_m_,
            telescopeBasisZ_x_, telescopeBasisZ_y_, telescopeBasisZ_z_, xl, yl, zl, opl,
            gpl, lambdaNm, taul, dlx, dly, dlz, 100.0);
        // WAVE 5 probe: the legacy route's own xi and perigee flag, from ITS arrival state.
        // Compared against the table route's on the same photon -- if they differ, that is the
        // whole remaining W2 discrepancy and it is a convention bug, not physics.
        if (okL && xiTable >= 0.0) {
          double const rL = std::sqrt(xl * xl + yl * yl + zl * zl);
          double const vRadArrL = (rL > 0.0) ? (xl * dlx + yl * dly + zl * dlz) / rL : 0.0;
          double const xiLegacy = std::sqrt(std::max(0.0, 1.0 - vRadArrL * vRadArrL));
          double const rEmitL = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
          double const vRadEmitL =
              (rEmitL > 0.0) ? (x0 * dx + y0 * dy + z0 * dz) / rEmitL : 0.0;
          bool const visL = (vRadEmitL < 0.0 && vRadArrL > 0.0);
          double const xiD = std::abs(xiTable - xiLegacy);
          absorptionWeights_.xiDiffSum += xiD;
          if (xiD > absorptionWeights_.xiDiffMax) absorptionWeights_.xiDiffMax = xiD;
          ++absorptionWeights_.xiCmpCount;
          if (visL != visPerigeeTable) ++absorptionWeights_.xiFlagMismatch;
        }
        TraceABRow row;
        row.hEmit = hEmit;
        row.mu = std::cos(r.theta_obs);
        row.lambdaNm = lambdaNm;
        row.statusTable = static_cast<int>(r.status);
        row.xTable = x_ground;
        row.yTable = y_ground;
        row.xLegacy = xl;
        row.yLegacy = yl;
        row.legacyOk = okL ? 1 : 0;
        row.Ltable = gp;
        row.Llegacy = gpl;
        row.dLrel = (gpl > 0.0) ? (gp - gpl) / gpl : 0.0;
        row.dImpact = okL ? std::sqrt((x_ground - xl) * (x_ground - xl) +
                                     (y_ground - yl) * (y_ground - yl) +
                                     (z_ground - zl) * (z_ground - zl))
                          : -1.0;
        row.tauTable = tauRoute;
        row.tauLegacy = taul;
        double const rT =
            std::sqrt(x_ground * x_ground + y_ground * y_ground + z_ground * z_ground);
        double const chT = std::sqrt((x_ground - x0) * (x_ground - x0) +
                                     (y_ground - y0) * (y_ground - y0) +
                                     (z_ground - z0) * (z_ground - z0));
        row.tauChordTable = absorptionTable_.getOpticalDepthChord(
            lambdaNm, hEmit, rT - earthRadiusM_, chT, earthRadiusM_);
        if (okL) {
          double const rL = std::sqrt(xl * xl + yl * yl + zl * zl);
          double const chL = std::sqrt((xl - x0) * (xl - x0) + (yl - y0) * (yl - y0) +
                                       (zl - z0) * (zl - z0));
          row.tauChordLegacy = absorptionTable_.getOpticalDepthChord(
              lambdaNm, hEmit, rL - earthRadiusM_, chL, earthRadiusM_);
        }
        // step A' referee: exact quadrature along the refracted ray to the observation
        // level (the same endpoints the tabulated tau covers)
        double const xiAB = std::fabs(std::sin(r.theta_obs));
        double const vRadEmitAB = (rEmit > 0.0) ? (x0 * dx + y0 * dy + z0 * dz) / rEmit : 0.0;
        row.tauRefQuad = absorptionTauQuadrature(
            lambdaNm, hEmit, groundZ_m_ - earthRadiusM_, xiAB, vRadEmitAB < 0.0);
        // WAVE 5: the same integral at 500 samples.  In u = sqrt(h - h_t) the integrand is
        // smooth, so midpoint convergence should be fast -- if 500 matches 20000, the
        // tangency tail can be integrated for ~1/40 of the cost, which is what makes the
        // widened geometry band affordable.
        row.tauRefQuadLo = absorptionTauQuadrature(
            lambdaNm, hEmit, groundZ_m_ - earthRadiusM_, xiAB, vRadEmitAB < 0.0, 500);
        writeTraceAB_(row);
        ++abRows_;
      }
      return cont;
    }
    if (r.status == CurvedRayTable::Status::Hit) {
      ++st.hits;
      x_ground = r.x;
      y_ground = r.y;
      z_ground = r.z;
      // optical_path_m is the GROUP optical path L + f_group*Sigma - the same
      // quantity rayTraceSpherical() returns - and the caller derives the
      // propagation time from it. The geometric+delay route is the `else` of
      // used_ray_trace, so the delay is applied EXACTLY ONCE.
      optical_path = r.optical_path_m;
      geometric_path = r.geometric_path_m;
      // tau along the same ray: CurvedAbsorptionTable, or the counted chord last
      // resort.  Zero only when absorption is disabled (then it is meaningless, not
      // missing) or when the absorption table failed to build - which is counted in
      // absFallbackTau and warned at build time, never silently absorbed.
      accumulated_tau = tauRoute;
      usedTable = true;
      return true;
    }
    ++st.rejected[static_cast<int>(r.status)];
    bool const ok = legacy();
    if (ok) {
      ++st.legacyFallbackSucceeded;
      // ---- WAVE 5 DECOUPLING: exact-tracer POSITION, table TAU ---------------------
      // The table refused the GEOMETRY (grazing / upgoing / beyond-critical), so the plane
      // tracer supplies the position -- 0.50 m class there, against the 4-66 m the table
      // would have produced.  It must NOT supply tau, though: its integral equals the CHORD's
      // (measured 0.998-1.000 in every mu bin), worth 2.4 pp of absorbed fraction.  The
      // absorption table is a separate object and is valid for the same ray, so tau comes
      // from it -- on the endpoints the tracer actually crossed, with xi from the ARRIVAL
      // state (both measured to agree with the table route: dxi 8e-5, 0 flag mismatches).
      // Measured cost of this decoupling: 0.06 pp (9.81 -> 9.87 %), i.e. free.
      if (st.absTable && extinctionCoeffFunc_ && refractiveIndexTable_) {
        double const rEm = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
        double const rIm =
            std::sqrt(x_ground * x_ground + y_ground * y_ground + z_ground * z_ground);
        double const vRadEmit = (rEm > 0.0) ? (x0 * dx + y0 * dy + z0 * dz) / rEm : 0.0;
        double const vRadArr =
            (rIm > 0.0)
                ? (x_ground * ground_dir_x + y_ground * ground_dir_y + z_ground * ground_dir_z) /
                      rIm
                : 0.0;
        double const xi = std::sqrt(std::max(0.0, 1.0 - vRadArr * vRadArr));
        bool const visitsPerigee = (vRadEmit < 0.0 && vRadArr > 0.0);
        auto const ta = st.absTable->tauBetween(rEm - earthRadiusM_, rIm - earthRadiusM_, xi,
                                                lambdaNm, visitsPerigee);
        if (ta.ok()) {
          accumulated_tau = ta.tau;
          ++st.absHits;
        } else {
          // The table refuses inside the tangency (T diverges there).  The quadrature handles
          // exactly that region by the u = sqrt(h - h_t) substitution, and 500 samples match
          // 20 000 to 0.039 % (measured, 2 000 photons) at 1/40 the cost -- so the accurate tau
          // costs ~0.05 s per event on the whole refused population.
          accumulated_tau = absorptionTauQuadrature(lambdaNm, rEm - earthRadiusM_,
                                                    rIm - earthRadiusM_, xi, vRadEmit < 0.0,
                                                    500);
          ++st.absFallbacks;
        }
      }
    } else {
      ++st.legacyFallbackFailed;
      ++st.finalDrop;
    }
    return ok;
  }

  inline CherenkovProcess::CurvedTableStats CherenkovProcess::curvedTableStats() const {
    CurvedTableStats s;
    s.gatePassed = gatePassed_.load();
    s.gateRejected = gateRejected_.load();
    s.upgoingPreCut = upgoingPreCut_.load();
    s.escapedAtmosphere = escapedAtmosphere_.load();
    s.maxStepsReached = maxStepsReached_.load();
    if (!curvedState_) {
      return s;
    }
    CurvedTableState const& st = *curvedState_;
    s.bypassed = st.bypassed.load();
    s.attempted = st.attemptedLookups.load();
    s.hits = st.hits.load();
    for (int k = 0; k < 9; ++k) s.rejected[k] = st.rejected[k].load();
    s.legacyFallbackSucceeded = st.legacyFallbackSucceeded.load();
    s.legacyFallbackFailed = st.legacyFallbackFailed.load();
    s.bypassedLegacySucceeded = st.bypassedLegacySucceeded.load();
    s.bypassedLegacyFailed = st.bypassedLegacyFailed.load();
    s.finalDrop = st.finalDrop.load();
    s.tableHitAbove = st.tableHitAbove.load();
    s.tableHitBelow = st.tableHitBelow.load();
    s.absTableTau = st.absHits.load();
    s.absSumRule = st.absSumRule.load();
    s.absFallbackTau = st.absFallbacks.load();
    s.tauMissing = st.absMissing.load();
    s.tableAvailable = st.available;
    s.tableAttempted = st.attempted;
    s.unavailableReason = st.unavailableReason;
    return s;
  }

  inline void CherenkovProcess::resetCurvedTableStats() const {
    if (!curvedState_) return;
    CurvedTableState& st = *curvedState_;
    st.bypassed = 0;
    st.attemptedLookups = 0;
    st.hits = 0;
    for (int k = 0; k < 9; ++k) st.rejected[k] = 0;
    st.legacyFallbackSucceeded = 0;
    st.legacyFallbackFailed = 0;
    st.bypassedLegacySucceeded = 0;
    st.bypassedLegacyFailed = 0;
    st.finalDrop = 0;
    st.gatePassed = 0;
    st.gateRejected = 0;
    st.absHits = 0;
    st.absSumRule = 0;
    st.absFallbacks = 0;
    st.absMissing = 0;
    st.tableHitAbove = 0;
    st.tableHitBelow = 0;
    gatePassed_ = 0;
    gateRejected_ = 0;
    upgoingPreCut_ = 0;
    escapedAtmosphere_ = 0;
    maxStepsReached_ = 0;
  }
  // ---- Phase 3: tracer counters -------------------------------------------------
  inline TracerStats CherenkovProcess::tracerStats() const {
    TracerStats s;
    s.accepted = trAccepted_.load();
    s.rejected = trRejected_.load();
    s.fEvaluations = trF_.load();
    s.minStepFailures = trMinStepFail_.load();
    s.cappedSteps = trCapped_.load();
    s.crossingIterations = trCrossIters_.load();
    s.dopriRays = trDopriRays_.load();
    s.rk4Rays = trRk4Rays_.load();
    s.tangentRays = trTangentRays_.load();
    s.tangentUncertified = trTangentUncertified_.load();
    s.tangentFallbacks = trTangentFallbacks_.load();
    s.impactRejections = trImpactRej_.load();
    s.impactAccepted = trImpactOk_.load();
    s.hMin = trHMin_.load();
    s.hMax = trHMax_.load();
    // Mean DOPRI5 step from the accumulated arc length would need another counter; the
    // per-ray flush already records min/max, and the mean is available from
    // accepted/arc length in the benchmark harness.
    s.hMean = (s.accepted > 0) ? s.hMax : 0.0;
    return s;
  }

  inline void CherenkovProcess::resetTracerStats() const {
    trAccepted_ = 0;
    trRejected_ = 0;
    trF_ = 0;
    trMinStepFail_ = 0;
    trCapped_ = 0;
    trCrossIters_ = 0;
    trDopriRays_ = 0;
    trRk4Rays_ = 0;
    trTangentRays_ = 0;
    trTangentUncertified_ = 0;
    trTangentFallbacks_ = 0;
    trImpactRej_ = 0;
    trImpactOk_ = 0;
    impactUncertified_ = 0;
    trHMin_ = 1.0e300;
    trHMax_ = 0.0;
  }



  inline void CherenkovProcess::logCurvedTableStats() const {
    CurvedTableStats const s = curvedTableStats();
    std::size_t const gated = s.gatePassed + s.gateRejected;
    CORSIKA_LOGGER_INFO(
        logger_,
        "Cherenkov curved route: gate {}/{} pass ({:.3}%), table hits {} (above {}, below {}), "
        "rejected[upgoing={} dispersion={} grazing={} belowTable={} aboveTable={}], "
        "legacy after rejection {} ok/{} fail, blind legacy {} ok/{} fail, dropped {}",
        s.gatePassed, gated, gated ? 100.0 * double(s.gatePassed) / double(gated) : 0.0, s.hits,
        s.tableHitAbove, s.tableHitBelow, s.rejected[5], s.rejected[7], s.rejected[8],
        s.rejected[3], s.rejected[4], s.legacyFallbackSucceeded, s.legacyFallbackFailed,
        s.bypassedLegacySucceeded, s.bypassedLegacyFailed, s.finalDrop);
    // tau provenance for the same population (wave 4b: the geometry table no longer
    // carries any).  These count every ray that reached the lookup, i.e. including the
    // few whose GEOMETRY the table then rejects (those go to the legacy tracer with a
    // tau that nobody samples) -- so absTableTau + absFallbackTau is the number of
    // traced rays, not the number of served ones.  tauMissing counts rays where
    // absorption is on but no absorption table exists, i.e. a build failure: they get
    // tau = 0 and are never silent.
    std::size_t const withTau = s.absTableTau + s.absFallbackTau;
    CORSIKA_LOGGER_INFO(logger_,
                        "Cherenkov curved route tau: absorption table {}/{} (dip/sum rule {}), "
                        "exact chord (last resort) {}/{}, missing (no abs table) {}",
                        s.absTableTau, withTau, s.absSumRule, s.absFallbackTau, withTau,
                        s.tauMissing);
    // Phase 2: cost accounting.  upgoingPreCut is a subset of gatePassed (Telescope
    // photons refused before tracing because they move outward from more than one
    // aperture above the dish); escaped/maxSteps are tracer failures that used to emit
    // one warning per photon.
    CORSIKA_LOGGER_INFO(
        logger_,
        "Cherenkov tracer cost: upgoing pre-cut {} (of {} gated), escaped atmosphere {}, "
        "max-steps {}",
        s.upgoingPreCut, s.gatePassed, s.escapedAtmosphere, s.maxStepsReached);
  }

  inline bool CherenkovProcess::curvedTableActive() const {
    return enableCurvedTable_ && curvedState_ && curvedState_->available;
  }

  inline std::vector<CurvedAbsorptionTable::Rejection>
  CherenkovProcess::absorptionRejections() const {
    if (!curvedState_ || !curvedState_->absTable) return {};
    return curvedState_->absTable->recentRejections();
  }

  inline bool CherenkovProcess::rayTraceSpherical(
    double x0, double y0, double z0_altitude,
    double dx, double dy, double dz,
    double& x_ground, double& y_ground, double& z_ground,
    double& optical_path, double& geometric_path,
    double lambdaNm,
    double& accumulated_tau,
    double maxRadialStep_m) const {

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

    // --- target classification (2g-h2k, step 2) -----------------------------
    // The target is the radius groundZ_m_.  The previous "already at or below
    // ground" short circuit returned tau = 0 whenever the target was at or above
    // the start, which silently zeroed every upward and same-altitude ray (17 of
    // 42 cells of the frozen set).  Classify the target explicitly instead:
    //   targetDir = -1 : target below the start -> first descending crossing
    //   targetDir = +1 : target above the start -> first ascending crossing
    //   targetDir =  0 : target at the start radius (|dh| < kSameAltitude_m)
    //                    -> the SECOND crossing, where the ray returns to that
    //                    radius after dipping below it
    constexpr double kSameAltitude_m = 1.0e-3;
    double const dh_target = h_ground - (r0 - R_E);
    int const targetDir =
        (dh_target > kSameAltitude_m) ? 1 : ((dh_target < -kSameAltitude_m) ? -1 : 0);

    // Step 3b (2g-h2k) quadrature: the extinction integral is sampled on
    // sub-intervals fixed in altitude (10 m of radius), with a floor set by the
    // transverse length (50 m) so that near-horizontal steps are not
    // under-sampled.  Both are independent of the geometric step size.
    constexpr double kTauSubStep_m = 10.0;
    constexpr double kTauSubStepTransverse_m = 50.0;

    // Radial unit vector at emission point (points outward/up from Earth center)
    double const r_hat_x = px / r0;
    double const r_hat_y = py / r0;
    double const r_hat_z = pz / r0;

    // cos(theta) = dot(direction, radial)
    double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
    double const sin_theta0 = std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

    // Phase 2: a photon emitted above the target sphere and moving outward has a
    // monotonically increasing radius (n*r monotonic, as the table's constructor
    // requires), so it can never reach it.  Provable, not a heuristic: this replaces
    // integrating the ray to the 100 km ceiling and returning "escaped".  Rays emitted
    // below the target (Earth-skimming) and the same-altitude family are untouched.
    if (r0 > r_ground + kSameAltitude_m && cos_theta0 >= 0.0) {
      ++upgoingPreCut_;
      return false;
    }

    // Refractive index at emission altitude
    double const h0 = r0 - R_E;
    double const n0 = refractiveIndexFunc_(h0, lambdaNm).first;  // phase index for Bouguer

    // ---- Phase 1: index and group path straight from the refractive-index table ----
    // Same contract as rayTraceSphericalToPlane(): with a separable index the table
    // supplies N, dN/dz and f_phase/f_group, so no refractiveIndexFunc_ call is made
    // per step and dN/dz is the exact log-slope instead of a +-1 m difference.
    bool const useIndexTable = separableIndex_ && refractiveIndexTable_ &&
                               refractiveIndexTable_->valid() && refractiveIndexFunc_;
    double fPhase = 1.0, fGroup = 1.0;
    if (useIndexTable) {
      auto const ng0 = refractiveIndexFunc_(h0, lambdaNm);
      double const Nr0 = refractiveIndexTable_->N(r0);
      if (Nr0 > 1.0e-12) {
        fPhase = (ng0.first - 1.0) / Nr0;
        fGroup = (ng0.second - 1.0) / Nr0;
      } else {
        fPhase = 0.0;
        fGroup = 0.0;
      }
    }

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

    // Plane frame convention: the propagation plane is kept CANONICAL - x_plane
    // decreases (alpha <= 0) for every azimuth, which is the branch the dr/dz_step
    // formulas were derived in. The plane -> physical mapping is the single fixed
    // mirror  pos = z_plane*r_hat - x_plane*t_hat  applied at both boundaries (the
    // straight-line basis below and the curved-crossing conversion). Since t_hat is
    // defined as the tangential part of the photon direction, the ray always moves
    // along +t_hat, so the physical transverse displacement -x_plane must be
    // positive. Deriving the sign from t_hat_y (the previous az_sign hack) put the
    // plane on the wrong branch for t_hat_y < 0 and gave an 8.1% position error.

    // 2D coordinates in the propagation plane (Earth-centered Cartesian)
    double x_plane = 0.0;
    double z_plane = r0;
    double r_cur = r0;

    // --- Step 3c (2g-h2k): unit direction in the propagation plane ------------
    // The direction is carried as state and advanced by RK4 together with the
    // position, instead of being reconstructed at the left end of every step from
    // the Bouguer invariant.  The plane -> physical mapping is
    // pos = z_plane*r_hat - x_plane*t_hat, and t_hat is built from the tangential
    // part of the photon direction (d.t_hat = sin(theta0) > 0), so
    //     d pos/ds = u_z*r_hat - u_x*t_hat  ==  cos(theta0)*r_hat + sin(theta0)*t_hat
    // gives u_z = cos(theta0) and u_x = -sin(theta0).
    double u_x_plane = -sin_theta0;
    double u_z_plane = cos_theta0;

    optical_path = 0.0;
    geometric_path = 0.0;
    accumulated_tau = 0.0;

    int const max_steps = 10000; 
    double const max_altitude = 100000.0;  // 100 km maximum

    double x_prev = x_plane;
    double z_prev = z_plane;
    double r_prev = r_cur;
    double up_x = u_x_plane, up_z = u_z_plane;  // direction at the step start (Hermite)
    double optical_prev = 0.0;
    double geometric_prev = 0.0;
    double tau_prev = 0.0;  // accumulated_tau at the start of the current step
    bool dipped_below = false; // targetDir == 0: has the ray been below the target?

    for (int step = 0; step < max_steps; ++step) {
      double const h_cur = r_cur - R_E;
      // Check if we've escaped atmosphere
      if (h_cur > max_altitude) {
        // Phase 2: counted always, logged only for the first three (see the plane tracer).
        std::size_t const nEsc = escapedAtmosphere_.fetch_add(1) + 1;
        if (nEsc <= 3) {
          double const photon_theta_rad = std::acos(std::max(-1.0, std::min(1.0, -dz)));
          double const photon_theta_deg = photon_theta_rad * 180.0 / M_PI;
          CORSIKA_LOGGER_WARN(logger_, "rayTraceSpherical: escaped atmosphere at step {}, h_cur={:.1f}m > max_altitude={:.1f}m, photon zenith=({:.2f} deg) [logged for the first 3 only]", 
                             step, h_cur, max_altitude, photon_theta_deg);
        }
        return false;
      }

      // Phase 1: only the group index for the non-table fallback is needed here; the
      // phase index and its derivative are evaluated PER RK4 STAGE (see kappaAt).
      double n_group = 1.0;
      if (!useIndexTable) {
        n_group = refractiveIndexFunc_(h_cur, lambdaNm).second;
      }
      // --- Step 2h (2g-h2k): LOCAL zenith angles from the state --------------
      // u_x_plane / u_z_plane are Cartesian components in the FIXED propagation
      // plane frame (pos = z_plane*r_hat0 - x_plane*t_hat0), not the local
      // (transverse, radial) pair.  Reading them as (sin, cos) of the zenith angle
      // is only correct at the emission azimuth: once the ray has swept an angle
      // phi, u_z = cos(chi) with chi = psi + phi, so a ray that is still climbing
      // shows u_z < 0 and was mistaken for a descending one.  On grazing rays
      // (phi up to ~1 deg over 100 km) that misfired the near-ground chord branch
      // and replaced tens of km of refracted path with a vacuum chord (measured:
      // -1.3 % to -8.6 % in the crossing path length, both signs).
      // Project the state direction onto the local radial/transverse basis of the
      // CURRENT point instead.  t_hat_plane = (z, -x)/r, so
      //   cos(psi) = u . r_hat_plane ,  sin(psi) = |u x r_hat_plane| (magnitude).
      double const rp_x = x_plane / r_cur;
      double const rp_z = z_plane / r_cur;
      double const cos_theta_signed = u_x_plane * rp_x + u_z_plane * rp_z;  // signed drift
      double const cos_theta_cur = std::min(1.0, std::fabs(cos_theta_signed));
      double const sin_theta_cur = std::min(1.0, std::fabs(u_x_plane * rp_z - u_z_plane * rp_x));


      // --- Straight-line shortcut for the FINAL approach to the target --------
      // The trigger must be the remaining PATH length, not the radial distance: a
      // ray 1 km above the target can still have tens to hundreds of km to go, and
      // the chord is only exact over the last kilometre (refraction over 1 km is
      // ~4e-5 m of sag).  The previous radial test fired over the whole path of
      // every near-horizontal ray (dist_to_ground is NEGATIVE for emitters below
      // the target), which bypassed refraction entirely for all six 3000 m cells.
      // The "genuinely above the target" margin keeps the same-altitude cells
      // (emitter on the sphere, initial drift negative) integrating the dip-and-
      // return branch instead of shortcutting it.
      constexpr double kAboveTarget_m = 1.0e-3;
      double const distRemaining_m =
          (cos_theta_signed < 0.0 && r_cur > r_ground + kAboveTarget_m)
              ? (r_cur - r_ground) / (-cos_theta_signed)
              : std::numeric_limits<double>::max();
      if (distRemaining_m <= 1000.0) {
        // Reconstruct the current 3D direction from the 2D propagation plane state
        // CONVENTION REPAIR (variant B): x_plane is a MIRRORED image of the
        // physical transverse coordinate t_hat. The curved loop below is built in
        // that mirrored frame - its dr/dz_step assume it, which is why its step
        // lengths and accumulated path are correct - so the conversions that leave
        // the frame must apply the mirror. Without it the curved part and the
        // straight-line part come out with opposite signs and partially cancel: the
        // impact point lands on the wrong side of the emission (measured 18.4 km at
        // h = 20 km, mu = 0.9, where the true displacement is 9.655 km).
        // Step 3c: direction and position come straight from the plane state,
        // in the fixed emission frame pos = z_plane*r_hat - x_plane*t_hat.
        double const dir3d_x = u_z_plane * r_hat_x - u_x_plane * t_hat_x;
        double const dir3d_y = u_z_plane * r_hat_y - u_x_plane * t_hat_y;
        double const dir3d_z = u_z_plane * r_hat_z - u_x_plane * t_hat_z;

        double const pos3d_x = z_plane * r_hat_x - x_plane * t_hat_x;
        double const pos3d_y = z_plane * r_hat_y - x_plane * t_hat_y;
        double const pos3d_z = z_plane * r_hat_z - x_plane * t_hat_z;

        // Analytical straight-line intersection with ground sphere: |pos + t*dir|^2 = r_ground^2
        double const pos_dot_dir = pos3d_x * dir3d_x + pos3d_y * dir3d_y + pos3d_z * dir3d_z;
        double const pos_dot_pos = pos3d_x * pos3d_x + pos3d_y * pos3d_y + pos3d_z * pos3d_z;
        double const discriminant = pos_dot_dir * pos_dot_dir - (pos_dot_pos - r_ground * r_ground);
        if (discriminant < 0.0) {
          // Ray does not intersect the ground sphere — miss
          return false;
        }
        double const sqrt_disc = std::sqrt(discriminant);
        // Near / far root per travel direction (Step 2).
        double t_hit = (targetDir < 0) ? (-pos_dot_dir - sqrt_disc)
                                       : (-pos_dot_dir + sqrt_disc);
        if (t_hit < 0.0) {
          t_hit = (targetDir < 0) ? (-pos_dot_dir + sqrt_disc) : (-pos_dot_dir - sqrt_disc);
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
      // One transverse controller for every zenith angle.  The old H_scale/100
      // branch for sin(theta) > 0.99 was a leftover from the pre-3b quadrature:
      // with RK4 on (position, direction) and an altitude-fixed tau quadrature the
      // path is step-insensitive (measured: L changes by 1e-6 relative between
      // H/100 = 80 m and H/20 = 400 m steps on 100-250 km grazing rays), while the
      // 80 m step quintupled the cost of every near-horizontal ray.
      double const H_scale = computeLocalScaleHeight(h_cur);
      double const dx = std::min(H_scale / 20.0, 1000.0);

      // --- Step 3c (2g-h2k): arc length from the same step policy -------------
      // The transverse increment dx is the controller input (H_scale/20 capped at
      // 1000 m).  With |dx_plane/ds| = sin(theta) the arc length is ds = dx/sin(theta),
      // and the radial change per step is |dh| = ds*cos(theta), so the Step-3 radial
      // cap reads ds <= maxRadialStep/cos(theta).
      double ds = (sin_theta_cur > 1.0e-9) ? dx / sin_theta_cur
                                          : (maxRadialStep_m > 0.0 ? maxRadialStep_m : 100.0);
      if (maxRadialStep_m > 0.0 && cos_theta_cur > 1.0e-9) {
        ds = std::min(ds, maxRadialStep_m / cos_theta_cur);
      }

      // Store previous state for interpolation
      x_prev = x_plane;
      z_prev = z_plane;
      r_prev = r_cur;
      up_x = u_x_plane;
      up_z = u_z_plane;
      optical_prev = optical_path;
      geometric_prev = geometric_path;
      tau_prev = accumulated_tau;

      // --- RK4 on (position, direction) in the propagation plane --------------
      // du/ds = (n'/n) (r_hat - (u.r_hat) u) is the ray equation of a spherically
      // symmetric medium.  Unlike the previous left-endpoint reconstruction it needs
      // no sign conventions and no turning-point test: the direction rotates
      // smoothly through a tangent point.  n'/n is taken once per step (it varies by
      // ~1 % over a step at H = 8 km): from the table's stored log-slope when the
      // index is separable, otherwise from a central difference of the index.
      double sigma_inc = 0.0;  // int N ds over this step (Phase 1: 4th-order group path)
      {
        // Phase 1b: stage-local kappa -- see the long note in rayTraceSphericalToPlane.
        // Freezing it at the step start was the dominant geometry error of the tracer.
        auto const kappaAt = [&](double const rr) -> double {
          if (useIndexTable) {
            double const Ncur = refractiveIndexTable_->N(rr);
            double const nph = 1.0 + Ncur * fPhase;
            return (nph > 1.0e-12) ? (fPhase * refractiveIndexTable_->dNdz(rr)) / nph : 0.0;
          }
          double const hh = rr - R_E;
          double const hLo = std::max(0.0, hh - 1.0);
          double const hHi = std::max(1.0e-6, hh + 1.0);
          double const nLo = refractiveIndexFunc_(hLo, lambdaNm).first;
          double const nHi = refractiveIndexFunc_(hHi, lambdaNm).first;
          if (!(nLo > 1.0e-12) || !(nHi > 1.0e-12)) return 0.0;
          return (std::log(nHi) - std::log(nLo)) / (hHi - hLo);
        };
        auto const rayDeriv = [&](double const xq, double const zq, double const vx,
                                  double const vz, double& dxt, double& dzt, double& dvxt,
                                  double& dvzt) {
          double const rq = std::sqrt(xq * xq + zq * zq);
          double const rx = xq / rq;
          double const rz = zq / rq;
          double const ur = vx * rx + vz * rz;
          double const kq = kappaAt(rq);
          dxt = vx;
          dzt = vz;
          dvxt = kq * (rx - ur * vx);
          dvzt = kq * (rz - ur * vz);
        };
        double a1x, a1z, b1x, b1z;
        rayDeriv(x_plane, z_plane, u_x_plane, u_z_plane, a1x, a1z, b1x, b1z);
        double a2x, a2z, b2x, b2z;
        rayDeriv(x_plane + 0.5 * ds * a1x, z_plane + 0.5 * ds * a1z,
                 u_x_plane + 0.5 * ds * b1x, u_z_plane + 0.5 * ds * b1z, a2x, a2z, b2x, b2z);
        double a3x, a3z, b3x, b3z;
        rayDeriv(x_plane + 0.5 * ds * a2x, z_plane + 0.5 * ds * a2z,
                 u_x_plane + 0.5 * ds * b2x, u_z_plane + 0.5 * ds * b2z, a3x, a3z, b3x, b3z);
        double a4x, a4z, b4x, b4z;
        rayDeriv(x_plane + ds * a3x, z_plane + ds * a3z, u_x_plane + ds * b3x,
                 u_z_plane + ds * b3z, a4x, a4z, b4x, b4z);
        x_plane += ds / 6.0 * (a1x + 2.0 * a2x + 2.0 * a3x + a4x);
        z_plane += ds / 6.0 * (a1z + 2.0 * a2z + 2.0 * a3z + a4z);
        u_x_plane += ds / 6.0 * (b1x + 2.0 * b2x + 2.0 * b3x + b4x);
        u_z_plane += ds / 6.0 * (b1z + 2.0 * b2z + 2.0 * b3z + b4z);
        double const un = std::sqrt(u_x_plane * u_x_plane + u_z_plane * u_z_plane);
        u_x_plane /= un;
        u_z_plane /= un;

        // Phase 1: dSigma/ds = N(r) on the same RK4 stages (see the plane tracer).
        if (useIndexTable) {
          auto const NofStage = [&](double const xq, double const zq) {
            return refractiveIndexTable_->N(std::sqrt(xq * xq + zq * zq));
          };
          double const sg1 = NofStage(x_prev, z_prev);
          double const sg2 = NofStage(x_prev + 0.5 * ds * a1x, z_prev + 0.5 * ds * a1z);
          double const sg3 = NofStage(x_prev + 0.5 * ds * a2x, z_prev + 0.5 * ds * a2z);
          double const sg4 = NofStage(x_prev + ds * a3x, z_prev + ds * a3z);
          sigma_inc = ds / 6.0 * (sg1 + 2.0 * sg2 + 2.0 * sg3 + sg4);
        }
      }

      geometric_path += ds;
      double const op_inc = useIndexTable ? (ds + fGroup * sigma_inc) : (n_group * ds);
      optical_path += op_inc;

      // Dense output for everything that happens inside this step: the tau quadrature
      // and the crossing of the target radius (Phase 1; see HermiteStep).
      HermiteStep const seg{x_prev, z_prev, up_x, up_z, x_plane, z_plane,
                            u_x_plane, u_z_plane, ds};

      // --- Accumulate extinction optical depth along the curved ray ---
      // Step 3b (2g-h2k): the sub-intervals are fixed in ALTITUDE (10 m of radius)
      // and bounded below by the transverse length, so the quadrature no longer
      // depends on the geometric step size.  alpha is evaluated at the exact
      // radius of the interpolated position: Phase 1 takes that position from the
      // Hermite dense output (on the arc) instead of the chord, which removes the
      // sagitta-sized radius bias.
      double const r_new = std::sqrt(x_plane * x_plane + z_plane * z_plane);
      if (extinctionCoeffFunc_) {
        int const nByAltitude =
            static_cast<int>(std::ceil(std::abs(r_new - r_prev) / kTauSubStep_m));
        int const nByLength = static_cast<int>(std::ceil(ds / kTauSubStepTransverse_m));
        int const n_tau = std::max(1, std::max(nByAltitude, nByLength));
        double const ds_tau = ds / static_cast<double>(n_tau);
        for (int its = 0; its < n_tau; ++its) {
          double const fr = (static_cast<double>(its) + 0.5) / static_cast<double>(n_tau);
          double sx = 0.0, sz = 0.0;
          seg.pos(fr, sx, sz);
          accumulated_tau += extinctionCoeffFunc_(std::sqrt(sx * sx + sz * sz) - R_E, lambdaNm) *
                             ds_tau;
        }
      }

      // The position (and r_cur) have already been advanced by the RK4 block.
      r_cur = r_new;

      // Detect the crossing of the target radius in the direction of travel.
      // targetDir == 0 requires the ray to dip below the target first and then to
      // return to it (the same-altitude / horizontal case), otherwise the start
      // point itself would be taken as the hit.
      bool hit_target = false;
      if (targetDir < 0) {
        hit_target = (r_cur <= r_ground);
      } else if (targetDir > 0) {
        hit_target = (r_cur >= r_ground);
      } else {
        if (r_cur < r_ground - kSameAltitude_m) dipped_below = true;
        hit_target = dipped_below && (r_cur >= r_ground - kSameAltitude_m);
      }
      if (hit_target) {
        // Sub-step fraction of the crossing.  r(s) is far from linear over a step when
        // the ray grazes the target sphere (a 1 km same-altitude chord climbs 0.03 m
        // per 400 m step), so use the radial drift dr/ds = cos(psi) at BOTH ends of the
        // step and solve the resulting quadratic instead of interpolating r linearly.
        // Measured: linear interpolation cost 5.9e-7 in the n = 1 straight-ray probe
        // and 7e-4 on the 1 km horizontal cell; the quadratic is at round-off.
        double alpha = 0.0;
        {
          double const c0 = cos_theta_signed;  // dr/ds at the step start
          double const c1 = (u_x_plane * x_plane + u_z_plane * z_plane) / r_cur;  // at the end
          double const dang = r_ground - r_prev;
          double const kq = (c1 - c0) / (2.0 * ds);
          bool solved = false;
          if (std::abs(kq) > 1.0e-20 && std::abs(c0) > 1.0e-15) {
            double const disc = c0 * c0 + 4.0 * kq * dang;
            if (disc >= 0.0) {
              double const sq = std::sqrt(disc);
              double const den = (c0 >= 0.0) ? (c0 + sq) : (c0 - sq);
              if (std::abs(den) > 1.0e-15) {
                alpha = 2.0 * dang / (den * ds);
                solved = true;
              }
            }
          }
          if (!solved && std::abs(r_prev - r_cur) > 1.0e-12) {
            alpha = (r_prev - r_ground) / (r_prev - r_cur);
          }
          alpha = std::max(0.0, std::min(1.0, alpha));
        }

        // Phase 1: refine the crossing on the Hermite curve.  The quadratic above
        // treats r(s) as a parabola in the CHORD frame; Newton on F(t) = |P(t)| -
        // r_ground solves on the arc, which is what the photon follows (1-2 mm at the
        // production step sizes, measured; < 0.001 mm after).  Started from the
        // quadratic fraction, so 2 iterations reach round-off.
        double tc = alpha;
        for (int it = 0; it < 6; ++it) {
          double xq = 0.0, zq = 0.0, dxq = 0.0, dzq = 0.0;
          seg.pos(tc, xq, zq);
          seg.dpos(tc, dxq, dzq);
          double const rq = std::sqrt(xq * xq + zq * zq);
          if (!(rq > 0.0)) break;
          double const fq = rq - r_ground;
          double const fpq = (xq * dxq + zq * dzq) / rq;
          if (std::fabs(fpq) < 1.0e-300) break;
          double const dt = fq / fpq;
          double const tn = tc - dt;
          if (!(tn >= 0.0 && tn <= 1.0)) break;
          tc = tn;
          if (std::fabs(dt) <= 1.0e-13) break;
        }
        if (!(tc >= 0.0 && tc <= 1.0)) tc = alpha;

        // Partial step to the target (same step-3b quadrature: fixed in altitude,
        // evaluated at the exact radius of the interpolated position).
        // Step 2i (2g-h2k): the main quadrature above already added the FULL step, so
        // it must be dropped here -- the closing partial segment is the only part that
        // counts.  Resetting to the step start keeps this exact; adding the partial
        // step on top of the full one counted the last step (1 + alpha) times, a
        // first-order-in-ds tau bias (measured: 4.5e-2 on a 1 km vertical cell with a
        // 100 m cap, and all 9 non-zero step-halving violations).
        accumulated_tau = tau_prev;
        double const ds_partial = tc * ds;
        if (extinctionCoeffFunc_) {
          int const nByAltitude = static_cast<int>(
              std::ceil(std::abs(r_ground - r_prev) / kTauSubStep_m));
          int const nByLength = static_cast<int>(std::ceil(ds_partial / kTauSubStepTransverse_m));
          int const n_tau = std::max(1, std::max(nByAltitude, nByLength));
          double const ds_tau = ds_partial / static_cast<double>(n_tau);
          for (int its = 0; its < n_tau; ++its) {
            double const fr = (static_cast<double>(its) + 0.5) / static_cast<double>(n_tau);
            double sx = 0.0, sz = 0.0;
            seg.pos(tc * fr, sx, sz);
            accumulated_tau +=
                extinctionCoeffFunc_(std::sqrt(sx * sx + sz * sz) - R_E, lambdaNm) * ds_tau;
          }
        }

        seg.pos(tc, x_plane, z_plane);
        r_cur = r_ground;
        optical_path = optical_prev + tc * op_inc;
        geometric_path = geometric_prev + tc * ds;

        // Convert plane coordinates back to 3D world coordinates (Step 9)
        // mirror x_plane back to the physical transverse coordinate (see the
        // convention-repair note above)
        double const gx = -x_plane * t_hat_x + z_plane * r_hat_x;
        double const gy = -x_plane * t_hat_y + z_plane * r_hat_y;
        double const gz = -x_plane * t_hat_z + z_plane * r_hat_z;

        x_ground = gx;
        y_ground = gy;
        z_ground = gz;  // absolute z in Earth-centered coordinates
        return true;
      }

      // Safety: reject only when the ray has GENUINELY passed its tangent point.
      // The previous form compared r_cur against r_prev at the 1e-6 m level, which
      // fires on rounding for grazing rays (the radial change per 400 m step is
      // 0.4 m at |cos(theta)| = 1e-3, but the third derivative noise is not) and
      // rejected 11 otherwise fine cells once the near-ground shortcut was gone.
      // A 1 m margin above the running minimum is unambiguous.
      // Exact criterion, no tolerance: a ray that is rising while still above the
      // target can never reach it.  (A radius-margin test cannot be used: the dip of
      // a 1 km same-altitude chord is only 2 cm, so any margin above ~1e-3 m would
      // reject ordinary horizontal cells.)
      // Exact criterion using the LOCAL radial drift, no tolerance: a ray that is
      // rising away from the target radius while still above it can never reach it.
      if (cos_theta_signed > 0.0 && r_cur > r_ground) {
        return false;  // turned away above the target: unreachable
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

  namespace {
    /// Minimal scope timer used by the optional Cherenkov profiler.
    struct CherenkovProfScope {
      bool on_;
      std::chrono::steady_clock::time_point t0_;
      double* acc_;
      CherenkovProfScope(bool on, double* acc) : on_(on), acc_(acc) {
        if (on_) t0_ = std::chrono::steady_clock::now();
      }
      ~CherenkovProfScope() {
        if (on_) {
          *acc_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0_)
                       .count();
        }
      }
    };
  } // namespace

  // --- apply continuous Cherenkov photon production ---
  template <typename TParticle>
  inline corsika::ProcessReturn
  CherenkovProcess::doContinuous(corsika::Step<TParticle>& step,
                                  bool const /*stepLimit*/) const {

    CherenkovProfScope profTotal(enableProfiling_, &profile_.total_s);

    auto const& projectile = step.getParticlePre();

    auto const E    = projectile.getEnergy();
    auto const mass = projectile.getMass();
    auto const q    = projectile.getChargeNumber();

    if (q == 0) {
      return corsika::ProcessReturn::Ok;
    }

    auto const posPre = step.getPositionPre();
    auto const dir = projectile.getDirection();

    double const pos_x_m = posPre.getX() / 1_m;
    double const pos_y_m = posPre.getY() / 1_m;
    double const pos_z_m = posPre.getZ() / 1_m;
    double const dir_x = dir.getX();
    double const dir_y = dir.getY();
    double const dir_z = dir.getZ();

    handleTelescopePositionCorrection(E / 1_GeV, pos_x_m, pos_y_m, pos_z_m,
                                      dir_x, dir_y, dir_z,
                                      static_cast<int>(projectile.getPID()));

    double const m_over_E = mass / E;
    double const beta2    = 1. - m_over_E * m_over_E;

    if (beta2 <= 0.) {
      return corsika::ProcessReturn::Ok;
    }

    double const beta          = std::sqrt(beta2);
    // Local refractive index at production start height (absolute z in meters)
    // Use reference wavelength 400nm for threshold check (phase index)
    double const n_pre = refractiveIndexFunc_ ?
      refractiveIndexFunc_(pos_z_m - earthRadiusM_, 400.0).first
      : refractiveIndex_;
    double const betaThreshold = 1.0 / n_pre;

    if (beta < betaThreshold) {
      return corsika::ProcessReturn::Ok;
    }

    // Calculate step length
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
    {
    CherenkovProfScope profProd(enableProfiling_, &profile_.production_s);
    if (refractiveIndexFunc_ && enableDispersion_) {
      // Full 2D integration over position and wavelength when dispersion is enabled
      int const nPosSamples = 10;      // spatial samples along step
      int const nLambdaSamples = 10;   // wavelength samples
      double const dLambda_cm = (lambda_max_cm - lambda_min_cm) / nLambdaSamples;
      
      for (int i = 0; i < nPosSamples; ++i) {
        double const frac = (i + 0.5) / nPosSamples;
        auto const posIntegrate = posPre + displacement * frac;
        double const z_integrate = posIntegrate.getZ() / 1_m;
        
        // Wavelength integral: ∫ sin²θ(λ) / λ² dλ
        double wavelengthIntegral = 0.0;
        for (int j = 0; j < nLambdaSamples; ++j) {
          double const lambda_cm = lambda_min_cm + (j + 0.5) * dLambda_cm;
          double const lambda_nm = lambda_cm * 1e7;  // convert cm to nm
          
          // Use phase index for Cherenkov angle calculation (sin²θ)
          double const n_integrate = refractiveIndexFunc_(z_integrate - earthRadiusM_, lambda_nm).first;
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
        double const z_integrate = posIntegrate.getZ() / 1_m;
        // When dispersion disabled, wavelength parameter is ignored; use phase index
        double const n_integrate = refractiveIndexFunc_(z_integrate - earthRadiusM_, 400.0).first;
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
    } // production scope

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

    double const px = dir.getX();
    double const py = dir.getY();
    double const pz = dir.getZ();
    
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
      double const sa_x = showerDir.getX();
      double const sa_y = showerDir.getY();
      double const sa_z = showerDir.getZ();
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

    if (enableProfiling_) { profile_.bunches += static_cast<std::size_t>(nBunches); }
    CherenkovProfScope profLoop(enableProfiling_, &profile_.bunchloop_s);
    for (int i = 0; i < nBunches; ++i) {
      double const frac = posDist(rng);
      auto const posEmit = posPre + displacement * frac;

      double const x0 = posEmit.getX() / 1_m;
      double const y0 = posEmit.getY() / 1_m;
      double const z0 = posEmit.getZ() / 1_m;

      // Sample wavelength from 1/λ² distribution early (needed for ray tracing)
      // p(λ) ∝ 1/λ², so λ = 1 / (1/λ_min - u*(1/λ_min - 1/λ_max))
      double const u_wavelength = posDist(rng);
      double const lambdaNm = 1.0 / (lambda_min_inv - u_wavelength * lambda_inv_range);

      // Get altitude and wavelength-dependent refractive index at emission point
      // Use phase index for Cherenkov angle calculation
      double const n_local = refractiveIndexFunc_ ?
        refractiveIndexFunc_(z0 - earthRadiusM_, lambdaNm).first : refractiveIndex_;
      
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
      AbsorptionState absorptionState = AbsorptionState::NotEvaluated;
      if (projectionMode_ == ProjectionMode::Ground) {
        if (enableCurvature_ && refractiveIndexFunc_) {
          // Table by default for eligible curved-ground rays; the UNCHANGED legacy
          // tracer runs for ineligible rays and after a table rejection (accounted
          // separately, see traceCurvedGround).
          bool used_table_route = false;
          bool hit = traceCurvedGround(x0, y0, z0, dx, dy, dz, lambdaNm,
                                       xg, yg, zg,
                                       ray_optical_path, ray_geom_path,
                                       ray_accumulated_tau, used_table_route);
          if (hit) {
            used_ray_trace = true;
            // The ground tracer integrates tau along the CURVED ray and does not
            // sample absorption; the caller attenuates the weight by exp(-tau).
            absorptionState = AbsorptionState::OpticalDepthAvailable;
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

              // Gate: within 10x the radius, PLUS the refraction displacement the
              // straight line does not know about.  A ray bends by L^2/(2 R_eff) relative
              // to its vacuum chord (R_eff ~ 2.4e7 m for n' = -3.27e-4/8 km), i.e. 2 m at
              // 10 km but 208 m at 100 km and 833 m at 200 km: with the bare 10R margin the
              // gate was DROPPING photons that the curved trace would have put inside the
              // dish.  Measured cost of the wider gate: the acceptance fraction grows as
              // (margin/10R)^2, so the per-photon cost stays at the ~20 ns gate itself
              // (the 82 us curved trace still runs for <1e-4 of photons).
              double const dist_check = std::sqrt(xg_check*xg_check + yg_check*yg_check);
              constexpr double kCurvRadius_m = 2.4e7;
              constexpr double kMarginSafety = 1.5;  // L^2/(2R) is the MEAN displacement
              double const lCheck = std::abs(t_check);
              double const margin = 10.0 * telescopeRadius_m_ +
                                    kMarginSafety * lCheck * lCheck / (2.0 * kCurvRadius_m);
              if (dist_check <= margin) {
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
          ++gateRejected_;
          missCounter++;
          continue;
        }
        ++gatePassed_;

        // ---- Phase 2: provably unreachable photons, refused before any tracing ------
        // In a spherically symmetric medium whose n*r increases monotonically (checked
        // by the table's constructor; true for every real atmosphere) a photon emitted
        // with a non-negative radial direction cosine has a monotonically increasing
        // radius.  Its crossing X of the (infinite) dish plane therefore satisfies
        //     |X| >= r_emit   =>   |X - C| >= r_emit - r_obs   =>   sqrt(x^2+y^2) >= r_emit - r_obs
        // because the dish coordinates are the plane components of X - C.  So when the
        // emission is more than one aperture above the dish, no crossing can be inside
        // the dish and the photon can be refused exactly -- no trace, no hit lost.
        // Measured on the 43 km / 2944 m / 90 deg point: 65 033 of 199 315 gated photons,
        // 61 078 of which used to be traced across the plane tens of km above the dish
        // only to be discarded by the aperture.
        {
          double const rEmitCut = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
          double const rObsCut =
              std::sqrt(telescopeX_m_ * telescopeX_m_ + telescopeY_m_ * telescopeY_m_ +
                        telescopeZ_m_ * telescopeZ_m_);
          // CHER_NO_UPGOING_CUT=1 is a DIAGNOSTIC switch (same pattern as
          // CHER_AB_PLANE_BACK): it disables the cut so the saving can be measured as
          // an A/B on the same binary instead of inferred from the counter.
          if (rEmitCut > 0.0 && getenv("CHER_NO_UPGOING_CUT") == nullptr) {
            double const cosRadial = (x0 * dx + y0 * dy + z0 * dz) / rEmitCut;
            if (cosRadial >= 0.0 && rEmitCut > rObsCut + telescopeRadius_m_) {
              ++upgoingPreCut_;
              ++missCounter;
              continue;
            }
          }
        }

        // Curved tracing: the TABLE first (eligible rays), else the plane tracer.
        // traceCurvedGround now serves Telescope mode as well and dispatches its own
        // fallback to rayTraceSphericalToPlane; it also reports whether the tabulated
        // tau is available, which decides who owns the absorption draw.
        if (enableCurvature_ && refractiveIndexFunc_) {
          bool used_table_route = false;
          bool const hit = traceCurvedGround(x0, y0, z0, dx, dy, dz, lambdaNm, det_x_world,
                                            det_y_world, det_z_world, ray_optical_path,
                                            ray_geom_path, ray_accumulated_tau, used_table_route);
          if (hit) {
            used_ray_trace = true;
            // The tracer always returns tau and NEVER samples (deterministic
            // absorption), whichever route answered -- table, continuation or legacy.
            absorptionState = AbsorptionState::OpticalDepthAvailable;
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
            double const H = computeLocalScaleHeight(z_avg - earthRadiusM_);
            double const n_avg = refractiveIndexFunc_ ? refractiveIndexFunc_(z_avg - earthRadiusM_, lambdaNm).first : refractiveIndex_;
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

        // --- closed-form atmospheric delay (P1) -------------------------------
        // For a straight path the group-delay integral factorises as
        //   int (n_group-1) ds = f_group(lam) * int N(z) ds ,
        // and int N ds = (A(z_end) - A(z_start)) / u_z  (u_z = dz/ds), with A the
        // precomputed vertical column integral. This replaces the 10-100 segment
        // loop exactly. The branch threshold is on Delta z (not u_z): below it the
        // midpoint form is used (its relative error ~ (Delta z/H)^2/24 ~ 6.5e-8 at
        // Delta z = 10 m with H ~ 8 km), which also avoids the Delta-A cancellation
        // (~1e-13 relative) that would otherwise occur.
        bool delayComputed = false;
        CherenkovProfScope profTiming(enableProfiling_, &profile_.timing_s);
        if (enableFastTiming_ && path_len > 0.0) {
          constexpr double kSlantIntegralMinDeltaZ_m = 10.0;
          double const delta_z_fast = std::abs(det_z_world - z0);
          double const N0_fast = refractiveIndexTable_->N(z0);
          // f_group(lam) from the emission-point group index: n_group(z0,lam)-1 =
          // N(z0)*f_group(lam).
          double const n_group0 = refractiveIndexFunc_(z0 - earthRadiusM_, lambdaNm).second;
          double const f_group = (N0_fast > 1e-12) ? (n_group0 - 1.0) / N0_fast : 1.0;
          double delay_m;
          if (delta_z_fast >= kSlantIntegralMinDeltaZ_m) {
            double const uz = (det_z_world - z0) / path_len;
            delay_m = f_group * (refractiveIndexTable_->A(det_z_world) -
                                 refractiveIndexTable_->A(z0)) / uz;
          } else {
            double const z_mid = 0.5 * (z0 + det_z_world);
            delay_m = f_group * refractiveIndexTable_->N(z_mid) * path_len;
          }
          delay_ns = (delay_m * path_correction / c_light) * 1e9;
          delayComputed = true;
        }

        if (!delayComputed) {
        if (path_len < 1e-2) {
          // Very short path: single evaluation at emission altitude using group index
          auto const [n_phase_emit, n_group_emit] = refractiveIndexFunc_ ? 
              refractiveIndexFunc_(z0 - earthRadiusM_, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
          double const N_group = n_group_emit - 1.0;
          delay_ns = (path_len * path_correction * N_group / c_light) * 1e9;
        } else {
          double const uz_path = dz_path / path_len;
          double const delta_z = std::abs(dz_path);

          int segments = 10;
          if (timingSegmentsOverride_ > 0) {
            segments = timingSegmentsOverride_; // high-resolution validation reference
          } else if (delta_z < 100.0) {
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
                refractiveIndexFunc_(z_avg - earthRadiusM_, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
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
                  refractiveIndexFunc_(z_mid - earthRadiusM_, lambdaNm) : std::make_pair(refractiveIndex_, refractiveIndex_);
              double const N_group = n_group_mid - 1.0;
              delay_m += N_group * ds;
            }
            // Apply curvature correction to delay path as well
            delay_ns = (delay_m * path_correction / c_light) * 1e9;
          }
        }
        } // if (!delayComputed)

        totalTime_ns = (timeEmit / 1_ns) + t_geom_ns + delay_ns;
      }

      // Convert production height from Earth-centered z to altitude above sea level
      // z0 is in Earth-centered coordinates (Earth radius + altitude)
      double const z0_above_sea_level = z0 - earthRadiusM_;

      // ---- atmospheric absorption: deterministic weight attenuation -----------------
      // The bunch weight is multiplied by exp(-tau); nothing is sampled.  This is
      // unbiased: exp(-tau) is exactly the expected fraction of the bunch's photons that
      // survive, and all photons of a bunch share one path (same lambda, same trajectory)
      // so a single factor is correct for the whole bunch.  What disappears is the
      // binomial noise of the old survival draw -- a variance reduction, not a bias.
      double survival = 1.0;
      if (enableAbsorption_) {
        double const h_det_m = det_z_world - earthRadiusM_;
        double const dx_path = det_x_world - x0;
        double const dy_path = det_y_world - y0;
        double const dz_path = det_z_world - z0;
        double const path_length =
            std::sqrt(dx_path * dx_path + dy_path * dy_path + dz_path * dz_path);
        double const cos_zenith = (path_length > 1e-6) ? std::abs(dz_path) / path_length : 1.0;
        survival = absorptionTransmission(absorptionState, ray_accumulated_tau,
                                          lambdaNm, z0_above_sea_level, h_det_m, path_length,
                                          cos_zenith);
        absorptionWeights_.emittedWeight += static_cast<long double>(weightPerBunch);
        ++absorptionWeights_.evaluatedBunches;
        // WAVE 5 instrument: record the tau that reached the weight and the survival.
        if (absorptionState == AbsorptionState::OpticalDepthAvailable) {
          double const t = std::max(0.0, ray_accumulated_tau);
          absorptionWeights_.tauSum += t;
          absorptionWeights_.tauSumSq += t * t;
          if (t > absorptionWeights_.tauMax) absorptionWeights_.tauMax = t;
          ++absorptionWeights_.tauCount;
        }
        absorptionWeights_.survSum += survival;
        if (survival < 0.1) {
          ++absorptionWeights_.bSurvLt01;
        } else if (survival < 0.5) {
          ++absorptionWeights_.bSurvLt05;
        } else if (survival < 0.9) {
          ++absorptionWeights_.bSurvLt09;
        } else {
          ++absorptionWeights_.bSurvHigh;
        }
        if (survival < kAbsorptionKeepFloor) {
          // Deeply absorbed bunch: dropped to keep the output small, and the weight it
          // would have carried is ADDED UP so the truncation bias is reported, never
          // silently lost.
          absorptionWeights_.cutWeight += static_cast<long double>(weightPerBunch) *
                                          static_cast<long double>(survival);
          ++absorptionWeights_.cutBunches;
          continue;
        }
      }
      double const weightOut = weightPerBunch * survival;
      absorptionWeights_.recordedWeight += static_cast<long double>(weightOut);

      // Record the bunch once with weight (number of photons it represents)
      // For Ground mode: directions are in world coordinates
      // For ShowerAxis/Telescope modes: directions are in shower coordinates
      // Pass optical_path for timing calculations (accounts for refractive index)
      double const pathLengthForRecording = used_ray_trace ? ray_optical_path : 
                                             (path_len + (delay_ns * c_light / 1e9));
      CherenkovProfScope profRec(enableProfiling_, &profile_.recorder_s);
      recorder_(xg, yg, zg, lambdaNm, z0_above_sea_level, dx_out, dy_out, dz_out,
            dx, dy, dz,
            totalTime_ns, weightOut, n_local, pathLengthForRecording);
      if (enableProfiling_) { ++profile_.recorded; }
    }

    return corsika::ProcessReturn::Ok;
  }

} // namespace corsika::cherenkov
