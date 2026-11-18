#pragma once

#include <corsika/framework/process/ContinuousProcess.hpp>
#include <corsika/framework/geometry/FourVector.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/ParticleProperties.hpp>
#include <corsika/framework/core/Step.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/modules/cherenkov/AtmosphericAbsorption.hpp>
#include <corsika/modules/cherenkov/CurvedRayTable.hpp>
#include <corsika/modules/cherenkov/CurvedAbsorptionTable.hpp>
#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>
#include <corsika/modules/cherenkov/Dopri5.hpp>
#include <corsika/modules/cherenkov/TangentPropagator.hpp>

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <fstream>   // step-A A/B CSV (see enableTraceAB)
#include <utility>
#include <vector>

namespace corsika::cherenkov {

  // Projection mode for Cherenkov photons
  enum class ProjectionMode {
    Ground,      // Project to fixed ground z-plane (fails for horizontal showers)
    ShowerAxis,  // Project to plane perpendicular to shower axis
    Telescope    // Project to telescope surface perpendicular to shower axis with finite radius
  };

  /// Which integrator the exact tracers use.  RK4 is the default, the reference and the
  /// fallback: it is the integrator every acceptance test and every production number in
  /// diagnostics/phase1_geometry.md was measured with.  DOPRI5 is an opt-in adaptive
  /// 5(4) pair (see Dopri5.hpp); it is only eligible when the index is separable and the
  /// refractive-index table is valid, because its error control scales the optical-path
  /// component by f_group, which needs the separable decomposition.
  /// TANGENT is the Phase-4 regularised mode: the same ray, integrated in the
  /// turning-point variable (Phi, s, Sigma, v) with r = r_t + v^2.  It is opt-in and
  /// independent of DOPRI5: both are selected by this one switch.
  /// AUTO (Phase 4 step 2) is the production setting: the tangent mode is attempted
  /// first for every ray, and a REFUSAL falls back to DOPRI5 for that same ray.  It is
  /// selection per ray, not a mid-ray switch, because the regularised coordinate covers
  /// the whole path (v runs from negative, through the turning point, to positive), so a
  /// ray never needs two integrations; the safeguard that hysteresis would provide is
  /// the counted fallback instead.
  enum class TracerIntegrator { RK4, DOPRI5, TANGENT, AUTO };

  /// Per-run counters for the exact tracers (Phase 3).  Atomic: the photon loop is
  /// parallel.  `cappedSteps` counts accepted steps whose size was clipped by the
  /// configured maximum, i.e. steps where the controller wanted to go larger -- the
  /// direct measure of a restrictive cap.
  struct TracerStats {
    std::size_t accepted{};
    std::size_t rejected{};
    std::size_t fEvaluations{};
    std::size_t minStepFailures{};
    std::size_t cappedSteps{};
    std::size_t crossingIterations{};
    std::size_t dopriRays{};
    std::size_t rk4Rays{};
    std::size_t tangentRays{};        ///< Phase 4: rays served by the tangent mode
    std::size_t tangentUncertified{}; ///< ... refused (degenerate turning point, step
                                      ///< guard, above the atmosphere): never silent
    std::size_t tangentFallbacks{};   ///< AUTO: refused and re-traced by DOPRI5
    std::size_t impactRejections{};   ///< Phase 4 step 3: steps shortened for the hit budget
    std::size_t impactAccepted{};     ///< ... steps accepted under the impact criterion
    double hMin{0.0};
    double hMean{0.0};
    double hMax{0.0};
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

    /// Lightweight per-step profiling accumulators (diagnostics; off by default).
    struct ProfileAccumulators {
      double total_s{0.0};      // whole doContinuous()
      double production_s{0.0}; // photon-production integral + Poisson draw
      double bunchloop_s{0.0};  // per-bunch projection + timing + absorption + record
      double timing_s{0.0};     // atmospheric-delay computation (fast or segment loop)
      double recorder_s{0.0};   // HitRecorder call (e.g. CSV/eventio write)
      std::size_t ri_calls{0};  // refractiveIndexFunc_ calls
      std::size_t bunches{0};   // bunches generated
      std::size_t recorded{0};  // bunches handed to the recorder
    };

    // ---- step C: the plane crossing that lies BEHIND the launch state ---------
    // rayTraceSphericalToPlane() accepts a crossing at tPlane <= 0, which is legitimate
    // when the state handed to it is not the emission (launch-and-continue from the
    // table's sphere impact: the photon may have crossed the dish plane just before it).
    // From the EMISSION it is a phantom hit, because a photon cannot cross a plane that
    // is behind it.  These counters measure how often each case happens before any
    // policy is chosen.
    struct PlaneBackStats {
      std::size_t calls{};        ///< calls that reached the crossing test
      std::size_t backward{};     ///< tPlane <= 0 entered
      std::size_t nearZero{};     ///< ... and |tPlane| <= 1 m (numerical tangency)
      std::size_t rejected{};     ///< ... and farther than the tolerance: NOT a hit
      double minTPlane{};         ///< most negative tPlane seen [m]
    };
    PlaneBackStats planeBackStats() const {
      PlaneBackStats s;
      s.calls = nPlaneCalls_.load();
      s.backward = nPlaneBack_.load();
      s.nearZero = nPlaneTangent_.load();
      s.minTPlane = planeMinT_.load();
      s.rejected = nPlaneRejected_.load();
      return s;
    }
    void resetPlaneBackStats() const {
      nPlaneCalls_ = 0; nPlaneBack_ = 0; nPlaneTangent_ = 0; nPlaneRejected_ = 0;
      planeMinT_ = 0.0;
    }

    /// Deterministic absorption bookkeeping (see absorptionWeightStats()).
    /// Absorption multiplies the bunch weight by exp(-tau); bunches whose survival falls
    /// below kAbsorptionKeepFloor are dropped for file size and the weight they would
    /// have carried is accumulated here, so the truncation bias is reported.
    struct AbsorptionWeightStats {
      long double emittedWeight{0.0};   ///< weight handed to the absorption step
      long double recordedWeight{0.0};  ///< weight actually recorded
      long double cutWeight{0.0};       ///< weight lost to the keep-floor cut
      std::size_t evaluatedBunches{0};
      std::size_t cutBunches{0};
      // WAVE 5 instrument: the tau / survival distribution that ACTUALLY reaches the
      // weight.  Built in-process (a CSV write throttled the run and was lost to an OOM),
      // printed once at the end.  It answers the one open W2 question: is the run-level
      // mean survival consistent with the per-bunch tau, or is the bookkeeping wrong?
      std::size_t tauCount{0};
      double tauSum{0.0};
      double tauSumSq{0.0};
      double tauMax{0.0};
      double survSum{0.0};
      std::size_t bSurvLt01{0};
      std::size_t bSurvLt05{0};
      std::size_t bSurvLt09{0};
      std::size_t bSurvHigh{0};
      // WAVE 5 probe: do the two routes agree on xi (sin theta_arrival) and on the perigee
      // flag for the SAME photon?  When both exist (A/B sample) this settles the convention.
      std::size_t xiCmpCount{0};
      double xiDiffSum{0.0};
      double xiDiffMax{0.0};
      std::size_t xiFlagMismatch{0};
      // WAVE 5: the gated substitution's own tau against the tracer's, on the same photon.
      std::size_t substCount{0};
      double substTracerSum{0.0};
      double substTableSum{0.0};
      double substDiffMaxAbs{0.0};
      std::size_t substRejected{0};
    };

  private:
    double refractiveIndex_;
    double fLambdaMin;
    double fLambdaMax;
    double fBunchSize;
    double groundZ_m_;
    double earthRadiusM_;  // Earth's radius in meters (for height above sea level calculation)
    ProjectionMode projectionMode_;
    // Fixed shower axis (primary particle direction) used for all projections
    std::shared_ptr<corsika::media::ShowerAxis const> primaryShowerAxis_;
    // Shower core coordinates (for ShowerAxis projection mode)
    double showerCoreX_m_;
    double showerCoreY_m_;
    double showerCoreZ_m_;
    // Telescope parameters (for Telescope projection mode)
    mutable double telescopeX_m_;
    mutable double telescopeY_m_;
    mutable double telescopeZ_m_;
    double telescopeRadius_m_;
    // Original telescope center (for per-event reset and correction offsets)
    double telescopeXOriginal_m_;
    double telescopeYOriginal_m_;
    double telescopeZOriginal_m_;
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
    mutable std::size_t absorptionDraws_{0};  ///< survival draws taken (always 0 now)
    /// Weight-attenuation bookkeeping (deterministic absorption, see AbsorptionWeightStats)
    mutable AbsorptionWeightStats absorptionWeights_;
    // ---- step A A/B diagnostics (mutable: traceCurvedGround is const) ----------
    mutable std::size_t abMaxRows_{0};
    mutable std::size_t abRows_{0};
    mutable std::string abPath_;
    mutable std::mutex abMutex_;  ///< serialises the CSV writes across threads
    /// WAVE 5: rows written to <abPath>.legacy.csv, the legacy route's own tau against
    /// the tabulated one on the SAME photon.  This is the instrument that replaced a
    /// six-photon sample: it answers whether the factor-2 bulk difference is uniform
    /// (a bug) or tail-only (a path effect), and whether the legacy integral is a chord.
    mutable std::size_t abLegacyRows_{0};
    mutable std::once_flag legacyAbsOnce_;
    mutable std::shared_ptr<CurvedAbsorptionTable const> legacyAbsTable_;
    // step-C counters (atomic: the photon loop is parallel)
    mutable std::atomic<std::size_t> nPlaneCalls_{0}, nPlaneBack_{0}, nPlaneTangent_{0},
        nPlaneRejected_{0};
    mutable std::atomic<double> planeMinT_{0.0};
    struct TraceABRow {
      double hEmit{}, mu{}, lambdaNm{};
      int statusTable{-1};
      double Ltable{}, Llegacy{}, dLrel{}, dImpact{};
      double xTable{}, yTable{}, xLegacy{}, yLegacy{};
      double tauTable{}, tauLegacy{}, tauChordTable{}, tauChordLegacy{};
      double tauRefQuad{};  ///< step A' referee: exact quadrature along the refracted ray
      double tauRefQuadLo{};  ///< WAVE 5: the same quadrature at 500 samples (convergence test)
      int legacyOk{0};
    };
    void writeTraceAB_(TraceABRow const& row) const;
    /// Per-photon gate classification (2g-h2n wave 1).  Atomic because the gate runs in
    /// the same parallel loop as the tracing and is not owned by CurvedTableState
    /// (which is created lazily by traceCurvedGround, i.e. after the gate).
    mutable std::atomic<std::size_t> gatePassed_{0}, gateRejected_{0};
    // ---- Phase 2: cost accounting (same reasoning as the gate counters: the cuts fire
    // before CurvedTableState exists, so they cannot live there) -------------------
    /// Photons refused before any tracing because they are provably unreachable:
    /// (a) Telescope: moving outward from a radius more than one aperture above the
    ///     dish, so no crossing of the dish plane can be inside the aperture;
    /// (b) rayTraceSpherical: emitted above the target sphere and moving outward.
    mutable std::atomic<std::size_t> upgoingPreCut_{0};
    /// Tracer failures by kind.  Counted ALWAYS, logged only for the first few: the
    /// "escaped atmosphere" warning used to fire once per photon (40 032 lines on one
    /// production event).
    mutable std::atomic<std::size_t> escapedAtmosphere_{0};
    mutable std::atomic<std::size_t> maxStepsReached_{0};

    // ---- Phase 3: integrator selection and its counters --------------------------
    TracerIntegrator tracerIntegrator_{TracerIntegrator::RK4};
    double dopriAtolPos_m_{1.0e-3}; ///< [m]      position, per component
    double dopriAtolDir_{1.0e-9};   ///< [-]      direction cosine
    double dopriAtolOpt_m_{1.0e-4}; ///< [m]      Sigma (group path)
    double dopriMinStep_m_{1.0e-3};
    double dopriMaxStep_m_{0.0};    ///< 0 => the RK4 heuristic cap
    double impactBudget_m_{0.0};    ///< Phase 4 step 3: 0 => impact control off
    bool crossingRelanding_{false}; ///< Phase 4 step 3: land the crossing on the plane
    mutable std::atomic<std::size_t> impactUncertified_{0};
    mutable std::atomic<std::size_t> trAccepted_{0}, trRejected_{0}, trF_{0},
        trMinStepFail_{0}, trCapped_{0}, trCrossIters_{0}, trDopriRays_{0}, trRk4Rays_{0},
        trTangentRays_{0}, trTangentUncertified_{0}, trTangentFallbacks_{0},
        trImpactRej_{0}, trImpactOk_{0};
    mutable std::atomic<double> trHMin_{1.0e300}, trHMax_{0.0};
    // Extinction coefficient function: returns extinction per meter at altitude and wavelength
    std::function<double(double, double)> extinctionCoeffFunc_;

    // Optional vertical refractivity table for the closed-form straight-line timing
    // path (see RefractiveIndexTable.hpp). When present (and the refractive-index
    // function is separable), the 10-100 segment timing loop is replaced by an
    // exact closed-form column integral.
    std::shared_ptr<RefractiveIndexTable const> refractiveIndexTable_;
    bool enableFastTiming_{false};  // set in the ctor from refractiveIndexTable_
    // Separability of n(z,lam)-1 = N(z)*f(lam), measured once in the ctor.  This is
    // what lets the tracers take BOTH the phase and the group index and their radial
    // derivative from the refractive-index table (N, dN/dz) instead of three
    // refractiveIndexFunc_ calls per step.  Deliberately independent of
    // enableFastTiming_, which is a TIMING switch: --disable-cherenkov-fast-timing
    // must not silently change how the trajectory is integrated.
    bool separableIndex_{false};

    // ---- curved-ray table: lazy, synchronized, process-lifetime state --------
    // Held through a shared_ptr so that copies of CherenkovProcess share ONE table
    // and ONE set of counters (documented semantics: counters are per-run and are
    // reset only with no tracing in flight; independent concurrent runs need
    // separate CherenkovProcess instances). Initialization runs exactly once under
    // a once_flag; the counters are atomic so concurrent lookups do not race.
    struct CurvedTableState {
      std::once_flag once;
      std::shared_ptr<CurvedRayTable const> table;
      /// Wave 2c: accumulated optical depth along the refracted ray, (h, xi, lambda).
      /// It replaces the geometry table's own tau for the curved route (19x smaller and
      /// measured 3.8e-4 against an independent quadrature vs 3.4e-3 for the 5-node
      /// in-table version).  Wave 4b removed that in-table version entirely: when this
      /// pointer is null the route falls back to the exact chord optical depth and
      /// counts it in absFallbackTau (mirrored as CurvedTableStats::absFallbackTau).
      std::shared_ptr<CurvedAbsorptionTable const> absTable;
      std::atomic<std::size_t> absHits{0};
      /// of those, the rays whose path visited the perigee (both legs accumulated)
      std::atomic<std::size_t> absSumRule{0};
      std::atomic<std::size_t> absFallbacks{0};
      /// Rays where absorption is on but no absorption table exists: tau = 0 and counted
      /// (absFallbackTau counts the chord last resort; this one is the build failure).
      std::atomic<std::size_t> absMissing{0};
      bool attempted{false};
      bool available{false};
      std::string unavailableReason;
      std::atomic<std::size_t> bypassed{0};
      std::atomic<std::size_t> attemptedLookups{0};
      std::atomic<std::size_t> hits{0};
      std::atomic<std::size_t> rejected[9]{};
      std::atomic<std::size_t> legacyFallbackSucceeded{0};
      std::atomic<std::size_t> legacyFallbackFailed{0};
      std::atomic<std::size_t> bypassedLegacySucceeded{0};
      std::atomic<std::size_t> bypassedLegacyFailed{0};
      std::atomic<std::size_t> finalDrop{0};
      // ---- per-photon classification (2g-h2n wave 1) -------------------------
      // These replace the Monte-Carlo estimates of the acceptance fractions with the
      // numbers a production run actually sees: how many photons reach the gate, how
      // many of those the table answers, and split by the branch that answered.
      std::atomic<std::size_t> gatePassed{0};   // pre-check gate accepted
      std::atomic<std::size_t> gateRejected{0}; // pre-check gate dropped (straight line far off)
      std::atomic<std::size_t> tableHitAbove{0}; // table Hit, emission above r_obs (outward block)
      std::atomic<std::size_t> tableHitBelow{0}; // table Hit, emission below r_obs (inward block)
    };
    mutable std::shared_ptr<CurvedTableState> curvedState_;
    bool enableCurvedTable_{true};
    /// On-disk table cache directory ("" = disabled); see setTableCacheDir().
    std::string tableCacheDir_{};
    /// WAVE 5: grazing rejection band [m] — lookup guard, not part of the table or its cache
    /// key. 2400 m puts the cut at mu ~= 0.0347, just above the last coarse bin of the impact
    /// cliff (0.244 m at mu 0.031-0.033); above 0.033 the table is millimetre-accurate
    /// (0.011 m at 0.033-0.040, 0.006 m above 0.040). Calibration: 164 m per 0.001 of mu,
    /// from the measured cuts at 400 m (mu 0.0225) and 1500 m (mu 0.0292).
    double grazingClearanceM_{2400.0};
    int timingSegmentsOverride_{0}; // 0 = built-in heuristic; >0 = fixed count
                                    // (high-resolution validation reference)

    // Lightweight per-step profiling (diagnostics; disabled by default)
    mutable ProfileAccumulators profile_;
    bool enableProfiling_{false};

    // Telescope basis vectors (local axes in world coordinates)
    double telescopeBasisX_x_{};
    double telescopeBasisX_y_{};
    double telescopeBasisX_z_{};
    double telescopeBasisY_x_{};
    double telescopeBasisY_y_{};
    double telescopeBasisY_z_{};
    double telescopeBasisZ_x_{};
    double telescopeBasisZ_y_{};
    double telescopeBasisZ_z_{};

    // Telescope position correction configuration
    bool enableTelescopePositionCorrection_{};
    double telescopePositionCorrectionFraction_{};
    std::size_t telescopePositionCorrectionInterval_{};

    // Telescope position correction state (mutable because doContinuous is const)
    mutable std::size_t telescopePositionCorrectionStepCounter_{};
    mutable bool telescopePositionCorrectionCandidateValid_{};
    mutable double telescopePositionCorrectionCandidateEnergyGeV_{};
    mutable int telescopePositionCorrectionCandidatePid_{};
    mutable double telescopePositionCorrectionCandidatePosX_m_{};
    mutable double telescopePositionCorrectionCandidatePosY_m_{};
    mutable double telescopePositionCorrectionCandidatePosZ_m_{};
    mutable double telescopePositionCorrectionCandidateDirX_{};
    mutable double telescopePositionCorrectionCandidateDirY_{};
    mutable double telescopePositionCorrectionCandidateDirZ_{};
    mutable bool telescopePositionCorrectionLastUpdateValid_{};
    mutable int telescopePositionCorrectionLastUpdatePid_{};
    mutable std::mutex telescopePositionCorrectionMutex_;

    // Helper function to compute local atmospheric scale height from refractive index gradient
    double computeLocalScaleHeight(double z) const;

    // Helper to compute telescope local basis vectors from pointing angles
    void computeTelescopeBasis() noexcept;

    // Track and apply telescope position correction using the most energetic charged particle
    void handleTelescopePositionCorrection(double energy_GeV,
                         double pos_x_m,
                         double pos_y_m,
                         double pos_z_m,
                         double dir_x,
                         double dir_y,
                         double dir_z,
                         int pid) const;

    // Apply the current correction candidate (mutex must be held)
    bool applyTelescopePositionCorrectionCandidate() const;

  public:
    // Spherical ray tracing to tilted plane using Bouguer's formula: n(r)*r*sin(theta) = const
    // Returns true if ray reaches the plane; outputs impact point and optical/geometric path lengths
    // Exposed for the Telescope-mode acceptance set (test-only use), on the same
    // terms as rayTraceSpherical() below; the production path is unchanged and
    // still goes through doContinuous().
    bool rayTraceSphericalToPlane(double x0, double y0, double z0_altitude,
           double dx, double dy, double dz,
           double plane_x, double plane_y, double plane_z,
           double norm_x, double norm_y, double norm_z,
           double& x_impact, double& y_impact, double& z_impact,
           double& optical_path, double& geometric_path,
           double lambdaNm,
           double& accumulated_tau,
           double& dx_impact, double& dy_impact, double& dz_impact,
           double maxRadialStep_m = 100.0) const;

    /// Adaptive DOPRI5 5(4) variant of rayTraceSphericalToPlane (Phase 3).  Same
    /// contract, same args; used only when tracerIntegrator() == DOPRI5 and the index
    /// is separable with a valid table.  The RK4 version above stays the reference and
    /// the fallback.
    bool rayTracePlaneDopri(double x0, double y0, double z0_altitude,
                            double dx, double dy, double dz,
                            double plane_x, double plane_y, double plane_z,
                            double norm_x, double norm_y, double norm_z,
                            double& x_impact, double& y_impact, double& z_impact,
                            double& optical_path, double& geometric_path,
                            double lambdaNm,
                            double& accumulated_tau,
                            double& dx_impact, double& dy_impact, double& dz_impact,
                            double maxRadialStep_m = 100.0) const;

    /// Phase 4: regularised (turning-point) variant, same contract.  Integrates in
    /// y = (Phi, s, Sigma, v) with r = r_t + v^2 and the conserved invariant
    /// p = n r sin(psi); selected by tracerIntegrator() == TANGENT.
    bool rayTracePlaneTangent(double x0, double y0, double z0_altitude,
                              double dx, double dy, double dz,
                              double plane_x, double plane_y, double plane_z,
                              double norm_x, double norm_y, double norm_z,
                              double& x_impact, double& y_impact, double& z_impact,
                              double& optical_path, double& geometric_path,
                              double lambdaNm,
                              double& accumulated_tau,
                              double& dx_impact, double& dy_impact, double& dz_impact,
                              double maxRadialStep_m = 100.0) const;

    /// Curved ground tracing with the table as the default for eligible rays.
    /// Returns true on a hit and fills the same outputs as rayTraceSpherical().
    /// Eligibility (one predicate, evaluated in one place):
    ///   enableCurvature_ && enableCurvedTable_ && Ground projection
    ///   && !(enableAbsorption_ && extinctionCoeffFunc_)   [absorption -> legacy]
    ///   && refractiveIndexTable_ valid && table built/valid && coverage holds
    /// On a table rejection or an ineligible ray the UNCHANGED legacy tracer runs
    /// and the outcome is accounted for separately.
    bool traceCurvedGround(double x0, double y0, double z0, double dx, double dy, double dz,
                           double lambdaNm, double& x_ground, double& y_ground, double& z_ground,
                           double& optical_path, double& geometric_path,
                           double& accumulated_tau, bool& usedTable) const;

    // Spherical ray tracing using Bouguer's formula: n(r)*r*sin(theta) = const
    // Returns true if ray reaches ground; outputs ground x/y (Cartesian), optical and geometric path lengths.
    // Exposed for the CurvedRayTable legacy comparison (test-only use); the
    // production path is unchanged and still goes through doContinuous().
    bool rayTraceSpherical(double x0, double y0, double z0_altitude,
                 double dx, double dy, double dz,
                 double& x_ground, double& y_ground, double& z_ground,
                 double& optical_path, double& geometric_path,
                 double lambdaNm,
                 double& accumulated_tau,
                 double maxRadialStep_m = 100.0) const;

  private:

  public:
    CherenkovProcess(double refractiveIndex,
                     double lambdaMin   = 300.0,
                     double lambdaMax   = 900.0,
                     double bunchSize   = 5.0,
                     double groundZ_m   = 0.0,
                     HitRecorder recorder = HitRecorder{},
                     ProjectionMode mode = ProjectionMode::Ground,
                     std::shared_ptr<corsika::media::ShowerAxis const> primaryShowerAxis = nullptr,
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
                     const std::string& absorptionFile = "",
                     double telescopePositionCorrectionFraction = 0.0,
                     std::shared_ptr<RefractiveIndexTable const> refractiveIndexTable = nullptr,
                     bool enableCurvedTable = true);

    // ---- curved-ray table accounting ----------------------------------------
    // The five classes below satisfy, by construction:
    //   total     = bypassed + attempted
    //   attempted = hits + sum(rejected)
    //   sum(rejected) = legacyFallbackSucceeded + legacyFallbackFailed
    //   finalDrop = bypassedLegacyFailed + legacyFallbackFailed
    // A table rejection rescued by the (unchanged) legacy tracer reproduces the
    // baseline acceptance for that ray; a legacy failure does not. Whether a table
    // hit would have failed in the legacy tracer is NOT observable from these
    // counters -- that needs a separate both-routes comparison run.
    struct CurvedTableStats {
      std::size_t bypassed{};                 // ineligible by the predicate
      std::size_t attempted{};                // table lookups actually made
      std::size_t hits{};                     // table hits used
      std::size_t rejected[9]{};              // by CurvedRayTable::Status value
      std::size_t legacyFallbackSucceeded{};  // after a table rejection
      std::size_t legacyFallbackFailed{};
      std::size_t bypassedLegacySucceeded{};  // legacy used because ineligible
      std::size_t bypassedLegacyFailed{};
      std::size_t finalDrop{};
      // per-photon classification (see CurvedTableState)
      std::size_t gatePassed{};
      std::size_t gateRejected{};
      std::size_t tableHitAbove{};
      std::size_t tableHitBelow{};
      /// Phase 2: Telescope photons refused before tracing because they move outward
      /// from a radius more than one aperture above the dish (provably unrecordable).
      std::size_t upgoingPreCut{};
      /// Phase 2: tracer failures by kind, counted always and logged only a few times.
      std::size_t escapedAtmosphere{};
      std::size_t maxStepsReached{};
      /// Wave 2c/4b: tau taken from the absorption table, and the rays answered by the
      /// exact chord instead (the table rejects them inside a tangency).  Their sum is
      /// the number of rays the curved route gave a tau to; tauMissing counts the rays
      /// where absorption is on but no absorption table exists (build failure), which
      /// get tau = 0 rather than an invented value.
      std::size_t absTableTau{};
      std::size_t absFallbackTau{};
      /// rays answered with the dip/sum rule (path visited the perigee)
      std::size_t absSumRule{};
      /// absorption on but no absorption table at all (build failure): tau = 0, counted
      std::size_t tauMissing{};
      bool tableAvailable{};                  // built and valid this run
      bool tableAttempted{};                  // initialization was run
      std::string unavailableReason;
    };
    CurvedTableStats curvedTableStats() const;
    /// One INFO block with the per-photon classification of the curved route (gate,
    /// table hits split by branch, rejections by status, legacy outcomes).  Call once
    /// at the end of a run: the module has no finalize hook, so this is pull-based.
    void logCurvedTableStats() const;
    /// The last few absorption-table rejections (hA, hB, xi, lambda, status), i.e. the
    /// rays whose tau had to keep the geometry table's own value.  Diagnostics only: the
    /// buffer lives in CurvedAbsorptionTable and holds the last 32 entries.  Empty when
    /// the absorption table does not exist.
    std::vector<CurvedAbsorptionTable::Rejection> absorptionRejections() const;
    /// Reset the counters. Only safe with no tracing in flight.
    void resetCurvedTableStats() const;
    /// True when the table is built, valid and enabled.
    bool curvedTableActive() const;
    void setEnableCurvedTable(bool v) { enableCurvedTable_ = v; }
    bool enableCurvedTable() const { return enableCurvedTable_; }
    /// Directory for the on-disk table cache ("" = disabled).  The geometry and
    /// absorption tables cost 31.3 s + 17.9 s and 249 MB to BUILD and are pure
    /// deterministic functions of the profile and the configuration, so a second
    /// run with the same parameters can read them in ~0.3 s instead.  A cached file
    /// is used ONLY if its embedded key matches this run's parameters and sampled
    /// profile; otherwise it is ignored and the table is rebuilt, so a stale cache
    /// can never affect a result silently.
    void setTableCacheDir(std::string dir) { tableCacheDir_ = std::move(dir); }
    std::string const& tableCacheDir() const { return tableCacheDir_; }
    /// WAVE 5: grazing rejection band [m] applied by the geometry table's lookup guard.
    /// 0 keeps only the structural dh_in_ guard (the pre-wave-5 behaviour); larger values
    /// send the horizon cliff to the exact tracer, which now keeps the table's tau thanks
    /// to the decoupling in the rejection path.
    void setGrazingClearance(double m) { grazingClearanceM_ = (m > 0.0) ? m : 0.0; }
    double grazingClearance() const { return grazingClearanceM_; }

    // ---- 2g-h2n step A: per-photon A/B against the plane integrator -----------
    /// Diagnostic switch.  When count > 0, the first `count` Telescope-mode photons that
    /// the TABLE answers are ALSO integrated with the plane tracer from the same emission
    /// state, and one CSV row per photon is appended to `path`:
    ///
    ///   id, h_emit_m, mu, lambda_nm, status_table,
    ///   L_table_m, L_legacy_m, dL_rel, d_impact_m, x_table, y_table, x_legacy, y_legacy,
    ///   tau_table, tau_legacy, tau_chord_table, tau_chord_legacy, legacy_ok
    ///
    /// The production result is untouched: the A/B trace is an extra call, discarded
    /// afterwards.  `tau_table` is the tabulated value BEFORE the launch-and-continue
    /// segment, `tau_legacy` the plane tracer's own quadrature, and the two `tau_chord`
    /// columns are the exact chord integrals (§24) to each route's own impact point, so a
    /// disagreement can be attributed to the optical depth or to the geometry.
    void enableTraceAB(std::size_t count, std::string path) const;

    /// Step A': the REFEREE for the step-A comparison.  Exact quadrature of
    /// alpha(h, lambda) along the refracted ray with invariant xi = sin(theta_arrival),
    /// using the SAME refractive-index table the absorption table is built from, and
    /// visiting the perigee twice when the photon descends to it and comes back up.
    ///
    /// This replaces the straight-chord referee used first: a chord is not the photon's
    /// path, and at grazing the refracted ray dips to a perigee far below the emission, so
    /// the chord referee made the (correct) table look 537 % high.
    double absorptionTauQuadrature(double lambdaNm, double hEmitAboveSea,
                                   double hObsAboveSea, double xi,
                                   bool descendsAtEmit, int nSamples = 20000) const;
    std::size_t traceABRows() const { return abRows_; }


    // ---- absorption ownership contract --------------------------------------
    // Absorption is DETERMINISTIC: every route produces an optical depth and the caller
    // multiplies the bunch weight by exp(-tau).  No RNG is involved anywhere, so there
    // is nothing to double-count; the state only says WHERE the tau came from:
    //   NotEvaluated          - no absorption integral was computed; the chord-based
    //                           model applies (straight line / curvature disabled)
    //   OpticalDepthAvailable - the tracer integrated tau along the ray (curved table,
    //                           launch-and-continue, or the sphere integrator)
    enum class AbsorptionState { NotEvaluated, OpticalDepthAvailable };

    /// The single seam for the absorption step, so it is testable in isolation.
    /// Returns the survival probability exp(-tau) -- NEVER samples, never consumes
    /// randomness.  tau is the state-appropriate optical depth: the tracer's value for
    /// OpticalDepthAvailable, the exact chord integral for NotEvaluated.
    double absorptionTransmission(AbsorptionState state, double accumulatedTau, double lambdaNm,
                                  double hEmitAboveSea, double hDetAboveSea, double pathLength,
                                  double cosZenith) const;

    /// Number of absorption draws taken.  Always 0 now (the model is deterministic);
    /// kept so tests can assert that no randomness is consumed on any route.
    std::size_t absorptionDrawCount() const { return absorptionDraws_; }
    void resetAbsorptionDrawCount() const { absorptionDraws_ = 0; }

    /// ---- deterministic absorption bookkeeping (weight attenuation) --------------
    /// exp(-tau) is the exact expectation of how many of a bunch's photons survive, so
    /// multiplying the recorded weight by it keeps the yield unbiased and removes only
    /// the binomial variance of the old survival sampling.
    AbsorptionWeightStats const& absorptionWeightStats() const { return absorptionWeights_; }
    void resetAbsorptionWeightStats() const { absorptionWeights_ = AbsorptionWeightStats{}; }
    /// Survival below this is dropped instead of recorded (see AbsorptionWeightStats).
    static constexpr double kAbsorptionKeepFloor = 1.0e-4;


    // Reset telescope position correction at the start of each event
    void resetTelescopePositionCorrectionForEvent();

    struct TelescopePositionCorrectionResult {
      bool enabled{};
      bool has_update{};
      double offset_local_x_m{};
      double offset_local_y_m{};
      double offset_local_z_m{};
      int last_update_pid{};
    };

    TelescopePositionCorrectionResult getTelescopePositionCorrectionResult() const;

    /// Enable/disable the lightweight per-step profiler (off by default).
    void setProfiling(bool const on) { enableProfiling_ = on; }
    /// Accumulated profile since the last resetProfile().
    ProfileAccumulators const& getProfile() const { return profile_; }
    /// Reset the accumulated profile.
    void resetProfile() const { profile_ = ProfileAccumulators{}; }
    /// Force a fixed segment count in the (non-fast) timing path; 0 = heuristic.
    /// Useful as a high-resolution validation reference.
    void setTimingSegmentsOverride(int const n) { timingSegmentsOverride_ = n; }
    /// True when the closed-form straight-line timing path is active.
    bool fastTimingEnabled() const { return enableFastTiming_; }
    /// Force the closed-form straight-line timing path on/off (validation toggle).
    /// Disabling it falls back to the segment loop.
    void setFastTiming(bool const on) { enableFastTiming_ = on; }

    // ---- Phase 3: tracer integrator selection (see TracerIntegrator) ------------
    void setTracerIntegrator(TracerIntegrator const t) { tracerIntegrator_ = t; }
    TracerIntegrator tracerIntegrator() const { return tracerIntegrator_; }
    /// DOPRI5 tolerances: absolute position [m], absolute+relative direction
    /// (dimensionless), absolute optical-path/group-path [m].  Deliberately NOT scaled
    /// by the Earth-centred radius: 1e-9 of 6.4e6 m would be 6 mm of impact error.
    void setDopriTolerances(double const atolPos_m, double const atolDir,
                            double const atolOpt_m) {
      dopriAtolPos_m_ = (atolPos_m > 0.0) ? atolPos_m : 1.0e-3;
      dopriAtolDir_ = (atolDir > 0.0) ? atolDir : 1.0e-9;
      dopriAtolOpt_m_ = (atolOpt_m > 0.0) ? atolOpt_m : 1.0e-4;
    }
    double dopriAtolPos() const { return dopriAtolPos_m_; }
    /// Minimum step (failure guard) and maximum step [m].  A maximum <= 0 means "use the
    /// RK4 heuristic cap" (H_local/20 capped at 1000 m, plus the radial cap), which is
    /// what keeps the DOPRI5 path inside the region where it is validated.
    void setDopriStepLimits(double const minStep_m, double const maxStep_m) {
      dopriMinStep_m_ = (minStep_m > 0.0) ? minStep_m : 1.0e-3;
      dopriMaxStep_m_ = (maxStep_m > 0.0) ? maxStep_m : 0.0;
    }
    double dopriMaxStep() const { return dopriMaxStep_m_; }
    double dopriMinStep() const { return dopriMinStep_m_; }

    // ---- Phase 4, step 3: impact-aware accuracy control --------------------------
    /// Impact budget [m].  When > 0, a step is shortened until the estimated error
    /// PROJECTED ONTO THE TARGET PLANE is within that ray's share of the budget.  The
    /// projection is what the plain per-component tolerance cannot express: a position
    /// error along the plane normal moves the hit by delta/(n.u), which for a shallow ray
    /// is tens of times the error itself (phase4_plan.md, step 3).
    /// 0 (the default) disables the control and leaves the tracer bit-identical to Phase 3.
    void setImpactBudget(double const b_m) { impactBudget_m_ = (b_m > 0.0) ? b_m : 0.0; }
    double impactBudget() const { return impactBudget_m_; }
    /// Crossing correction: re-land the crossing step ON the plane so the 5th-order state,
    /// rather than the quartic interpolant's interior value, is the hit.
    void setCrossingRelanding(bool const on) { crossingRelanding_ = on; }
    bool crossingRelanding() const { return crossingRelanding_; }
    /// Rays whose impact budget could not be met even at the minimum step, i.e. traced but
    /// NOT certified.  Reported, never silent.
    std::size_t impactUncertified() const { return impactUncertified_.load(); }

    /// Aggregated tracer counters for the run (see TracerStats).
    TracerStats tracerStats() const;
    void resetTracerStats() const;

    /// Set the Cherenkov bunch size (photons per bunch). Larger bunches reduce the
    /// number of bunches (hence output size and CPU cost) at the price of coarser
    /// statistical granularity. Intended to be called per event with an
    /// energy-dependent value (CORSIKA 7 CERSIZ-like). Ignored if <= 0.
    void setBunchSize(double const b) { if (b > 0.0) { fBunchSize = b; } }
    double getBunchSize() const { return fBunchSize; }

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
