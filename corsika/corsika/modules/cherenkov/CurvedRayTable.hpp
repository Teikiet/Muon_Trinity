#pragma once

/*
 * CurvedRayTable.hpp
 *
 * Lookup table for curved photon paths in a spherically symmetric atmosphere
 * (CherenkovProcess, --enable-curvature). Header-only; no framework deps.  The ONLY
 * I/O is the opt-in on-disk build cache (saveToFile / loadFromFile), which is off
 * unless a cache directory is configured.
 *
 * Physics
 * -------
 * Bouguer's invariant p = n(r) r sin(theta) fixes the whole ray. With r_obs the
 * observation level: p_crit = n(r_obs) r_obs; a downgoing ray reaches r_obs iff
 * p <= p_crit; upgoing rays never return if n*r is monotonic; rays with p > p_crit
 * miss the observation level.
 *
 * Storage (vacuum residuals, so the interpolated quantities are smooth)
 * -------------------------------------------------------------------
 * Indexed by h (emission height above r_obs, two-level grid) x mu_obs
 * (|cos theta_obs| = sqrt(1-(p/p_crit)^2)). The mu axis is uniform in a monotone
 * coordinate s(mu) -- log in mu below kMuB = 0.9, uniform in u = sqrt(1-mu) above
 * -- so the grid is fine at BOTH ends of mu. Per node we store the residuals
 * against an ANALYTIC vacuum reference ray with the same mu_obs:
 *
 *   Dd = d - d_vac(h,mu)          DL = L - L_vac(h,mu)
 *   dDd/dh, dDL/dh                (exact, from the ODE)
 *   dDd/dmu, dDL/dmu              (differenced across mu columns)
 *   d2Dd/dhdmu, d2DL/dhdmu        (differenced from the exact d/dh)
 *   Dd, DL at f = 1-2eps, 1-eps, 1, 1+eps, 1+2eps
 *                                 (stored as DELTAS from the f = 1 node, quartic
 *                                  Lagrange in f; the refractivity scales the
 *                                  whole ray, so a low-order response in f is not
 *                                  accurate enough at the 1 cm / 0.1 ns level)
 *   Sigma = int N ds              (bilinear; justified in the unit test)
 *   dSigma/df                     (central difference of the f = 1 -/+ eps traces)
 *
 * The vacuum reference carries the sqrt(mu^2 + 2h/r_obs) cusp and the steep
 * near-ground h-slope; subtracting it analytically leaves only the smooth
 * refraction part, so bicubic Hermite (in h and s, with the cross term) is
 * accurate to ~mm. The residual still has a sqrt(1-mu^2) structure at the
 * near-vertical end (leading term sin(theta)*(L - L_vac)); the s axis makes that
 * LINEAR in s, which is what the u = sqrt(1-mu) node placement is for.
 *
 * Timing is the group optical path along the PHASE ray:
 *   t * c = L(f_phase) + f_group * Sigma(f_phase) .
 * Both f-dependences are carried: L moves by ~4.4e4 m per unit f (the DL deltas)
 * and Sigma by ~8 m per unit f (dSigma/df). Dropping the Sigma term leaves ~1 ns
 * at grazing mu. Fermat's principle does NOT let these be folded into one
 * "optical path" quantity: stationarity requires fixed endpoints, but the build
 * fixes (r_obs, mu_obs, h), so the arrival point on the r_obs + h sphere moves
 * with f. mu is computed from the photon's own f_phase in p and p_crit.
 *
 * Photons whose f_phase lies outside [f_lo_, f_hi_] = [1-2eps, 1+2eps] are
 * REJECTED (never extrapolated); the caller falls back to the exact tracer.
 *
 * Build frame: the ray is traced OUTWARD from r_obs with theta_start = arccos(mu)
 * in [0,pi/2] and ds > 0, so cos(theta) > 0 and tan(theta) > 0 throughout, hence
 * dDd/dh > 0 and dDL/dh > 0.
 *
 * Inward block (2g-h2n): photons emitted BELOW the observation level
 * ------------------------------------------------------------------
 * Earth-skimming nu_tau photons are all emitted below the telescope and RISE to
 * it.  That family is still 2-D: an emission at h < 0 that arrives at r_obs with
 * radial direction cosine +mu is the time reverse of a ray started at r_obs going
 * inward at cos(theta) = -mu, so azimuth drops out again.  The block is therefore
 * the same n_h x n_mu array extended to h < 0 (rows [0, i0), ascending, uniform on
 * dh_lo_), down to the ground r = R.  Ground mode (r_obs == R) has no room below
 * and gets an EMPTY block, so that path keeps its memory and behaviour exactly.
 *
 * The reference has to change: the vacuum chord with the arrival impact parameter
 * p = n(r_obs) r_obs sin(theta_obs) only reaches down to r = p, but a ray that
 * skims the Earth is emitted BELOW its own vacuum perigee (refraction holds it up;
 * measured on the frozen acceptance set: 500 m / 200 km arrives with mu = 0.021,
 * i.e. 1.2 deg above the local horizontal, and p is 1.4 km above the ground while
 * the emission is 2.5 km below it).  There is no vacuum chord for those rays at
 * all, which is why plain bilinear interpolation of the residual gave 33 % (and
 * why "just delete the h < 0 rejection" is not an option).
 *
 * Instead the inward branch uses the ray in a UNIFORM medium of the emission's
 * index n_e = n(r_obs + h):
 *
 *   b = p / n_e ,  L_ref = sqrt(r_obs^2 - b^2) - sqrt(r_e^2 - b^2)
 *   D_ref = r_obs * (asin(b/r_e) - asin(b/r_obs))
 *   dL_ref/dh = -r_e/sqrt(r_e^2-b^2) + b (db/dh)(1/sqrt(r_e^2-b^2) - 1/sqrt(r_obs^2-b^2))
 *   dPhi_ref/dh = (db/dh - b/r_e)/sqrt(r_e^2-b^2) - (db/dh)/sqrt(r_obs^2-b^2)
 *
 * with db/dh = -p n'(r_e) f / n_e^2.  By construction sin(psi_emit) = b/r_e equals
 * the PHYSICAL emission angle p/(n_e r_e), so this reference leaves the emission
 * point along the exact photon direction and exists everywhere the photon does
 * (b <= r_e always, with equality exactly at the apparent horizon).  Only the
 * variation of n along the path is left in the residual, so the residual is smooth
 * and carries the same sqrt(mu - mu_horizon) cusp as the physical ray (whose
 * coefficient differs only by ~delta n/n).  The h = 0 row is shared with the
 * outward block (both give d = L = 0 there); its stored h-derivative is the
 * OUTWARD one, so the single h cell that touches it ([ -dh_lo, 0 ]) is rejected --
 * the table never needs the emission to be within one node of r_obs.
 *
 * Lookups are rejected (caller falls back to the exact tracer) inside the apparent
 * horizon, where the reference degenerates: the guard is r_e - b >= dh_in/2, which
 * is the physical reachability condition (a photon cannot leave its emission point
 * tangentially) plus one half node of margin for the sqrt cusp.
 */

#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <functional>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

namespace corsika::cherenkov {

  // NOTE (gap 4: the DESCENDING-emission branch).  A photon emitted below r_obs whose
  // direction points DOWN dips through its perigee and rises back; it shares (h, mu)
  // with the rising branch, so it would need a second residual set.  It is NOT served:
  // in its only physical window (mu <~ 0.031, because the perigee must clear the
  // ground) the uniform-index chord is a 10-20 % reference whose sqrt cusp sits inside
  // the interpolated path (measured |dL|/L = 1.7e-3 .. 3.7e-2, |dP| up to 2.5 km,
  // worse than the 82 us fallback).  It needs a constant-bending (arc) reference --
  // see diagnostics/b_table/README.md, section "Gap 4".
  //
  // Those photons keep falling back through Status::Upgoing, which is correct.

  class CurvedRayTable {
  public:
    /// Why a lookup produced (or failed to produce) a result. Callers must switch on
    /// this rather than inferring a rejection from the shared counters, which are
    /// process-wide and racy under concurrent use.
    enum class Status {
      Hit,                  ///< traced to the observation level
      InvalidTable,         ///< table not built / not valid
      InvalidInput,         ///< emission radius not positive
      EmissionBelowObs,     ///< emission below the tabulated inward block (i.e. under ground)
      AboveTable,           ///< emission above the tabulated altitude range
      Upgoing,              ///< moving outward: never returns to the observation level
      BeyondCritical,       ///< p >= p_crit: misses the observation level
      DispersionOutOfRange, ///< f_phase outside [fLo(), fHi()]: no extrapolation
      Grazing,              ///< mu < muMin(), or (inward block) inside the apparent horizon
    };
    static char const* statusName(Status s) {
      switch (s) {
        case Status::Hit: return "hit";
        case Status::InvalidTable: return "invalid-table";
        case Status::InvalidInput: return "invalid-input";
        case Status::EmissionBelowObs: return "emission-below-obs";
        case Status::AboveTable: return "above-table";
        case Status::Upgoing: return "upgoing";
        case Status::BeyondCritical: return "beyond-critical";
        case Status::DispersionOutOfRange: return "dispersion-out-of-range";
        case Status::Grazing: return "grazing";
      }
      return "unknown";
    }

    struct Result {
      bool hit{false};
      Status status{Status::InvalidTable};
      double x{0.0}, y{0.0}, z{0.0};    ///< impact position (absolute, root CS) [m]
      double dx{0.0}, dy{0.0}, dz{0.0}; ///< arrival direction (unit)
      double time_s{0.0};               ///< propagation time emission -> impact [s]
      double theta_obs{0.0};            ///< arrival zenith angle [rad]
      /// Geometric path length L(f_phase) along the traced ray [m].
      double geometric_path_m{0.0};
      /// GROUP optical path length int(n_group) ds = L(f_phase) + f_group*Sigma(f_phase) [m].
      /// This is the quantity the propagation time is derived from; it is NOT the
      /// phase optical path.
      double optical_path_m{0.0};
      /// Sigma(f_phase) = int N ds along the ray [m].
      double sigma_m{0.0};
    };

    CurvedRayTable() = default;

    CurvedRayTable(RefractiveIndexTable const& ri, double R, double rObs,
                   double hTop = 120000.0, int nMu = 800, double dhLo = 10.0,
                   double dhHi = 100.0, double zBreak = 30000.0, double fEps = 0.02);

    bool valid() const { return valid_; }
    std::size_t hNodes() const { return z_nodes_.size(); }
    std::size_t muNodes() const { return n_mu_; }
    double memoryMB() const;

    // ---- on-disk build cache (wave 5) --------------------------------------
    // The build is 31 s and 229 MB of pure deterministic product of (refractive-
    // index profile, R, r_obs, exactly geometry).  Nothing about a run's physics
    // enters it, so it can be built once and reused across runs and across
    // configurations that share those parameters.
    //
    // Safety: the file carries `key`, a fingerprint of every build input including
    // 32 sampled N(z)/dNdz(z) pairs, plus a checksum of the payload.  A file whose
    // key does not match THIS configuration - or that is truncated, corrupt, or
    // written by an older layout - is REFUSED, and the caller rebuilds.  Silently
    // reusing a stale table is therefore impossible.

    /// Fingerprint of everything the build depends on (parameters + profile).
    static std::uint64_t buildKey(RefractiveIndexTable const& ri, double R, double rObs,
                                  double hTop, int nMu, double dhLo, double dhHi, double zBreak,
                                  double fEps);
    /// Write the table atomically (tmp file + rename).  False on any I/O failure;
    /// the in-memory table is never modified.
    bool saveToFile(std::string const& path, std::uint64_t key) const;
    /// Load a table written by saveToFile() and bound to `ri`.  On any mismatch the
    /// object is left INVALID and `why` (if given) says which check failed.
    bool loadFromFile(std::string const& path, RefractiveIndexTable const& ri, std::uint64_t key,
                      std::string* why = nullptr);

    /// Rows of the inward block (emission below the observation level): 0 when
    /// r_obs == R (Ground mode), ceil((r_obs-R)/dhLo) otherwise.  Stored at rows
    /// [0, zeroRow()) in ascending h, i.e. row 0 is the ground.
    std::size_t inwardNodes() const { return n_in_; }
    /// Row index of h = 0 (the first row of the outward block).
    std::size_t zeroRow() const { return i0_; }
    /// Uniform spacing of the inward block [m]; 0 when the block is empty.
    double inwardSpacing() const { return dh_in_; }

    /// Trace a photon emitted at (x0,y0,z0) with unit direction (dx,dy,dz), using
    /// the photon's own dispersion factors (fPhase for geometry, fGroup for time).
    /// No wavelength: this table carries geometry only -- the optical depth comes from
    /// CurvedAbsorptionTable (or, for the rays that table refuses, the exact chord).
    Result trace(double x0, double y0, double z0, double dx, double dy, double dz,
                 double fPhase, double fGroup) const;

    std::size_t rejectedUpgoing() const { return n_upgoing_; }
    std::size_t rejectedBeyondCrit() const { return n_beyond_; }
    /// photons whose f_phase lies outside the bracketed dispersion range
    std::size_t rejectedDispersion() const { return n_disp_; }
    /// photons more grazing than the tabulated mu range (mu < muMin())
    std::size_t rejectedGrazing() const { return n_murej_; }
    void resetCounters() const { n_upgoing_ = 0; n_beyond_ = 0; n_disp_ = 0; n_murej_ = 0; }
    /// The table is validated for mu >= muMin(); below that the grazing rays are
    /// not resolved (the near-horizontal integrand is too stiff) and the caller
    /// must fall back to the exact tracer.
    double muMin() const { return mu_min_; }
    /// WAVE 5: set the grazing rejection band [m]; 0 keeps only the structural dh_in_ guard.
    void setGrazingClearance(double m) { grazingClearance_m_ = (m > 0.0) ? m : 0.0; }
    double grazingClearance() const { return grazingClearance_m_; }
    /// solved break value and the number of nodes above it (diagnostics)
    double muBreak() const { return mu_b_; }

    /// Index of the h cell for altitude h_above_obs [m] (two-level grid, O(1)).
    std::size_t hCell(double h) const;

    // ---- diagnostics for the unit tests ------------------------------------
    double muNode(std::size_t j) const { return mu_[j]; }
    double hNode(std::size_t i) const { return z_nodes_[i]; }
    double storedDd(std::size_t i, std::size_t j) const { return Dd_[i * n_mu_ + j]; }
    double storedDdh(std::size_t i, std::size_t j) const { return Ddh_[i * n_mu_ + j]; }
    double storedDdmu(std::size_t i, std::size_t j) const { return Ddmu_[i * n_mu_ + j]; }
    double storedDL(std::size_t i, std::size_t j) const { return DL_[i * n_mu_ + j]; }
    double storedDLh(std::size_t i, std::size_t j) const { return DLh_[i * n_mu_ + j]; }
    double storedSigma(std::size_t i, std::size_t j) const { return S_[i * n_mu_ + j]; }
    double storedDSigma(std::size_t i, std::size_t j) const { return DSig_[i * n_mu_ + j]; }
    /// dispersion-node path-length deltas (t = -2,-1,+1,+2)
    double storedDeltaL(std::size_t i, std::size_t j, int t) const {
      std::size_t const k = i * n_mu_ + j;
      switch (t) {
        case -2: return DLm2_[k];
        case -1: return DLm1_[k];
        case 1: return DLp1_[k];
        default: return DLp2_[k];
      }
    }
    /// dispersion-node geometry deltas (t = -2,-1,+1,+2)
    double storedDelta(std::size_t i, std::size_t j, int t) const {
      std::size_t const k = i * n_mu_ + j;
      switch (t) {
        case -2: return Dfm2_[k];
        case -1: return Dfm1_[k];
        case 1: return Dfp1_[k];
        default: return Dfp2_[k];
      }
    }
    /// grid-coordinate diagnostics: s = muAxis(mu), the node s values, and the
    /// (cell index, local position) the lookup would use. Test-only hooks.
    double muAxisOf(double mu) const { return muAxis(mu); }
    double nodeSOf(std::size_t j) const { return nodeS(j); }
    double deltaS() const { return s_tot_ / double(n_mu_ - 1); }
    double sTot() const { return s_tot_; }
    void muCell(double mu, std::size_t& j, double& t) const {
      double const ss = muAxis(mu);
      std::size_t jj = static_cast<std::size_t>(ss / deltaS());
      if (jj + 1 >= n_mu_) jj = n_mu_ - 2;
      while (jj > 0 && nodeS(jj) > ss) --jj;
      while (jj + 2 < n_mu_ && nodeS(jj + 1) <= ss) ++jj;
      j = jj;
      t = (ss - nodeS(jj)) / (nodeS(jj + 1) - nodeS(jj));
    }
    double fLo() const { return f_lo_; }
    double fHi() const { return f_hi_; }
    /// diagnostics: which mu column failed the build and why (0 = ok)
    int buildFailCol() const { return buildFailCol_; }
    int buildFailReason() const { return buildFailReason_; }
    /// bilinear Sigma lookup (diagnostics only; mirrors the trace() path).  NOTE: this
    /// hook covers the OUTWARD block only (h is clamped to [0, hTop]); the inward rows
    /// are reachable through trace() alone.
    double sigmaBilinear(double h, double mu) const;
    /// analytic vacuum reference (same convention as the build)
    double vacuumD(double h, double mu) const;
    double vacuumL(double h, double mu) const;

  private:
    void buildColumn(double mu, double f, std::vector<double>& dOut,
                     std::vector<double>& lOut, std::vector<double>& sOut,
                     std::vector<double>& thOut) const;

    /// Inward pass: fills rows [0, i0_) with the state of the ray started at r_obs
    /// heading INWARD at cos(theta) = -mu (the time reverse of a rising photon
    /// arriving with +mu).  Rows beyond the turning point are left at the analytic
    /// reference (residual 0) and are never read: trace() rejects the h cell that
    /// touches them (see the "Inward block" section in the file header).
    void buildColumnInward(double mu, double f, std::vector<double>& dOut,
                           std::vector<double>& lOut, std::vector<double>& sOut,
                           std::vector<double>& thOut) const;

    /// Analytic reference for the inward branch: the ray in a uniform medium of the
    /// emission's index n(r_obs + h).  Exists wherever the photon does.
    struct InwardRef {
      double L{0.0};  ///< path length emission -> r_obs [m], positive
      double D{0.0};  ///< arc distance r_obs*(psi_emit - psi_obs) [m], positive
      double dL{0.0}; ///< dL/dh at fixed (mu, f)
      double dD{0.0}; ///< dD/dh at fixed (mu, f)
    };
    InwardRef inwardRef(double h, double mu, double f) const;

    /// Substep count per 10 m cell for the column integration (scales with 1/mu so
    /// the near-horizontal columns stay 4th order).
    int nSubBaseFor(double mu, double f) const;

    double N(double zAbs) const { return ri_ ? ri_->N(zAbs) : 0.0; }
    double dNdz(double zAbs) const { return ri_ ? ri_->dNdz(zAbs) : 0.0; }

    RefractiveIndexTable const* ri_ = nullptr;
    double R_ = 0.0, r_obs_ = 0.0;
    double z_break_ = 0.0, dh_lo_ = 0.0, dh_hi_ = 0.0;
    std::size_t n_lo_ = 0;
    std::vector<double> z_nodes_;
    std::size_t i0_ = 0;      ///< row index of h = 0 (start of the outward block)
    std::size_t n_in_ = 0;    ///< rows of the inward block (h < 0), at [0, i0_)
    double dh_in_ = 0.0;      ///< uniform inward spacing [m] (0 when the block is empty)
    std::size_t n_mu_ = 0;
    std::vector<double> mu_;
    // Monotone mu-axis coordinate s(mu) in which the node grid is UNIFORM:
    //   mu <= kMuB :  s = ln(mu/mu_min_)                (fine at small mu)
    //   mu >  kMuB :  s = s_b + K*(u_b - sqrt(1-mu))    (fine at mu -> 1)
    // The residual's sqrt(1-mu^2) cusp near mu = 1 is LINEAR in s there, so the
    // bicubic Hermite reproduces it.
    //
    // WAVE 5: a sqrt(mu - mu_min_) axis was TRIED AND REFUTED here.  The impact error
    // is a cliff at the apparent horizon (median 66.6 m below mu 0.02, 5.5 m in
    // 0.02-0.03, 12 mm above 0.03), which looked like a cusp in mu.  Tabulating in
    // sqrt(mu - mu_min) with mu_min = 0.018 -- 1.08 % -> 0.5-0.9 % spacing, twice the
    // useful resolution -- changed NOTHING above mu = 0.03 (0.0118 m to 4 digits) and
    // was 4-10 % WORSE below it.  So on THIS axis the error is not interpolation.
    // The likely reason it could not work is that the cusp sits at the PER-COLUMN
    // horizon mu_h(h) = sqrt(2*(r_obs - r_emit)/r_obs), which runs from ~0.018 to
    // ~0.028 over the emission heights of this event: a single global axis cannot
    // linearise a singularity whose position moves with h.  Testing the idea properly
    // needs a per-column axis, which the shared-node grid does not support -- see the
    // wave-5 note in diagnostics/b_table/README.md before trying it again.
    // ALL stored mu-derivatives are d/ds, and the one-sided difference at the last
    // node is exact because the function is linear in s at that end.
    // Nominal break; the constructor nudges it slightly so the break lands on a
    // node (the coordinate is C1 but not C2 there).
    static constexpr double kMuBNominal = 0.9;
    double mu_b_ = kMuBNominal;
    double u_b_ = 0.0, s_b_ = 0.0, K_ = 0.0, s_tot_ = 0.0;
    double muAxis(double mu) const {
      if (mu <= mu_b_) return std::log(mu / mu_min_);
      return s_b_ + K_ * (u_b_ - std::sqrt(std::max(0.0, 1.0 - mu)));
    }
    double axisMu(double s) const {
      if (s <= s_b_) return mu_min_ * std::exp(s);
      double const u = u_b_ - (s - s_b_) / K_;
      return 1.0 - u * u;
    }
    double nodeS(std::size_t j) const { return s_tot_ * double(j) / double(n_mu_ - 1); }
    // node values (size n_h*n_mu_, row-major i*n_mu_+j)
    // geometry residual d - d_vac and its bicubic (h,mu) derivatives
    std::vector<float> Dd_, Ddh_, Ddmu_, Ddhm_;
    // path-length residual L - L_vac and its bicubic (h,mu) derivatives
    std::vector<float> DL_, DLh_, DLmu_, DLhm_;
    // DELTAS at the four extra dispersion nodes f = 1-2eps, 1-eps, 1+eps, 1+2eps
    // (base node is f = 1) for the geometry (Df*) and the path length (DL*)
    std::vector<float> Dfm2_, Dfm1_, Dfp1_, Dfp2_;
    std::vector<float> DLm2_, DLm1_, DLp1_, DLp2_;
    // Sigma = int N ds along the f=1 ray, plus dSigma/df (central difference of
    // the f = 1 -/+ eps traces). The group optical path is L(f_p) + f_g*Sigma(f_p),
    // so BOTH the f-dependence of L (via the DL* deltas) and of Sigma (via DSig_)
    // have to be carried: Sigma alone moves ~0.4% over +-2%, which is ~1 ns.
    std::vector<float> S_, DSig_;
    double f_lo_ = 0.0, f_hi_ = 0.0; ///< valid f_phase range (no extrapolation)
    double f_eps_ = 0.02;
    bool valid_ = false;
    mutable std::size_t n_upgoing_ = 0, n_beyond_ = 0, n_disp_ = 0, n_murej_ = 0;
    // Lowest mu_obs the table is validated for, and the start of the s axis.
    // Measured (log-spaced mu scan, worst h and lambda): with the s grid the d and
    // t errors stay inside 1 cm / 0.1 ns for mu >= 2e-4; below it they break down
    // (1.2e-2 m at mu = 1e-4, 3.9e-2 m at 2.5e-5, 2.5 m at 1e-6, worst at small h),
    // because the near-horizontal integrand is too stiff for the resolution.
    // Photons below mu_min fall back to the exact tracer. The physical 1 PeV
    // Cherenkov range is mu ~ 0.002 - 0.06 (near-horizontal showers at
    // mu ~ 0.012 - 0.058), i.e. 10-300x above this cutoff.
    double mu_min_ = 2.0e-4;
    /// WAVE 5: grazing rejection band, as a clearance from the apparent horizon [m].  The
    /// impact error is a cliff there (66.6 m below mu 0.02) because d diverges as the ray
    /// approaches tangency, so rays within this clearance are refused and the exact tracer
    /// supplies the position.  Runtime-settable so the band can be swept without rebuilding,
    /// and it does NOT enter the table contents or the cache key (it is a lookup guard).
    /// 2400 m is the measured default. The cliff runs to mu 0.033 (1.05 m at 0.029-0.030,
    /// 0.66 m at 0.030-0.031, 0.244 m at 0.031-0.033) and stops dead above it (0.011 m at
    /// 0.033-0.040, 0.006 m above 0.040). The band is calibrated at 164 m of clearance per
    /// 0.001 of mu from the measured cuts at 400 m (mu 0.0225) and 1500 m (mu 0.0292), so
    /// 2400 m cuts at mu ~= 0.0347 and every remaining table answer is millimetre-accurate.
    double grazingClearance_m_ = 2400.0;
    int buildFailCol_ = -1, buildFailReason_ = 0;
  };


  inline double CurvedRayTable::vacuumL(double h, double mu) const {
    double const a = r_obs_ * mu;
    return -a + std::sqrt(a * a + 2.0 * r_obs_ * h + h * h);
  }

  inline double CurvedRayTable::vacuumD(double h, double mu) const {
    double const Lv = vacuumL(h, mu);
    double const st0 = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    double const phi = std::atan2(Lv * st0, r_obs_ + Lv * mu);
    return r_obs_ * phi;
  }

  // The two functions above are the OUTWARD (h >= 0) reference.  The inward block
  // uses inwardRef() below instead: a vacuum chord with the arrival impact
  // parameter does not reach an emission that lies below its own perigee.

  inline CurvedRayTable::InwardRef CurvedRayTable::inwardRef(double h, double mu, double f) const {
    InwardRef r;
    double const rE = r_obs_ + h;
    if (!(rE > 0.0)) return r;
    double const nObs = 1.0 + N(r_obs_) * f;
    double const nE = 1.0 + N(rE) * f;
    double const st0 = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    double const p = nObs * r_obs_ * st0; // Bouguer invariant of the arrival
    double const nEpos = std::max(nE, 1.0e-12);
    // Geometric impact parameter of the uniform-index reference.  b <= r_emit by
    // construction (b/r_emit is the PHYSICAL emission sine), equality exactly at
    // the apparent horizon; the clamp only keeps the formula real round-off away
    // from it.
    double const b = std::min(p / nEpos, std::min(r_obs_, rE) * (1.0 - 1.0e-12));
    double const A = std::sqrt(std::max(0.0, r_obs_ * r_obs_ - b * b));
    double const B = std::sqrt(std::max(0.0, rE * rE - b * b));
    r.L = A - B;
    // Phi = asin(b/r_emit) - asin(b/r_obs): the arc between the emission and the
    // impact point.  Evaluated as one atan2 of the numerators of sin(Phi) and
    // cos(Phi) -- sin(Phi) = b(A-B)/(r_obs rE), cos(Phi) = (AB+b^2)/(r_obs rE) --
    // which avoids both the two asin calls (~100 ns of a 500 ns lookup) and their
    // cancellation for small Phi.
    r.D = r_obs_ * std::atan2(b * (A - B), A * B + b * b);
    // h-derivatives at fixed (mu, f).  b depends on h through n(r_obs + h).
    double const dbdh = -p * dNdz(rE) * f / (nEpos * nEpos);
    double const floor = 1.0e-9 * r_obs_; // keeps 1/B finite at the horizon
    double const Bs = std::max(B, floor), As = std::max(A, floor);
    r.dL = -rE / Bs + b * dbdh * (1.0 / Bs - 1.0 / As);
    r.dD = r_obs_ * ((dbdh - b / rE) / Bs - dbdh / As);
    return r;
  }

  inline int CurvedRayTable::nSubBaseFor(double mu, double f) const {
    int nSubBase = 4;
    if (mu < 1.0e-3) {
      // Scale with the f-dependent stiffness as well as mu: the bending rate is
      // a(f) = 1/r + f*N'/(1+f*N), and the near-ground integrand stiffness ~ 1/a.
      // The effect is only ~2%, so this stays cheap.
      double const r0 = r_obs_;
      double const Nv = N(r0);
      double const dN = ri_->logSlope(ri_->cellIndex(r0)) * Nv;
      double const a1 = 1.0 / r0 + dN / (1.0 + Nv);
      double const af = 1.0 / r0 + (f * dN) / (1.0 + f * Nv);
      double const stiff = (af > 0.0 && a1 > 0.0) ? std::max(1.0, a1 / af) : 1.0;
      double const k = std::min(256.0, std::ceil(1.0e-3 / mu) * stiff);
      nSubBase = static_cast<int>(4.0 * k);
    }
    return nSubBase;
  }

  inline CurvedRayTable::CurvedRayTable(RefractiveIndexTable const& ri, double R, double rObs,
                                        double hTop, int nMu, double dhLo, double dhHi,
                                        double zBreak, double fEps)
      : ri_(&ri), R_(R), r_obs_(rObs), z_break_(zBreak), dh_lo_(dhLo), dh_hi_(dhHi),
        f_eps_(fEps) {
    if (!ri.valid()) return;
    if (!(r_obs_ >= R_) || !(hTop > 0.0) || nMu < 8) return;
    if (!(dh_lo_ > 0.0) || !(dh_hi_ >= dh_lo_)) return;
    if (!(z_break_ > 0.0) || !(f_eps_ > 0.0 && f_eps_ < 0.5)) return;

    // n*r must be monotonic increasing over the table's altitude range: only then
    // do upgoing rays never return and is the p > p_crit rejection safe.
    for (int m = 0; m <= 1000; ++m) {
      double const hh = hTop * double(m) / 1000.0;
      double const rr = r_obs_ + hh;
      if (!(1.0 + N(rr) + rr * dNdz(rr) > 0.0)) return;
    }

    // ---- h grid: inward block (h < 0), then two levels above -----------------
    // The inward block serves photons emitted below the observation level (see the
    // file header).  It is uniform on dh_lo_ down to the ground r = R_ and EMPTY
    // when r_obs == R_ (Ground mode: nothing below the ground), so Ground mode keeps
    // its node count, memory and numbers bit-for-bit.
    z_nodes_.clear();
    n_in_ = 0;
    dh_in_ = 0.0;
    double const depth = r_obs_ - R_;
    if (depth > 2.0 * dh_lo_) {
      // Adjust the spacing so the ground lands exactly on a node: the hCell index
      // arithmetic assumes a uniform block.
      int const nIn = static_cast<int>(std::ceil(depth / dh_lo_));
      dh_in_ = depth / double(nIn);
      n_in_ = static_cast<std::size_t>(nIn);
      z_nodes_.reserve(n_in_ + 2);
      for (int k = nIn; k >= 1; --k) z_nodes_.push_back(-dh_in_ * double(k));
    }
    i0_ = z_nodes_.size();
    z_nodes_.push_back(0.0);
    for (double h = 0.0;;) {
      double const step = (h < z_break_) ? dh_lo_ : dh_hi_;
      double nh = h + step;
      if (nh > hTop) nh = hTop;
      z_nodes_.push_back(nh);
      if (nh >= hTop - 1e-9) break;
      h = nh;
    }
    // n_lo_ counts the nodes of the OUTWARD block up to z_break_ (the hCell index
    // arithmetic adds i0_ to it).
    n_lo_ = 0;
    for (std::size_t i = i0_; i < z_nodes_.size(); ++i)
      if (z_nodes_[i] <= z_break_ + 1e-9) ++n_lo_;
    if (n_lo_ < 2) return;

    // ---- mu grid: LOG-spaced (uniform relative resolution) ----------------
    // The residual carries a ~1/mu^k structure at small h (there d_vac ~ h/mu),
    // so a uniform RELATIVE spacing is what keeps the bicubic interpolation error
    // bounded; the previous (j/N)^2 grading left 1e-2 m errors at small h.
    n_mu_ = static_cast<std::size_t>(nMu);
    mu_.resize(n_mu_);
    // K is FIXED by the C1 condition 1/mu_b == K/(2*u_b), so ds/dmu is continuous
    // at the break. mu_b is then nudged so that s_b/Delta_s is an integer, i.e. the
    // break lands exactly on a node (d2mu/ds2 still jumps there, so a cell that
    // straddled it would drop to O(Delta_s^2)). No free node-count parameters.
    auto setup = [&](double mb) {
      mu_b_ = mb;
      u_b_ = std::sqrt(std::max(1e-12, 1.0 - mb));
      s_b_ = std::log(mb / mu_min_);
      K_ = 2.0 * u_b_ / mb;
      s_tot_ = s_b_ + K_ * u_b_;
    };
    double bestMb = kMuBNominal, bestDev = 1e30;
    for (int k = 0; k <= 40000; ++k) {
      double const mb = 0.88 + 1.0e-6 * double(k);
      setup(mb);
      double const ratio = s_b_ / (s_tot_ / double(n_mu_ - 1));
      double const dev = std::abs(ratio - std::round(ratio));
      if (dev < bestDev) { bestDev = dev; bestMb = mb; }
    }
    setup(bestMb);
    for (std::size_t j = 0; j < n_mu_; ++j) mu_[j] = axisMu(nodeS(j));

    std::size_t const n_h = z_nodes_.size();
    std::size_t const n = n_h * n_mu_;
    Dd_.assign(n, 0.f);
    Ddh_.assign(n, 0.f);
    Ddmu_.assign(n, 0.f);
    Ddhm_.assign(n, 0.f);
    DL_.assign(n, 0.f);
    DLh_.assign(n, 0.f);
    DLmu_.assign(n, 0.f);
    DLhm_.assign(n, 0.f);
    Dfm2_.assign(n, 0.f);
    Dfm1_.assign(n, 0.f);
    Dfp1_.assign(n, 0.f);
    Dfp2_.assign(n, 0.f);
    DLm2_.assign(n, 0.f);
    DLm1_.assign(n, 0.f);
    DLp1_.assign(n, 0.f);
    DLp2_.assign(n, 0.f);
    S_.assign(n, 0.f);
    DSig_.assign(n, 0.f);
    // valid f_phase range: the dispersion nodes bracket it exactly, so the lookup
    // never extrapolates (out-of-range photons are rejected and fall back).
    f_lo_ = 1.0 - 2.0 * f_eps_;
    f_hi_ = 1.0 + 2.0 * f_eps_;

    std::vector<double> d0(n_h), l0(n_h), s0(n_h), t0(n_h);
    std::vector<double> d2(n_h), l2(n_h), s2(n_h), t2(n_h);
    std::vector<double> dA(n_h), lA(n_h), sA(n_h), tA(n_h);
    std::vector<double> dP(n_h), lP(n_h), sP(n_h), tP(n_h);
    std::vector<double> dQ(n_h), lQ(n_h), sQ(n_h), tQ(n_h);

    for (std::size_t j = 0; j < n_mu_; ++j) {
      double const mu = mu_[j];
      // Geometry at the base node f = 1 and at the four dispersion nodes.  The inward
      // pass fills rows [0, i0_) and the outward pass rows [i0_, n_h), so the two
      // share the scratch arrays; the inward pass must run first (it seeds every row
      // with the analytic reference, i.e. residual 0, and the outward pass then
      // overwrites its own rows).
      if (n_in_ > 0) buildColumnInward(mu, 1.0, d0, l0, s0, t0);
      buildColumn(mu, 1.0, d0, l0, s0, t0);
      if (n_in_ > 0) buildColumnInward(mu, 1.0 - 2.0 * f_eps_, d2, l2, s2, t2);
      buildColumn(mu, 1.0 - 2.0 * f_eps_, d2, l2, s2, t2);
      if (n_in_ > 0) buildColumnInward(mu, 1.0 - f_eps_, dA, lA, sA, tA);
      buildColumn(mu, 1.0 - f_eps_, dA, lA, sA, tA);
      if (n_in_ > 0) buildColumnInward(mu, 1.0 + f_eps_, dP, lP, sP, tP);
      buildColumn(mu, 1.0 + f_eps_, dP, lP, sP, tP);
      if (n_in_ > 0) buildColumnInward(mu, 1.0 + 2.0 * f_eps_, dQ, lQ, sQ, tQ);
      buildColumn(mu, 1.0 + 2.0 * f_eps_, dQ, lQ, sQ, tQ);
      double const fs[5] = {1.0, 1.0 - 2.0 * f_eps_, 1.0 - f_eps_, 1.0 + f_eps_,
                            1.0 + 2.0 * f_eps_};
      for (std::size_t i = 0; i < n_h; ++i) {
        double const h = z_nodes_[i];
        double const r = r_obs_ + h;
        std::size_t const k = i * n_mu_ + j;
        // Reference: the vacuum chord above the observation level, the uniform-index
        // chord below it (see the file header).  The vacuum chord is f-independent, so
        // the outward block keeps the plain "deltas from the base" dispersion form;
        // the inward reference does depend on f, so its deltas are built from
        // RESIDUALS (traced minus reference at the same f).
        double Lref[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double Dref[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double dDv = 0.0, dLv = 0.0;
        if (i < i0_) {
          for (int q = 0; q < 5; ++q) {
            InwardRef const ref = inwardRef(h, mu, fs[q]);
            Lref[q] = ref.L;
            Dref[q] = ref.D;
            if (q == 0) { dDv = ref.dD; dLv = ref.dL; }
          }
        } else {
          double const Lvq = vacuumL(h, mu);
          // vacuum h-derivatives (same closed forms as the traced ray, vacuum theta)
          double const ctv = std::min(1.0, (r_obs_ * mu + Lvq) / r);
          double const stv = std::sqrt(std::max(0.0, 1.0 - ctv * ctv));
          Lref[0] = Lvq;
          Dref[0] = vacuumD(h, mu);
          dDv = (ctv > 0.0) ? r_obs_ * (stv / ctv) / r : 0.0;
          dLv = (ctv > 0.0) ? 1.0 / ctv : 0.0;
        }
        // traced h-derivatives = ODE RHS at the node (exact, not differenced).
        // |cta| covers both branches: the inward build runs at cos(theta) < 0.
        double const cta = std::cos(t0[i]);
        double const sta = std::sin(t0[i]);
        double const dDa = (std::abs(cta) > 1e-15) ? r_obs_ * (sta / cta) / r : 0.0;
        double const dLa = (std::abs(cta) > 1e-15) ? 1.0 / cta : 0.0;
        Dd_[k] = float(d0[i] - Dref[0]);
        // optical path P = L + Sigma at f = 1 (Fermat makes dP/df = Sigma exactly)
        DL_[k] = float(l0[i] - Lref[0]);
        Ddh_[k] = float(dDa - dDv);
        DLh_[k] = float(dLa - dLv);
        // Dispersion nodes stored as DELTAS from the base residual: the base uses
        // the bicubic (h,mu) Hermite while these use the cheaper bilinear form, so
        // storing deltas keeps the bilinear error proportional to a small number.
        Dfm2_[k] = float((d2[i] - d0[i]) - (Dref[1] - Dref[0]));
        Dfm1_[k] = float((dA[i] - d0[i]) - (Dref[2] - Dref[0]));
        Dfp1_[k] = float((dP[i] - d0[i]) - (Dref[3] - Dref[0]));
        Dfp2_[k] = float((dQ[i] - d0[i]) - (Dref[4] - Dref[0]));
        DLm2_[k] = float((l2[i] - l0[i]) - (Lref[1] - Lref[0]));
        DLm1_[k] = float((lA[i] - l0[i]) - (Lref[2] - Lref[0]));
        DLp1_[k] = float((lP[i] - l0[i]) - (Lref[3] - Lref[0]));
        DLp2_[k] = float((lQ[i] - l0[i]) - (Lref[4] - Lref[0]));
        S_[k] = float(s0[i]);
        // dSigma/df at f = 1 from the f = 1 -/+ eps traces (Sigma is an integral
        // along the traced ray; the second-order term is measured in the test).
        DSig_[k] = float((sP[i] - sA[i]) / (2.0 * f_eps_));
      }
      // Per-column sanity check: a single NaN or a stalled column (e.g. f = 0.96 or
      // 1.04 at grazing mu) must fail the build loudly instead of leaking a bogus
      // value into the lookups. Strict node-by-node monotonicity is too strong
      // (the residual is not required to be monotone at every 10 m node), so we
      // require finiteness and a strictly positive total column growth.
      for (std::size_t i = 0; i < n_h; ++i) {
        std::size_t const k = i * n_mu_ + j;
        if (!std::isfinite(Dd_[k]) || !std::isfinite(DL_[k]) || !std::isfinite(Dfm2_[k]) ||
            !std::isfinite(Dfm1_[k]) || !std::isfinite(Dfp1_[k]) ||
            !std::isfinite(Dfp2_[k]) || !std::isfinite(DLm2_[k]) ||
            !std::isfinite(DLm1_[k]) || !std::isfinite(DLp1_[k]) ||
            !std::isfinite(DLp2_[k]) || !std::isfinite(S_[k]) || !std::isfinite(DSig_[k]) ||
            !std::isfinite(Ddh_[k])) {
          buildFailCol_ = int(j);
          buildFailReason_ = 3; // non-finite value in this column
          return;
        }
      }
      // h = 0 residuals must vanish (it is the row the two blocks share)
      {
        std::size_t const k0 = i0_ * n_mu_ + j;
        if (!(std::abs(double(Dd_[k0])) < 1.0) || !(std::abs(double(DL_[k0])) < 1.0)) {
          buildFailCol_ = int(j);
          buildFailReason_ = 1; // h = 0 residuals must vanish
          return;
        }
      }
      // Stalled-column check on the OUTWARD block: a constant NON-ZERO column means
      // the trace failed.  (A constant ZERO column is legitimate: the vertical mu = 1
      // ray has no refraction at all, so both residuals vanish identically.)
      // The inward block is excluded on purpose: past the turning point it is filled
      // with the analytic reference, i.e. a constant residual of zero by design.
      std::size_t const kStart = i0_ * n_mu_ + j;
      double mnD = Dd_[kStart], mxD = Dd_[kStart], mnP = DL_[kStart], mxP = DL_[kStart];
      for (std::size_t i = i0_ + 1; i < n_h; ++i) {
        std::size_t const k = i * n_mu_ + j;
        mnD = std::min(mnD, double(Dd_[k]));
        mxD = std::max(mxD, double(Dd_[k]));
        mnP = std::min(mnP, double(DL_[k]));
        mxP = std::max(mxP, double(DL_[k]));
      }
      if (mnD < -1.0 || mnP < -1.0) {
        buildFailCol_ = int(j);
        buildFailReason_ = 4; // residual went negative
        return;
      }
      if (mxD == mnD && mnD != 0.0) {
        buildFailCol_ = int(j);
        buildFailReason_ = 2; // stalled d column
        return;
      }
      if (mxP == mnP && mnP != 0.0) {
        buildFailCol_ = int(j);
        buildFailReason_ = 2; // stalled P column
        return;
      }
      // Inward block: the residual must stay a bounded perturbation of the reference
      // (no sign requirement: refraction can push the arc either way there).  A
      // runaway value would mean the inward integration diverged.
      for (std::size_t i = 0; i < i0_; ++i) {
        std::size_t const k = i * n_mu_ + j;
        if (!(std::abs(double(Dd_[k])) < 1.0e5) || !(std::abs(double(DL_[k])) < 1.0e5) ||
            !(std::abs(double(Ddh_[k])) < 1.0e12) || !(std::abs(double(DLh_[k])) < 1.0e12)) {
          buildFailCol_ = int(j);
          buildFailReason_ = 6; // inward residual out of bounds
          return;
        }
      }
    }

    // ---- mu-derivatives (non-uniform 3-point central differences) ----------
    auto dmu_at = [&](std::vector<float> const& V, std::size_t i, std::size_t j) -> double {
      if (n_mu_ < 3) return 0.0;
      std::size_t const k = i * n_mu_;
      if (j == 0) {
        double const h1 = nodeS(1) - nodeS(0);
        return (double(V[k + 1]) - double(V[k + 0])) / h1;
      }
      if (j + 1 == n_mu_) {
        double const h1 = nodeS(n_mu_ - 1) - nodeS(n_mu_ - 2);
        return (double(V[k + n_mu_ - 1]) - double(V[k + n_mu_ - 2])) / h1;
      }
      double const h1 = nodeS(j) - nodeS(j - 1);
      double const h2 = nodeS(j + 1) - nodeS(j);
      double const c1 = -h2 / (h1 * (h1 + h2));
      double const c2 = (h2 - h1) / (h1 * h2);
      double const c3 = h1 / (h2 * (h1 + h2));
      return c1 * double(V[k + j - 1]) + c2 * double(V[k + j]) + c3 * double(V[k + j + 1]);
    };

    for (std::size_t i = 0; i < n_h; ++i) {
      for (std::size_t j = 0; j < n_mu_; ++j) {
        std::size_t const k = i * n_mu_ + j;
        Ddmu_[k] = float(dmu_at(Dd_, i, j));
        DLmu_[k] = float(dmu_at(DL_, i, j));
        Ddhm_[k] = float(dmu_at(Ddh_, i, j)); // d/dmu of the exact d/dh
        DLhm_[k] = float(dmu_at(DLh_, i, j));
      }
    }

    // ---- validity: reject NaN/inf anywhere ---------------------------------
    valid_ = true;
    for (std::size_t k = 0; k < n && valid_; ++k) {
      if (!std::isfinite(double(Dd_[k])) || !std::isfinite(double(DL_[k])) ||
          !std::isfinite(double(Ddh_[k])) || !std::isfinite(double(DLh_[k])) ||
          !std::isfinite(double(Ddmu_[k])) || !std::isfinite(double(DLmu_[k])) ||
          !std::isfinite(double(Ddhm_[k])) || !std::isfinite(double(DLhm_[k])) ||
          !std::isfinite(double(Dfm2_[k])) || !std::isfinite(double(Dfm1_[k])) ||
          !std::isfinite(double(Dfp1_[k])) || !std::isfinite(double(Dfp2_[k])) ||
          !std::isfinite(double(DLm2_[k])) || !std::isfinite(double(DLm1_[k])) ||
          !std::isfinite(double(DLp1_[k])) || !std::isfinite(double(DLp2_[k])) ||
          !std::isfinite(double(DSig_[k])) || !std::isfinite(S_[k])) {
        valid_ = false;
      }
    }
  }


  inline double CurvedRayTable::memoryMB() const {
    double const n = double(Dd_.size());
    // 18 float arrays per node (17 residuals + Sigma).  Wave 4b removed the tau
    // block (5 lambda nodes, 64.1 MiB): the curved route now takes tau from
    // CurvedAbsorptionTable, so no extinction is stored here at all.
    double bytes = 18.0 * 4.0 * n;
    bytes += 8.0 * double(mu_.size()) + 8.0 * double(z_nodes_.size());
    return bytes / 1048576.0;
  }

  inline std::size_t CurvedRayTable::hCell(double h) const {
    if (z_nodes_.size() < 2) return 0;
    if (h < 0.0) {
      // Inward block: uniform on dh_in_, rows [0, i0_) with row 0 at the ground, so
      // the index counts UP from z_nodes_[0] (not down from 0).
      if (n_in_ == 0) return 0;
      std::size_t i = static_cast<std::size_t>((h - z_nodes_.front()) / dh_in_);
      if (i >= n_in_) i = n_in_ - 1;
      // guard against off-by-one from the index arithmetic; note the STRICT
      // comparison, so a query exactly on a node stays in the cell below it (the
      // cell that touches h = 0 is rejected by trace(), and this keeps its
      // neighbour usable)
      while (i > 0 && z_nodes_[i] > h) --i;
      while (i + 1 < i0_ && z_nodes_[i + 1] < h) ++i;
      return i;
    }
    std::size_t i = i0_;
    if (h < z_break_) {
      i = i0_ + static_cast<std::size_t>(h / dh_lo_);
    } else {
      i = i0_ + n_lo_ - 1 + static_cast<std::size_t>((h - z_break_) / dh_hi_);
    }
    if (i + 1 >= z_nodes_.size()) i = z_nodes_.size() - 2;
    // guard against off-by-one from the two-level index arithmetic
    while (i > i0_ && z_nodes_[i] > h) --i;
    while (i + 2 < z_nodes_.size() && z_nodes_[i + 1] <= h) ++i;
    return i;
  }


  inline void CurvedRayTable::buildColumn(double mu, double f, std::vector<double>& dOut,
                                          std::vector<double>& lOut, std::vector<double>& sOut,
                                          std::vector<double>& thOut) const {
    std::size_t const n_h = z_nodes_.size();
    dOut.resize(n_h);
    lOut.resize(n_h);
    sOut.resize(n_h);
    thOut.resize(n_h);

    double th = std::acos(std::min(1.0, std::max(0.0, mu))); // build-frame zenith at r_obs
    double phi = 0.0, L = 0.0, Sig = 0.0;
    // Row i0_ is h = 0 (the shared row with the inward block); this pass owns
    // [i0_, n_h).
    dOut[i0_] = 0.0;
    lOut[i0_] = 0.0;
    sOut[i0_] = 0.0;
    thOut[i0_] = th;

    // ODE with h as the independent variable:
    //   dtheta/dh = -tan(theta)*(1/r + n'/n)
    //   dPhi/dh   = tan(theta)/r
    //   dL/dh     = 1/cos(theta)
    //   dSigma/dh = N/cos(theta)
    // The last three are exactly the stored h-derivatives (no differencing).
    // The refractive-index log-slope is fixed per SUBSTEP, to the cell the substep
    // starts in, and used by all four RK4 stages.  Taking it from each evaluation
    // point instead (the previous behaviour) let the stage at a step's end read the
    // NEXT cell's slope, because cellIndex() truncates; that made the column
    // integration first order (measured: p = 0.98-1.00, 2.17 m on the Stage 3 ray)
    // instead of converging at the scheme's order.
    double bStep = 0.0; // log-slope of the substep's start cell
    auto rhs = [&](double hh, double t, double& dth, double& dph, double& dl, double& ds) {
      // NOTE: RefractiveIndexTable is indexed by the ABSOLUTE z, i.e. by the
      // radius r = R + height (its z_min_ is the sea-level radius), so r is the
      // right argument here -- never r - R.
      double const r = r_obs_ + hh;
      double const Nz = N(r);
      double const n = 1.0 + Nz * f;
      double const dN = bStep * Nz; // dN/dz at r, from this substep's own cell
      double const npn = (f * dN) / n;                         // n'/n
      double const k = 1.0 / r + npn;
      double const c = std::cos(t), s = std::sin(t);
      double const cc = (c > 1e-12) ? c : 1e-12;
      dth = -s * k / cc;
      dph = s / (r * cc);
      dl = 1.0 / cc;
      ds = Nz / cc;
    };

    // Near-horizontal rays (small mu) have a sqrt-stiff integrand near h = 0
    // (dPhi/dh ~ tan(theta)/r with tan(theta) ~ 1/mu). Scale the substep count
    // with 1/mu so the build stays 4th-order accurate there; the extra cost only
    // applies to the handful of small-mu columns.
    int const nSubBase = nSubBaseFor(mu, f);
    for (std::size_t i = i0_; i + 1 < n_h; ++i) {
      double const hA = z_nodes_[i];
      double const hB = z_nodes_[i + 1];
      int const nSub = (i == i0_) ? 4 * nSubBase : nSubBase; // extra at the start
      double const dh = (hB - hA) / double(nSub);
      double hh = hA;
      for (int m = 0; m < nSub; ++m) {
        bStep = ri_->logSlope(ri_->cellIndex(r_obs_ + hh)); // this substep's own cell
        double k1t, k1p, k1l, k1s, k2t, k2p, k2l, k2s;
        double k3t, k3p, k3l, k3s, k4t, k4p, k4l, k4s;
        rhs(hh, th, k1t, k1p, k1l, k1s);
        rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2p, k2l, k2s);
        rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3p, k3l, k3s);
        rhs(hh + dh, th + dh * k3t, k4t, k4p, k4l, k4s);
        th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
        phi += dh / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
        Sig += dh / 6.0 * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
        hh += dh;
      }
      dOut[i + 1] = r_obs_ * phi; // impact distance d = r_obs * Phi
      lOut[i + 1] = L;
      sOut[i + 1] = Sig;
      thOut[i + 1] = th;
    }
  }


  inline void CurvedRayTable::buildColumnInward(
      double mu, double f, std::vector<double>& dOut, std::vector<double>& lOut,
      std::vector<double>& sOut, std::vector<double>& thOut) const {
    if (n_in_ == 0) return;
    // Seed every inward row with its own analytic reference (residual 0).  Rows past
    // the turning point stay at that seed and are never read: trace() rejects any h
    // cell that reaches inside the apparent horizon.  Seeding (rather than zeroing)
    // keeps the array finite and, should a guard ever slip, gives the reference
    // answer instead of a wild one.
    double const th0 = std::acos(std::min(1.0, std::max(-1.0, -mu)));
    for (std::size_t i = 0; i < i0_; ++i) {
      InwardRef const r1 = inwardRef(z_nodes_[i], mu, f);
      dOut[i] = r1.D;
      lOut[i] = r1.L;
      sOut[i] = 0.0;
      thOut[i] = th0; // keeps the stored h-derivative finite in the filled region
    }

    // ODE with h as the independent variable, advanced from h = 0 DOWNWARD (dh < 0).
    // The build direction is the time reverse of the rising photon: cos(theta) = -mu
    // at r_obs, i.e. theta in (pi/2, pi), and theta decreases towards pi/2 as the ray
    // goes down until it turns (the apparent horizon of that arrival angle).
    //   dtheta/dh = -tan(theta)*(1/r + n'/n)
    //   dPhi/dh   = tan(theta)/r          (Phi grows: dPhi/ds = sin(theta)/r > 0)
    //   dL/dh     = 1/cos(theta)          (both negative: L grows downwards)
    //   dSigma/dh = N/cos(theta)
    double th = th0;
    double phi = 0.0, L = 0.0, Sig = 0.0;
    double hCur = 0.0;
    int const nSubBase = nSubBaseFor(mu, f);
    double bStep = 0.0;
    auto rhs = [&](double hh, double t, double& dth, double& dph, double& dl, double& ds) {
      double const r = r_obs_ + hh;
      double const Nz = N(r);
      double const n = 1.0 + Nz * f;
      double const dN = bStep * Nz; // dN/dr at r, from this substep's own cell
      double const k = 1.0 / r + (f * dN) / n;
      double const c = std::cos(t), s = std::sin(t);
      // cos(theta) < 0 here; the guard keeps its SIGN (the outward pass clamps to
      // +1e-12, which would flip the whole column).
      double const cc = (std::abs(c) > 1e-12) ? c : (c < 0.0 ? -1e-12 : 1e-12);
      dth = -s * k / cc;
      dph = s / (r * cc);
      dl = 1.0 / cc;
      ds = Nz / cc;
    };
    constexpr double kPi = 3.14159265358979323846;
    // Advance the state from hCur to hTarget (either direction).  Returns false when
    // the ray reached its turning point on the way down; hCur then holds where it
    // stopped so the caller can hand over to the upward leg.
    auto advance = [&](double hTarget) -> bool {
      bool const descending = (hTarget < hCur);
      int nSub = nSubBase;
      // Step control near the apparent horizon.  With eps = |theta - pi/2|,
      // d(eps)/dh = c/eps (c = 1/r + n'/n ~ 1.27e-7 m^-1), so as eps falls the fixed
      // 2.5 m substep eventually exceeds eps itself and the RK4 JUMPS THROUGH the
      // tangency (cos(theta) flips sign and the column integrates away to nonsense --
      // measured: |Dd| = 1.7e5 m, build reason 6).  Refine so that d(eps)/eps stays
      // under 5%: dh <= 0.05 eps^2 / c.
      double const epsNow = std::abs(th - 0.5 * kPi);
      if (epsNow > 0.0 && epsNow < 5.0e-2) {
        double const rCur = r_obs_ + hCur;
        double const Nz = N(rCur);
        double const nCur = 1.0 + Nz * f;
        double const cA = 1.0 / rCur + (f * ri_->logSlope(ri_->cellIndex(rCur)) * Nz) / nCur;
        double const dhMax = 0.05 * epsNow * epsNow / std::max(1.0e-30, cA);
        int const nNeed = static_cast<int>(
            std::ceil(std::abs(hTarget - hCur) / std::max(1.0e-3, dhMax)));
        nSub = std::max(nSub, std::min(20000, nNeed));
      }
      double const dh = (hTarget - hCur) / double(nSub);
      double hh = hCur;
      for (int m = 0; m < nSub; ++m) {
        bStep = ri_->logSlope(ri_->cellIndex(r_obs_ + hh)); // this substep's own cell
        double k1t, k1p, k1l, k1s, k2t, k2p, k2l, k2s;
        double k3t, k3p, k3l, k3s, k4t, k4p, k4l, k4s;
        rhs(hh, th, k1t, k1p, k1l, k1s);
        rhs(hh + 0.5 * dh, th + 0.5 * dh * k1t, k2t, k2p, k2l, k2s);
        rhs(hh + 0.5 * dh, th + 0.5 * dh * k2t, k3t, k3p, k3l, k3s);
        rhs(hh + dh, th + dh * k3t, k4t, k4p, k4l, k4s);
        th += dh / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
        phi += dh / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        L += dh / 6.0 * (k1l + 2.0 * k2l + 2.0 * k3l + k4l);
        Sig += dh / 6.0 * (k1s + 2.0 * k2s + 2.0 * k3s + k4s);
        hh += dh;
        // Stop a hair BEFORE the turning point: there cos(theta) -> 0 so the stored
        // dL/dh = 1/cos(theta) = -1/eps diverges.  eps = 1.6e-3 corresponds to a
        // reference clearance r_e - b = r_e eps^2/2 = 8 m, inside the 10 m the lookup
        // guard requires, so no usable cell is ever truncated.
        if (descending && th <= 0.5 * kPi + 1.6e-3) {
          hCur = hh;
          return false;
        }
      }
      hCur = hTarget;
      return true;
    };

    auto storeRow = [&](std::size_t i, std::vector<double>& dO, std::vector<double>& lO,
                        std::vector<double>& sO, std::vector<double>& tO) {
      dO[i] = r_obs_ * phi; // impact distance d = r_obs * Phi (Phi > 0 downwards)
      lO[i] = L;
      sO[i] = Sig;
      tO[i] = th;
    };

    // ---- first leg: down from r_obs, rows [0, i0_) in decreasing index order ------
    // Down from r_obs, rows [0, i0_) in decreasing index order.  A turn means the
    // node is inside the apparent horizon: it keeps its reference seed and trace()
    // rejects every cell that touches it.
    for (std::size_t i = i0_; i-- > 0;) {
      if (!advance(z_nodes_[i])) break;
      storeRow(i, dOut, lOut, sOut, thOut);
    }
  }


  inline CurvedRayTable::Result CurvedRayTable::trace(double x0, double y0, double z0, double dx,
                                                      double dy, double dz, double fPhase,
                                                      double fGroup) const {
    Result res;
    if (!valid_) { res.status = Status::InvalidTable; return res; }

    double const rEmit = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
    if (!(rEmit > 0.0)) { res.status = Status::InvalidInput; return res; }
    double const h = rEmit - r_obs_;
    // Two branches: h >= 0 (emission ABOVE r_obs, the photon must come down) and
    // h < 0 (emission BELOW r_obs, Earth-skimming nu_tau: the photon must rise).
    bool const inward = (h < 0.0);
    if (inward && (n_in_ == 0 || h < z_nodes_.front() - 1.0e-9)) {
      ++n_beyond_; // no inward block (Ground mode) or under the tabulated range
      res.status = Status::EmissionBelowObs;
      return res;
    }
    if (h > z_nodes_.back()) {
      ++n_beyond_; // above the table
      res.status = Status::AboveTable;
      return res;
    }

    // radial unit vector and the photon's radial direction cosine
    double const irr = 1.0 / rEmit;
    double const rx = x0 * irr, ry = y0 * irr, rz = z0 * irr;
    double const cRad = dx * rx + dy * ry + dz * rz;
    // Above r_obs the photon has to move inward (cRad < 0) to reach it.  Below r_obs
    // it has to move outward (cRad > 0) to rise to it; a photon emitted below and
    // going DOWN dips through its perigee and rises, which is a DIFFERENT ray family
    // for the same (h, mu) -- it is not tabulated (see the note at the top of the
    // file) and falls back here.
    if (inward ? (cRad <= 0.0) : (cRad >= 0.0)) {
      ++n_upgoing_; // moving away from the observation level / descending emission
      res.status = Status::Upgoing;
      return res;
    }

    // Bouguer invariant, using the photon's own phase index.
    // RefractiveIndexTable is indexed by the ABSOLUTE z = radius (R + height).
    double const nObs = 1.0 + N(r_obs_) * fPhase;
    double const pCrit = nObs * r_obs_;
    double const nEmit = 1.0 + N(rEmit) * fPhase;
    double const sinTh = std::sqrt(std::max(0.0, 1.0 - cRad * cRad));
    double const p = nEmit * rEmit * sinTh;
    if (p >= pCrit) {
      ++n_beyond_; // misses the observation level
      res.status = Status::BeyondCritical;
      return res;
    }
    // The dispersion (f_phase) response is only tabulated between f_lo_ and f_hi_;
    // outside that range we must NOT extrapolate -- fall back to the exact tracer.
    if (fPhase < f_lo_ || fPhase > f_hi_) {
      ++n_disp_;
      res.status = Status::DispersionOutOfRange;
      return res;
    }

    double mu = std::sqrt(std::max(0.0, 1.0 - (p / pCrit) * (p / pCrit)));
    // Grazing cutoff: below mu_min_ the near-horizontal integrand is too stiff for
    // the tabulated resolution, so reject and let the caller fall back to the exact
    // tracer (never return a low-accuracy value).
    if (mu < mu_min_) {
      ++n_murej_;
      res.status = Status::Grazing;
      return res;
    }
    mu = std::min(1.0, mu);

    // ---- cell location -----------------------------------------------------
    std::size_t const i = hCell(h);
    double const s = muAxis(mu); // uniform grid coordinate
    std::size_t j = static_cast<std::size_t>(s / (s_tot_ / double(n_mu_ - 1)));
    if (j + 1 >= n_mu_) j = n_mu_ - 2;
    while (j > 0 && nodeS(j) > s) --j;
    while (j + 2 < n_mu_ && nodeS(j + 1) <= s) ++j;
    double const dH = z_nodes_[i + 1] - z_nodes_[i];
    double const dM = nodeS(j + 1) - nodeS(j);
    double const th = (h - z_nodes_[i]) / dH;
    double const tm = (s - nodeS(j)) / dM;

    // ---- inward-block validity guard ---------------------------------------
    // Two conditions, both of which must hold for the whole 2x2 stencil:
    //  * the cell must not touch h = 0: there the stored h-derivative is the OUTWARD
    //    one (the two branches have opposite dL/dh at the shared row), so the cell
    //    [ -dh_in, 0 ] has no usable h-derivative;
    //  * the cell must clear the apparent horizon.  b = p/n(r_emit) is the geometric
    //    impact parameter of the reference chord, whose clearance is r_emit - b =
    //    r_emit*eps^2/2 with eps = theta - pi/2.  The build stops integrating at
    //    eps = 1.6e-3 (8 m of clearance) and every deeper row is left at the bare
    //    reference, so the guard demands dh_in = 10 m: a deeply conservative margin
    //    (eps >= 1.8e-3) that no truncated row can satisfy.
    //    mu_[j-1] is used because b grows as mu falls, and the bicubic also reads the
    //    mu-derivatives of column j, which are differenced against column j-1.
    if (i < i0_) {
      std::size_t const jlo = (j > 0) ? (j - 1) : 0;
      double const nEmit = 1.0 + N(rEmit) * fPhase;
      double const st0lo = std::sqrt(std::max(0.0, 1.0 - mu_[jlo] * mu_[jlo]));
      double const bRef = (nEmit > 1.0e-12) ? (pCrit * st0lo / nEmit) : 0.0;
      // WAVE 5 EXPERIMENT: the impact error is a CLIFF at the apparent horizon (66.6 m below
      // mu 0.02, 5.5 m in 0.02-0.03, 0.012 m above 0.03) and the physical cause is that d
      // diverges as the ray approaches tangency -- no grid can represent that, which is why
      // doubling the resolution changed nothing.  The table already refuses rays within
      // dh_in_ (10 m) of the horizon; this widens the band to kGrazingClearance_m, sending
      // the whole cliff to the exact tracer (0.50 m class geometry) at the cost of its own
      // tau on those photons.
      if (i + 1 >= i0_ || (rEmit - bRef) < std::max(dh_in_, grazingClearance_m_)) {
        ++n_murej_;
        res.status = Status::Grazing;
        return res;
      }
    }

    // ---- bicubic Hermite in (h, mu), including the cross term --------------
    auto b00 = [](double t) { return ((2.0 * t - 3.0) * t) * t + 1.0; };
    auto b10 = [](double t) { return ((t - 2.0) * t + 1.0) * t; };
    auto b01 = [](double t) { return (3.0 - 2.0 * t) * t * t; };
    auto b11 = [](double t) { return (t - 1.0) * t * t; };
    double const Ah[2] = {b00(th), b01(th)};
    double const Bh[2] = {dH * b10(th), dH * b11(th)};
    double const Am[2] = {b00(tm), b01(tm)};
    double const Bm[2] = {dM * b10(tm), dM * b11(tm)};
    double const Lh[2] = {1.0 - th, th}; // bilinear basis (Sigma, f-sensitivity)
    double const Lm[2] = {1.0 - tm, tm};

    double Dd = 0.0, DL = 0.0, Sig = 0.0, DSig = 0.0;
    double Sm2 = 0.0, Sm1 = 0.0, Sp1 = 0.0, Sp2 = 0.0;
    double Tm2 = 0.0, Tm1 = 0.0, Tp1 = 0.0, Tp2 = 0.0;
    for (int p = 0; p < 2; ++p) {
      for (int q = 0; q < 2; ++q) {
        std::size_t const k = (i + p) * n_mu_ + (j + q);
        Dd += Ah[p] * Am[q] * double(Dd_[k]) + Bh[p] * Am[q] * double(Ddh_[k]) +
              Ah[p] * Bm[q] * double(Ddmu_[k]) + Bh[p] * Bm[q] * double(Ddhm_[k]);
        DL += Ah[p] * Am[q] * double(DL_[k]) + Bh[p] * Am[q] * double(DLh_[k]) +
              Ah[p] * Bm[q] * double(DLmu_[k]) + Bh[p] * Bm[q] * double(DLhm_[k]);
        Sm2 += Lh[p] * Lm[q] * double(Dfm2_[k]);
        Sm1 += Lh[p] * Lm[q] * double(Dfm1_[k]);
        Sp1 += Lh[p] * Lm[q] * double(Dfp1_[k]);
        Sp2 += Lh[p] * Lm[q] * double(Dfp2_[k]);
        Tm2 += Lh[p] * Lm[q] * double(DLm2_[k]);
        Tm1 += Lh[p] * Lm[q] * double(DLm1_[k]);
        Tp1 += Lh[p] * Lm[q] * double(DLp1_[k]);
        Tp2 += Lh[p] * Lm[q] * double(DLp2_[k]);
        Sig += Lh[p] * Lm[q] * S_[k];
        DSig += Lh[p] * Lm[q] * double(DSig_[k]);
      }
    }
    // Dispersion: quartic Lagrange in t = (f-1)/eps over the five traced nodes
    // t = -2,-1,0,+1,+2 (f = 1-2eps .. 1+2eps). The node values are stored as
    // DELTAS from the base, and the weights sum to 1, so the base keeps weight 1:
    //     V(f) = V(1) + sum_k w_k * delta_k .
    // Multiplying the base by the centre weight (as if the node values were
    // absolute) collapses the result to the bare delta at a node.
    double const t = (fPhase - 1.0) / f_eps_; // inside [f_lo_, f_hi_] by construction
    double const tp1 = t + 1.0, tm1 = t - 1.0, tm2 = t - 2.0, tp2 = t + 2.0;
    double const wm2 = tp1 * t * tm1 * tm2 / 24.0;  // t = -2
    double const wm1 = -tp2 * t * tm1 * tm2 / 6.0;  // t = -1
    double const wp1 = -tp2 * tp1 * t * tm2 / 6.0;  // t = +1
    double const wp2 = tp2 * tp1 * t * tm1 / 24.0;  // t = +2
    Dd += wm2 * Sm2 + wm1 * Sm1 + wp1 * Sp1 + wp2 * Sp2;
    DL += wm2 * Tm2 + wm1 * Tm1 + wp1 * Tp1 + wp2 * Tp2;

    // Reference at the query point: the vacuum chord above r_obs, the uniform-index
    // chord below it -- the same functions the build subtracted, evaluated with the
    // photon's own f_phase.
    double Lref = 0.0, Dref = 0.0;
    if (inward) {
      InwardRef const ref = inwardRef(h, mu, fPhase);
      Lref = ref.L;
      Dref = ref.D;
    } else {
      Lref = vacuumL(h, mu);
      Dref = vacuumD(h, mu);
    }
    double const d = Dref + Dd;
    // Transverse component of the emission direction.  It is exactly zero for a radial
    // emission (mu = 1), where the transverse direction is undefined and tn is a pure
    // round-off residual.  That is a VALID trace, not a table failure: the arrival point
    // is the observation sphere along the emission radial (Phi = 0, so d = 0) and the
    // arrival direction is the radial itself.  The tolerance is dimensionless, since tn
    // is sin(theta) of the emission direction: the round-off floor is ~1e-16, while the
    // most nearly vertical ray in the eligible domain is mu = 0.6, i.e. sin(theta) = 0.8,
    // so 1e-9 sits orders of magnitude above the floor and below any real ray.  A ray
    // with sin(theta) < 1e-9 would be displaced by under R*1e-9 ~ 6 mm anyway, so the
    // cut is safe.
    double tx = dx - cRad * rx, ty = dy - cRad * ry, tz = dz - cRad * rz;
    double const tn = std::sqrt(tx * tx + ty * ty + tz * tz);
    bool const degenerate = !(tn > 1e-9);
    // Time is the group optical path along the PHASE ray:
    //     t*c = L(f_phase) + f_group * Sigma(f_phase) .
    // Both f-dependences matter. L moves by ~4.4e4 m per unit f (carried by the
    // DL* node deltas), and Sigma moves by ~8 m per unit f (carried by DSig_);
    // dropping the Sigma term leaves ~1 ns at grazing mu. Note that Fermat does
    // NOT let us fold this into a single "optical path" quantity: stationarity
    // needs fixed endpoints, but the table fixes (r_obs, mu_obs, h), so the
    // arrival point on the r_obs + h sphere moves with f.
    double const SigF = Sig + (fPhase - 1.0) * DSig;
    double const Lpath = Lref + DL + fGroup * SigF;
    double const Phi = degenerate ? 0.0 : d / r_obs_;

    // ---- impact point and arrival direction --------------------------------
    // No early return here: a degenerate (radial) emission is a valid trace.  The
    // previous `if (!(tn > 0.0)) return res;` left res.status at its default
    // Status::InvalidTable, and where round-off made tn non-zero it built the arrival
    // direction from a meaningless transverse vector (measured: a spurious 93 m
    // displacement at h = 46000).  Zeroing the transverse vector makes the formulas
    // below produce the exact radial answer.
    if (!degenerate) {
      tx /= tn;
      ty /= tn;
      tz /= tn;
    }
    // build-frame transverse direction at r_obs: t0 = -(r_emit_hat x t_hat) x r_obs_hat
    double const nx = ry * tz - rz * ty, ny = rz * tx - rx * tz, nz = rx * ty - ry * tx;
    double const cp = std::cos(Phi), sp = std::sin(Phi);
    double const ox = cp * rx + sp * tx, oy = cp * ry + sp * ty, oz = cp * rz + sp * tz;
    double const t0x = -(ny * oz - nz * oy);
    double const t0y = -(nz * ox - nx * oz);
    double const t0z = -(nx * oy - ny * ox);
    double const st0 = std::sqrt(std::max(0.0, 1.0 - mu * mu));

    // Arrival direction.  The radial term carries the branch sign: above r_obs the
    // photon arrives going DOWN (cos = -mu), below it arrives RISING (cos = +mu).
    // The transverse term is the SAME in both branches: t0 is the time-reversed
    // build transverse direction, which the rotation by Phi above already flipped.
    double const sgn = inward ? 1.0 : -1.0;
    res.x = r_obs_ * ox;
    res.y = r_obs_ * oy;
    res.z = r_obs_ * oz;
    res.dx = sgn * mu * ox - st0 * t0x;
    res.dy = sgn * mu * oy - st0 * t0y;
    res.dz = sgn * mu * oz - st0 * t0z;
    double const dn = std::sqrt(res.dx * res.dx + res.dy * res.dy + res.dz * res.dz);
    if (dn > 0.0) {
      res.dx /= dn;
      res.dy /= dn;
      res.dz /= dn;
    }
    res.theta_obs = std::acos(std::min(1.0, std::max(-1.0, -sgn * mu)));
    res.time_s = Lpath / 299792458.0;
    res.geometric_path_m = Lpath - fGroup * SigF; // = vacuumL + DL
    res.optical_path_m = Lpath;                   // GROUP optical path
    res.sigma_m = SigF;
    res.hit = true;
    res.status = Status::Hit;
    return res;
  }

  inline double CurvedRayTable::sigmaBilinear(double h, double mu) const {
    if (!valid_ || z_nodes_.size() < 2) return 0.0;
    h = std::min(std::max(h, 0.0), z_nodes_.back());
    mu = std::min(1.0, std::max(mu_[0], mu));
    std::size_t const i = hCell(h);
    double const ss = muAxis(mu); // uniform grid coordinate
    std::size_t j = static_cast<std::size_t>(ss / (s_tot_ / double(n_mu_ - 1)));
    if (j + 1 >= n_mu_) j = n_mu_ - 2;
    while (j > 0 && nodeS(j) > ss) --j;
    while (j + 2 < n_mu_ && nodeS(j + 1) <= ss) ++j;
    double const th = (h - z_nodes_[i]) / (z_nodes_[i + 1] - z_nodes_[i]);
    double const tm = (ss - nodeS(j)) / (nodeS(j + 1) - nodeS(j));
    double const Lh[2] = {1.0 - th, th};
    double const Lm[2] = {1.0 - tm, tm};
    double s = 0.0;
    for (int p = 0; p < 2; ++p)
      for (int q = 0; q < 2; ++q) s += Lh[p] * Lm[q] * S_[(i + p) * n_mu_ + (j + q)];
    return s;
  }

  // ==========================================================================
  // wave 5: on-disk build cache
  // ==========================================================================
  inline std::uint64_t CurvedRayTable::buildKey(RefractiveIndexTable const& ri, double R,
                                                double rObs, double hTop, int nMu, double dhLo,
                                                double dhHi, double zBreak, double fEps) {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&h](double v) {
      std::uint64_t b = 0;
      std::memcpy(&b, &v, sizeof b);
      h = (h ^ b) * 1099511628211ull;
    };
    mix(R); mix(rObs); mix(hTop); mix(double(nMu));
    mix(dhLo); mix(dhHi); mix(zBreak); mix(fEps);
    mix(double(ri.size())); mix(ri.zMin()); mix(ri.zMax());
    // Sample the PROFILE too: the parameters alone are not enough, because two runs
    // can share every number here and still differ in N(z) (a different atmosphere
    // file).  Sampling N and dNdz makes that difference change the key.
    for (int i = 0; i < 32; ++i) {
      double const t = (double(i) + 0.5) / 32.0;
      double const z = ri.zMin() + t * (ri.zMax() - ri.zMin());
      mix(ri.N(z));
      mix(ri.dNdz(z));
    }
    return h;
  }

  inline bool CurvedRayTable::saveToFile(std::string const& path, std::uint64_t key) const {
    if (!valid_) return false;
    std::string const tmp = path + ".tmp." + std::to_string(key);
    std::ofstream f(tmp, std::ios::binary | std::ios::trunc);
    if (!f) return false;
    std::uint64_t sum = 1469598103934665603ull;
    auto put = [&](void const* p, std::size_t n) {
      f.write(static_cast<char const*>(p), static_cast<std::streamsize>(n));
      auto const* c = static_cast<unsigned char const*>(p);
      for (std::size_t i = 0; i < n; ++i) sum = (sum ^ c[i]) * 1099511628211ull;
    };
    auto putD = [&](double v) { put(&v, sizeof v); };
    auto putU = [&](std::uint64_t v) { put(&v, sizeof v); };
    char magic[8] = {'C', 'R', 'V', 'R', 'T', 'B', '4', '\0'};
    std::uint32_t const ver = 4u;
    put(magic, sizeof magic);
    put(&ver, sizeof ver);
    put(&key, sizeof key);
    putD(R_); putD(r_obs_); putD(z_break_); putD(dh_lo_); putD(dh_hi_); putD(dh_in_);
    putD(mu_b_); putD(u_b_); putD(s_b_); putD(K_); putD(s_tot_);
    putD(f_lo_); putD(f_hi_); putD(f_eps_); putD(mu_min_);
    putU(std::uint64_t(n_lo_)); putU(std::uint64_t(i0_)); putU(std::uint64_t(n_in_));
    putU(std::uint64_t(n_mu_));
    putU(std::uint64_t(mu_.size())); putU(std::uint64_t(z_nodes_.size()));
    std::vector<float> const* arr[18] = {&Dd_,  &Ddh_,  &Ddmu_, &Ddhm_, &DL_,   &DLh_,
                                         &DLmu_, &DLhm_, &Dfm2_, &Dfm1_, &Dfp1_, &Dfp2_,
                                         &DLm2_, &DLm1_, &DLp1_, &DLp2_, &S_,   &DSig_};
    for (auto const* v : arr) {
      putU(std::uint64_t(v->size()));
      put(v->data(), v->size() * sizeof(float));
    }
    put(mu_.data(), mu_.size() * sizeof(double));
    put(z_nodes_.data(), z_nodes_.size() * sizeof(double));
    f.write(reinterpret_cast<char const*>(&sum), sizeof sum);  // trailing checksum
    f.flush();
    if (!f) { std::remove(tmp.c_str()); return false; }
    f.close();
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {  // atomic on POSIX
      std::remove(tmp.c_str());
      return false;
    }
    return true;
  }

  inline bool CurvedRayTable::loadFromFile(std::string const& path, RefractiveIndexTable const& ri,
                                           std::uint64_t key, std::string* why) {
    auto fail = [why](char const* m) {
      if (why) *why = m;
      return false;
    };
    std::ifstream f(path, std::ios::binary);
    if (!f) return fail("cannot open");
    std::uint64_t sum = 1469598103934665603ull;
    auto get = [&](void* p, std::size_t n) -> bool {
      f.read(static_cast<char*>(p), static_cast<std::streamsize>(n));
      if (!f) return false;
      auto const* c = static_cast<unsigned char const*>(p);
      for (std::size_t i = 0; i < n; ++i) sum = (sum ^ c[i]) * 1099511628211ull;
      return true;
    };
    char magic[8] = {};
    std::uint32_t ver = 0;
    std::uint64_t fileKey = 0;
    if (!get(magic, sizeof magic)) return fail("truncated header");
    if (std::memcmp(magic, "CRVRTB4", 7) != 0) return fail("wrong magic (older cache format?)");
    if (!get(&ver, sizeof ver)) return fail("truncated version");
    if (ver != 4u) return fail("version mismatch");
    if (!get(&fileKey, sizeof fileKey)) return fail("truncated key");
    if (fileKey != key) return fail("build key mismatch (profile or parameters changed)");

    // Fill a scratch instance so a rejected file can never leave a half-built object
    // reachable through this one.
    CurvedRayTable t;
    auto getD = [&](double& v) { return get(&v, sizeof v); };
    auto getU = [&](std::uint64_t& v) { return get(&v, sizeof v); };
    if (!getD(t.R_) || !getD(t.r_obs_) || !getD(t.z_break_) || !getD(t.dh_lo_) ||
        !getD(t.dh_hi_) || !getD(t.dh_in_) || !getD(t.mu_b_) || !getD(t.u_b_) || !getD(t.s_b_) ||
        !getD(t.K_) || !getD(t.s_tot_) || !getD(t.f_lo_) || !getD(t.f_hi_) || !getD(t.f_eps_) ||
        !getD(t.mu_min_))
      return fail("truncated scalars");
    std::uint64_t nlo = 0, i0 = 0, nin = 0, nmu = 0, nmus = 0, nzs = 0;
    if (!getU(nlo) || !getU(i0) || !getU(nin) || !getU(nmu) || !getU(nmus) || !getU(nzs))
      return fail("truncated sizes");
    t.n_lo_ = std::size_t(nlo);
    t.i0_ = std::size_t(i0);
    t.n_in_ = std::size_t(nin);
    t.n_mu_ = std::size_t(nmu);
    if (t.n_mu_ < 2 || nmus != t.n_mu_) return fail("inconsistent mu axis size");
    if (nzs < 2) return fail("inconsistent h axis size");
    std::vector<float>* arr[18] = {&t.Dd_,   &t.Ddh_,  &t.Ddmu_, &t.Ddhm_, &t.DL_,   &t.DLh_,
                                   &t.DLmu_, &t.DLhm_, &t.Dfm2_, &t.Dfm1_, &t.Dfp1_, &t.Dfp2_,
                                   &t.DLm2_, &t.DLm1_, &t.DLp1_, &t.DLp2_, &t.S_,    &t.DSig_};
    for (auto* v : arr) {
      std::uint64_t n = 0;
      if (!getU(n)) return fail("truncated array size");
      if (n != nmu * nzs) return fail("inconsistent array size");
      v->resize(std::size_t(n));
      if (!get(v->data(), v->size() * sizeof(float))) return fail("truncated array payload");
    }
    t.mu_.resize(std::size_t(nmus));
    if (!get(t.mu_.data(), t.mu_.size() * sizeof(double))) return fail("truncated mu axis");
    t.z_nodes_.resize(std::size_t(nzs));
    if (!get(t.z_nodes_.data(), t.z_nodes_.size() * sizeof(double)))
      return fail("truncated h axis");
    std::uint64_t fileSum = 0;
    if (!f.read(reinterpret_cast<char*>(&fileSum), sizeof fileSum))
      return fail("truncated checksum");
    if (fileSum != sum) return fail("checksum mismatch (corrupt or partially written file)");

    t.ri_ = &ri;
    t.valid_ = true;
    t.buildFailCol_ = -1;
    t.buildFailReason_ = 0;
    *this = std::move(t);
    return true;
  }

} // namespace corsika::cherenkov
