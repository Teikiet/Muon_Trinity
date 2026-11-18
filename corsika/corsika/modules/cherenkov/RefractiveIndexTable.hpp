#pragma once

/*
 * RefractiveIndexTable.hpp
 *
 * Vertical refractive-index table used by the "closed-form" straight-line timing
 * path of CherenkovProcess.
 *
 * The application's runtime lookup (refrIndexFunc) evaluates the refractivity
 * N(z) = (n-1)(z) by linearly interpolating ln(N) on a uniform altitude grid,
 * clamping below the grid and exponentially extrapolating above it (with a hard
 * floor at N = 1e-10 -> N = 0). The full index factorises as
 *
 *      n_phase(z,lam) = 1 + N(z) * f_phase(lam)
 *      n_group(z,lam) = 1 + N(z) * f_group(lam)
 *
 * For a straight photon path dz/ds is constant (= u_z), so the group-delay
 * integral collapses to
 *
 *      int N(z(s)) ds = (A(z_end) - A(z_start)) / u_z ,
 *
 * with the vertical column integral A(z) = int_{z_min}^{z} N dz'. This table
 * stores N and A on the SAME uniform grid and integrates the SAME
 * piecewise-exponential form that N(z) uses, so N and A are mutually consistent
 * and the closed form is exact (to the accuracy of the interpolation).
 *
 * The in-cell A lookup and the cumulative cell sum share the same small-|b|
 * branch, so they agree in the degenerate limit.
 *
 * Header-only helper; no I/O, no framework dependencies.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace corsika::cherenkov {

  class RefractiveIndexTable {
  public:
    RefractiveIndexTable() = default;

    /**
     * Build the table.
     *
     * @param alt_km      altitudes [km] on a uniform grid (the fine grid used by
     *                    the runtime lookup)
     * @param lnN         ln(n-1) at those altitudes (the values the runtime lookup
     *                    interpolates)
     * @param seaLevelZ_m absolute z [m] of the first grid node, i.e.
     *                    R_earth + alt_km.front()*1000
     */
    RefractiveIndexTable(std::vector<double> const& alt_km,
                         std::vector<double> const& lnN, double seaLevelZ_m);

    bool valid() const { return valid_; }
    std::size_t size() const { return N_.size(); }

    /// Refractivity N(z) = (n-1) at absolute z [m]; same clamp/extrapolation rule
    /// as the runtime lookup.
    double N(double abs_z_m) const;

    /// Vertical column integral A(z) = int_{z_min}^{z} N dz' [m]; A(z_min) = 0.
    double A(double abs_z_m) const;

    double zMin() const { return z_min_; }
    double zMax() const { return z_min_ + double(N_.size() - 1) * dz_; }

    /// Cell index used by N()/A() for absolute z [m], clamped to [0, size-2].
    /// Same indexing as the runtime lookup (direct, O(1), uniform grid).
    std::size_t cellIndex(double abs_z_m) const {
      if (!valid_ || N_.size() < 2) return 0;
      double const z = abs_z_m - z_min_;
      if (z <= 0.0) return 0;
      std::size_t const i = static_cast<std::size_t>(z / dz_);
      return (i >= N_.size() - 1) ? (N_.size() - 2) : i;
    }

    /// Log-slope b_i = ln(N_{i+1}/N_i)/dz of cell i (0 if a node is non-positive).
    double logSlope(std::size_t const i) const {
      return (valid_ && i < B_.size()) ? B_[i] : 0.0;
    }

    /// dN/dz at absolute z [m] using the stored piecewise-exponential log-slope
    /// (no finite differences): N(z) = N_i*exp(b_i*(z-z_i))  =>  N'(z) = b_i*N(z).
    double dNdz(double abs_z_m) const { return logSlope(cellIndex(abs_z_m)) * N(abs_z_m); }

  private:
    /// b_i = ln(N_{i+1}/N_i) / dz (0 when either node is non-positive)
    double cellB(std::size_t i) const {
      if (N_[i] <= 0.0 || N_[i + 1] <= 0.0) return 0.0;
      return std::log(N_[i + 1] / N_[i]) / dz_;
    }

    double z_min_ = 0.0;    ///< absolute z [m] of node 0
    double dz_ = 0.0;       ///< uniform node spacing [m]
    std::vector<double> N_; ///< refractivity at the nodes
    std::vector<double> A_; ///< cumulative vertical integral at the nodes (A_[0]=0)
    /// Pre-computed cell log-slopes b_i.  The tracers evaluate N and dN/dz at four RK4
    /// stages per step, so recomputing ln(N_{i+1}/N_i) per lookup dominated their index
    /// cost (one log each).  Storing it costs 8 bytes per node (96 kB for the 12000-node
    /// production grid) and turns dN/dz into a multiply, which is what makes the
    /// stage-local kappa of the Phase-1 geometry fix cheaper than the frozen one it
    /// replaced.  Empty when the table is invalid.
    std::vector<double> B_;
    double N_min_ = 0.0;    ///< N at node 0 (clamp value below the grid)
    double N_max_ = 0.0;    ///< N at the last node
    double H_ext_ = 0.0;    ///< exponential scale height above the grid [m]
    bool valid_ = false;
  };

  // ---------------------------------------------------------------------------
  // implementation
  // ---------------------------------------------------------------------------

  inline RefractiveIndexTable::RefractiveIndexTable(std::vector<double> const& alt_km,
                                                    std::vector<double> const& lnN,
                                                    double seaLevelZ_m) {
    if (alt_km.size() < 3 || alt_km.size() != lnN.size()) return;

    std::size_t const n = alt_km.size();
    z_min_ = seaLevelZ_m + alt_km.front() * 1000.0;
    dz_ = (alt_km[1] - alt_km[0]) * 1000.0;
    if (!(dz_ > 0.0)) return;

    // Require a (near-)uniform grid so the O(1) indexing is valid.
    for (std::size_t i = 1; i < n; ++i) {
      double const d = (alt_km[i] - alt_km[i - 1]) * 1000.0;
      if (std::abs(d - dz_) > 1e-6 * dz_) return;
    }

    N_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      N_[i] = (lnN[i] > -690.0) ? std::exp(lnN[i]) : 0.0;
    }
    N_min_ = N_.front();
    N_max_ = N_.back();

    // Exponential scale height above the grid, from the last two nodes (identical
    // to the runtime extrapolation).
    double const ln_ratio = lnN[n - 1] - lnN[n - 2];
    H_ext_ = (ln_ratio < -1e-12) ? (dz_ / (-ln_ratio)) : dz_;

    // Cumulative vertical integral, integrating the piecewise-exponential N
    // exactly. The cell sum and the in-cell lookup share the same small-|b| branch.
    // The cell log-slopes are cached first: N(), A(), dNdz() and the tracers all use
    // them, and recomputing ln(N_{i+1}/N_i) per call was several logs per evaluation.
    B_.resize(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) B_[i] = cellB(i);
    A_.resize(n);
    A_[0] = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      double const b = B_[i];
      double dA;
      if (std::abs(b) < 1e-12) {
        dA = dz_ * 0.5 * (N_[i] + N_[i + 1]); // trapezoid (N ~ constant)
      } else {
        dA = (N_[i + 1] - N_[i]) / b; // exact: (N_{i+1}-N_i)/b
      }
      A_[i + 1] = A_[i] + dA;
    }
    valid_ = true;
  }

  inline double RefractiveIndexTable::N(double abs_z_m) const {
    if (!valid_) return 0.0;
    double const z = abs_z_m - z_min_;
    if (z <= 0.0) return N_min_; // clamp below the grid
    double const idx = z / dz_;
    std::size_t const i = static_cast<std::size_t>(idx);
    if (i >= N_.size() - 1) {
      // Exponential extrapolation above the grid, with the 1e-10 floor.
      if (H_ext_ <= 0.0) return N_max_;
      double const dz_up = abs_z_m - zMax();
      double const lnN = std::log(N_max_ > 1e-300 ? N_max_ : 1e-300) - dz_up / H_ext_;
      return (lnN < std::log(1e-10)) ? 0.0 : std::exp(lnN);
    }
    double const t = idx - double(i); // fractional cell coordinate in [0,1)
    double const b = B_[i];
    if (std::abs(b) < 1e-12) return N_[i] + t * (N_[i + 1] - N_[i]);
    return N_[i] * std::exp(b * t * dz_); // == exp(lerp(lnN))
  }

  inline double RefractiveIndexTable::A(double abs_z_m) const {
    if (!valid_) return 0.0;
    double const z = abs_z_m - z_min_;
    if (z <= 0.0) {
      // Below the grid N is clamped to N_min (constant) -> linear extension.
      return N_min_ * z;
    }
    double const idx = z / dz_;
    std::size_t const i = static_cast<std::size_t>(idx);
    if (i >= A_.size() - 1) {
      // Above the grid: integrate N_max*exp(-(z'-z_max)/H_ext), saturating once the
      // 1e-10 floor sets N = 0.
      if (H_ext_ <= 0.0) return A_.back();
      double const dz_up = abs_z_m - zMax();
      double const dz_floor =
          (N_max_ > 1e-10) ? H_ext_ * (std::log(N_max_) - std::log(1e-10)) : 0.0;
      double const dz_eff = std::min(dz_up, dz_floor);
      return A_.back() + N_max_ * H_ext_ * (1.0 - std::exp(-dz_eff / H_ext_));
    }
    // Inside cell i.
    double const t = z - double(i) * dz_; // local offset [m]
    double const b = B_[i];
    if (std::abs(b) < 1e-12) {
      // N linear in z -> trapezoid-consistent quadratic integral.
      return A_[i] + N_[i] * t + 0.5 * (N_[i + 1] - N_[i]) * (t * t) / dz_;
    }
    // Exact integral of N_i*exp(b*t'): A_i + (N(z)-N_i)/b.
    double const Nz = N_[i] * std::exp(b * t);
    return A_[i] + (Nz - N_[i]) / b;
  }

} // namespace corsika::cherenkov
