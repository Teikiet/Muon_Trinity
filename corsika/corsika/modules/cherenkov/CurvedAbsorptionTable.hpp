#pragma once

/*
 * CurvedAbsorptionTable.hpp
 *
 * Wave 2c: the absorption of the curved-ray track as a TABLE LOOKUP instead of an
 * integral, for BOTH routes (straight chord and refracted ray).
 *
 * Physics and storage
 * -------------------
 * Along a refracted ray the Bouguer invariant p = n(r) r sin(psi) is conserved, so a
 * single number labels the whole ray.  Normalising it by p_crit = n(r_obs) r_obs gives
 *
 *     xi = p / p_crit = sin(theta_arrival)   in (0, 1]
 *
 * and we tabulate the accumulated optical depth of that ray from the GROUND up to
 * height h:
 *
 *     T(h, xi, lambda) = int_0^h  alpha(h', lambda) ds(h'),   ds = dh' / cos(psi(h'))
 *
 * Two heights on the same ray then give the path integral between them exactly:
 *
 *     tau(A -> B) = |T(hA, xi, lambda) - T(hB, xi, lambda)|,
 *
 * because p (hence xi) is the same at both ends.  That is the whole content of the
 * table: no air mass, no Chapman function, no per-photon quadrature.
 *
 * This is candidate A of diagnostics/b_table (README sections 9, 22, 24), promoted to a
 * production class.  The axes are the ones that diagnostic validated (56 h x 256 xi),
 * but with FIVE wavelength nodes instead of 121: alpha(lambda) is smooth and Rayleigh
 * dominated, and log-log interpolation in lambda between 300/350/400/470/600 nm was
 * measured at 3.4e-3 on the tau of the geometry table (the 400 nm node is mandatory).
 * That takes the table from 6.6 MB to 0.29 MB.
 *
 * Domain, and why a lookup can legitimately fail
 * ----------------------------------------------
 * T diverges as h approaches the tangency of that ray (cos psi -> 0, the air mass
 * diverges like 1/sqrt(h_t - h)).  Cells inside the tangency are therefore NOT filled
 * (NaN) and tauAt() reports InsideTangency; the caller must fall back.  Likewise xi
 * below the tabulated minimum is rejected rather than clamped.  Above the profile top
 * (50 km) alpha is identically zero, so clamping there is exact, not an approximation.
 *
 * The table depends on {N(z) content, dispersion model, R, r_obs} -- the invalidation
 * contract of README section 8.2 -- and cannot be shared across configurations that
 * differ in any of them.  It is built once per process and const/stateless afterwards,
 * so it is safe to share between threads.
 */

#include <corsika/modules/cherenkov/RefractiveIndexTable.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace corsika::cherenkov {

  class CurvedAbsorptionTable {
  public:
    using ExtinctionFn = std::function<double(double, double)>; // alpha(h_asl, lambda_nm)

    enum class Status {
      Ok,
      InvalidTable,    ///< not built / not valid
      BelowXiMin,      ///< grazing ray the table does not resolve
      InsideTangency,  ///< h is at or below the tangency of that ray: T diverges
    };
    static char const* statusName(Status s) {
      switch (s) {
        case Status::Ok: return "ok";
        case Status::InvalidTable: return "invalid-table";
        case Status::BelowXiMin: return "below-xi-min";
        case Status::InsideTangency: return "inside-tangency";
      }
      return "unknown";
    }

    struct Result {
      double tau{0.0};
      Status status{Status::InvalidTable};
      bool ok() const { return status == Status::Ok; }
    };

    static constexpr double kXiMin = 2.0e-4; // smallest mu = sqrt(1-xi^2) resolved
    static constexpr double kMuMin = 2.0e-4;
    static constexpr double kHTop = 50000.0;

    /// f_phase(lambda): the dispersion model of the configuration (README 8.2).  It only
    /// enters through p_crit(lam) = n(r_obs, f(lam)) r_obs; with xi defined as
    /// sin(theta_arrival) the photon's own f cancels exactly, so this is what keeps the
    /// table exact at grazing angles (a 2e-4 error in p moves the tangency by ~1.5 km).
    using DispersionFn = std::function<double(double)>;

    CurvedAbsorptionTable() = default;
    CurvedAbsorptionTable(RefractiveIndexTable const& ri, double R, double rObs,
                          ExtinctionFn alpha, std::vector<double> lambdas,
                          DispersionFn fOfLambda = DispersionFn{}, int nXi = 800);

    bool valid() const { return valid_; }
    std::size_t hNodes() const { return h_.size(); }
    std::size_t xiNodes() const { return xi_.size(); }
    std::size_t lambdaNodes() const { return lam_.size(); }

    // ---- on-disk build cache (wave 5) --------------------------------------
    // 17.9 s and 19.45 MB of deterministic product of (profile, extinction function,
    // dispersion model, R, r_obs, lambda nodes).  Cached for the same reason as the
    // geometry table, and refused on the same terms: the key covers the lambda grid
    // AND sampled values of BOTH callbacks, so a different atmosphere file or a
    // different dispersion model cannot silently reuse this table.
    static std::uint64_t buildKey(double R, double rObs, std::vector<double> const& lambdas,
                                  ExtinctionFn const& alpha, DispersionFn const& fOfLambda);
    bool saveToFile(std::string const& path, std::uint64_t key) const;
    bool loadFromFile(std::string const& path, RefractiveIndexTable const& ri, std::uint64_t key,
                      std::string* why = nullptr);
    double memoryMB() const;
    double xiMin() const { return xi_.empty() ? 0.0 : xi_.front(); }
    double hTop() const { return h_.empty() ? 0.0 : h_.back(); }
    /// Tangency height of the ray xi (hTop if it does not turn below the top).
    double tangencyHeight(double xi) const;
    /// p_crit = n(r_obs) r_obs of the configuration this table was built for.
    double pCrit() const { return pCrit_; }

    /// T(h, xi, lambda): accumulated tau from the ground up to h along the ray xi.
    Result tauAt(double h_asl, double xi, double lambdaNm) const;
    /// tau between two heights of the same ray.
    /// The cumulative is referenced at the TRUE PERIGEE of the column (T(h_t) = 0), so
    ///   pathVisitsPerigee == false : tau = |T(hA) - T(hB)|  (monotone segment)
    ///   pathVisitsPerigee == true  : tau =  T(hA) + T(hB)   (the photon descended to the
    ///                                perigee and came back up: both legs count)
    /// The caller knows the radial velocity at both ends and therefore knows which case
    /// it is; two heights alone cannot tell the table.
    Result tauBetween(double hA, double hB, double xi, double lambdaNm,
                      bool pathVisitsPerigee = false) const;

    /// One rejected lookup, kept for attribution (debugging rule: never let a ray be
    /// rejected without saying which one and why).
    struct Rejection {
      double hA{0.0}, hB{0.0}, xi{0.0}, lambdaNm{0.0};
      Status status{Status::InvalidTable};
    };
    /// The last 32 rejections (process-wide, diagnostics only).
    std::vector<Rejection> const& recentRejections() const { return recent_; }
    std::size_t nLookups() const { return nLookups_; }
    std::size_t nRejectedXi() const { return nRejXi_; }
    std::size_t nRejectedTangency() const { return nRejTan_; }
    void resetCounters() const { nLookups_ = nRejXi_ = nRejTan_ = 0; }
    std::size_t filledNodes() const { return nFilled_; }
    int buildFailCol() const { return buildFailCol_; }

  private:
    double N(double zAbs) const { return ri_ ? ri_->N(zAbs) : 0.0; }
    double nOf(double zAbs, double f) const { return 1.0 + N(zAbs) * f; }
    /// Cell of a monotone grid containing x, clamped (same helper as the other tables).
    static std::size_t gridCell(std::vector<double> const& g, double x);
    static double axisFrac(std::vector<double> const& g, std::size_t i, double x);

    RefractiveIndexTable const* ri_ = nullptr;
    double R_ = 0.0, r_obs_ = 0.0, pCrit_ = 0.0;
    ExtinctionFn alpha_;
    DispersionFn fOfLam_;
    std::vector<double> h_, xi_, lam_;
    std::size_t n_h_{0}, n_xi_{0}, n_lam_{0};
    /// cumulative tau, row-major [lambda][h][xi]
    std::vector<float> T_;
    /// per (lambda, xi): TRUE perigee height of that column [m ASL] (0 when it lies at or
    /// below the ground).  p_crit carries f(lambda), so the perigee is wavelength dependent.
    std::vector<double> hTang_;
    /// per (lambda, xi): first h node at or above the perigee (index into h_)
    std::vector<std::size_t> startIdx_;
    bool valid_ = false;
    std::size_t nFilled_ = 0;
    int buildFailCol_ = -1;
    mutable std::size_t nLookups_ = 0, nRejXi_ = 0, nRejTan_ = 0;
    mutable std::vector<Rejection> recent_;
    void noteRejection(double hA, double hB, double xi, double lam, Status st) const {
      if (recent_.size() >= 32) recent_.erase(recent_.begin());
      recent_.push_back({hA, hB, xi, lam, st});
    }
  };

  // ==========================================================================
  // implementation
  // ==========================================================================

  inline CurvedAbsorptionTable::CurvedAbsorptionTable(RefractiveIndexTable const& ri, double R,
                                                      double rObs, ExtinctionFn alpha,
                                                      std::vector<double> lambdas,
                                                      DispersionFn fOfLambda, int nXi)
      : ri_(&ri), R_(R), r_obs_(rObs), alpha_(std::move(alpha)), fOfLam_(std::move(fOfLambda)) {
    if (!ri.valid() || !(alpha_) || lambdas.empty()) return;
    if (!(rObs >= R) || !(R > 0.0)) return;
    lam_ = std::move(lambdas);
    n_lam_ = lam_.size();

    // ---- h axis: 25-500 m geometric, then 500 m steps to 5 km, then 1 km ----------
    // The step size matters more than it looks: at grazing incidence dT/dh = alpha/mu
    // changes fast with h, so T(h) is convex and a 1 km cell costs 3.8 % on a thin
    // grazing pair (measured: 2500 -> 3000 m at xi = 0.9995).  500 m cells below 5 km
    // bring that under 1 % at the same memory as the coarser axis.
    h_.clear();
    for (double v : {0.0, 25.0, 50.0, 100.0, 250.0}) h_.push_back(v);
    // 100 m cells below 5 km, not 500 m: the FIRST cell above a column's perigee is
    // where T(h) curves hardest, and its height sets the interpolation error for every
    // near-horizon photon (measured at mu = 0.019: 500 m cells gave -7.4 % against a
    // direct quadrature, the quadratic fit -3.3 %; 100 m cells are the fix).
    for (int i = 5; i <= 50; ++i) h_.push_back(100.0 * i);         // 500 m .. 5 km
    for (int i = 6; i * 1000.0 <= kHTop; ++i) h_.push_back(1000.0 * i); // 6 .. 50 km
    n_h_ = h_.size();

    // ---- xi axis: the SAME axis as the geometry table's mu, because xi and mu carry the
    // same information: mu = cos(theta_obs) = sqrt(1-xi^2), and the tau cusp is at
    // mu -> 0 (grazing), where the air mass diverges.  Log-spaced in mu below the break,
    // then uniform in u = sqrt(1-mu).  Using a uniform-in-sqrt(1-xi) axis instead put the
    // acceptance set's rays (mu ~ 0.03, i.e. xi ~ 0.9995) 20-30 % apart per cell and cost
    // 3.3 % on tau.
    xi_.clear();
    n_xi_ = nXi;
    double const muMin = kMuMin;
    double const muB = 0.9;
    double const uB = std::sqrt(std::max(1.0e-12, 1.0 - muB));
    double const sB = std::log(muB / muMin);
    double const K = 2.0 * uB / muB;
    double const sTot = sB + K * uB;
    for (std::size_t j = 0; j < n_xi_; ++j) {
      double const s = sTot * double(j) / double(n_xi_ - 1);
      double mu = (s <= sB) ? muMin * std::exp(s) : (1.0 - std::pow(std::max(0.0, uB - (s - sB) / K), 2));
      mu = std::max(1.0e-9, std::min(1.0, mu));
      xi_.push_back(std::sqrt(std::max(0.0, 1.0 - mu * mu)));
    }
    // ascending in xi, i.e. descending in mu
    std::reverse(xi_.begin(), xi_.end());

    T_.assign(n_lam_ * n_h_ * n_xi_, std::numeric_limits<float>::quiet_NaN());
    hTang_.assign(n_lam_ * n_xi_, 0.0);
    startIdx_.assign(n_lam_ * n_xi_, n_h_);

    // ---- build: one pass per (lambda, xi) column ----------------------------
    // The cumulative is referenced at the TRUE PERIGEE of the column, T(h_t) = 0.  The
    // previous build referenced it at the first NODE where cos psi > 1e-3, which sits
    // 32-428 m ABOVE the real perigee, so the path between the perigee and that node was
    // dropped from every difference and the 500 m cell that brackets the tangency had to
    // represent a sqrt singularity.  A per-photon A/B against the exact chord measured
    // the consequence (1e6 GeV sample, 3000 photons): +537 % tau residual at mu < 0.03,
    // +18 % at 0.03-0.1 and -70 % at 0.1-0.3, with the SIGN flipping (diagnostics
    // README 25).
    //
    // Perigee: A(h) = n(R+h)(R+h) increases monotonically, and the ray exists where
    // A(h) >= p, so h_t solves A(h_t) = p.  When A(0) >= p the perigee is at or below the
    // ground and h_t = 0 (the ground is then the reference, as before).
    std::size_t nFilled = 0;
    for (std::size_t il = 0; il < n_lam_; ++il) {
      double const lam = lam_[il];
      double const f = fOfLam_ ? fOfLam_(lam) : 1.0;
      double const pCrit = nOf(r_obs_, f) * r_obs_;
      for (std::size_t ix = 0; ix < n_xi_; ++ix) {
        double const p = xi_[ix] * pCrit;
        auto Aof = [&](double hh) { double const z = R_ + hh; return nOf(z, f) * z; };
        auto cosPsi = [&](double hh) {
          double const s = p / Aof(hh);
          return (s < 1.0) ? std::sqrt(std::max(0.0, 1.0 - s * s)) : 0.0;
        };
        double h_t = 0.0;
        if (Aof(0.0) < p) {
          double lo = 0.0, hi = kHTop;
          for (int it = 0; it < 80; ++it) {
            double const mid = 0.5 * (lo + hi);
            (Aof(mid) < p) ? (lo = mid) : (hi = mid);
          }
          h_t = 0.5 * (lo + hi);
        }
        hTang_[il * n_xi_ + ix] = h_t;
        std::size_t iStart = 0;
        while (iStart < n_h_ && h_[iStart] < h_t) ++iStart;
        startIdx_[il * n_xi_ + ix] = iStart;
        double tau = 0.0;
        for (std::size_t ih = iStart; ih < n_h_; ++ih) {
          std::size_t const k = il * n_h_ * n_xi_ + ih * n_xi_ + ix;
          double const hHi = h_[ih];
          double const hLo = (ih == iStart) ? h_t : h_[ih - 1];
          if (hHi > hLo) {
            double const dh = hHi - hLo;
            double const cMid = cosPsi(0.5 * (hLo + hHi));
            if (ih == iStart && hLo > 0.0) {
              // First cell: h - h_t = (hHi - h_t) s^2 makes the integrand smooth at the
              // perigee (cos psi ~ s, dh = 2 (hHi-h_t) s ds, so alpha dh / cos psi is
              // finite).  Without it the sqrt is sampled by a uniform midpoint rule and
              // the cell value is meaningless.
              int const nSub = 128;
              double const H = hHi - hLo;
              for (int m = 0; m < nSub; ++m) {
                double const s = (double(m) + 0.5) / double(nSub);
                double const hm = hLo + H * s * s;
                double const c = cosPsi(hm);
                if (!(c > 0.0)) break;
                tau += alpha_(hm, lam) * (2.0 * H * s) / c / double(nSub);
              }
            } else {
              int const nSub = (cMid < 1.0e-1) ? 64 : 4;
              double const ds = dh / double(nSub);
              for (int m = 0; m < nSub; ++m) {
                double const hm = hLo + (double(m) + 0.5) * ds;
                double const c = cosPsi(hm);
                if (!(c > 1.0e-6)) break;
                tau += alpha_(hm, lam) * ds / c;
              }
            }
          }
          T_[k] = static_cast<float>(tau);
          ++nFilled;
        }
      }
      if (nFilled == 0) { buildFailCol_ = int(il); return; }
    }
    nFilled_ = nFilled;
    valid_ = true;
  }

  inline std::size_t CurvedAbsorptionTable::gridCell(std::vector<double> const& g, double x) {
    std::size_t const n = g.size();
    if (n < 2) return 0;
    double const d0 = g[1] - g[0];
    bool const uniform = (n > 2) && (d0 > 0.0) &&
                         (std::abs((g[n - 1] - g[0]) / d0 - double(n - 1)) < 1.0e-6);
    if (uniform) {
      double const t = std::ceil((x - g[0]) / d0) - 1.0;
      if (t <= 0.0) return 0;
      if (t >= double(n - 2)) return n - 2;
      return static_cast<std::size_t>(t);
    }
    std::size_t best = 0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      best = i;
      if (g[i + 1] >= x) break;
    }
    return best;
  }

  inline double CurvedAbsorptionTable::axisFrac(std::vector<double> const& g, std::size_t i,
                                                double x) {
    if (i + 1 >= g.size()) return 0.0;
    double const d = g[i + 1] - g[i];
    if (!(d > 0.0)) return 0.0;
    return std::max(0.0, std::min(1.0, (x - g[i]) / d));
  }

  inline double CurvedAbsorptionTable::memoryMB() const {
    double bytes = 4.0 * double(T_.size());
    bytes += 8.0 * double(h_.size() + xi_.size() + lam_.size() + hTang_.size());
    return bytes / 1048576.0;
  }

  inline double CurvedAbsorptionTable::tangencyHeight(double xi) const {
    if (!valid_) return 0.0;
    std::size_t const ix = gridCell(xi_, std::max(xi, xi_.front()));
    return hTang_[ix];  // first wavelength node; reporting only
  }

  // (h, xi) interpolation at one wavelength node.  u = sqrt(h - h_t) is formed PER
  // COLUMN because each column has its own perigee (p_crit carries f(lambda)); T is
  // proportional to u for constant extinction, so linear interpolation in u removes the
  // sqrt curvature that a 500 m cell cannot resolve.  The two columns are then blended
  // in xi with the same fraction the old bilinear used.
  inline CurvedAbsorptionTable::Result CurvedAbsorptionTable::tauAt(double h_asl, double xi,
                                                                   double lambdaNm) const {
    Result res;
    if (!valid_) return res;
    ++nLookups_;
    if (!(xi >= xi_.front())) {
      ++nRejXi_;
      res.status = Status::BelowXiMin;
      return res;
    }
    double const hh = std::max(0.0, std::min(h_.back(), h_asl));
    double const xv = std::min(1.0, xi);
    std::size_t const ix = gridCell(xi_, xv);
    double const fx = axisFrac(xi_, ix, xv);  // weight of column ix+1

    auto nodeInterp = [&](std::size_t il) -> double {
      double t[2] = {0.0, 0.0};
      for (int b = 0; b < 2; ++b) {
        std::size_t const c = std::min(n_xi_ - 1, ix + std::size_t(b));
        double const ht = hTang_[il * n_xi_ + c];
        if (!(hh > ht)) return std::numeric_limits<double>::quiet_NaN();  // inside tangency
        double const u = std::sqrt(hh - ht);
        std::size_t const i0c = startIdx_[il * n_xi_ + c];
        std::size_t j = i0c;
        while (j + 2 < n_h_ && h_[j + 1] < hh) ++j;
        double const uLo = std::sqrt(std::max(0.0, h_[j] - ht));
        double const uHi = std::sqrt(std::max(0.0, h_[j + 1] - ht));
        double const tLo = double(T_[il * n_h_ * n_xi_ + j * n_xi_ + c]);
        double const tHi = double(T_[il * n_h_ * n_xi_ + (j + 1) * n_xi_ + c]);
        double v = std::numeric_limits<double>::quiet_NaN();
        if (j == i0c && uLo > 0.0 && uHi > uLo) {
          // FIRST cell: T(h_t) = 0 is known EXACTLY, so fit T(u) = a u + b u^2 through
          // (0,0) and the two node values instead of the straight chord.  This is the
          // cell that brackets the tangency, where the sqrt curvature is worst and where
          // the emission of a near-horizon photon usually sits (measured: a linear fit
          // was 7.4 % low there, and it is what pushed |d tau|/tau to 8.8e-2 on the
          // acceptance set).
          double const det = uLo * uHi * (uHi - uLo);
          if (det > 0.0) {
            double const a = (tLo * uHi * uHi - tHi * uLo * uLo) / det;
            double const bb = (tHi * uLo - tLo * uHi) / det;
            /* a*uLo + b*uLo^2 = tLo and a*uHi + b*uHi^2 = tHi  (solved above) */
            v = a * u + bb * u * u;
          }
        }
        if (!std::isfinite(v) && uHi > uLo) {
          double const w = (u - uLo) / (uHi - uLo);
          v = tLo * (1.0 - w) + tHi * w;
        }
        if (!std::isfinite(v)) return std::numeric_limits<double>::quiet_NaN();
        t[b] = v;
      }
      return t[0] * (1.0 - fx) + t[1] * fx;
    };

    if (n_lam_ == 1 || !(lambdaNm > lam_.front())) {
      res.tau = nodeInterp(0);
    } else if (lambdaNm >= lam_.back()) {
      res.tau = nodeInterp(n_lam_ - 1);
    } else {
      std::size_t const il = gridCell(lam_, lambdaNm);
      double const t0 = nodeInterp(il);
      double const t1 = nodeInterp(il + 1);
      double const w = (std::log(lambdaNm) - std::log(lam_[il])) /
                       (std::log(lam_[il + 1]) - std::log(lam_[il]));
      if (t0 > 0.0 && t1 > 0.0) {
        res.tau = std::exp((1.0 - w) * std::log(t0) + w * std::log(t1));
      } else {
        res.tau = (1.0 - w) * t0 + w * t1;
      }
    }
    if (!std::isfinite(res.tau)) {
      ++nRejTan_;
      res.tau = 0.0;
      res.status = Status::InsideTangency;
      return res;
    }
    res.status = Status::Ok;
    return res;
  }

  inline CurvedAbsorptionTable::Result CurvedAbsorptionTable::tauBetween(double hA, double hB,
                                                                        double xi,
                                                                        double lambdaNm,
                                                                        bool visitsPerigee) const {
    Result const rA = tauAt(hA, xi, lambdaNm);
    if (!rA.ok()) {
      noteRejection(hA, hB, xi, lambdaNm, rA.status);
      return rA;
    }
    Result const rB = tauAt(hB, xi, lambdaNm);
    if (!rB.ok()) {
      noteRejection(hA, hB, xi, lambdaNm, rB.status);
      return rB;
    }
    Result out;
    // T is referenced at the perigee (T(h_t) = 0): a path that descends to the perigee
    // and climbs back accumulates BOTH legs, a monotone segment only the difference.
    out.tau = visitsPerigee ? (rA.tau + rB.tau) : std::fabs(rA.tau - rB.tau);
    out.status = Status::Ok;
    return out;
  }

  // ==========================================================================
  // wave 5: on-disk build cache
  // ==========================================================================
  inline std::uint64_t CurvedAbsorptionTable::buildKey(double R, double rObs,
                                                       std::vector<double> const& lambdas,
                                                       ExtinctionFn const& alpha,
                                                       DispersionFn const& fOfLambda) {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&h](double v) {
      std::uint64_t b = 0;
      std::memcpy(&b, &v, sizeof b);
      h = (h ^ b) * 1099511628211ull;
    };
    mix(R);
    mix(rObs);
    mix(double(lambdas.size()));
    for (double l : lambdas) mix(l);
    // Fingerprint the two callbacks by SAMPLING them: that is what the table is a
    // function of, and it is independent of which file or functor produced them.
    if (alpha) {
      for (int i = 0; i <= 6; ++i) {
        double const hh = 50000.0 * double(i) / 6.0;
        for (int j = 0; j <= 6; ++j) {
          double const lam = 300.0 + 300.0 * double(j) / 6.0;
          mix(alpha(hh, lam));
        }
      }
    }
    if (fOfLambda) {
      for (int j = 0; j <= 7; ++j) mix(fOfLambda(300.0 + 320.0 * double(j) / 7.0));
    }
    return h;
  }

  inline bool CurvedAbsorptionTable::saveToFile(std::string const& path,
                                                std::uint64_t key) const {
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
    char magic[8] = {'C', 'R', 'V', 'A', 'B', 'S', '2', '\0'};
    std::uint32_t const ver = 2u;
    put(magic, sizeof magic);
    put(&ver, sizeof ver);
    put(&key, sizeof key);
    putD(R_); putD(r_obs_); putD(pCrit_);
    putU(std::uint64_t(n_h_)); putU(std::uint64_t(n_xi_)); putU(std::uint64_t(n_lam_));
    putU(std::uint64_t(nFilled_));
    putU(std::uint64_t(h_.size())); putU(std::uint64_t(xi_.size())); putU(std::uint64_t(lam_.size()));
    put(h_.data(), h_.size() * sizeof(double));
    put(xi_.data(), xi_.size() * sizeof(double));
    put(lam_.data(), lam_.size() * sizeof(double));
    putU(std::uint64_t(T_.size()));
    put(T_.data(), T_.size() * sizeof(float));
    putU(std::uint64_t(hTang_.size()));
    put(hTang_.data(), hTang_.size() * sizeof(double));
    putU(std::uint64_t(startIdx_.size()));
    for (std::size_t v : startIdx_) {
      std::uint64_t const u = std::uint64_t(v);
      put(&u, sizeof u);
    }
    f.write(reinterpret_cast<char const*>(&sum), sizeof sum);
    f.flush();
    if (!f) { std::remove(tmp.c_str()); return false; }
    f.close();
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
      std::remove(tmp.c_str());
      return false;
    }
    return true;
  }

  inline bool CurvedAbsorptionTable::loadFromFile(std::string const& path,
                                                  RefractiveIndexTable const& ri,
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
    if (std::memcmp(magic, "CRVABS2", 7) != 0) return fail("wrong magic (older cache format?)");
    if (!get(&ver, sizeof ver)) return fail("truncated version");
    if (ver != 2u) return fail("version mismatch");
    if (!get(&fileKey, sizeof fileKey)) return fail("truncated key");
    if (fileKey != key) return fail("build key mismatch (atmosphere file or dispersion changed)");

    CurvedAbsorptionTable t;
    auto getD = [&](double& v) { return get(&v, sizeof v); };
    auto getU = [&](std::uint64_t& v) { return get(&v, sizeof v); };
    if (!getD(t.R_) || !getD(t.r_obs_) || !getD(t.pCrit_)) return fail("truncated scalars");
    std::uint64_t nh = 0, nx = 0, nl = 0, nf = 0, sh = 0, sx = 0, sl = 0;
    if (!getU(nh) || !getU(nx) || !getU(nl) || !getU(nf) || !getU(sh) || !getU(sx) || !getU(sl))
      return fail("truncated sizes");
    t.n_h_ = std::size_t(nh);
    t.n_xi_ = std::size_t(nx);
    t.n_lam_ = std::size_t(nl);
    t.nFilled_ = std::size_t(nf);
    if (t.n_h_ < 2 || t.n_xi_ < 2 || t.n_lam_ < 1) return fail("inconsistent axis sizes");
    if (sh != nh || sx != nx || sl != nl) return fail("inconsistent axis storage");
    t.h_.resize(std::size_t(sh));
    t.xi_.resize(std::size_t(sx));
    t.lam_.resize(std::size_t(sl));
    if (!get(t.h_.data(), t.h_.size() * sizeof(double))) return fail("truncated h axis");
    if (!get(t.xi_.data(), t.xi_.size() * sizeof(double))) return fail("truncated xi axis");
    if (!get(t.lam_.data(), t.lam_.size() * sizeof(double))) return fail("truncated lambda axis");
    std::uint64_t nT = 0;
    if (!getU(nT)) return fail("truncated T size");
    if (nT != nh * nx * nl) return fail("inconsistent T size");
    t.T_.resize(std::size_t(nT));
    if (!get(t.T_.data(), t.T_.size() * sizeof(float))) return fail("truncated T payload");
    std::uint64_t nHT = 0;
    if (!getU(nHT)) return fail("truncated hTang size");
    if (nHT != nl * nx) return fail("inconsistent hTang size");
    t.hTang_.resize(std::size_t(nHT));
    if (!get(t.hTang_.data(), t.hTang_.size() * sizeof(double))) return fail("truncated hTang");
    std::uint64_t nSI = 0;
    if (!getU(nSI)) return fail("truncated startIdx size");
    if (nSI != nHT) return fail("inconsistent startIdx size");
    t.startIdx_.resize(std::size_t(nSI));
    for (std::size_t i = 0; i < t.startIdx_.size(); ++i) {
      std::uint64_t u = 0;
      if (!get(&u, sizeof u)) return fail("truncated startIdx payload");
      t.startIdx_[i] = std::size_t(u);
    }
    std::uint64_t fileSum = 0;
    if (!f.read(reinterpret_cast<char*>(&fileSum), sizeof fileSum))
      return fail("truncated checksum");
    if (fileSum != sum) return fail("checksum mismatch (corrupt or partially written file)");

    t.ri_ = &ri;
    t.valid_ = true;
    t.buildFailCol_ = -1;
    *this = std::move(t);
    return true;
  }

} // namespace corsika::cherenkov
