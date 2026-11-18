#pragma once

/*
 * Dopri5.hpp
 *
 * Allocation-free Dormand-Prince 5(4) stepper with FSAL, quartic dense output and a
 * bounded step controller.  Header-only, no framework dependencies.
 *
 * Why: the tracers' fixed-step RK4 is 4th order and its step size is a heuristic
 * (H_local/20, capped).  A 5(4) pair costs 7 evaluations per step but carries a fifth
 * order solution plus an embedded fourth order one, so every step comes with a LOCAL
 * ERROR ESTIMATE.  On the long, smoothly varying paths of a horizontal shower that is
 * expected to buy many fewer steps for the same accuracy; whether it does is a
 * measurement, not an assumption (diagnostics/phase1_geometry.md, Phase 3).
 *
 * Deliberate properties
 * --------------------
 * * Fixed-size storage: the state has N components and the seven stage derivatives are
 *   `double k[7][N]` inside the object.  Nothing is allocated inside the stepping loop.
 * * FSAL: the seventh stage IS the derivative at the new state, so an ACCEPTED step
 *   hands k[6] to the next step's first stage.  invalidate() must be called whenever the
 *   state, the step or the model changes discontinuously (a crossing solve, a rejected
 *   step with a changed model, a kappa/model switch) -- the caller owns that contract.
 * * The dense output is the standard quartic interpolant (Hairer's P matrix, the one
 *   scipy's RK45 uses).  It is what LOCATES a surface crossing inside a large step: it is
 *   smooth and monotone there, so the root find is well posed.  It is NOT accurate enough
 *   to BE the crossing state: measured against a converged invariant quadrature, the
 *   interpolant of a 15-unit step is off by +2.08 mm at t = 0.25 while both step endpoints
 *   are inside 0.02 mm (diagnostics/phase4_plan.md, "Step 3, task 1").  Locate on the
 *   interpolant, then re-land the step on the surface and read the STATE -- re-landing
 *   drops that residual to 0.018 mm.
 * * Error control per component, with a floor: the position components are scaled by an
 *   absolute tolerance in METRES (never by the Earth-centred radius, which would let a
 *   1e-9 relative tolerance hide 6 mm of impact error), the direction by its own scale,
 *   and the integral components (optical path, optical depth) by their own.
 * * Controller: h_new = h * clip(0.2, 5, 0.9 * e^(-1/5)), with e == 0 handled
 *   explicitly (grow by the 5x cap, never divide by zero).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace corsika::cherenkov {

  /// One try's outcome.
  enum class Dopri5Result {
    Accepted,     ///< error <= 1, the state was advanced
    Rejected,     ///< error > 1, h was reduced and the state is untouched
    MinStepFailed ///< rejected at h <= hMin: the caller must fail the ray
  };

  template <int N>
  class Dopri5 {
  public:
    /// atol/rtol per component.  `rtol` may be zero for a purely absolute component.
    void setTolerances(double const (&atol)[N], double const (&rtol)[N]) {
      for (int i = 0; i < N; ++i) {
        atol_[i] = atol[i];
        rtol_[i] = rtol[i];
      }
    }

    /// Seeds the state and invalidates FSAL.
    void setState(double const (&y)[N]) {
      for (int i = 0; i < N; ++i) y_[i] = y[i];
      fsal_ = false;
    }

    double const (&state() const)[N] { return y_; }

    /// Invalidates FSAL reuse (call after any state/model discontinuity).
    void invalidate() { fsal_ = false; }

    void setMinStep(double const m) { hMin_ = (m > 0.0) ? m : 0.0; }
    double minStep() const { return hMin_; }

    /// Weighted max error of the last attempt, and the component that produced it.
    double lastError() const { return lastErr_; }
    int worstComponent() const { return worstComp_; }

    /// The UNWEIGHTED local error estimate of the last attempt, per component.  Exposed
    /// because a caller may need to weight it by something the scalar cannot express --
    /// e.g. the projection of a position error onto a target plane, which is the whole
    /// point of the impact-weighted acceptance (diagnostics/phase4_plan.md, step 3).
    double const (&lastErrorVector() const)[N] { return dErr_; }

    /// Size of the last ACCEPTED step (the one the dense output belongs to).
    double lastStepSize() const { return hAccepted_; }

    // ---- counters (per instance; the caller aggregates them) --------------------
    std::size_t nAccepted{0};
    std::size_t nRejected{0};

    /**
     * One attempted step of size `h`.  On Accepted the state holds the 5th-order
     * solution, k_[0] is the derivative there (FSAL armed) and the dense interpolation
     * of THIS step is available at dense().  On Rejected the state is untouched and
     * `hSuggested` holds the reduced step.  Never allocates.
     */
    template <typename Rhs>
    Dopri5Result tryStep(Rhs const& f, double const h, double& hSuggested) {
      hSuggested = h;
      if (!fsal_) {
        f(y_, k_[0]);
        ++nF;
      }

      // --- stages 2..6 (stage 1 is k_[0], from FSAL or from the evaluation above) ---
      // The ray ODE is AUTONOMOUS (kappa depends on position only), so the c-nodes of
      // the tableau do not enter the stage state: yt = y + h * sum_s a[s] k[s].
      // Passing c*h here as well used to multiply the increment by the node twice,
      // which made the scheme low order while keeping its error estimate small.
      stage(f, h, k2a_, 1);
      stage(f, h, k3a_, 2);
      stage(f, h, k4a_, 3);
      stage(f, h, k5a_, 4);
      stage(f, h, k6a_, 5);

      // --- 5th-order solution (also the 7th stage's state) -------------------------
      double yNew[N];
      for (int i = 0; i < N; ++i) {
        double acc = 0.0;
        for (int j = 0; j < 6; ++j) acc += b5_[j] * k_[j][i];
        yNew[i] = y_[i] + h * acc;
      }
      f(yNew, k_[6]); // FSAL derivative at the new state
      ++nF;

      // --- embedded 4th-order error: err = h * sum (b5 - b4) * k -------------------
      double err = 0.0;
      int worst = 0;
      for (int i = 0; i < N; ++i) {
        double d = 0.0;
        for (int j = 0; j < 7; ++j) d += e5_[j] * k_[j][i];
        d *= h;
        dErr_[i] = d;
        double const scale =
            atol_[i] + rtol_[i] * std::max(std::fabs(y_[i]), std::fabs(yNew[i]));
        double const e = (scale > 0.0) ? std::fabs(d) / scale : 0.0;
        if (e > err) {
          err = e;
          worst = i;
        }
      }
      lastErr_ = err;
      worstComp_ = worst;

      // --- accept / reject --------------------------------------------------------
      if (err <= 1.0) {
        // snapshot for the dense output BEFORE the state moves
        for (int i = 0; i < N; ++i) {
          yStart_[i] = y_[i];
          for (int j = 0; j < 7; ++j) kAcc_[j][i] = k_[j][i];
        }
        hAccepted_ = h;
        for (int i = 0; i < N; ++i) y_[i] = yNew[i];
        for (int i = 0; i < N; ++i) k_[0][i] = k_[6][i]; // FSAL hand-off
        fsal_ = true;
        ++nAccepted;
        hMinSeen = std::min(hMinSeen, h);
        hMaxSeen = std::max(hMaxSeen, h);
        hSum += h;
        hSuggested = std::min(h * growth(err, true), 5.0 * h);
        return Dopri5Result::Accepted;
      }

      hSuggested = h * std::max(0.2, growth(err, false));
      if (!(hSuggested > hMin_)) {
        ++nMinStepFailed;
        return Dopri5Result::MinStepFailed;
      }
      ++nRejected;
      return Dopri5Result::Rejected;
    }

    /**
     * Quartic dense output of the LAST ACCEPTED step at `t` in [0,1] (t = 1 is the new
     * state).  Valid until the next accepted step.
     */
    void dense(double const t, double (&out)[N]) const {
      double const th = std::min(1.0, std::max(0.0, t));
      double const th2 = th * th;
      double const th3 = th2 * th;
      double const th4 = th3 * th;
      double p[7];
      for (int j = 0; j < 7; ++j)
        p[j] = denseP_[j][0] * th + denseP_[j][1] * th2 + denseP_[j][2] * th3 +
               denseP_[j][3] * th4;
      for (int i = 0; i < N; ++i) {
        double acc = yStart_[i];
        for (int j = 0; j < 7; ++j) acc += hAccepted_ * p[j] * kAcc_[j][i];
        out[i] = acc;
      }
    }

    // (remaining public counters; the first two are declared above the stepping code)
    std::size_t nMinStepFailed{0};
    std::size_t nF{0}; ///< RHS evaluations (an FSAL reuse is not counted again)
    double hMinSeen{1.0e300};
    double hMaxSeen{0.0};
    double hSum{0.0};

  private:
    // ---- one stage: yt = y + h * sum_s a[s] * k[s], then k[j] = f(yt) -------------
    template <typename Rhs, int M>
    void stage(Rhs const& f, double const h, double const (&a)[M], int const j) {
      double yt[N];
      for (int i = 0; i < N; ++i) {
        double acc = 0.0;
        for (int s = 0; s < M; ++s) acc += a[s] * k_[s][i];
        yt[i] = y_[i] + h * acc;
      }
      f(yt, k_[j]);
      ++nF;
    }

    /// Controller gain 0.9 * e^(-1/5); e == 0 is handled explicitly (never divided).
    static double growth(double const err, bool const accepted) {
      if (!(err > 0.0)) return accepted ? 5.0 : 0.2;
      double const g = 0.9 * std::pow(err, -0.2);
      return accepted ? g : std::max(g, 0.2);
    }

    // ---- Dormand-Prince (4)5 tableau -------------------------------------------
    inline static constexpr double k2a_[1] = {1.0 / 5.0};
    inline static constexpr double k3a_[2] = {3.0 / 40.0, 9.0 / 40.0};
    inline static constexpr double k4a_[3] = {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0};
    inline static constexpr double k5a_[4] = {19372.0 / 6561.0, -25360.0 / 2187.0,
                                              64448.0 / 6561.0, -212.0 / 729.0};
    inline static constexpr double k6a_[5] = {9017.0 / 3168.0, -355.0 / 33.0,
                                              46732.0 / 5247.0, 49.0 / 176.0,
                                              -5103.0 / 18656.0};
    inline static constexpr double b5_[6] = {35.0 / 384.0, 0.0, 500.0 / 1113.0,
                                             125.0 / 192.0, -2187.0 / 6784.0,
                                             11.0 / 84.0};
    inline static constexpr double e5_[7] = {35.0 / 384.0 - 5179.0 / 57600.0, 0.0,
                                             500.0 / 1113.0 - 7571.0 / 16695.0,
                                             125.0 / 192.0 - 393.0 / 640.0,
                                             -2187.0 / 6784.0 + 92097.0 / 339200.0,
                                             11.0 / 84.0 - 187.0 / 2100.0, -1.0 / 40.0};

    /// Hairer's dense-output P matrix (the one scipy's RK45 uses): [stage][power-1].
    inline static constexpr double denseP_[7][4] = {
        {1.0, -8048581381.0 / 2820520608.0, 8663915743.0 / 2820520608.0,
         -12715105075.0 / 11282082432.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0,
         87487479700.0 / 32700410799.0},
        {0.0, -1754552775.0 / 470086768.0, 14199869525.0 / 1410260304.0,
         -10690763975.0 / 1880347072.0},
        {0.0, 127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0,
         701980252875.0 / 199316789632.0},
        {0.0, -282668133.0 / 205662961.0, 2019193451.0 / 616988883.0,
         -1453857185.0 / 822651844.0},
        {0.0, 40617522.0 / 29380423.0, -110615467.0 / 29380423.0,
         69997945.0 / 29380423.0}};

    double y_[N]{};
    double k_[7][N]{};
    double dErr_[N]{};    ///< unweighted local error estimate of the last attempt
    double kAcc_[7][N]{}; ///< stages of the last accepted step (dense output)
    double yStart_[N]{};  ///< state at the start of the last accepted step
    double hAccepted_{0.0};
    double atol_[N]{};
    double rtol_[N]{};
    double hMin_{0.0};
    double lastErr_{0.0};
    int worstComp_{0};
    bool fsal_{false};
  };

} // namespace corsika::cherenkov
