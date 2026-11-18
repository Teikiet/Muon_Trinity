#pragma once
// ============================================================================
// ref_tracer_b.hpp -- STANDALONE reference copy of production's curved-ray
// tracer (diagnostic only; NOT part of the build).
//
// Source of truth: corsika/modules/cherenkov/Cherenkov.inl
//   CherenkovProcess::rayTraceSphericalToPlane   lines 358-719
//   CherenkovProcess::computeLocalScaleHeight    lines 199-226
//
// This file preserves expression order, constants, branch conditions and
// floating-point types. Apparent defects are preserved deliberately; they are
// NOT fixed here. Every mechanical adaptation is labelled [ADAPTED].
//
// DEPENDENCY MANIFEST
//   refractiveIndexFunc_(h,lam) -> (n_phase,n_group)  LIFTED as injected
//       std::function (production Cherenkov.hpp member, set by the caller).
//       Signature and call sites unchanged.
//   extinctionCoeffFunc_(h,lam)                       LIFTED as injected
//       std::function. Returns 1/m. Same two call sites as production.
//   earthRadiusM_                                     LIFTED as a member,
//       value supplied by the caller (production member default).
//   CORSIKA_LOGGER_*                                  REMOVED (logging only;
//       no numerical effect).
//   std::mt19937 / uniform_real_distribution draw      [ADAPTED] replaced by
//       an explicit survivalThreshold argument with the same semantics.
//       Pass 0.0 to reproduce the "no early termination" branch exactly.
//   accumulated_tau entry value                        [ADAPTED] explicit
//       tau_start argument (production relies on caller-initialised ref).
// ============================================================================

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>

namespace refdiag {

// ============================================================================
// [2e-c AUDIT] Explicit termination outcome.
// Exactly one value is assigned before EVERY return of trace(), so an
// unpopulated field can never be mistaken for a physical result.
// ============================================================================
enum class Termination {
  None = 0,             // trace() never reached a return site (cannot occur)
  NoRefractiveIndex,    // Cherenkov.inl:369-372  no refractiveIndexFunc_
  Escape,               // Cherenkov.inl:482-488  h_cur > max_altitude
  AbsorptionMain,       // Cherenkov.inl:653-656  exp(-tau) < threshold, main step
  AbsorptionNearPlane,  // Cherenkov.inl:597-602  same test, straight segment
  NearPlane,            // Cherenkov.inl:604      hit via the straight segment
  NearPlaneParallel,    // Cherenkov.inl:558-561  |den| < 1e-10
  NearPlaneBehind,      // Cherenkov.inl:563-566  t_hit < 0
  Crossing,             // Cherenkov.inl:679-702  crossing interpolation
  Diverging,            // Cherenkov.inl:706-714  ray diverging from the plane
  StepLimit             // Cherenkov.inl:717-718  max_steps reached
};

inline char const* terminationName(Termination t) {
  switch (t) {
    case Termination::None:
      return "none";
    case Termination::NoRefractiveIndex:
      return "no_refractive_index";
    case Termination::Escape:
      return "escape";
    case Termination::AbsorptionMain:
      return "absorption_main";
    case Termination::AbsorptionNearPlane:
      return "absorption_near_plane";
    case Termination::NearPlane:
      return "near_plane";
    case Termination::NearPlaneParallel:
      return "near_plane_parallel";
    case Termination::NearPlaneBehind:
      return "near_plane_behind";
    case Termination::Crossing:
      return "crossing";
    case Termination::Diverging:
      return "diverging";
    case Termination::StepLimit:
      return "step_limit";
  }
  return "unknown";
}

// Legacy 2e-b numeric code, now DERIVED from the termination outcome. In 2e-b it
// was written by hand at a few return sites only, so the near-plane, near-parallel
// and near-behind returns left it at its default 0. -1 = not assigned.
inline int branchCodeOf(Termination t) {
  switch (t) {
    case Termination::Crossing:
      return 1;
    case Termination::NearPlane:
      return 2;
    case Termination::Diverging:
      return 3;
    case Termination::StepLimit:
      return 4;
    case Termination::Escape:
      return 5;
    case Termination::AbsorptionMain:
    case Termination::AbsorptionNearPlane:
      return 6;
    case Termination::NearPlaneParallel:
      return 7;
    case Termination::NearPlaneBehind:
      return 8;
    case Termination::NoRefractiveIndex:
      return 9;
    case Termination::None:
      return -1;
  }
  return -1;
}

struct RefTracer {
  // production: std::function<std::pair<double,double>(double,double)>
  //             -> (phase index, group index)
  using RefractiveIndexFn =
      std::function<std::pair<double, double>(double z_m, double lambda_nm)>;
  // production: std::function<double(double,double)> -> extinction, 1/m
  using ExtinctionFn = std::function<double(double h_m, double lambda_nm)>;

  double earthRadiusM_ = 0.0;
  RefractiveIndexFn refractiveIndexFunc_;
  ExtinctionFn extinctionCoeffFunc_;
  bool enableAbsorption_ = true;

  struct Result {
    bool hit = false;
    double x_impact = 0.0;
    double y_impact = 0.0;
    double z_impact = 0.0;
    double optical_path = 0.0;    // metres (group-index weighted)
    double geometric_path = 0.0;  // metres
    double accumulated_tau = 0.0; // dimensionless
    int steps = 0;                // iterations consumed in the loop

    // --- [2e-c AUDIT OBSERVATION ONLY] instrumentation -------------------------
    // These fields are written but never read by any numerical expression, so
    // they cannot alter the arithmetic or the evaluation order.
    //
    // COUNTER SEMANTICS (explicit):
    //   steps          loop iterations consumed (includes turning-point retries
    //                  and the iteration in which the trace terminates).
    //   main_advances  COMPLETED main-loop position advances only: incremented
    //                  after x_plane/z_plane/r_cur have been updated to the new
    //                  position (Cherenkov.inl:659-663). A loop entry, a
    //                  turning-point revert/retry, and the near-plane straight
    //                  segment do NOT count as an advance.
    //   turning        turning-point events (Cherenkov.inl:497-514).
    //   tau_main       net tau retained in the RETURNED accumulated_tau that was
    //                  added by main-loop steps: accumulated_tau at the last
    //                  main-loop tau checkpoint minus tau_start. Production does
    //                  NOT restore accumulated_tau at a turning point, so such a
    //                  step's tau stays inside this value; nothing is rolled back.
    //   tau_near       tau added by the near-plane straight segment, i.e.
    //                  accumulated_tau - (tau_main + tau_start). Meaningful only
    //                  when near_applicable; otherwise it is exactly 0 because no
    //                  near-plane substep ran.
    //   near_seg       near-plane segment length t_hit [m]; valid only when
    //                  near_applicable.
    //   near_sub       near-plane tau substep count; valid only when
    //                  near_applicable.
    //   cross_*        crossing-branch quantities; valid only when
    //                  cross_applicable.
    double tau_main = 0.0;
    double tau_near = 0.0;
    int main_advances = 0;
    int turning = 0;          // turning-point events
    int branch = -1;          // legacy 2e-b code, derived from termination
    Termination termination = Termination::None;
    bool near_applicable = false;  // near_seg / near_sub / tau_near meaningful
    double near_seg = 0.0;
    int near_sub = 0;
    bool cross_applicable = false;
    double cross_alpha = 0.0; // crossing fraction alpha
    double cross_ds = 0.0;    // final full-step ds at the crossing
    double cross_kext = 0.0;  // extinction used on that final step

  };

  // --- Cherenkov.inl:199-226, lifted verbatim (arithmetic unchanged) -------
  double computeLocalScaleHeight(double z) const {
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

    // Guard against N2 being too close to N1 (where log(N2/N1) -> 0)
    double const ratio = N2 / N1;
    double const log_ratio = std::log(ratio);

    // If magnitude of logarithm is too small, the gradient is too flat
    if (std::abs(log_ratio) < 1e-6) {
      return 8000.0;  // fallback when gradient is effectively zero
    }

    // H = -dz / ln(N2/N1)
    return -dz / log_ratio;
  }

  Result trace(double x0, double y0, double z0_altitude,
               double dx, double dy, double dz,
               double plane_x, double plane_y, double plane_z,
               double norm_x, double norm_y, double norm_z,
               double lambdaNm,
               double survivalThreshold_in = 0.0,
               double tau_start = 0.0) const;
};

// ============================================================================
// Definition -- Cherenkov.inl:358-719 transcribed statement for statement.
// [ADAPTED] out-parameters -> Result fields; logging removed; RNG draw ->
// survivalThreshold_in; accumulated_tau reference -> local seeded by tau_start.
// Nothing else is reordered.
// ============================================================================
inline RefTracer::Result RefTracer::trace(double x0, double y0,
                                          double z0_altitude, double dx,
                                          double dy, double dz, double plane_x,
                                          double plane_y, double plane_z,
                                          double norm_x, double norm_y,
                                          double norm_z, double lambdaNm,
                                          double survivalThreshold_in,
                                          double tau_start) const {
  Result res;

  // --- Cherenkov.inl:369-372 ---
  if (!refractiveIndexFunc_) {
    res.termination = Termination::NoRefractiveIndex;  // [2e-c obs]
    res.branch = branchCodeOf(res.termination);        // [2e-c obs]
    return res;  // production returns false; res.hit stays false
  }

  double const R_E = earthRadiusM_;

  // --- Cherenkov.inl:381-383 [ADAPTED] threshold supplied by caller ---
  double const survivalThreshold =
      (enableAbsorption_ && extinctionCoeffFunc_) ? survivalThreshold_in : 0.0;

  // --- Cherenkov.inl:385-406 ---
  double const px = x0;
  double const py = y0;
  double const pz = z0_altitude;
  double const r0 = std::sqrt(px * px + py * py + pz * pz);

  double const r_hat_x = px / r0;
  double const r_hat_y = py / r0;
  double const r_hat_z = pz / r0;

  double const cos_theta0 = dx * r_hat_x + dy * r_hat_y + dz * r_hat_z;
  double const sin_theta0 =
      std::sqrt(std::max(0.0, 1.0 - cos_theta0 * cos_theta0));

  double const h0 = r0 - R_E;
  double const n0 = refractiveIndexFunc_(h0, lambdaNm).first;  // phase

  double const ray_inv = n0 * r0 * sin_theta0;

  // --- Cherenkov.inl:410-412 ---
  auto distanceToPlane = [&](double x, double y, double z) -> double {
    return (x - plane_x) * norm_x + (y - plane_y) * norm_y +
           (z - plane_z) * norm_z;
  };

  // --- Cherenkov.inl:415-445 ---
  double tang_x = dx - cos_theta0 * r_hat_x;
  double tang_y = dy - cos_theta0 * r_hat_y;
  double tang_z = dz - cos_theta0 * r_hat_z;
  double const tang_mag =
      std::sqrt(tang_x * tang_x + tang_y * tang_y + tang_z * tang_z);

  double t_hat_x, t_hat_y, t_hat_z;
  if (tang_mag > 1e-12) {
    t_hat_x = tang_x / tang_mag;
    t_hat_y = tang_y / tang_mag;
    t_hat_z = tang_z / tang_mag;
  } else {
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

    double const t_norm =
        std::sqrt(t_hat_x * t_hat_x + t_hat_y * t_hat_y + t_hat_z * t_hat_z);
    t_hat_x /= t_norm;
    t_hat_y /= t_norm;
    t_hat_z /= t_norm;
  }

  // --- Cherenkov.inl:449 ---
  double const az_sign = (t_hat_y >= 0.0) ? 1.0 : -1.0;

  // --- Cherenkov.inl:452-476 ---
  double x_plane = 0.0;
  double z_plane = r0;
  double r_cur = r0;

  double x_cur = x0;
  double y_cur = y0;
  double z_cur = z0_altitude;

  int direction = (cos_theta0 > 0.0) ? 1 : -1;

  double optical_path = 0.0;
  double geometric_path = 0.0;
  double accumulated_tau = tau_start;  // [ADAPTED] seeded by caller

  int const max_steps = 10000;
  double const max_altitude = 100000.0;

  double r_prev = r_cur;
  double x_prev = x_cur;
  double y_prev = y_cur;
  double z_prev = z_cur;
  double optical_prev = 0.0;
  double geometric_prev = 0.0;
  double dist_prev = distanceToPlane(x_cur, y_cur, z_cur);

  int step = 0;

  // [2e-c OBSERVATION ONLY] hoisted observations; written but never read by any
  // numerical expression, so arithmetic and evaluation order are unchanged.
  double obs_tau_before_step = accumulated_tau;
  double obs_last_ds = 0.0;
  double obs_last_kext = 0.0;
  Termination obs_term = Termination::None;  // set at every return site
  int obs_main_advances = 0;                 // completed main-loop advances
  double obs_tau_main = tau_start;           // last main-loop tau checkpoint
  bool obs_near_applicable = false;
  double obs_near_seg = 0.0;
  int obs_near_sub = 0;
  bool obs_cross_applicable = false;

  // [ADAPTED] production returns from many points; state is published before
  // every return so partial accumulations stay visible, as with the original
  // reference out-parameters. publish() also derives the audit bookkeeping from
  // the state it publishes, so the counters always describe the RETURNED result
  // rather than merely attempted additions.
  auto publish = [&]() {
    res.x_impact = x_cur;
    res.y_impact = y_cur;
    res.z_impact = z_cur;
    res.optical_path = optical_path;
    res.geometric_path = geometric_path;
    res.accumulated_tau = accumulated_tau;
    res.steps = step;
    // --- [2e-c obs] derived bookkeeping -------------------------------------
    res.termination = obs_term;
    res.branch = branchCodeOf(obs_term);
    res.main_advances = obs_main_advances;
    res.tau_main = obs_tau_main - tau_start;
    res.tau_near = accumulated_tau - obs_tau_main;
    res.near_applicable = obs_near_applicable;
    res.near_seg = obs_near_seg;
    res.near_sub = obs_near_sub;
    res.cross_applicable = obs_cross_applicable;
  };


  for (; step < max_steps; ++step) {
    double const h_cur = r_cur - R_E;

    // --- Cherenkov.inl:482-488: escaped atmosphere ---
    if (h_cur > max_altitude) {
      // [ADAPTED] the acos/deg expression in production only feeds the log
      // message; it does not affect state and is therefore omitted.
      res.branch = 5;  // [2e-b obs] escaped atmosphere
      obs_term = Termination::Escape;  // [2e-c obs]
      publish();
      return res;  // production: return false
    }

    // --- Cherenkov.inl:491-494 ---
    auto const [n_phase, n_group] = refractiveIndexFunc_(h_cur, lambdaNm);
    double sin_theta_cur = ray_inv / (n_phase * r_cur);

    // --- Cherenkov.inl:497-514: turning point ---
    if (sin_theta_cur > 1.0) {
      ++res.turning;  // [2e-b obs] turning-point branch taken
      r_cur = r_prev;
      x_plane = x_prev;
      z_plane = z_prev;
      x_cur = x_prev;
      y_cur = y_prev;
      z_cur = z_prev;
      optical_path = optical_prev;
      geometric_path = geometric_prev;
      // NOTE: accumulated_tau is NOT restored here. Production behaviour,
      // preserved deliberately.
      direction = -direction;
      continue;
    }

    // --- Cherenkov.inl:516-517 ---
    double const cos_theta_cur =
        std::sqrt(std::max(0.0, 1.0 - sin_theta_cur * sin_theta_cur));
    double const tan_theta_cur =
        (cos_theta_cur > 1e-12) ? sin_theta_cur / cos_theta_cur : 1e12;

    // --- Cherenkov.inl:520-527: signed distance to plane ---
    double const cos_alpha_d = z_plane / r_cur;
    double const sin_alpha_d = x_plane / r_cur;
    double const cur3d_x =
        r_cur * (cos_alpha_d * r_hat_x + sin_alpha_d * t_hat_x);
    double const cur3d_y =
        r_cur * (cos_alpha_d * r_hat_y + sin_alpha_d * t_hat_y);
    double const cur3d_z =
        r_cur * (cos_alpha_d * r_hat_z + sin_alpha_d * t_hat_z);
    double const dist_to_plane = (cur3d_x - plane_x) * norm_x +
                                 (cur3d_y - plane_y) * norm_y +
                                 (cur3d_z - plane_z) * norm_z;
    double const abs_dist_to_plane = std::abs(dist_to_plane);

    // --- Cherenkov.inl:530: straight-line segment within 1000 m ---
    if (abs_dist_to_plane <= 1000.0 && direction == -1) {
      double const cos_alpha_c = z_plane / r_cur;
      double const sin_alpha_c = x_plane / r_cur;

      double const r3d_x = cos_alpha_c * r_hat_x + sin_alpha_c * t_hat_x;
      double const r3d_y = cos_alpha_c * r_hat_y + sin_alpha_c * t_hat_y;
      double const r3d_z = cos_alpha_c * r_hat_z + sin_alpha_c * t_hat_z;

      double const t3d_x = -sin_alpha_c * r_hat_x + cos_alpha_c * t_hat_x;
      double const t3d_y = -sin_alpha_c * r_hat_y + cos_alpha_c * t_hat_y;
      double const t3d_z = -sin_alpha_c * r_hat_z + cos_alpha_c * t_hat_z;

      double const dir3d_x = -cos_theta_cur * r3d_x + sin_theta_cur * t3d_x;
      double const dir3d_y = -cos_theta_cur * r3d_y + sin_theta_cur * t3d_y;
      double const dir3d_z = -cos_theta_cur * r3d_z + sin_theta_cur * t3d_z;

      double const pos3d_x =
          r_cur * (cos_alpha_c * r_hat_x + sin_alpha_c * t_hat_x);
      double const pos3d_y =
          r_cur * (cos_alpha_c * r_hat_y + sin_alpha_c * t_hat_y);
      double const pos3d_z =
          r_cur * (cos_alpha_c * r_hat_z + sin_alpha_c * t_hat_z);

      // --- Cherenkov.inl:555-571 ---
      double const num = (plane_x - pos3d_x) * norm_x +
                         (plane_y - pos3d_y) * norm_y +
                         (plane_z - pos3d_z) * norm_z;
      double const den = dir3d_x * norm_x + dir3d_y * norm_y + dir3d_z * norm_z;
      if (std::abs(den) < 1e-10) {
        obs_term = Termination::NearPlaneParallel;  // [2e-c obs]
        publish();
        return res;  // ray parallel to plane
      }
      double const t_hit = num / den;
      if (t_hit < 0.0) {
        obs_term = Termination::NearPlaneBehind;  // [2e-c obs]
        publish();
        return res;  // intersection behind current position
      }

      // [ADAPTED] locals instead of writing the out-params early; identical
      // arithmetic, same use in h_impact.
      double const xi = pos3d_x + t_hit * dir3d_x;
      double const yi = pos3d_y + t_hit * dir3d_y;
      double const zi = pos3d_z + t_hit * dir3d_z;

      // --- Cherenkov.inl:574-580 ---
      double const h_cur_plane = r_cur - R_E;
      double const h_impact = std::sqrt(xi * xi + yi * yi + zi * zi) - R_E;
      double const h_mid = 0.5 * (h_cur_plane + h_impact);
      auto const [n_phase_mid, n_group_mid] =
          refractiveIndexFunc_(h_mid, lambdaNm);
      (void)n_phase_mid;  // unused in production as well; kept for parity
      double const seg_len = t_hit;
      geometric_path += seg_len;
      optical_path += n_group_mid * seg_len;

      // --- Cherenkov.inl:582-595 ---
      int const n_tau_steps = std::max(1, static_cast<int>(seg_len / 50.0));
      double const ds_tau = seg_len / static_cast<double>(n_tau_steps);
      // [2e-c obs] the segment and its substep count exist from here on
      obs_near_applicable = true;
      obs_near_seg = seg_len;
      obs_near_sub = n_tau_steps;
      for (int ts = 0; ts < n_tau_steps; ++ts) {
        double const s_mid = (static_cast<double>(ts) + 0.5) * ds_tau;
        double const tau_x = pos3d_x + s_mid * dir3d_x;
        double const tau_y = pos3d_y + s_mid * dir3d_y;
        double const tau_z = pos3d_z + s_mid * dir3d_z;
        double const tau_r =
            std::sqrt(tau_x * tau_x + tau_y * tau_y + tau_z * tau_z);
        double const tau_h = tau_r - R_E;
        if (enableAbsorption_ && extinctionCoeffFunc_) {
          accumulated_tau += extinctionCoeffFunc_(tau_h, lambdaNm) * ds_tau;
        }
      }

      // --- Cherenkov.inl:597-602 ---
      if (enableAbsorption_ && extinctionCoeffFunc_) {
        if (std::exp(-accumulated_tau) < survivalThreshold) {
          obs_term = Termination::AbsorptionNearPlane;  // [2e-c obs]
          publish();
          return res;  // photon already dead
        }
      }

      // --- Cherenkov.inl:604 ---
      obs_term = Termination::NearPlane;  // [2e-c obs]
      publish();
      res.hit = true;
      res.x_impact = xi;
      res.y_impact = yi;
      res.z_impact = zi;
      return res;
    }

    // --- Cherenkov.inl:608-615 ---
    double const H_scale = computeLocalScaleHeight(h_cur);
    double const dx_step_base = H_scale / 20.0;
    double dx;
    if (sin_theta_cur > 0.99) {
      dx = H_scale / 100.0;
    } else {
      dx = std::min(dx_step_base, 1000.0);
    }

    // --- Cherenkov.inl:617-624 ---
    r_prev = r_cur;
    x_prev = x_cur;
    y_prev = y_cur;
    z_prev = z_cur;
    optical_prev = optical_path;
    geometric_prev = geometric_path;
    double const dist_at_step_start = dist_prev;

    // --- Cherenkov.inl:627-639 ---
    double const cos_alpha = z_plane / r_cur;
    double const sin_alpha = x_plane / r_cur;

    double dx_sign = -1.0;
    double dz_sign = -1.0;
    if (static_cast<double>(direction) > 0) {
      dz_sign = 1.0;
    }
    double const dr = std::abs(dx * cos_theta_cur /
                               (sin_theta_cur * cos_alpha +
                                cos_theta_cur * sin_alpha));
    double const dx_step = dx_sign * dx;
    double const dz_step =
        dz_sign * std::abs(dr * (cos_alpha - tan_theta_cur * sin_alpha));

    // --- Cherenkov.inl:642-657 ---
    double const ds = std::sqrt(dx_step * dx_step + dz_step * dz_step);
    geometric_path += ds;
    optical_path += n_group * ds;

    if (enableAbsorption_ && extinctionCoeffFunc_) {
      // NOTE: r_cur is still the PREVIOUS radius here (it is updated at
      // Cherenkov.inl:663), so r_cur == r_prev and h_mid collapses to the
      // step's START altitude, not its midpoint. Preserved deliberately.
      double const h_mid = 0.5 * (r_prev - R_E + r_cur - R_E);
      double const kext = extinctionCoeffFunc_(h_mid, lambdaNm);
      obs_tau_before_step = accumulated_tau;  // [2e-b obs]
      obs_last_ds = ds;                       // [2e-b obs]
      obs_last_kext = kext;                   // [2e-b obs]
      accumulated_tau += kext * ds;
      obs_tau_main = accumulated_tau;         // [2e-c obs] main-loop checkpoint
      if (std::exp(-accumulated_tau) < survivalThreshold) {
        obs_term = Termination::AbsorptionMain;  // [2e-c obs]
        publish();
        return res;  // photon already dead
      }
    }

    // --- Cherenkov.inl:661-674 ---
    x_plane += az_sign * dx_step;
    z_plane += dz_step;
    r_cur = std::sqrt(x_plane * x_plane + z_plane * z_plane);
    ++obs_main_advances;  // [2e-c obs] a COMPLETED position advance

    double const cos_alpha_step = z_plane / r_cur;
    double const sin_alpha_step = x_plane / r_cur;
    double const gx =
        r_cur * (cos_alpha_step * r_hat_x + sin_alpha_step * t_hat_x);
    double const gy =
        r_cur * (cos_alpha_step * r_hat_y + sin_alpha_step * t_hat_y);
    double const gz =
        r_cur * (cos_alpha_step * r_hat_z + sin_alpha_step * t_hat_z);

    x_cur = gx;
    y_cur = gy;
    z_cur = gz;

    // --- Cherenkov.inl:677-703: plane crossing ---
    double const dist_cur = distanceToPlane(x_cur, y_cur, z_cur);

    if (dist_at_step_start * dist_cur <= 0.0) {
      double alpha;
      if (std::abs(dist_at_step_start - dist_cur) > 1e-12) {
        alpha = dist_at_step_start / (dist_at_step_start - dist_cur);
      } else {
        alpha = 0.5;
      }

      alpha = std::max(0.0, std::min(1.0, alpha));

      double const xi = x_prev + alpha * (x_cur - x_prev);
      double const yi = y_prev + alpha * (y_cur - y_prev);
      double const zi = z_prev + alpha * (z_cur - z_prev);

      // NOTE: optical_path and geometric_path are interpolated to the crossing
      // point; accumulated_tau is NOT interpolated and retains the FULL step's
      // contribution. Preserved deliberately.
      optical_path = optical_prev + alpha * (optical_path - optical_prev);
      geometric_path =
          geometric_prev + alpha * (geometric_path - geometric_prev);

      obs_cross_applicable = true;         // [2e-c obs]
      res.cross_alpha = alpha;             // [2e-c obs]
      res.cross_ds = obs_last_ds;          // [2e-c obs]
      res.cross_kext = obs_last_kext;      // [2e-c obs]
      // accumulated_tau keeps the FULL step here (production does not
      // interpolate it), so it is also the main-loop tau checkpoint.
      obs_tau_main = accumulated_tau;      // [2e-c obs]
      obs_term = Termination::Crossing;    // [2e-c obs]

      publish();
      res.hit = true;
      res.x_impact = xi;
      res.y_impact = yi;
      res.z_impact = zi;
      return res;
    }
    dist_prev = dist_cur;

    // --- Cherenkov.inl:706-714 ---
    if (direction < 0 && r_cur > r_prev + 1e-6) {
      if (dist_prev * dist_cur > 0.0 &&
          std::abs(dist_cur) > std::abs(dist_prev)) {
        res.branch = 3;  // [2e-b obs] diverging
        obs_term = Termination::Diverging;  // [2e-c obs]
        publish();
        return res;  // ray diverging from plane
      }
    }
  }

  obs_term = Termination::StepLimit;  // [2e-c obs] loop exit
  publish();
  return res;
}

}  // namespace refdiag
