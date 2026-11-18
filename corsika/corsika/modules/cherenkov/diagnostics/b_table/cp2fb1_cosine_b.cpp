// ============================================================================
// cp2fb1_cosine_b.cpp -- Checkpoint 2f-b1: isolate the failed vacuum-cosine
// check of T2.  DIAGNOSTIC ONLY.  Separate from cp2fb_domain_b.cpp, which is
// preserved unchanged.
//
//   ./cp2fb1_cosine_b
//
// Reproduces T2's exact double-precision chain for the vacuum case (f = 0), then
// separates (i) evaluation error of the double formulas from (ii) inconsistent
// rounded representations of the intended ray, using 100-decimal-digit
// arithmetic (boost::multiprecision::cpp_dec_float_100 >= 80 digits).
// ============================================================================
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

using mp = boost::multiprecision::cpp_dec_float_100;
static constexpr double R_E = 6371000.0;

struct Row {
  double zObs, h, muObs;
  bool isWorst;
};

struct Res {
  double p{}, rObs{}, rE{}, muE{}, L{}, sObs{}, rs{}, lineT2{}, invT2{}, factT2{}, resid{};
  int kAtMax{-1};
  double maxResid{0.0};
};

// T2's vacuum chain, verbatim in structure (Cherenkov-free; n = 1, f = 0)
static Res chain(double zObs, double hE, double muObs) {
  Res r;
  r.rObs = R_E + zObs;
  r.rE = r.rObs + hE;
  r.p = r.rObs * std::sqrt(std::max(0.0, 1.0 - muObs * muObs));
  double const pr = r.p / r.rE;
  r.muE = std::sqrt(std::max(0.0, 1.0 - pr * pr));
  r.sObs = std::sqrt(std::max(0.0, (r.rObs - r.p) * (r.rObs + r.p)));
  r.L = hE * (r.rE + r.rObs) / (r.rE * r.muE + r.sObs);
  // T2 samples s = L*k/4 and takes the max |mu_line - mu_inv|
  for (int k = 0; k <= 4; ++k) {
    double const s = r.L * double(k) / 4.0;
    double const rs = std::hypot(s - r.rE * r.muE, r.p);
    double const muLine = (r.rE * r.muE - s) / rs;
    double const q = r.p / rs;
    double const muInv = std::sqrt(std::max(0.0, 1.0 - q * q));
    double const d = std::fabs(muLine - muInv);
    if (d > r.maxResid) {
      r.maxResid = d;
      r.kAtMax = k;
    }
    if (k == 4) {
      r.rs = rs;
      r.lineT2 = muLine;
      r.invT2 = muInv;
      r.factT2 = std::sqrt(std::max(0.0, (rs - r.p) * (rs + r.p))) / rs;
      r.resid = d;
    }
  }
  return r;
}

static double ulpOf(double x) {
  return std::nextafter(x, std::numeric_limits<double>::infinity()) - x;
}

int main() {
  Row const rows[6] = {{0.0, 50000.0, 1.0, false},      {0.0, 50000.0, 0.5, false},
                       {0.0, 50000.0, 1.0e-3, false},   {0.0, 50000.0, 2.0e-4, false},
                       {0.0, 50000.0, 1.0e-6, false},   {0.0, 1.0, 1.0e-6, true}};
  std::printf("=== Checkpoint 2f-b1: isolate the failed vacuum-cosine check ===\n");
  std::printf("arithmetic: boost::multiprecision::cpp_dec_float_100 (100 decimal digits)\n");
  std::printf("case: vacuum (n = 1, f = 0), z_obs = 0, h = 50000 m (plus the T2 worst row)\n\n");

  // ---------------- 2. sanity table --------------------------------------
  std::printf("--- 2. sanity table (before any explanation) ---\n");
  std::printf("  %-10s %-24s %-24s %-13s %-13s %-13s\n", "mu_obs", "g_obs - p [m]",
              "g_obs mu^2/2 [m]", "ulp(p) [m]", "observed resid", "ulp(p)/(g_obs mu)");
  for (int i = 0; i < 6; ++i) {
    Res const r = chain(rows[i].zObs, rows[i].h, rows[i].muObs);
    double const gObs = (1.0 + 0.0 * 0.0) * r.rObs;  // f = 0 -> g = r
    double const gap = gObs - r.p;
    double const est = gObs * rows[i].muObs * rows[i].muObs / 2.0;
    double const sen = ulpOf(r.p) / (gObs * rows[i].muObs);
    std::printf("  %-10.17g %-24.17g %-24.17g %-13.6e %-13.6e %-13.6e%s\n", rows[i].muObs, gap,
                est, ulpOf(r.p), r.maxResid, sen, rows[i].isWorst ? "   <- T2 worst row" : "");
  }
  std::printf("\n");

  // ---------------- 1b. the worst row's stored data -----------------------
  Res const w = chain(0.0, 1.0, 1.0e-6);
  std::printf("--- 1b. worst row recovered by recomputation (not inferred) ---\n");
  std::printf("  z_obs = 0 m, h = 1 m, mu_obs = 1e-6 ; max residual at sample k = %d (s = L)\n",
              w.kAtMax);
  std::printf("  stored: p = %.17g m ; g_obs = %.17g m ; rE = %.17g m\n", w.p, w.rObs, w.rE);
  std::printf("          mu_E = %.17g ; L = %.17g m ; r(s=L) = %.17g m\n", w.muE, w.L, w.rs);
  std::printf("  direction at emission: (sin(theta_E), 0, -mu_E) = (%.17g, 0, %.17g)\n",
              std::sqrt(std::max(0.0, 1.0 - w.muE * w.muE)), -w.muE);
  std::printf("  mu_line  = (rE mu_E - L)/r(s)      = %.17g\n", w.lineT2);
  std::printf("  mu_inv   = sqrt(1 - (p/r(s))^2)    = %.17g   <- T2's invariant form\n", w.invT2);
  std::printf("  mu_fact  = sqrt((r-p)(r+p))/r      = %.17g   <- factored form\n", w.factT2);
  std::printf("  |mu_line - mu_inv| = %.17g ; matches the T2 gate residual 4.274971e-11: %s\n",
              w.resid, (std::fabs(w.resid - 4.274971e-11) < 1e-15) ? "YES" : "NO");
  std::printf("  DEPENDENCE: mu_E = sqrt(1-(p/rE)^2) and L both use the SAME stored p,\n");
  std::printf("              so mu_line also descends from that p. Neither quantity is an\n");
  std::printf("              independent measurement of the intended ray.\n\n");

  // ---------------- 3. multiprecision separation --------------------------
  std::printf("--- 3. multiprecision separation (stored doubles imported exactly) ---\n");
  std::printf("  A = intended launch cosine ; B = factored invariant from the STORED p\n");
  std::printf("  %-10s %-14s %-14s %-14s %-14s %-14s %-14s %-12s\n", "mu_obs", "B - A",
              "invT2 - B", "fact - B", "line - lineMP", "closure", "|inv-line|", "k_max");
  double wInvB = 0.0, wFactB = 0.0, wLine = 0.0, wClos = 0.0, wBA = 0.0, wResid = 0.0;
  for (int i = 0; i < 6; ++i) {
    Res const r = chain(rows[i].zObs, rows[i].h, rows[i].muObs);
    mp const P(r.p), RE(r.rE), RR(r.rs), GO(r.rObs);
    mp const Brs = sqrt((RR - P) * (RR + P)) / RR;        // same radius as T2's sample
    mp const invMP = sqrt(mp(1) - (P / RR) * (P / RR));   // MP of T2's unfactored form
    mp const factMP = sqrt((RR - P) * (RR + P)) / RR;     // MP of the factored form
    mp const lineMP = (RE * mp(r.muE) - mp(r.L)) / RR;    // MP of T2's line formula
    mp const A = mp(rows[i].muObs);                       // intended launch parameter
    double const dBA = std::fabs((Brs - A).convert_to<double>());
    double const dInvB = std::fabs((mp(r.invT2) - Brs).convert_to<double>());
    double const dFactB = std::fabs((mp(r.factT2) - Brs).convert_to<double>());
    double const dLine = std::fabs((mp(r.lineT2) - lineMP).convert_to<double>());
    double const dBL = std::fabs((Brs - lineMP).convert_to<double>());
    double const closure = r.resid - (dInvB + dBL + dLine);
    std::printf("  %-10.17g %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-14.3e %-12d\n",
                rows[i].muObs, dBA, dInvB, dFactB, dLine, closure, r.resid, r.kAtMax);
    wBA = std::max(wBA, dBA);
    wInvB = std::max(wInvB, dInvB);
    wFactB = std::max(wFactB, dFactB);
    wLine = std::max(wLine, dLine);
    wClos = std::max(wClos, closure);
    wResid = std::max(wResid, r.resid);
    (void)invMP; (void)factMP;
  }
  std::printf("\n");

  // ---------------- 4. gates and classification ---------------------------
  std::printf("--- 4. gates (the brief's tolerances, applied as stated) ---\n");
  auto gate = [](char const* what, double v, double tol) {
    bool const ok = (v <= tol);
    std::printf("  %-58s value=%.6e  tol=%.1e  %s\n", what, v, tol, ok ? "PASS" : "FAIL");
    return ok;
  };
  bool const g1 = gate("invT2 (T2's unfactored double form) vs B", wInvB, 1.0e-12);
  bool const g2 = gate("lineT2 (double) vs its multiprecision result", wLine, 1.0e-12);
  bool const g3 = gate("fact (factored double form) vs B  [demonstration]", wFactB, 1.0e-12);
  bool const g4 = gate("unexplained remainder of the original discrepancy", wClos, 1.0e-12);
  std::printf("\n  input rounding |B - A| (max, not a gate): %.6e\n", wBA);
  std::printf("  original T2 cosine residual (max): %.6e\n\n", wResid);

  std::printf("--- classification ---\n");
  std::printf("  (1) numerical evaluation error : %s\n",
              (!g1 && g3 && g2) ? "PRESENT and OWNS the failed T2 check" : "see the gate rows");
  std::printf("      - the unfactored double form sqrt(1-(p/r)^2) exceeds 1e-12 (%.3e);\n", wInvB);
  std::printf("        the factored double form sqrt((r-p)(r+p))/r stays at %.3e and the\n", wFactB);
  std::printf("        line form at %.3e.\n", wLine);
  std::printf("  (2) inconsistent rounded representations : %s\n",
              (wBA > 1.0e-12) ? "PRESENT (|B - A| > 1e-12)" : "not above 1e-12");
  std::printf("      - the stored p implies a cosine %.3e away from the intended launch value;\n",
              wBA);
  std::printf("        both T2 quantities descend from that same p, so T2 could not see it.\n");
  std::printf("  (4) unresolved defect : %s\n",
              (g4 && g3) ? "NONE within 1e-12" : "REMAINDER ABOVE 1e-12 -- investigate");
  std::printf("=== end of 2f-b1 ===\n");
  return (g3 && g2 && g4) ? 0 : 1;
}
