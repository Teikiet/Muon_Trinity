// cp2fh_window_identity_b.cpp -- SYNTHETIC sanity check of the +-25 m window identity
// Diagnostic ONLY. No production profile, no boundary semantics, no table.
//   alpha_window(h) = (A(h+25) - A(h-25)) / 50
// with A piecewise-linear (constant slope per segment).
#include <algorithm>
#include <cmath>
#include <cstdio>

int main() {
  double const sL = 1.0e-4;   // slope below the knot   [m^-1]
  double const sR = 3.0e-4;   // slope above the knot   [m^-1]
  double const hk = 5000.0;   // interior knot [m]
  // Cumulative optical depth A(h) = int_0^h alpha, piecewise linear:
  //   h <= 1000        : A = 0
  //   1000 < h <= 5000 : slope sL
  //   h > 5000         : slope sR
  auto A = [&](double h) -> double {
    if (h <= 1000.0) return 0.0;
    if (h <= hk) return sL * (h - 1000.0);
    return sL * (hk - 1000.0) + sR * (h - hk);
  };
  auto aw = [&](double h) -> double { return (A(h + 25.0) - A(h - 25.0)) / 50.0; };

  double const gate = 1.0e-14 * std::max(std::fabs(sL), std::fabs(sR)) + 1.0e-16;
  int bad = 0;

  struct Case {
    double h, expect;
    char const* what;
  } cases[] = {
      {3000.0, sL, "constant-slope segment"},
      {hk, 0.5 * (sL + sR), "interior knot -> (sL+sR)/2"},
  };
  std::printf("gate = 1e-14*max(|sL|,|sR|) + 1e-16 = %.3e m^-1\n", gate);
  for (Case const& c : cases) {
    double const got = aw(c.h);
    double const d = std::fabs(got - c.expect);
    bool const ok = (d <= gate);
    if (!ok) ++bad;
    std::printf("  %-28s h=%.1f expect=%.17g got=%.17g |d|=%.3e  %s\n", c.what, c.h, c.expect,
                got, d, ok ? "PASS" : "FAIL");
  }
  // Point derivative at the knot: one-sided limits sL and sR; not generally equal
  // to the window average, and not defined as a point value.
  double const winAtKnot = aw(hk);
  bool const derivEqualsWindow = (sL == winAtKnot) || (sR == winAtKnot);
  std::printf("  point derivative at knot: left=%.17g right=%.17g ; window=%.17g ; "
              "one-sided derivative equals window? %s\n",
              sL, sR, winAtKnot, derivEqualsWindow ? "YES" : "NO");
  // Guard: with the production max(0, .) floor, a negative window average is
  // clipped; here slopes are positive so it is inactive (reported, not asserted).
  std::printf("  floor max(0,·) active in these cases? %s\n",
              (winAtKnot < 0.0 || aw(3000.0) < 0.0) ? "YES" : "NO");
  std::printf("RESULT: %s\n", bad == 0 ? "SYNTHETIC IDENTITIES HOLD" : "SYNTHETIC IDENTITIES FAIL");
  return bad == 0 ? 0 : 1;
}
