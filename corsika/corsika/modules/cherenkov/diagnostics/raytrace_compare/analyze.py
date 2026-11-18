#!/usr/bin/env python3
"""analyze.py -- summarise the six-configuration ray-tracing comparison."""
import glob
import os
import re
import statistics as st

OUT = "/scratch/general/vast/u1520754/chercmp"
NAMES = {"A": "straight, no abs", "B": "straight, abs on",
         "C": "curved legacy, no abs", "D": "curved legacy, abs on",
         "E": "curved + table, no abs", "F": "straight, no abs, legacy timing"}
SEEDS = [600, 601, 602, 603, 604]

def parse(path):
    txt = open(path, errors="replace").read()
    g = lambda pat: (re.search(pat, txt).group(1) if re.search(pat, txt) else None)
    prof = re.search(r"Cherenkov profile: total=([\d.]+) s, production=([\d.]+) s, "
                     r"bunchloop=([\d.]+) s, timing=([\d.]+) s, recorder=([\d.]+) s, "
                     r"bunches=(\d+), recorded=(\d+)", txt)
    return {
        "hits": int(g(r"Total Cherenkov hits recorded: (\d+)") or 0),
        "run": float(g(r"Run time: ([\d.]+) s") or 0.0),
        "total": float(prof.group(1)) if prof else 0.0,
        "timing": float(prof.group(4)) if prof else 0.0,
        "bunchloop": float(prof.group(3)) if prof else 0.0,
        "bunches": int(prof.group(6)) if prof else 0,
        "abs": "on" if re.search(r"absorption=on", txt) else "off",
        "fast": "on" if re.search(r"Cherenkov fast timing: on", txt) else "off",
        "curv": "on" if re.search(r"enable-curvature|curved", txt) else "off",
    }

print(f"{'cfg':<3} {'geometry/absorption':<28} {'hits mean+-SE':<20} {'timing_s mean':<14} "
      f"{'run_s mean':<11} {'abs':<4} {'fast':<5} {'n'}")
for name, label in NAMES.items():
    rows = []
    for s in SEEDS:
        p = f"{OUT}/m_{name}_s{s}.log"
        if os.path.exists(p):
            r = parse(p)
            if r["hits"] is not None:
                rows.append(r)
    if not rows:
        print(f"{name:<3} {label:<28} (no runs yet)")
        continue
    hits = [r["hits"] for r in rows]
    se = st.stdev(hits) / len(hits) ** 0.5 if len(hits) > 1 else 0.0
    print(f"{name:<3} {label:<28} {st.mean(hits):8.1f} +- {se:6.1f}   "
          f"{st.mean([r['timing'] for r in rows]):8.3f}      "
          f"{st.mean([r['run'] for r in rows]):7.2f}     {rows[0]['abs']:<4} "
          f"{rows[0]['fast']:<5} {len(rows)}")

# yield ratios vs the straight-line no-absorption baseline
base = [parse(f"{OUT}/m_A_s{s}.log")["hits"] for s in SEEDS
        if os.path.exists(f"{OUT}/m_A_s{s}.log")]
if base:
    mb = st.mean(base)
    print(f"\nyield ratios relative to A (straight line, no absorption), A = {mb:.1f} hits:")
    for name in "BCDEF":
        rows = [parse(f"{OUT}/m_{name}_s{s}.log")["hits"] for s in SEEDS
                if os.path.exists(f"{OUT}/m_{name}_s{s}.log")]
        if rows:
            print(f"  {name}: {st.mean(rows)/mb:6.4f}   ({st.mean(rows):.1f} hits)")
