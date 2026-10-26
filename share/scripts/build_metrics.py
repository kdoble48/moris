#!/usr/bin/env python3
"""Snapshot MORIS build metrics for before/after comparison of build-system changes.

Policy: every build-system change is justified by a measured build-time or size
reduction. This script produces the measurement.

Usage:
    python3 build_metrics.py <build_dir> [--label baseline] [--out <dir>] [--before <log>] [--wall-s <s>]
    python3 build_metrics.py --diff A.json B.json

Metrics captured:
  - .ninja_log: total edge CPU time, link-vs-compile split, top-20 most
    expensive edges (latest entry per output wins). By default the whole tree:
    every entry in the log. To measure one build (an incremental rebuild),
    copy the log first and pass the copy as --before; only the entries ninja
    wrote after it count:
        cp build_opt/.ninja_log /tmp/before.log; ninja -C build_opt moris
        python3 build_metrics.py build_opt --label incr --before /tmp/before.log
  - --wall-s: the build's wall time as measured by the caller (ninja's log has
    no invocation start/end, and its edge times restart when ninja re-runs cmake)
  - sizes: every executable (bin/, *.exe, projects/**), every static lib (*.a),
    every shared object (lib/*.so), and the whole build tree
Writes one JSON per invocation: <out>/<label>_<builddir>_<date>.json and prints
a human summary. --diff compares edge times only between snapshots of the same
scope (whole_tree or since_before); sizes are always compared.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import date
from pathlib import Path

LINK_SUFFIXES = (".exe", ".so", ".a")


def log_entries(log: Path):
    """The entry lines of a .ninja_log (start_ms, end_ms, mtime, output, hash), header and blanks dropped."""
    if not log.exists():
        return []
    return [line for line in log.read_text(errors="replace").splitlines() if line and not line.startswith("#")]


def parse_ninja_log(build_dir: Path, before: Path | None = None):
    log = build_dir / ".ninja_log"
    if not log.exists():
        return None
    entries = log_entries(log)
    if before is not None:
        # Entries already in the copy are older builds. Ninja appends one line per finished edge and, when it
        # recompacts the log, rewrites the surviving entries unchanged, so the set difference is exactly what was
        # written after the copy: robust to log order, a cmake re-run inside the build, and no-op builds (empty).
        seen = set(log_entries(before))
        entries = [line for line in entries if line not in seen]
    latest = {}  # output -> (start_ms, end_ms)
    for line in entries:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            start, end = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        latest[parts[3]] = (start, end)

    edges = []
    for out, (start, end) in latest.items():
        dur = (end - start) / 1000.0
        is_link = out.endswith(LINK_SUFFIXES) or (
            "/" in out and not out.endswith((".o", ".stamp", ".cmake", ".h", ".hpp", ".cpp"))
            and Path(out).suffix == ""
        )
        kind = "link" if is_link else ("compile" if out.endswith(".o") else "other")
        edges.append({"output": out, "seconds": round(dur, 1), "kind": kind})

    total = sum(e["seconds"] for e in edges)
    by_kind = {}
    for e in edges:
        k = by_kind.setdefault(e["kind"], {"count": 0, "seconds": 0.0})
        k["count"] += 1
        k["seconds"] = round(k["seconds"] + e["seconds"], 1)
    top = sorted(edges, key=lambda e: -e["seconds"])[:20]
    return {
        "scope": "whole_tree" if before is None else "since_before",
        "edges": len(edges),
        "total_edge_cpu_min": round(total / 60.0, 1),
        "by_kind": by_kind,
        "link_share_pct": round(100.0 * by_kind.get("link", {}).get("seconds", 0) / total, 1) if total else 0.0,
        "top20": top,
    }


def collect_sizes(build_dir: Path):
    exes, libs_a, libs_so = [], [], []
    for p in build_dir.rglob("*"):
        if not p.is_file() or p.is_symlink():
            continue
        n, sz = str(p.relative_to(build_dir)), p.stat().st_size
        if p.suffix == ".a":
            libs_a.append((n, sz))
        elif ".so" in p.name and "CMakeFiles" not in n:
            libs_so.append((n, sz))
        elif (p.suffix == ".exe" or (os.access(p, os.X_OK) and p.suffix == "" and "CMakeFiles" not in n
                                     and ("/bin/" in n or n.startswith("bin/") or "mains" in n))):
            exes.append((n, sz))

    def summarize(items):
        items = sorted(items, key=lambda x: -x[1])
        return {
            "count": len(items),
            "total_mb": round(sum(s for _, s in items) / 1e6, 1),
            "top10": [{"file": n, "mb": round(s / 1e6, 1)} for n, s in items[:10]],
        }

    du = subprocess.run(["du", "-sb", str(build_dir)], capture_output=True, text=True)
    tree_gb = round(int(du.stdout.split()[0]) / 1e9, 2) if du.returncode == 0 else None
    return {"executables": summarize(exes), "static_libs": summarize(libs_a),
            "shared_libs": summarize(libs_so), "tree_gb": tree_gb}


def diff(a_path: str, b_path: str):
    a, b = json.loads(Path(a_path).read_text()), json.loads(Path(b_path).read_text())
    print(f"delta ({Path(b_path).name} - {Path(a_path).name}):")
    an, bn = a.get("ninja") or {}, b.get("ninja") or {}
    if an and bn:
        a_scope, b_scope = an.get("scope", "whole_tree"), bn.get("scope", "whole_tree")  # pre-scope JSONs: whole tree
        if a_scope != b_scope:
            print(f"  edge times NOT compared: scope {a_scope} vs {b_scope} (one build vs the whole tree)")
        else:
            print(f"  edge CPU:  {an['total_edge_cpu_min']} -> {bn['total_edge_cpu_min']} min "
                  f"({bn['total_edge_cpu_min'] - an['total_edge_cpu_min']:+.1f}), scope {a_scope}")
            print(f"  link share: {an['link_share_pct']}% -> {bn['link_share_pct']}%")
    aw, bw = a.get("wall_s"), b.get("wall_s")
    if aw is not None and bw is not None:
        print(f"  wall:      {aw} -> {bw} s ({bw - aw:+.1f})")
    for key in ("executables", "static_libs", "shared_libs"):
        at, bt = a["sizes"][key]["total_mb"], b["sizes"][key]["total_mb"]
        print(f"  {key}: {at} -> {bt} MB ({bt - at:+.1f}), "
              f"count {a['sizes'][key]['count']} -> {b['sizes'][key]['count']}")
    print(f"  tree: {a['sizes']['tree_gb']} -> {b['sizes']['tree_gb']} GB")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("build_dir", nargs="?", type=Path)
    p.add_argument("--label", default="snapshot")
    p.add_argument("--out", type=Path, help="output dir (default $MORIS_RUNS_DIR/benchmarks/build_metrics)")
    p.add_argument("--before", type=Path, help="copy of .ninja_log taken before the build to measure")
    p.add_argument("--wall-s", type=float, help="the build's wall time in seconds, measured by the caller")
    p.add_argument("--diff", nargs=2, metavar=("A.json", "B.json"))
    a = p.parse_args()
    if a.diff:
        diff(*a.diff)
        return
    if a.build_dir is None:
        p.error("build_dir is required (or --diff A.json B.json)")
    build_dir = a.build_dir.resolve()
    if not build_dir.is_dir():
        p.error(f"not a directory: {build_dir}")
    if a.before is not None and not a.before.is_file():
        p.error(f"--before: no such file: {a.before}")
    out_dir = a.out or Path(os.environ.get("MORIS_RUNS_DIR", ".")) / "benchmarks" / "build_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    snap = {"label": a.label, "build_dir": str(build_dir), "date": date.today().isoformat(),
            "ninja": parse_ninja_log(build_dir, a.before), "sizes": collect_sizes(build_dir)}
    if a.wall_s is not None:
        snap["wall_s"] = round(a.wall_s, 1)
    out = out_dir / f"{a.label}_{build_dir.name}_{snap['date']}.json"
    out.write_text(json.dumps(snap, indent=2))

    n, s = snap["ninja"], snap["sizes"]
    if n:
        wall = f", wall {a.wall_s:.1f} s" if a.wall_s is not None else ""
        print(f"[{a.label}] {build_dir.name} ({n['scope'].replace('_', ' ')}): {n['edges']} edges, "
              f"{n['total_edge_cpu_min']} CPU-min{wall}, link share {n['link_share_pct']}%")
    print(f"  exes: {s['executables']['count']} = {s['executables']['total_mb']} MB | "
          f".a: {s['static_libs']['count']} = {s['static_libs']['total_mb']} MB | "
          f".so: {s['shared_libs']['count']} = {s['shared_libs']['total_mb']} MB | "
          f"tree {s['tree_gb']} GB")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
