#!/usr/bin/env python3
"""Snapshot MORIS build metrics for before/after comparison of build-system changes.

Policy: every build-system change is justified by a measured build-time or size
reduction. This script produces the measurement.

Usage:
    python3 build_metrics.py <build_dir> [--label baseline] [--out <dir>] [--all-sessions]

Metrics captured:
  - .ninja_log: total edge CPU time, wall span, link-vs-compile split, top-20
    most expensive edges. By default only the LAST ninja session that did work
    (the incremental build just run); --all-sessions keeps every session, latest
    entry per output (the whole tree's cost, as before).
  - sizes: every executable (bin/, *.exe, projects/**), every static lib (*.a),
    every shared object (lib/*.so), and the whole build tree
Writes one JSON per invocation: <out>/<label>_<builddir>_<date>.json and prints
a human summary. Compare two snapshots with:  python3 build_metrics.py --diff A.json B.json
"""

import json
import os
import subprocess
import sys
from datetime import date
from pathlib import Path

LINK_SUFFIXES = (".exe", ".so", ".a")


def read_sessions(log: Path):
    """Split .ninja_log into ninja sessions: lists of (output, start_ms, end_ms).

    Times are relative to the start of each ninja invocation, and ninja appends an entry when an edge finishes,
    so end times never decrease within a session; a decrease starts a new one. A no-op invocation writes nothing.
    After ninja recompacts the log, the rewritten old entries are unordered and split into spurious short
    sessions; the entries appended afterwards still form the last session.
    """
    sessions, prev_end = [[]], -1
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        try:
            start, end = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        if end < prev_end:
            sessions.append([])
        prev_end = end
        sessions[-1].append((parts[3], start, end))
    return [s for s in sessions if s]


def parse_ninja_log(build_dir: Path, all_sessions: bool = False):
    log = build_dir / ".ninja_log"
    if not log.exists():
        return None
    sessions = read_sessions(log)
    rows = [r for s in sessions for r in s] if all_sessions else (sessions[-1] if sessions else [])
    latest = {out: (start, end) for out, start, end in rows}  # output -> (start_ms, end_ms), latest wins
    # Wall span is only meaningful within one session (times restart at 0 in each).
    wall_s = None if all_sessions or not rows else (max(e for _, _, e in rows) - min(s for _, s, _ in rows)) / 1000.0

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
        "scope": "all_sessions" if all_sessions else "last_session",
        "wall_s": round(wall_s, 1) if wall_s is not None else None,
        "edges": len(edges),
        "total_edge_cpu_min": round(total / 60.0, 1),
        "by_kind": by_kind,
        "link_share_pct": round(100.0 * by_kind.get("link", {}).get("seconds", 0) / total, 1) if total else 0,
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
        print(f"  edge CPU:  {an['total_edge_cpu_min']} -> {bn['total_edge_cpu_min']} min "
              f"({bn['total_edge_cpu_min'] - an['total_edge_cpu_min']:+.1f})")
        print(f"  link share: {an['link_share_pct']}% -> {bn['link_share_pct']}%")
        if an.get("wall_s") is not None and bn.get("wall_s") is not None:
            print(f"  wall:      {an['wall_s']} -> {bn['wall_s']} s ({bn['wall_s'] - an['wall_s']:+.1f})")
        if an.get("scope", "all_sessions") != bn.get("scope", "all_sessions"):
            print(f"  warning: comparing scope {an.get('scope', 'all_sessions')} with {bn.get('scope', 'all_sessions')}")
    for key in ("executables", "static_libs", "shared_libs"):
        at, bt = a["sizes"][key]["total_mb"], b["sizes"][key]["total_mb"]
        print(f"  {key}: {at} -> {bt} MB ({bt - at:+.1f}), "
              f"count {a['sizes'][key]['count']} -> {b['sizes'][key]['count']}")
    print(f"  tree: {a['sizes']['tree_gb']} -> {b['sizes']['tree_gb']} GB")


def main():
    if "--diff" in sys.argv:
        i = sys.argv.index("--diff")
        diff(sys.argv[i + 1], sys.argv[i + 2])
        return
    build_dir = Path(sys.argv[1]).resolve()
    label = sys.argv[sys.argv.index("--label") + 1] if "--label" in sys.argv else "snapshot"
    out_dir = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else Path(
        os.environ.get("MORIS_RUNS_DIR", ".")) / "benchmarks" / "build_metrics"
    out_dir.mkdir(parents=True, exist_ok=True)

    snap = {"label": label, "build_dir": str(build_dir), "date": date.today().isoformat(),
            "ninja": parse_ninja_log(build_dir, "--all-sessions" in sys.argv), "sizes": collect_sizes(build_dir)}
    out = out_dir / f"{label}_{build_dir.name}_{snap['date']}.json"
    out.write_text(json.dumps(snap, indent=2))

    n, s = snap["ninja"], snap["sizes"]
    if n:
        wall = "" if n["wall_s"] is None else (f", wall {n['wall_s']:.1f} s" if n["wall_s"] < 120
                                                else f", wall {n['wall_s'] / 60:.1f} min")
        print(f"[{label}] {build_dir.name} ({n['scope'].replace('_', ' ')}): {n['edges']} edges, "
              f"{n['total_edge_cpu_min']} CPU-min{wall}, link share {n['link_share_pct']}%")
    print(f"  exes: {s['executables']['count']} = {s['executables']['total_mb']} MB | "
          f".a: {s['static_libs']['count']} = {s['static_libs']['total_mb']} MB | "
          f".so: {s['shared_libs']['count']} = {s['shared_libs']['total_mb']} MB | "
          f"tree {s['tree_gb']} GB")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
