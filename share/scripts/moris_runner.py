#!/usr/bin/env python3

"""
Thin Python orchestrator for MORIS examples.

Commands:
  - compile: Build a shared object (.so) from a C++ input file using MORIS's script
  - run    : Run the MORIS main with mpirun on a given .so/.xml
  - level-set-beam: Convenience wrapper for the Level_Set_Beam_SIMP_Hole_Seeding example

Requires:
  - MORISROOT environment variable pointing to the MORIS repo root
  - A configured build directory (e.g., build_dbg) generated with CMake (prefer Ninja)
"""

import argparse
import ast
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def require_env(name: str) -> str:
    val = os.environ.get(name)
    if not val:
        sys.stderr.write(f"Error: environment variable {name} is not set\n")
        sys.exit(2)
    return val


def check_file(path: Path, kind: str = "file"):
    if kind == "file" and not path.is_file():
        sys.stderr.write(f"Error: {path} does not exist or is not a file\n")
        sys.exit(2)
    if kind == "dir" and not path.is_dir():
        sys.stderr.write(f"Error: {path} does not exist or is not a directory\n")
        sys.exit(2)


def run_cmd(cmd, cwd=None, env=None, capture=False):
    if capture:
        proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        return proc.returncode, proc.stdout
    else:
        proc = subprocess.run(cmd, cwd=cwd, env=env)
        return proc.returncode, None


def _deep_merge_dict(into: dict, update: dict) -> dict:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(into.get(key), dict):
            _deep_merge_dict(into[key], value)
        else:
            into[key] = value
    return into


def _parse_override_text(text: str, *, format_hint: str | None, source: str) -> dict:
    errors: list[str] = []

    def _ensure_dict(obj: object) -> dict:
        if obj is None:
            return {}
        if not isinstance(obj, dict):
            raise ValueError(f"Override payload from {source} must be a mapping, got {type(obj).__name__}")
        return obj

    if format_hint in (None, "json"):
        try:
            return _ensure_dict(json.loads(text))
        except json.JSONDecodeError as exc:
            if format_hint == "json":
                raise ValueError(f"Failed to parse JSON overrides from {source}: {exc}") from exc
            errors.append(f"json: {exc}")

    if format_hint in (None, "yaml", "yml"):
        try:
            import yaml  # type: ignore

        except ImportError as exc:
            if format_hint in ("yaml", "yml"):
                raise RuntimeError("PyYAML is required to parse YAML overrides but is not installed") from exc
        else:
            try:
                return _ensure_dict(yaml.safe_load(text))
            except yaml.YAMLError as exc:  # type: ignore[name-defined]
                if format_hint in ("yaml", "yml"):
                    raise ValueError(f"Failed to parse YAML overrides from {source}: {exc}") from exc
                errors.append(f"yaml: {exc}")

    if format_hint in (None, "python"):
        try:
            return _ensure_dict(ast.literal_eval(text))
        except (ValueError, SyntaxError) as exc:
            if format_hint == "python":
                raise ValueError(f"Failed to parse Python literal overrides from {source}: {exc}") from exc
            errors.append(f"python: {exc}")

    raise ValueError(f"Could not parse overrides from {source}; attempted formats: {', '.join(errors) or 'none'}")


def _load_override_source(value: str, format_hint: str | None) -> dict:
    candidate = Path(value)
    text: str
    hint = format_hint
    if candidate.exists():
        text = candidate.read_text()
        if hint is None:
            suffix = candidate.suffix.lower()
            if suffix == ".json":
                hint = "json"
            elif suffix in (".yaml", ".yml"):
                hint = "yaml"
            elif suffix == ".py":
                hint = "python"
        source = str(candidate)
    else:
        text = value
        source = "<inline override>"

    return _parse_override_text(text, format_hint=hint, source=source)


def load_overrides(args) -> dict:
    merged: dict = {}
    sources: list[tuple[str, str | None]] = []
    for item in getattr(args, "overrides", []) or []:
        sources.append((item, None))
    for item in getattr(args, "overrides_json", []) or []:
        sources.append((item, "json"))
    for item in getattr(args, "overrides_yaml", []) or []:
        sources.append((item, "yaml"))
    for item in getattr(args, "overrides_python", []) or []:
        sources.append((item, "python"))

    for value, hint in sources:
        override = _load_override_source(value, hint)
        _deep_merge_dict(merged, override)

    return merged


def compile_shared_object(morisroot: Path, build: str, source: Path, workdir: Path, keep_link: bool = False) -> Path:
    create_so = morisroot / "share" / "scripts" / "create_shared_object.sh"
    check_file(create_so, "file")

    workdir.mkdir(parents=True, exist_ok=True)
    base = source.stem
    dest_cpp = workdir / f"{base}.cpp"
    if source.resolve() != dest_cpp.resolve():
        shutil.copy2(source, dest_cpp)

    args = [str(create_so), str(workdir), build, base]
    if keep_link:
        args.append("dbg")

    # Run the script in the workdir because it checks for <base>.cpp in CWD
    code, out = run_cmd(args, cwd=str(workdir), capture=True)
    if code != 0:
        sys.stderr.write(out or "\n")
        sys.stderr.write("Compilation script failed\n")
        sys.exit(code)

    so_path = workdir / f"{base}.so"
    if not so_path.is_file():
        sys.stderr.write(out or "\n")
        sys.stderr.write(f".so not produced: {so_path}\n")
        sys.exit(1)
    return so_path


def build_logger_args(outputlog: str | None, severity: int | None, directoutput: int | None):
    args = []
    if outputlog:
        args += ["--outputlog", outputlog]
    if severity is not None:
        args += ["--severity", str(severity)]
    if directoutput is not None:
        args += ["--directoutput", str(directoutput)]
    return args


def run_moris(morisroot: Path, build: str, np: int, input_path: Path,
              outputlog: str | None = None, severity: int | None = None, directoutput: int | None = None,
              mpirun: str = "mpirun", extra: list[str] | None = None, cwd: Path | None = None) -> int:
    moris_bin = morisroot / build / "projects" / "mains" / "moris"
    check_file(moris_bin, "file")
    logger_args = build_logger_args(outputlog, severity, directoutput)
    cmd = [mpirun, "-np", str(np), str(moris_bin)] + logger_args + [str(input_path)]
    if extra:
        cmd.extend(extra)
    code, _ = run_cmd(cmd, cwd=str(cwd) if cwd else None)
    return code


def cmd_compile(args):
    morisroot = Path(require_env("MORISROOT"))
    source = Path(args.source)
    check_file(source, "file")
    workdir = Path(args.workdir) if args.workdir else Path.cwd()
    so = compile_shared_object(morisroot, args.build, source, workdir, keep_link=args.dbg_link)
    print(str(so))


def cmd_run(args):
    morisroot = Path(require_env("MORISROOT"))
    input_path = Path(args.input)
    check_file(input_path, "file")
    rc = run_moris(morisroot, args.build, args.np, input_path,
                   outputlog=args.outputlog, severity=args.severity, directoutput=args.directoutput,
                   mpirun=args.mpirun)
    result = {
        "command": "run",
        "np": args.np,
        "build": args.build,
        "input": str(input_path),
        "outputlog": args.outputlog,
        "severity": args.severity,
        "directoutput": args.directoutput,
        "exit_code": rc,
        "parameter_receipt": str(Path.cwd() / "Parameter_Receipt.xml"),
    }
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2))
    print(json.dumps(result))
    sys.exit(rc)


def cmd_level_set_beam(args):
    morisroot = Path(require_env("MORISROOT"))
    example_cpp = morisroot / "projects" / "EXA" / "optimization" / "Level_Set_Beam_SIMP_Hole_Seeding" / "Level_Set_Beam_SIMP_Hole_Seeding.cpp"
    if args.source:
        example_cpp = Path(args.source)
    check_file(example_cpp, "file")

    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        workdir = Path(tempfile.mkdtemp(prefix="moris_level_set_"))

    so = compile_shared_object(morisroot, args.build, example_cpp, workdir, keep_link=args.dbg_link)
    rc = run_moris(morisroot, args.build, args.np, so,
                   outputlog=args.outputlog, severity=args.severity, directoutput=args.directoutput,
                   mpirun=args.mpirun, cwd=workdir)
    result = {
        "command": "level-set-beam",
        "np": args.np,
        "build": args.build,
        "workdir": str(workdir),
        "source": str(example_cpp),
        "so": str(so),
        "exit_code": rc,
        "parameter_receipt": str(workdir / "Parameter_Receipt.xml"),
    }
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2))
    print(json.dumps(result))
    sys.exit(rc)


def cmd_py_run(args):
    try:
        original_flags = None
        if hasattr(sys, "getdlopenflags") and hasattr(os, "RTLD_GLOBAL") and hasattr(os, "RTLD_NOW"):
            original_flags = sys.getdlopenflags()
            sys.setdlopenflags(os.RTLD_GLOBAL | os.RTLD_NOW)
        import moris_py  # local import to avoid dependency for other commands
    except ImportError as exc:
        sys.stderr.write("Error: moris_py module not found. Rebuild MORIS with BUILD_PYBINDINGS=ON.\n")
        sys.exit(2)
    finally:
        if hasattr(sys, "setdlopenflags") and original_flags is not None:
            sys.setdlopenflags(original_flags)

    input_path = Path(args.input)
    check_file(input_path, "file")

    overrides = load_overrides(args)
    logger_args = build_logger_args(args.outputlog, args.severity, args.directoutput)

    parameter_receipt = args.parameter_receipt or str(Path.cwd() / "Parameter_Receipt.xml")

    run_result: moris_py.RunResult
    suffix = input_path.suffix.lower()
    if suffix == ".so":
        run_result = moris_py.run_so(str(input_path), logger_args, args.meshgen, parameter_receipt, overrides)
    elif suffix == ".xml":
        run_result = moris_py.run_xml(str(input_path), logger_args, args.meshgen, parameter_receipt, overrides)
    else:
        sys.stderr.write(f"Error: Unsupported input extension '{input_path.suffix}'. Use .so or .xml.\n")
        sys.exit(2)

    result = {
        "command": "py-run",
        "input": str(input_path),
        "meshgen": args.meshgen,
        "overrides_applied": bool(overrides),
        "parameter_receipt": run_result.parameter_receipt,
        "iqi": run_result.iqi,
        "rank": run_result.rank,
        "size": run_result.size,
        "exit_code": run_result.ret,
    }

    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2))

    print(json.dumps(result))
    sys.exit(run_result.ret)


def add_logger_flags(p: argparse.ArgumentParser):
    p.add_argument("--outputlog", help="Logger output file prefix (per rank)")
    p.add_argument("--severity", type=int, help="Logger severity level (0..3)")
    p.add_argument("--directoutput", type=int, help="Logger direct output format (1..3)")


def main():
    ap = argparse.ArgumentParser(description="MORIS Python runner")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # compile
    ap_compile = sub.add_parser("compile", help="Compile a .cpp input to a .so with MORIS script")
    ap_compile.add_argument("--build", required=True, help="Build directory name (e.g., build_dbg)")
    ap_compile.add_argument("--source", required=True, help="Path to .cpp input file")
    ap_compile.add_argument("--workdir", help="Working directory to build in (default: CWD)")
    ap_compile.add_argument("--dbg-link", action="store_true", help="Keep soft link to input_file.cpp (script dbg mode)")
    ap_compile.set_defaults(func=cmd_compile)

    # run
    ap_run = sub.add_parser("run", help="Run moris on a .so or .xml input")
    ap_run.add_argument("--np", type=int, required=True, help="Number of ranks")
    ap_run.add_argument("--build", required=True, help="Build directory name (e.g., build_dbg)")
    ap_run.add_argument("--input", required=True, help="Path to input (.so or .xml)")
    ap_run.add_argument("--mpirun", default="mpirun", help="MPI runner command (default: mpirun)")
    ap_run.add_argument("--json", help="Write a JSON summary to this file")
    add_logger_flags(ap_run)
    ap_run.set_defaults(func=cmd_run)

    # py-run (in-process via bindings)
    ap_py = sub.add_parser("py-run", help="Run MORIS in-process via moris_py with optional overrides")
    ap_py.add_argument("--input", required=True, help="Path to input (.so or .xml)")
    ap_py.add_argument("--meshgen", action="store_true", help="Use mesh generation library variant")
    ap_py.add_argument("--parameter-receipt", help="Path for Parameter_Receipt.xml (default: ./Parameter_Receipt.xml)")
    ap_py.add_argument("--json", help="Write a JSON summary to this file")
    ap_py.add_argument("--overrides", action="append", help="Override payload (file path or inline JSON/YAML/Python dict)")
    ap_py.add_argument("--overrides-json", action="append", help="Override payload in JSON (file path or inline)")
    ap_py.add_argument("--overrides-yaml", action="append", help="Override payload in YAML (file path or inline)")
    ap_py.add_argument("--overrides-python", action="append", help="Override payload as Python literal (file path or inline)")
    add_logger_flags(ap_py)
    ap_py.set_defaults(func=cmd_py_run)

    # level-set-beam convenience
    ap_lvl = sub.add_parser("level-set-beam", help="Compile and run the Level_Set_Beam_SIMP_Hole_Seeding example")
    ap_lvl.add_argument("--np", type=int, default=1, help="Number of ranks (default: 1)")
    ap_lvl.add_argument("--build", required=True, help="Build directory name (e.g., build_dbg)")
    ap_lvl.add_argument("--source", help="Override path to example .cpp (defaults to repo example)")
    ap_lvl.add_argument("--workdir", help="Working directory (default: temp dir)")
    ap_lvl.add_argument("--mpirun", default="mpirun", help="MPI runner command (default: mpirun)")
    ap_lvl.add_argument("--dbg-link", action="store_true", help="Keep soft link to input_file.cpp (script dbg mode)")
    ap_lvl.add_argument("--json", help="Write a JSON summary to this file")
    add_logger_flags(ap_lvl)
    ap_lvl.set_defaults(func=cmd_level_set_beam)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
