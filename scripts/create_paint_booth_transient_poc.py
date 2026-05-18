#!/usr/bin/env python3
"""Create a minimal pimpleFoam transient PoC from a steady paint-booth case.

The script reuses a converged steady OpenFOAM case as the initial condition for a
transient step-response smoke test.  It copies the source case's constant/system
configuration, copies the latest steady fields into destination 0/, switches the
solver dictionaries to pimpleFoam, and changes the supply velocity boundary.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Iterable


FIELD_NAMES = ("U", "p", "phi", "k", "omega", "nut")


def numeric_time_dirs(case_dir: Path) -> list[Path]:
    out: list[Path] = []
    for path in case_dir.iterdir():
        if not path.is_dir():
            continue
        try:
            float(path.name)
        except ValueError:
            continue
        out.append(path)
    return sorted(out, key=lambda p: float(p.name))


def latest_time_dir(case_dir: Path) -> Path:
    times = [p for p in numeric_time_dirs(case_dir) if p.name != "0"]
    if not times:
        raise FileNotFoundError(f"No latest steady result time found in {case_dir}")
    return times[-1]


def copytree_clean(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def copy_initial_fields(source_latest: Path, source_zero: Path, dst_zero: Path) -> None:
    if dst_zero.exists():
        shutil.rmtree(dst_zero)
    dst_zero.mkdir(parents=True, exist_ok=True)
    for name in FIELD_NAMES:
        src = source_latest / name
        if not src.exists():
            src = source_zero / name
        if src.exists():
            shutil.copy2(src, dst_zero / name)
    # snappy layer auxiliary fields can exist in 0/; copy only if present and useful.
    for name in ("cellLevel", "pointLevel"):
        src = source_zero / name
        if src.exists():
            shutil.copy2(src, dst_zero / name)


def replace_supply_velocity(u_path: Path, target_velocity: float) -> None:
    text = u_path.read_text(encoding="utf-8")
    replacement = f"supplyInlet {{ type fixedValue; value uniform (0 0 -{target_velocity:g}); }}"
    text, n = re.subn(
        r"supplyInlet\s*\{[^{}]*type\s+fixedValue;[^{}]*value\s+uniform\s*\([^)]*\);[^{}]*\}",
        replacement,
        text,
        flags=re.S,
    )
    if n != 1:
        raise RuntimeError(f"Expected to replace one supplyInlet block in {u_path}, replaced {n}")
    u_path.write_text(text, encoding="utf-8")


def write_control_dict(path: Path, *, end_time: float, delta_t: float, max_delta_t: float, write_interval: float) -> None:
    path.write_text(
        f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}}
application     pimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time:g};
deltaT          {delta_t:g};
writeControl    adjustableRunTime;
writeInterval   {write_interval:g};
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable true;
adjustTimeStep  yes;
maxCo           1.0;
maxDeltaT       {max_delta_t:g};
""",
        encoding="utf-8",
    )


def write_fv_schemes(path: Path) -> None:
    path.write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss upwind;
    div(phi,k) bounded Gauss upwind;
    div(phi,omega) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
""",
        encoding="utf-8",
    )


def write_fv_solution(path: Path) -> None:
    path.write_text(
        """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.05;
        smoother        GaussSeidel;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }

    UFinal
    {
        $U;
        relTol          0;
    }

    "(k|omega|nut)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }

    "(k|omega|nut)Final"
    {
        $k;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor       yes;
    nOuterCorrectors        1;
    nCorrectors             2;
    nNonOrthogonalCorrectors 0;
    pRefCell                0;
    pRefValue               0;
}

relaxationFactors
{
    equations
    {
        U       0.7;
        k       0.7;
        omega   0.7;
    }
}
""",
        encoding="utf-8",
    )


def remove_generated_dirs(case_dir: Path, names: Iterable[str]) -> None:
    for name in names:
        path = case_dir / name
        if path.exists():
            shutil.rmtree(path)
    for path in case_dir.glob("log.*"):
        path.unlink()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-case", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--target-supply-velocity", type=float, default=5.45)
    ap.add_argument("--end-time", type=float, default=5.0)
    ap.add_argument("--delta-t", type=float, default=0.02)
    ap.add_argument("--max-delta-t", type=float, default=0.1)
    ap.add_argument("--write-interval", type=float, default=1.0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    repo = Path.cwd().resolve()
    source_case = args.source_case.resolve() if args.source_case.is_absolute() else (repo / args.source_case).resolve()
    case_dir = args.case_dir.resolve() if args.case_dir.is_absolute() else (repo / args.case_dir).resolve()
    if not source_case.exists():
        raise FileNotFoundError(source_case)
    if case_dir.exists() and args.force:
        shutil.rmtree(case_dir)
    if case_dir.exists():
        raise FileExistsError(f"Destination exists; use --force to replace: {case_dir}")

    source_latest = latest_time_dir(source_case)
    case_dir.mkdir(parents=True)
    copytree_clean(source_case / "constant", case_dir / "constant")
    copytree_clean(source_case / "system", case_dir / "system")
    copy_initial_fields(source_latest, source_case / "0", case_dir / "0")
    remove_generated_dirs(case_dir, ["VTK", "postProcessing"])

    write_control_dict(
        case_dir / "system" / "controlDict",
        end_time=args.end_time,
        delta_t=args.delta_t,
        max_delta_t=args.max_delta_t,
        write_interval=args.write_interval,
    )
    write_fv_schemes(case_dir / "system" / "fvSchemes")
    write_fv_solution(case_dir / "system" / "fvSolution")
    replace_supply_velocity(case_dir / "0" / "U", args.target_supply_velocity)

    metadata = {
        "case_type": "paint_booth_transient_step_poc",
        "source_case": str(source_case.relative_to(repo) if source_case.is_relative_to(repo) else source_case),
        "source_latest_time": source_latest.name,
        "case_dir": str(case_dir.relative_to(repo) if case_dir.is_relative_to(repo) else case_dir),
        "solver": "pimpleFoam",
        "initial_condition": "latest steady fields copied to 0/",
        "target_supply_velocity_mps": args.target_supply_velocity,
        "end_time_s": args.end_time,
        "delta_t_s": args.delta_t,
        "max_delta_t_s": args.max_delta_t,
        "write_interval_s": args.write_interval,
    }
    (case_dir / "transient_poc_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
