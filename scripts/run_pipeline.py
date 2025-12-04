#!/usr/bin/env python3
"""
OpenFOAM CFD Pipeline Runner
Complete pipeline: Geometry -> Mesh -> Solve -> Visualize
"""

import os
import sys
import subprocess
import argparse
import shutil
import time

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from geometry.heatsink import create_heatsink_geometry
from mesh.convert_mesh import run_openfoam_command
from solver.cht_case_setup import setup_cht_case
from visualization.visualize_cht import visualize_cht_results


def run_pipeline(
    case_name: str = "heatsink_cht",
    output_base: str = None,
    num_fins: int = 5,
    fin_height: float = 25e-3,
    inlet_velocity: float = 1.0,
    inlet_temp: float = 300,
    heat_flux: float = 5000,
    end_time: float = 50,
    mesh_size: float = 2.0e-3,
    skip_geometry: bool = False,
    skip_mesh: bool = False,
    skip_solve: bool = False,
    skip_visualize: bool = False
):
    """
    Run complete CFD pipeline.

    Args:
        case_name: Name of the case
        output_base: Base directory for output (default: cases/)
        num_fins: Number of heatsink fins
        fin_height: Height of fins in meters
        inlet_velocity: Air inlet velocity (m/s)
        inlet_temp: Air inlet temperature (K)
        heat_flux: Heat flux at heatsink bottom (W/m²)
        end_time: Simulation end time (s)
        mesh_size: Base mesh size in meters
        skip_*: Skip specific pipeline stages
    """
    start_time = time.time()

    # Setup directories
    if output_base is None:
        output_base = os.path.join(os.path.dirname(__file__), '..', 'cases')

    case_dir = os.path.join(output_base, case_name)
    os.makedirs(case_dir, exist_ok=True)

    results_dir = os.path.join(case_dir, 'results')
    os.makedirs(results_dir, exist_ok=True)

    msh_file = os.path.join(case_dir, f"{case_name}.msh")

    print("=" * 70)
    print("OpenFOAM CFD Pipeline")
    print("=" * 70)
    print(f"Case: {case_name}")
    print(f"Directory: {case_dir}")
    print(f"Parameters:")
    print(f"  Fins: {num_fins}, Height: {fin_height*1000:.1f} mm")
    print(f"  Inlet: {inlet_velocity} m/s, {inlet_temp} K")
    print(f"  Heat flux: {heat_flux} W/m²")
    print(f"  End time: {end_time} s")
    print("=" * 70)

    # Stage 1: Geometry Generation
    if not skip_geometry:
        print("\n[Stage 1/4] Generating Geometry with Gmsh...")
        print("-" * 50)

        create_heatsink_geometry(
            num_fins=num_fins,
            fin_height=fin_height,
            mesh_size_solid=mesh_size,
            mesh_size_fluid=mesh_size * 1.5,
            mesh_size_interface=mesh_size * 0.7,
            output_file=msh_file
        )
    else:
        print("\n[Stage 1/4] Skipping geometry generation")

    # Stage 2: Mesh Conversion
    if not skip_mesh:
        print("\n[Stage 2/4] Converting Mesh to OpenFOAM...")
        print("-" * 50)

        # Create minimal system files for mesh conversion
        system_dir = os.path.join(case_dir, "system")
        os.makedirs(system_dir, exist_ok=True)

        # Minimal controlDict for mesh tools
        control_dict = """FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application chtMultiRegionFoam;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 1;
deltaT 0.01;
writeControl timeStep;
writeInterval 100;
"""
        with open(os.path.join(system_dir, "controlDict"), 'w') as f:
            f.write(control_dict)

        # fvSchemes
        fv_schemes = """FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""
        with open(os.path.join(system_dir, "fvSchemes"), 'w') as f:
            f.write(fv_schemes)

        # fvSolution
        fv_solution = """FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers { }
"""
        with open(os.path.join(system_dir, "fvSolution"), 'w') as f:
            f.write(fv_solution)

        # Convert mesh
        run_openfoam_command(
            f"gmshToFoam {case_name}.msh",
            case_dir,
            "Converting Gmsh mesh to OpenFOAM"
        )

        # Create cellZones from cellSets
        toposet_dict = """FoamFile { version 2.0; format ascii; class dictionary; object topoSetDict; }
actions
(
    { name fluid; type cellZoneSet; action new; source setToCellZone; set fluid; }
    { name solid; type cellZoneSet; action new; source setToCellZone; set solid; }
);
"""
        with open(os.path.join(system_dir, "topoSetDict"), 'w') as f:
            f.write(toposet_dict)

        run_openfoam_command("topoSet", case_dir, "Creating cellZones")

        # Split mesh into regions
        run_openfoam_command(
            "splitMeshRegions -cellZones -overwrite",
            case_dir,
            "Splitting mesh into regions"
        )

        print("Mesh conversion complete!")
    else:
        print("\n[Stage 2/4] Skipping mesh conversion")

    # Stage 3: Setup and Solve
    if not skip_solve:
        print("\n[Stage 3/4] Setting up and running CHT simulation...")
        print("-" * 50)

        # Setup case files
        setup_cht_case(
            case_dir,
            inlet_velocity=inlet_velocity,
            inlet_temp=inlet_temp,
            heat_flux=heat_flux,
            end_time=end_time
        )

        # Create region fvSchemes and fvSolution
        for region in ['fluid', 'solid']:
            region_system = os.path.join(case_dir, 'system', region)
            os.makedirs(region_system, exist_ok=True)

            # Copy main fvSchemes/fvSolution
            main_system = os.path.join(case_dir, 'system')
            for fname in ['fvSchemes', 'fvSolution']:
                src = os.path.join(main_system, fname)
                dst = os.path.join(region_system, fname)
                if os.path.exists(src) and not os.path.exists(dst):
                    shutil.copy(src, dst)

        # Run solver
        print("Running chtMultiRegionFoam solver...")
        try:
            run_openfoam_command(
                "chtMultiRegionFoam",
                case_dir,
                "Running CHT simulation"
            )
        except RuntimeError as e:
            print(f"Solver encountered an issue: {e}")
            print("Continuing with available results...")

        # Convert to VTK
        print("Converting results to VTK format...")
        for region in ['fluid', 'solid']:
            try:
                run_openfoam_command(
                    f"foamToVTK -region {region}",
                    case_dir,
                    f"Converting {region} region to VTK"
                )
            except Exception as e:
                print(f"Warning: VTK conversion for {region} failed: {e}")

    else:
        print("\n[Stage 3/4] Skipping solve")

    # Stage 4: Visualization
    if not skip_visualize:
        print("\n[Stage 4/4] Creating visualizations...")
        print("-" * 50)

        success = visualize_cht_results(case_dir, results_dir)

        if not success:
            print("Visualization had issues, but pipeline continues...")
    else:
        print("\n[Stage 4/4] Skipping visualization")

    # Summary
    elapsed = time.time() - start_time
    print("\n" + "=" * 70)
    print("Pipeline Complete!")
    print("=" * 70)
    print(f"Total time: {elapsed:.1f} seconds")
    print(f"Case directory: {case_dir}")
    print(f"Results directory: {results_dir}")
    print("=" * 70)

    return case_dir


def main():
    parser = argparse.ArgumentParser(
        description='Run OpenFOAM CFD Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline with defaults
  python run_pipeline.py

  # Custom heatsink parameters
  python run_pipeline.py --num-fins 7 --fin-height 30 --heat-flux 10000

  # Skip geometry (use existing mesh)
  python run_pipeline.py --skip-geometry

  # Only visualize existing results
  python run_pipeline.py --skip-geometry --skip-mesh --skip-solve
        """
    )

    parser.add_argument('--case-name', default='heatsink_cht', help='Case name')
    parser.add_argument('--output-dir', default=None, help='Output base directory')

    # Geometry parameters
    parser.add_argument('--num-fins', type=int, default=5, help='Number of fins')
    parser.add_argument('--fin-height', type=float, default=25, help='Fin height (mm)')
    parser.add_argument('--mesh-size', type=float, default=2.0, help='Mesh size (mm)')

    # Simulation parameters
    parser.add_argument('--velocity', type=float, default=1.0, help='Inlet velocity (m/s)')
    parser.add_argument('--temperature', type=float, default=300, help='Inlet temperature (K)')
    parser.add_argument('--heat-flux', type=float, default=5000, help='Heat flux (W/m²)')
    parser.add_argument('--end-time', type=float, default=50, help='End time (s)')

    # Skip flags
    parser.add_argument('--skip-geometry', action='store_true', help='Skip geometry generation')
    parser.add_argument('--skip-mesh', action='store_true', help='Skip mesh conversion')
    parser.add_argument('--skip-solve', action='store_true', help='Skip solver execution')
    parser.add_argument('--skip-visualize', action='store_true', help='Skip visualization')

    args = parser.parse_args()

    run_pipeline(
        case_name=args.case_name,
        output_base=args.output_dir,
        num_fins=args.num_fins,
        fin_height=args.fin_height * 1e-3,
        inlet_velocity=args.velocity,
        inlet_temp=args.temperature,
        heat_flux=args.heat_flux,
        end_time=args.end_time,
        mesh_size=args.mesh_size * 1e-3,
        skip_geometry=args.skip_geometry,
        skip_mesh=args.skip_mesh,
        skip_solve=args.skip_solve,
        skip_visualize=args.skip_visualize
    )


if __name__ == "__main__":
    main()
