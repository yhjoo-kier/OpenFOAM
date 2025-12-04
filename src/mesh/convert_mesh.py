#!/usr/bin/env python3
"""
Mesh conversion utilities for OpenFOAM
Converts Gmsh mesh to OpenFOAM format and prepares for multi-region CHT
"""

import subprocess
import os
import shutil


def run_openfoam_command(cmd, case_dir, description=""):
    """Run OpenFOAM command with proper environment."""
    env = os.environ.copy()
    env['WM_PROJECT_DIR'] = '/usr/share/openfoam'
    env['FOAM_ETC'] = '/usr/share/openfoam/etc'
    env['WM_PROJECT'] = 'OpenFOAM'

    print(f"Running: {description or cmd}")
    result = subprocess.run(
        cmd,
        shell=True,
        cwd=case_dir,
        env=env,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")

    return result.stdout


def convert_gmsh_to_openfoam(msh_file: str, case_dir: str):
    """
    Convert Gmsh mesh to OpenFOAM polyMesh format.

    Args:
        msh_file: Path to Gmsh .msh file
        case_dir: OpenFOAM case directory
    """
    # Ensure case directory exists
    os.makedirs(case_dir, exist_ok=True)
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)

    # Copy mesh file to case directory if not already there
    msh_basename = os.path.basename(msh_file)
    msh_dest = os.path.join(case_dir, msh_basename)
    if os.path.abspath(msh_file) != os.path.abspath(msh_dest):
        shutil.copy(msh_file, msh_dest)

    # Run gmshToFoam
    run_openfoam_command(
        f"gmshToFoam {msh_basename}",
        case_dir,
        "Converting Gmsh mesh to OpenFOAM format"
    )

    print(f"Mesh converted successfully to: {case_dir}/constant/polyMesh")


def create_cell_zones_from_physical_groups(case_dir: str):
    """
    Create cellZones from Gmsh physical groups.
    This is needed for splitMeshRegions to work properly.
    """
    # Create topoSetDict to convert cellSets to cellZones
    toposet_dict = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    {
        name    fluid;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     fluid;
    }
    {
        name    solid;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        set     solid;
    }
);

// ************************************************************************* //
"""

    system_dir = os.path.join(case_dir, "system")
    with open(os.path.join(system_dir, "topoSetDict"), 'w') as f:
        f.write(toposet_dict)

    # Run topoSet to create cellZones
    run_openfoam_command("topoSet", case_dir, "Creating cellZones from cellSets")


def split_mesh_regions(case_dir: str):
    """
    Split mesh into separate regions for multi-region simulation.
    """
    run_openfoam_command(
        "splitMeshRegions -cellZones -overwrite",
        case_dir,
        "Splitting mesh into fluid and solid regions"
    )

    print("Mesh split into regions: fluid, solid")


def check_mesh(case_dir: str, region: str = None):
    """Check mesh quality."""
    cmd = "checkMesh"
    if region:
        cmd += f" -region {region}"

    run_openfoam_command(cmd, case_dir, f"Checking mesh{' for ' + region if region else ''}")


def main():
    """Main function for mesh conversion."""
    import argparse

    parser = argparse.ArgumentParser(description='Convert Gmsh mesh to OpenFOAM')
    parser.add_argument('msh_file', help='Input Gmsh mesh file')
    parser.add_argument('case_dir', help='OpenFOAM case directory')
    parser.add_argument('--split', action='store_true', help='Split mesh into regions')

    args = parser.parse_args()

    convert_gmsh_to_openfoam(args.msh_file, args.case_dir)

    if args.split:
        create_cell_zones_from_physical_groups(args.case_dir)
        split_mesh_regions(args.case_dir)
        check_mesh(args.case_dir, "fluid")
        check_mesh(args.case_dir, "solid")


if __name__ == "__main__":
    main()
