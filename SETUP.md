# VM Setup Guide for OpenFOAM CFD Pipeline

This document provides instructions for setting up a VM environment to run the OpenFOAM CFD pipeline examples.

## Quick Setup (One-liner)

For Ubuntu 24.04 LTS, run:

```bash
sudo bash scripts/setup_vm.sh
```

## Manual Setup

If you prefer to set up manually or are on a different system:

### 1. Install System Packages

```bash
sudo apt-get update
sudo apt-get install -y openfoam gmsh build-essential libosmesa6-dev libglx-mesa0 python3 python3-pip
```

### 2. Install Python Packages

```bash
pip3 install gmsh pyvista matplotlib numpy vtk
```

### 3. Set OpenFOAM Environment Variables

Add to your `~/.bashrc` or run before each session:

```bash
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM
```

Or source the provided environment script:

```bash
source scripts/env.sh
```

## Verification

Run the verification script to check all components:

```bash
./scripts/verify_installation.sh
```

Expected output:
```
========================================
OpenFOAM CFD Pipeline - Verification
========================================

System Tools:
Checking blockMesh...                   OK
Checking laplacianFoam...               OK
Checking icoFoam...                     OK
Checking foamToVTK...                   OK
Checking gmsh...                        OK
Checking python3...                     OK

Python Packages:
Checking gmsh (Python)...               OK
Checking pyvista...                     OK
Checking matplotlib...                  OK
Checking numpy...                       OK
Checking vtk...                         OK

Project Files:
Checking simple_heatsink.py...          OK
Checking run_pipeline.py...             OK
Checking channelFlow case...            OK

========================================
Results: 14 passed, 0 failed
========================================
All checks passed! The environment is ready.
```

## Running Examples

### Example 1: Heatsink Heat Conduction

```bash
# Source environment
source scripts/env.sh

# Create case
python3 src/geometry/simple_heatsink.py cases/my_heatsink --heat-flux 5000

# Run simulation
cd cases/my_heatsink
blockMesh
laplacianFoam
foamToVTK
```

### Example 2: Channel Flow (Incompressible)

```bash
# Source environment
source scripts/env.sh

# Run simulation
cd channelFlow
blockMesh
icoFoam
foamToVTK
```

## Installed Components

| Component | Version | Description |
|-----------|---------|-------------|
| OpenFOAM | v1912 | CFD solver suite |
| Gmsh | v4.12 | Mesh generation |
| Python | 3.11+ | Scripting environment |
| pyvista | 0.43+ | 3D visualization |
| matplotlib | 3.8+ | 2D plotting |
| VTK | 9.3+ | Visualization toolkit |

## Troubleshooting

### OpenFOAM commands not found

Make sure environment variables are set:
```bash
source scripts/env.sh
```

### Python import errors

Reinstall Python packages:
```bash
pip3 install --upgrade gmsh pyvista matplotlib numpy vtk
```

### Permission denied on scripts

Make scripts executable:
```bash
chmod +x scripts/*.sh
```
