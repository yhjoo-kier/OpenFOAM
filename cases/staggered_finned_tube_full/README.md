# Staggered Finned-Tube Heat Exchanger - Full Domain (MPI Parallel)

This case simulates conjugate heat transfer through a **full staggered finned-tube heat exchanger** with N REV (Representative Elementary Volume) repetitions in the streamwise direction.

Unlike the periodic (cyclicAMI) case, this setup:
- Uses **real inlet/outlet boundary conditions**
- Captures **actual flow development** through the heat exchanger
- Enables **MPI parallel computation**

## Geometry

### Domain Configuration

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| REV length | 2×S_L | 0.12 m | Single REV streamwise length |
| Total length | N×2×S_L | 1.2 m (N=10) | Full domain streamwise |
| Transverse | S_T | 0.06 m | Channel height |
| Axial | P_fin | 0.05 m | Fin spacing |

### Tube Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Tube outer radius | 0.012 m | |
| Tube inner radius | 0.010 m | |
| Fin radius | 0.022 m | |
| Fin thickness | 0.001 m | |

### Staggered Arrangement

Each REV unit contains:
- **Row 1** (x = S_L/2): Tubes at y = 0 and y = S_T (half tubes at boundaries)
- **Row 2** (x = 3S_L/2): Full tube at y = S_T/2 (offset)

With N=10 repetitions, there are **30 tube positions** total.

## Boundary Conditions

| Boundary | Type | Condition |
|----------|------|-----------|
| inlet | fixedValue | U = (1, 0, 0) m/s, T = 300 K |
| outlet | zeroGradient | Natural outflow |
| top, bottom | symmetry | |
| front, back | symmetry | |
| fluid-solid interface | coupled | CHT coupling |

## Running the Case

### Prerequisites
- OpenFOAM v1912+
- Gmsh with Python API
- Python 3 with: `gmsh`, `numpy`, `pyvista`, `matplotlib`

### 1) Mesh Generation

```bash
cd cases/staggered_finned_tube_full

# Generate mesh with N=10 REV repetitions
python3 scripts/generate_mesh.py --output mesh.msh -N 10

# Convert to OpenFOAM format
gmshToFoam mesh.msh

# Split into fluid/solid regions
splitMeshRegions -cellZones -overwrite
```

### 2) Setup Boundary Conditions

```bash
python3 scripts/setup_openfoam.py --post-split --Ubar 1.0 --Tinlet 300 --Tsolid 350 -N 10 --nprocs 4
```

### 3) Domain Decomposition (MPI)

```bash
# Decompose for parallel run
decomposePar -allRegions
```

### 4) Parallel Execution

```bash
# Run with 4 MPI processes
mpirun -np 4 chtMultiRegionSimpleFoam -parallel > log.solver 2>&1

# Monitor progress
tail -f log.solver
```

### 5) Reconstruct and Post-process

```bash
# Reconstruct from parallel data
reconstructPar -allRegions -latestTime

# Convert to VTK
foamToVTK -latestTime -region fluid
foamToVTK -latestTime -region solid
```

## Expected Results

### Flow Development
- Entrance region with developing velocity profile
- Periodic fully-developed pattern after ~3-4 REV units
- Wake interactions between tube rows

### Heat Transfer
- Temperature increase along flow direction
- Local heat transfer enhancement at tube surfaces
- Fin effectiveness demonstration

## MPI Parallelization

This case uses the **scotch** decomposition method for automatic load balancing:

```
numberOfSubdomains 4;
method scotch;
```

For large runs, increase `numberOfSubdomains` based on available cores.

## Files

| Item | Path |
|------|------|
| Mesh script | `scripts/generate_mesh.py` |
| Setup script | `scripts/setup_openfoam.py` |
| Solver log | `log.solver` |
| Results | `{time}/fluid/`, `{time}/solid/` |
| VTK files | `VTK/fluid/`, `VTK/solid/` |

