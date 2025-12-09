# Staggered Finned-Tube Heat Exchanger - Full Domain (MPI Parallel)

This case simulates conjugate heat transfer through a **full staggered finned-tube heat exchanger** with N REV (Representative Elementary Volume) repetitions in the streamwise direction.

## Key Features

- **Real inlet/outlet boundary conditions** (not cyclic)
- **MPI parallel computation** enabled
- **Inlet/outlet extension regions** for numerical stability
- Captures **actual flow development** through the heat exchanger

## Domain Configuration

```
|←─ inlet ─→|←──────── HX Core (10 REV) ─────────→|←─── outlet ───→|
|  extension |  ●  ○  ●  ○  ●  ○  ●  ○  ●  ○     |   extension    |
|   (1×S_L)  |    ○  ●  ○  ●  ○  ●  ○  ●  ○  ●   |    (5×S_L)     |
|   60 mm    |              1200 mm                |    300 mm      |
|            |                                      |                |
└────────────┴──────────────────────────────────────┴────────────────┘
                        Total: 1560 mm
```

### Why Extension Regions?

| Region | Length | Purpose |
|--------|--------|---------|
| **Inlet** | 1×S_L (60mm) | Uniform flow before first tube row |
| **Outlet** | 5×S_L (300mm) | Wake recovery, prevents backflow at outlet |

The outlet extension is critical for numerical stability - without it, the simulation diverged at T=226.

## Geometry Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| S_L | 60 mm | Longitudinal pitch |
| S_T | 60 mm | Transverse pitch |
| P_fin | 50 mm | Fin spacing (z-direction) |
| Tube outer radius | 12 mm | |
| Tube inner radius | 10 mm | |
| Fin radius | 22 mm | |
| Fin thickness | 1 mm | |

## Boundary Conditions

| Boundary | Type | Value |
|----------|------|-------|
| inlet | fixedValue | U = (1, 0, 0) m/s, T = 300 K |
| outlet | zeroGradient/fixedValue | dU/dx = 0, p = 0 |
| top, bottom | symmetry | |
| front, back | symmetry | |
| fluid-solid | coupled | CHT (turbulentTemperatureCoupledBaffleMixed) |

## Running the Case

### Prerequisites
- OpenFOAM v1912+
- Gmsh with Python API
- Python 3: `gmsh`, `numpy`, `pyvista`, `matplotlib`

### 1) Mesh Generation

```bash
cd cases/staggered_finned_tube_full

# Generate mesh with extensions (default: inlet=1×S_L, outlet=5×S_L)
python3 scripts/generate_mesh.py --output mesh.msh -N 10 --inlet-ext 1 --outlet-ext 5

# Convert to OpenFOAM format
gmshToFoam mesh.msh

# Split into fluid/solid regions
splitMeshRegions -cellZones -overwrite
```

### 2) Setup Boundary Conditions

```bash
python3 scripts/setup_openfoam.py --post-split --Ubar 1.0 --Tinlet 300 --Tsolid 350 -N 10 --nprocs 4
```

### 3) Domain Decomposition & Parallel Run

```bash
# Setup decomposeParDict for all regions
# (done automatically by setup script)

# Decompose
decomposePar -allRegions

# Run with MPI
mpirun -np 4 chtMultiRegionSimpleFoam -parallel > log.solver 2>&1
```

### 4) Post-processing

```bash
# Reconstruct results
reconstructPar -allRegions -latestTime

# Convert to VTK
foamToVTK -latestTime -region fluid
foamToVTK -latestTime -region solid

# Generate visualization images
python3 scripts/visualize_results.py
```

## Results

### Temperature Development

| Location | Temperature |
|----------|-------------|
| Inlet | 300 K |
| Outlet | 315.5 K |
| Temperature rise | +15.5 K |

### Performance

| Metric | Value |
|--------|-------|
| Mesh cells (fluid) | 639,867 |
| Solid regions | 30 |
| MPI processors | 4 |
| Iterations | 1000 |
| Runtime | ~16 minutes |

## Generated Files

| File | Description |
|------|-------------|
| `images/temperature_field.png` | Temperature contours (z=0 plane) |
| `images/velocity_field.png` | Velocity magnitude (z=0 plane) |
| `images/temperature_profile.png` | Centerline temperature |
| `images/geometry_3d.png` | 3D geometry view |
| `results_summary.txt` | Numerical results summary |
| `log.solver` | Solver output log |

## Comparison with Cyclic Case

| Feature | Cyclic (cyclicAMI) | Full Domain |
|---------|-------------------|-------------|
| MPI parallel | ❌ Not supported | ✅ Supported |
| Flow development | Fully developed only | Includes entrance region |
| Domain size | 1 REV | N REV + extensions |
| Numerical stability | Good | Requires outlet extension |
| Computational cost | Low | ~N× higher |
