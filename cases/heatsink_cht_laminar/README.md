# Heatsink CHT Laminar Case

3D conjugate heat transfer (CHT) simulation of a pin-fin heatsink with laminar water cooling.

## Geometry

### Main Heatsink Domain
| Parameter | Value | Description |
|-----------|-------|-------------|
| L_dom | 80 mm | x-direction length |
| W_dom | 83 mm | y-direction length |
| H_fluid | 15 mm | Fluid/fin height |
| H_base | 2 mm | Base plate thickness |

### Pin Fins (8x8 Array)
| Parameter | Value | Description |
|-----------|-------|-------------|
| D_fin | 6 mm | Fin diameter |
| Pitch (x, y) | 10 mm | Center-to-center spacing |
| Origin | (5, 6.5) mm | First fin center position |
| Height | 15 mm | Fin height (same as fluid) |

### Inlet/Outlet Tubes
| Parameter | Value | Description |
|-----------|-------|-------------|
| D_tube | 5 mm | Tube diameter |
| L_tube | 25 mm | Tube length |
| Inlet position | y=5mm, z=9.5mm | At x=0, extending to x=-25mm |
| Outlet position | y=78mm, z=9.5mm | At x=80mm, extending to x=105mm |

## Boundary Conditions

### Fluid Region
| Boundary | Type | Value |
|----------|------|-------|
| inlet | fixedValue (pressure) | 10 Pa |
| outlet | fixedValue (pressure) | 0 Pa |
| inlet | fixedValue (temperature) | 293.15 K (20°C) |
| walls | noSlip | - |
| interface | coupled | turbulentTemperatureCoupledBaffleMixed |

### Solid Region
| Boundary | Type | Value |
|----------|------|-------|
| base_bottom | externalWallHeatFluxTemperature | 60000 W/m² |
| walls_solid | zeroGradient | - |
| interface | coupled | turbulentTemperatureCoupledBaffleMixed |

## Material Properties

### Fluid (Air at 20°C)
| Property | Value |
|----------|-------|
| Density | 1.18 kg/m³ |
| Specific heat | 1005 J/kg·K |
| Dynamic viscosity | 1.8e-5 Pa·s |
| Prandtl number | 0.71 |
| Hf (reference enthalpy) | 300000 J/kg |

### Solid (Aluminum)
| Property | Value |
|----------|-------|
| Density | 2700 kg/m³ |
| Specific heat | 900 J/kg·K |
| Thermal conductivity | 200 W/m·K |

## Running the Case

### Prerequisites
- OpenFOAM v1912 or later
- Gmsh with Python API
- Python 3 with: gmsh, numpy

### Quick Start
```bash
./Allrun
```

### Step-by-Step Execution

1. **Generate mesh**
```bash
python3 scripts/generate_mesh.py --output mesh.msh
```

2. **Convert to OpenFOAM format**
```bash
gmshToFoam mesh.msh
```

3. **Split into regions**
```bash
splitMeshRegions -cellZones -overwrite
```

4. **Setup boundary conditions**
```bash
python3 scripts/setup_openfoam.py --post-split
```

5. **Run solver**
```bash
chtMultiRegionSimpleFoam
```

### Mesh Options
```bash
# Default mesh (fine)
python3 scripts/generate_mesh.py --output mesh.msh

# Coarse mesh for testing
python3 scripts/generate_mesh.py --output mesh.msh --coarse

# Custom mesh size
python3 scripts/generate_mesh.py --output mesh.msh --mesh-size 0.002
```

### Custom Parameters
```bash
python3 scripts/setup_openfoam.py --post-split \
    --Tinlet 293.15 \
    --heatFlux 60000 \
    --p-inlet 10 \
    --p-outlet 0
```

## Post-Processing

### Convert to VTK
```bash
foamToVTK -latestTime -region fluid
foamToVTK -latestTime -region solid
```

### View in ParaView
```bash
paraview VTK/fluid/fluid_*.vtk
```

### Check convergence
```bash
grep "Time =" log.chtMultiRegionSimpleFoam | tail -10
```

## Output Files

| Path | Description |
|------|-------------|
| `log.chtMultiRegionSimpleFoam` | Solver log |
| `VTK/fluid/` | Fluid region VTK files |
| `VTK/solid/` | Solid region VTK files |
| `postProcessing/` | Function object output |

## Notes

- This is a **laminar** flow case. For turbulent flow, modify `turbulenceProperties` and add appropriate turbulence fields.
- The solver uses `chtMultiRegionSimpleFoam` for steady-state CHT analysis.
- Currently configured with velocity inlet (0.1 m/s) instead of pressure-driven flow.
- Heat flux applied at base plate bottom surface (currently reduced to 1000 W/m² for stability testing).

## Known Issues and Status

**Status**: Case setup complete, but numerical stability issues remain unresolved.

### Current Issues

1. **Temperature divergence**: The enthalpy equation (h) produces extreme values after a few iterations, causing the Newton-Raphson T(h) calculation to fail with "Maximum iterations exceeded" or "Negative initial temperature T0" errors.

2. **Root cause analysis (updated Jan 2026)**:
   - **NOT a mesh quality issue**: Testing with simple blockMesh hexahedral meshes shows identical instability
   - **Fundamental solver limitation**: chtMultiRegionSimpleFoam and buoyantSimpleFoam with heRhoThermo have numerical instability in pressure-velocity-energy coupling for low-speed, low-temperature-difference flows
   - The continuity errors are extremely large (40000+) from the first iteration
   - The h equation "converges" (residual ~1e-07) but produces extreme enthalpy values
   - This causes T to become negative, which crashes the density calculation (rho → 0 or ∞)

### Attempted Fixes

1. Reduced relaxation factors (h: 0.01, 0.001, 0.0001 → all fail)
2. Changed from water to air (lower Cp, less sensitive)
3. Added Hf offset (300000 J/kg) to keep h positive
4. Used bounded upwind scheme for h equation
5. Tried transient solver (chtMultiRegionFoam) - marginally more stable but still crashes
6. Tested with perfectGas equation of state - same issue
7. Tested with rhoConst equation of state - same issue
8. Tested with simple hexahedral mesh - same issue (confirms NOT mesh related)
9. Tested with zero initial velocity - still diverges due to induced velocities
10. Tested with zero gravity - no improvement

### Working Alternative: buoyantBoussinesqSimpleFoam

**buoyantBoussinesqSimpleFoam** (incompressible Boussinesq approximation) works correctly:
- Uses temperature T directly instead of enthalpy h
- Incompressible formulation avoids density-temperature coupling instability
- Converges reliably for the same geometry and boundary conditions
- **Limitation**: Single-region solver, cannot directly do CHT

### Recommended Approaches

1. **Boussinesq + Manual Coupling**:
   - Use buoyantBoussinesqSimpleFoam for fluid region
   - Use laplacianFoam for solid region
   - Iteratively exchange temperature at interface

2. **Use Established CHT Cases**:
   - Start from OpenFOAM tutorial cases (e.g., multiRegionHeater)
   - Modify geometry incrementally while maintaining working configuration

3. **Alternative CFD Packages**:
   - Commercial solvers (ANSYS Fluent, STAR-CCM+) may handle low-Ma CHT more gracefully
   - Or use specialized conjugate heat transfer codes

4. **Higher Fidelity Approach**:
   - Use finer time stepping with chtMultiRegionFoam (transient)
   - Very small time steps (1e-6 s) with maxCo < 0.1
   - Much longer computation time required

## Cleanup
```bash
./Allclean
```
