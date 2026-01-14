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

1. **Temperature divergence**: The enthalpy equation (h) produces extreme values after a few iterations, causing the Newton-Raphson T(h) calculation to fail with "Maximum iterations exceeded" error.

2. **Root cause hypothesis**:
   - The tetrahedral mesh combined with small inlet/outlet tubes creates high velocity gradients
   - The convection term in the energy equation produces unrealistic enthalpy values at certain cells
   - Large continuity errors indicate pressure-velocity coupling issues

### Attempted Fixes

1. Reduced relaxation factors (h: 0.0001 → still fails)
2. Changed from water to air (lower Cp, less sensitive)
3. Added Hf offset (300000 J/kg) to keep h positive
4. Used bounded upwind scheme for h equation
5. Tried transient solver (chtMultiRegionFoam) - same issue

### Recommended Next Steps

1. **Improve mesh quality**: Use hexahedral mesh (blockMesh) or snappyHexMesh
2. **Simplify geometry**: Remove inlet/outlet tubes, use flat patches
3. **Alternative solver**: Try buoyantSimpleFoam with coupled thermal BC (if single-region approximation acceptable)
4. **Validate on tutorial**: Test settings on multiRegionHeater tutorial first

## Cleanup
```bash
./Allclean
```
