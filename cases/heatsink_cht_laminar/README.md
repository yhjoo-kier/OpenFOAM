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

### Fluid (Water at 20°C)
| Property | Value |
|----------|-------|
| Density | 998 kg/m³ |
| Specific heat | 4182 J/kg·K |
| Thermal conductivity | 0.6 W/m·K |
| Dynamic viscosity | 0.001 Pa·s |
| Prandtl number | 7.0 |

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
- Pressure-driven flow: ΔP = 10 Pa between inlet and outlet.
- Heat flux applied at base plate bottom surface (60 kW/m²).

## Cleanup
```bash
./Allclean
```
