# Staggered circular finned-tube CHT case

This case provides a streamwise-periodic conjugate heat transfer setup for a finned-tube heat exchanger representative elementary volume (REV).

## Steps
1. Generate the conformal periodic mesh with Gmsh:
   ```bash
   python generate_mesh.py --output mesh.msh
   ```
2. Convert the mesh and build the multi-region topology inside OpenFOAM:
   ```bash
   gmshToFoam mesh.msh
   changeDictionary
   splitMeshRegions -cellZones -overwrite
   ```
3. Create the region dictionaries (already populated by `setup_openfoam.py` if committed) or regenerate with custom parameters:
   ```bash
   python setup_openfoam.py --Ubar 1.0 --Tinlet 300 --Tsolid 320
   ```
4. Run the steady CHT simulation:
   ```bash
   chtMultiRegionSimpleFoam
   ```
5. Post-process for ΔP, LMTD, Nu, and friction factor:
   ```bash
   python post_process.py . --Umean 1.0
   ```

### Notes
- Periodicity is enforced for inlet/outlet (x), top/bottom (y), and front/back (z) faces using `gmsh.model.mesh.setPeriodic` to guarantee node matching.
- `fvOptions` uses `meanVelocityForce` to drive the specified bulk velocity instead of fan pressure jumps; a commented energy source block is available for thermal mean-value control when periodic temperature boundaries are used.
- Function objects in `controlDict` record patch-averaged pressure/temperature and integrate wall heat flux for the post-processing script.
