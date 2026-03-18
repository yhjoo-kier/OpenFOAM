# Indoor ventilation OpenFOAM case

Generated from:
- Mesh: `/home/yhjoo/projects/OpenFOAM/generated/bench_a3_04.msh`
- Solver: `simpleFoam`
- Turbulence model: `kOmegaSST`

Expected mesh boundary patches:
- `inlet`
- `outlet`
- `roomWalls`
- `obstacleWalls`

Typical run sequence when OpenFOAM is available:
```bash
source /usr/share/openfoam/etc/bashrc  # or project env
cd /home/yhjoo/projects/OpenFOAM/cases/bench_a3_04
gmshToFoam bench_a3_04.msh
checkMesh
simpleFoam
```
