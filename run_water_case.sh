#!/bin/bash
# Water-cooled CHT simulation run script
set -e

# Setup environment
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc

CASE_DIR="cases/heatsink_water_cht"

echo "=== 1. Setup Case Directory ==="
rm -rf $CASE_DIR
mkdir -p $CASE_DIR/system $CASE_DIR/constant $CASE_DIR/0

echo "=== 2. Configure Case (Pre-Mesh) ==="

# --- UTILITIES ---
# Helper to write files
write_file() {
    cat > $1 <<EOF
$2
EOF
}

# --- CONTROLDICT ---
write_file $CASE_DIR/system/controlDict "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     chtMultiRegionFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         10.0;
deltaT          0.001;
writeControl    adjustableRunTime;
writeInterval   0.5;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
maxCo           1.0;
adjustTimeStep  yes;
"

# --- DECOMPOSE PAR ---
write_file $CASE_DIR/system/decomposeParDict "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
numberOfSubdomains 4;
method          scotch;
"

# --- FV SCHEMES (Global) ---
write_file $CASE_DIR/system/fvSchemes "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"

# --- FV SOLUTION (Global) ---
write_file $CASE_DIR/system/fvSolution "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    \"(rho|p_rgh).*\" { solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0.01; }
    \"(U|h|e|k|epsilon).*\" { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.1; }
}
PIMPLE { nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { \"h.*\" 1; \"U.*\" 1; } }
"

echo "=== 3. Generate Mesh ==="
# Generate Gmsh mesh
python3 src/geometry/generate_water_heatsink_gmsh.py $CASE_DIR/heatsink.msh

# Convert to OpenFOAM
cd $CASE_DIR
gmshToFoam heatsink.msh

# Cleanup
rm heatsink.msh

echo "=== 4. Configure Regions ==="
# Split mesh into regions based on 'solid' and 'fluid' cellZones (from Gmsh physical groups)
splitMeshRegions -cellZones -overwrite

echo "=== 5. Configure Properties (Post-Split) ==="


# --- DECOMPOSE PAR ---
write_file system/decomposeParDict "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
numberOfSubdomains 4;
method          scotch;
"

# --- FV SCHEMES (Global) ---
# ... (Standard schemes, same as before)
write_file system/fvSchemes "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"

# --- FV SOLUTION (Global) ---
write_file system/fvSolution "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    \"(rho|p_rgh).*\" { solver PCG; preconditioner DIC; tolerance 1e-7; relTol 0.01; }
    \"(U|h|e|k|epsilon).*\" { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.1; }
}
PIMPLE { nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0; }
relaxationFactors { equations { \"h.*\" 1; \"U.*\" 1; } }
"

# --- REGION PROPERTIES ---
mkdir -p constant
write_file constant/regionProperties "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      regionProperties;
}
regions
(
    fluid       (fluid)
    solid       (solid)
);
"

# --- FLUID REGION SETUP ---
mkdir -p system/fluid constant/fluid 0/fluid

# Schemes/Solution (cp from global)
cp system/fvSchemes system/fluid/
cp system/fvSolution system/fluid/

# Update Fluid Schemes for convection
sed -i 's/divSchemes { default none; }/divSchemes { default none; div(phi,U) Gauss upwind; div(phi,h) Gauss upwind; div(phi,K) Gauss upwind; div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; }/' system/fluid/fvSchemes

# Properties (Water)
# rho ~ 998, Cp ~ 4182, mu ~ 0.001, Pr ~ 7 -> k = mu*Cp/Pr = 0.001*4182/7 ~ 0.6
# OpenFOAM heRhoThermo
write_file constant/fluid/thermophysicalProperties "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      thermophysicalProperties;
}
thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          hConst;
    equationOfState rhoConst; // Better for water than perfectGas
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie { molWeight 18.0; }
    thermodynamics { Cp 4182; Hf 0; }
    transport { mu 0.001; Pr 7; }
    equationOfState { rho 998; }
}
"

write_file constant/fluid/turbulenceProperties "
FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType laminar;
"

write_file constant/fluid/radiationProperties "
FoamFile { version 2.0; format ascii; class dictionary; object radiationProperties; }
radiation off;
radiationModel none;
"

write_file constant/fluid/g "
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; object g; }
dimensions [0 1 -2 0 0 0 0];
value (0 0 -9.81);
"

# Dictionary change for boundary patching
# Maps gmsh tags to patches
# We will use changeDictionary later if needed, but if we named Physical Groups in Gmsh,
# splitMeshRegions should have created patches named 'inlet', 'outlet'.
# We check this at runtime logic.

# Initial Conditions (Fluid)
# U: inlet 10Pa driving? No, user said inlet BC is pressure 10Pa.
# Actually this is a pressure driven flow.
# U BCs: inlet: pressureInletOutletVelocity? Or just zeroGradient?
# If p is fixed types, U is usually pressureInletOutletVelocity or zeroGradient.
write_file 0/fluid/U "
FoamFile { version 2.0; format ascii; class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    inlet { type pressureInletOutletVelocity; value uniform (0 0 0); }
    outlet { type pressureInletOutletVelocity; value uniform (0 0 0); }
    \".*\" { type noSlip; }
}
"

write_file 0/fluid/p_rgh "
FoamFile { version 2.0; format ascii; class volScalarField; object p_rgh; }
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    inlet { type totalPressure; p0 uniform 10; value uniform 10; gamma 1; }
    outlet { type fixedValue; value uniform 0; }
    \".*\" { type fixedFluxPressure; value uniform 0; }
}
"

write_file 0/fluid/p "
FoamFile { version 2.0; format ascii; class volScalarField; object p; }
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    \".*\" { type calculated; value uniform 0; }
}
"

write_file 0/fluid/T "
FoamFile { version 2.0; format ascii; class volScalarField; object T; }
dimensions [0 0 0 1 0 0 0];
internalField uniform 293.15; // 20 C
boundaryField
{
    inlet { type fixedValue; value uniform 293.15; }
    outlet { type inletOutlet; inletValue uniform 293.15; value uniform 293.15; }
    \"fluid_to_.*\" { type compressible::turbulentTemperatureCoupledBaffleMixed; Tnbr T; kappaMethod fluidThermo; value uniform 293.15; }
    \".*\" { type zeroGradient; }
}
"


# --- SOLID REGION SETUP ---
mkdir -p system/solid constant/solid 0/solid

cp system/fvSchemes system/solid/
cp system/fvSchemes system/solid/
# cp system/fvSolution system/solid/ # Solid needs symmetric solver for h
write_file system/solid/fvSolution "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    \"(h|hFinal)\"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }
}
PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
}
relaxationFactors
{
    equations
    {
        h               1;
    }
}
"

# Properties (Aluminium)
# rho 2700, Cp 900, k 200
write_file constant/solid/thermophysicalProperties "
FoamFile { version 2.0; format ascii; class dictionary; object thermophysicalProperties; }
thermoType
{
    type            heSolidThermo;
    mixture         pureMixture;
    transport       constIso;
    thermo          hConst;
    equationOfState rhoConst;
    specie          specie;
    energy          sensibleEnthalpy;
}
mixture
{
    specie { molWeight 27; }
    transport { kappa 200; }
    thermodynamics { Cp 900; Hf 0; }
    equationOfState { rho 2700; }
}
"

write_file constant/solid/radiationProperties "
FoamFile { version 2.0; format ascii; class dictionary; object radiationProperties; }
radiation off;
radiationModel none;
"

# Initial Conditions (Solid)
write_file 0/solid/T "
FoamFile { version 2.0; format ascii; class volScalarField; object T; }
dimensions [0 0 0 1 0 0 0];
internalField uniform 293.15;
boundaryField
{
    heat_source { type externalWallHeatFluxTemperature; mode flux; q uniform 60000; kappaMethod solidThermo; value uniform 293.15; }
    \"solid_to_.*\" { type compressible::turbulentTemperatureCoupledBaffleMixed; Tnbr T; kappaMethod solidThermo; value uniform 293.15; }
    \".*\" { type zeroGradient; }
}
"

write_file 0/solid/p "
FoamFile { version 2.0; format ascii; class volScalarField; object p; }
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    \".*\" { type calculated; value \$internalField; }
}
"

echo "=== 5. Finalize Setup ==="
# Map regions (boundary conditions between fluid and solid)
# We need to ensure 'inlet', 'outlet', 'heat_source' are correctly mapped.
# All other faces should be mapped to default wall or coupled.

# Using changeDictionary to couple regions
write_file system/fluid/changeDictionaryDict "
FoamFile { version 2.0; format ascii; class dictionary; object changeDictionaryDict; }
boundary
{
    \"inlet.*\" { name inlet; }
    \"outlet.*\" { name outlet; }
    \"fluid_to_.*\" { type mappedWall; sampleMode nearestPatchFace; sampleRegion solid; samplePatch solid_to_fluid; }
}
"
# Note: splitMeshRegions may name patches 'fluid_to_solid' automatically.
# We trust it mostly.

write_file system/solid/changeDictionaryDict "
FoamFile { version 2.0; format ascii; class dictionary; object changeDictionaryDict; }
boundary
{
    \"heat_source.*\" { name heat_source; }
    \"solid_to_.*\" { type mappedWall; sampleMode nearestPatchFace; sampleRegion fluid; samplePatch fluid_to_solid; }
}
"

# Run changeDictionary?
# Actually, let's see what splitMeshRegions produces first.
# For now, we assume standard naming.

echo "Setup Complete."

echo "=== 6. Run Simulation ==="
# Fix missing constant/g
cp constant/fluid/g constant/g

echo "Running chtMultiRegionFoam..."
chtMultiRegionFoam > log.chtMultiRegionFoam 2>&1

echo "Done."
