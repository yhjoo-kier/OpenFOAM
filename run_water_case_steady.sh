#!/bin/bash
# Steady-State Water-cooled CHT simulation run script
# Uses chtMultiRegionSimpleFoam (SIMPLE algorithm)
set -e

# Setup environment
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc

# Source case (Transient)
SRC_CASE="cases/heatsink_water_cht"
# New case (Steady)
CASE_DIR="cases/heatsink_water_cht_steady"

echo "=== 1. Setup Steady Case Directory ==="
if [ ! -d "$SRC_CASE" ]; then
    echo "Error: Source case $SRC_CASE does not exist. Run run_water_case.sh first."
    exit 1
fi

rm -rf $CASE_DIR
cp -r $SRC_CASE $CASE_DIR

# Clean up results from source case (keep 0 and constant and system)
find $CASE_DIR -maxdepth 1 -name '[0-9]*' -not -name '0' -exec rm -rf {} +
rm -rf $CASE_DIR/VTK
rm -f $CASE_DIR/log.*

echo "=== 2. Modify Configuration for Steady State ==="

# --- UTILITIES ---
write_file() {
    cat > $1 <<EOF
$2
EOF
}

# 1. Update controlDict
# Switch application to chtMultiRegionSimpleFoam
# endTime becomes iteration count (e.g., 2000)
write_file $CASE_DIR/system/controlDict "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     chtMultiRegionSimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         2000;
deltaT          1;
writeControl    timeStep;
writeInterval   100;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
"

# 2. Update fvSchemes (Global)
# ddtSchemes must be steadyState
write_file $CASE_DIR/system/fvSchemes "
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default bounded Gauss upwind; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"
# Copy global schemes to regions
cp $CASE_DIR/system/fvSchemes $CASE_DIR/system/fluid/fvSchemes
cp $CASE_DIR/system/fvSchemes $CASE_DIR/system/solid/fvSchemes

# Update Fluid divSchemes specific
# Add bounded upwind for stability in steady RANS
write_file $CASE_DIR/system/fluid/fvSchemes "
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U)      bounded Gauss upwind;
    div(phi,h)      bounded Gauss upwind;
    div(phi,K)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,k)      bounded Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"

# 3. Update fvSolution (Global & Regions)
# Must use SIMPLE control instead of PIMPLE
# Add relaxationFactors for steady convergence
SIMPLE_SOLUTION="
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    \"(p|p_rgh|rho).*\"
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
    }
    \"(U|h|e|k|epsilon).*\"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }
    \"(h|hFinal)\" 
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
}
relaxationFactors
{
    fields
    {
        p_rgh           0.3;
        rho             1.0;
    }
    equations
    {
        U               0.7;
        h               0.7;
        k               0.7;
        epsilon         0.7;
    }
}
"


# 4. Fluid Solution (Asymmetric h)
FLUID_SOLUTION="
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
solvers
{
    \"(p|p_rgh|rho).*\"
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
    }
    \"(U|h|e|k|epsilon).*\"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-7;
        relTol          0.1;
    }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
}
relaxationFactors
{
    fields
    {
        p_rgh           0.3;
        rho             1.0;
    }
    equations
    {
        U               0.7;
        h               0.7;
        k               0.7;
        epsilon         0.7;
    }
}
"
write_file $CASE_DIR/system/fluid/fvSolution "$FLUID_SOLUTION"

# 5. Solid Solution (Symmetric h)
SOLID_SOLUTION="
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
SIMPLE
{
    nNonOrthogonalCorrectors 0;
}
relaxationFactors
{
    equations
    {
        h               0.7;
    }
}
"
write_file $CASE_DIR/system/solid/fvSolution "$SOLID_SOLUTION"


echo "=== 3. Finalize Steady Setup ==="

# Run Script
echo "=== 4. Run Steady Simulation ==="
cd $CASE_DIR
echo "Running chtMultiRegionSimpleFoam..."
# Ensure g exists
if [ ! -f constant/g ]; then
    cp constant/fluid/g constant/g
fi

chtMultiRegionSimpleFoam > log.chtMultiRegionSimpleFoam 2>&1

echo "Done. Converting results..."
foamToVTK -region fluid
foamToVTK -region solid

echo "Steady state analysis complete."
