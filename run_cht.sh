#!/bin/bash
# CHT simulation run script

set -e

export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc

CASE_DIR="cases/heatsink_cht_full"

echo "=== Setting up CHT case ==="
cd /home/user/OpenFOAM
rm -rf $CASE_DIR
python src/geometry/heatsink_cht_mesh.py $CASE_DIR --velocity 0.5 --heat-flux 5000

cd $CASE_DIR

# Fix controlDict header
cat > system/controlDict << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     chtMultiRegionFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.5;
deltaT          0.001;
writeControl    adjustableRunTime;
writeInterval   0.1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
maxCo           0.5;
maxDi           10;
adjustTimeStep  yes;

// ************************************************************************* //
EOF

echo "=== Running blockMesh ==="
blockMesh

echo "=== Running topoSet ==="
topoSet

echo "=== Running splitMeshRegions ==="
splitMeshRegions -cellZones -overwrite

# Copy g file
cp constant/fluid/g constant/g

# Create solid p file
mkdir -p 0/solid
cat > 0/solid/p << 'EOF'
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
dimensions      [1 -1 -2 0 0 0 0];
internalField   uniform 1e5;
boundaryField
{
    ".*"
    {
        type            calculated;
        value           $internalField;
    }
}
EOF

echo "=== Running changeDictionary ==="
changeDictionary -region fluid
changeDictionary -region solid

echo "=== Running chtMultiRegionFoam ==="
chtMultiRegionFoam

echo "=== Converting to VTK ==="
foamToVTK -region fluid
foamToVTK -region solid

echo "=== Done! ==="
