#!/bin/bash
#
# OpenFOAM CFD Pipeline - Installation Verification Script
# Run this script to verify that all components are properly installed
#
# Usage:
#   ./scripts/verify_installation.sh
#

# Don't use set -e as arithmetic operations return non-zero on 0

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PASS=0
FAIL=0

check() {
    local name="$1"
    local cmd="$2"
    printf "%-40s" "Checking $name..."
    if eval "$cmd" > /dev/null 2>&1; then
        echo -e "${GREEN}OK${NC}"
        PASS=$((PASS + 1))
    else
        echo -e "${RED}FAIL${NC}"
        FAIL=$((FAIL + 1))
    fi
}

echo "========================================"
echo "OpenFOAM CFD Pipeline - Verification"
echo "========================================"
echo ""

# Source OpenFOAM environment
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

echo "System Tools:"
check "blockMesh" "command -v blockMesh"
check "laplacianFoam" "command -v laplacianFoam"
check "icoFoam" "command -v icoFoam"
check "foamToVTK" "command -v foamToVTK"
check "gmsh" "command -v gmsh"
check "python3" "command -v python3"

echo ""
echo "Python Packages:"
check "gmsh (Python)" "python3 -c 'import gmsh'"
check "pyvista" "python3 -c 'import pyvista'"
check "matplotlib" "python3 -c 'import matplotlib'"
check "numpy" "python3 -c 'import numpy'"
check "vtk" "python3 -c 'import vtk'"

echo ""
echo "Project Files:"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
check "simple_heatsink.py" "test -f $PROJECT_DIR/src/geometry/simple_heatsink.py"
check "run_pipeline.py" "test -f $PROJECT_DIR/scripts/run_pipeline.py"
check "channelFlow case" "test -d $PROJECT_DIR/cases/channelFlow"

echo ""
echo "========================================"
echo "Results: ${PASS} passed, ${FAIL} failed"
echo "========================================"

if [ $FAIL -gt 0 ]; then
    echo -e "${RED}Some checks failed. Please run setup_vm.sh to install missing components.${NC}"
    exit 1
else
    echo -e "${GREEN}All checks passed! The environment is ready.${NC}"
    exit 0
fi
