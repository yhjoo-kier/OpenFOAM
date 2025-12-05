#!/bin/bash
#
# OpenFOAM CFD Pipeline - VM Setup Script
# For Ubuntu 24.04 LTS (Noble Numbat)
#
# This script sets up the complete development environment for the OpenFOAM CFD pipeline.
# It installs all required dependencies including OpenFOAM, Gmsh, and Python packages.
#
# Usage:
#   chmod +x setup_vm.sh
#   sudo ./setup_vm.sh
#
# After running this script, you can:
#   1. Run the heatsink example:
#      python3 src/geometry/simple_heatsink.py cases/my_heatsink
#      cd cases/my_heatsink && blockMesh && laplacianFoam && foamToVTK
#
#   2. Run the channel flow example:
#      cd channelFlow && blockMesh && icoFoam && foamToVTK
#

set -e  # Exit on error

echo "========================================"
echo "OpenFOAM CFD Pipeline - VM Setup"
echo "========================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if running as root
if [ "$EUID" -ne 0 ]; then
    print_error "Please run as root (use sudo)"
    exit 1
fi

# Get the original user if running with sudo
ORIG_USER=${SUDO_USER:-$USER}
ORIG_HOME=$(eval echo ~$ORIG_USER)

print_status "Setting up for user: $ORIG_USER"
print_status "Home directory: $ORIG_HOME"

# ============================================
# Step 1: Update system packages
# ============================================
print_status "Step 1/5: Updating system packages..."
apt-get update -qq

# ============================================
# Step 2: Install OpenFOAM and Gmsh
# ============================================
print_status "Step 2/5: Installing OpenFOAM and Gmsh..."
apt-get install -y -qq \
    openfoam \
    gmsh \
    build-essential \
    libosmesa6-dev \
    libglx-mesa0 \
    python3 \
    python3-pip \
    python3-venv

# Verify OpenFOAM installation
if ! command -v blockMesh &> /dev/null; then
    print_error "OpenFOAM installation failed - blockMesh not found"
    exit 1
fi

if ! command -v gmsh &> /dev/null; then
    print_error "Gmsh installation failed"
    exit 1
fi

print_status "OpenFOAM and Gmsh installed successfully"

# ============================================
# Step 3: Set up OpenFOAM environment variables
# ============================================
print_status "Step 3/5: Configuring OpenFOAM environment variables..."

# Create environment file
OPENFOAM_ENV_FILE="/etc/profile.d/openfoam.sh"
cat > "$OPENFOAM_ENV_FILE" << 'EOF'
# OpenFOAM environment variables
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM
EOF

chmod +x "$OPENFOAM_ENV_FILE"

# Source for current session
export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

print_status "OpenFOAM environment configured"

# ============================================
# Step 4: Install Python packages
# ============================================
print_status "Step 4/5: Installing Python packages..."

pip3 install --quiet \
    gmsh \
    pyvista \
    matplotlib \
    numpy \
    vtk

# Verify Python packages
python3 -c "import gmsh, pyvista, matplotlib, numpy, vtk; print('Python packages verified')" || {
    print_error "Python package verification failed"
    exit 1
}

print_status "Python packages installed successfully"

# ============================================
# Step 5: Verify installation with test run
# ============================================
print_status "Step 5/5: Verifying installation with test run..."

# Get the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Check if we're in the right directory
if [ -f "$PROJECT_DIR/src/geometry/simple_heatsink.py" ]; then
    print_status "Running verification test..."

    # Create test case
    TEST_DIR=$(mktemp -d)
    python3 "$PROJECT_DIR/src/geometry/simple_heatsink.py" "$TEST_DIR/verify_test" --heat-flux 1000 2>/dev/null

    # Run blockMesh
    cd "$TEST_DIR/verify_test"
    if blockMesh > /dev/null 2>&1; then
        print_status "blockMesh: OK"
    else
        print_warning "blockMesh test failed"
    fi

    # Run laplacianFoam (quick test with endTime=10)
    sed -i 's/endTime         1000/endTime         10/' system/controlDict
    if laplacianFoam > /dev/null 2>&1; then
        print_status "laplacianFoam: OK"
    else
        print_warning "laplacianFoam test failed"
    fi

    # Cleanup test directory
    rm -rf "$TEST_DIR"

    cd "$PROJECT_DIR"
else
    print_warning "Project files not found, skipping verification test"
fi

# ============================================
# Summary
# ============================================
echo ""
echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""
echo "Installed components:"
echo "  - OpenFOAM (v1912)"
echo "  - Gmsh (v4.12)"
echo "  - Python packages: gmsh, pyvista, matplotlib, numpy, vtk"
echo ""
echo "OpenFOAM environment variables are set in: $OPENFOAM_ENV_FILE"
echo ""
echo "Quick start:"
echo "  1. Heatsink heat conduction example:"
echo "     python3 src/geometry/simple_heatsink.py cases/my_heatsink"
echo "     cd cases/my_heatsink && blockMesh && laplacianFoam && foamToVTK"
echo ""
echo "  2. Channel flow example:"
echo "     cd channelFlow && blockMesh && icoFoam && foamToVTK"
echo ""
echo "For new terminal sessions, run: source /etc/profile.d/openfoam.sh"
echo ""
