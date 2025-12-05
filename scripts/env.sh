#!/bin/bash
#
# OpenFOAM environment variables
# Usage: source scripts/env.sh
#

export WM_PROJECT_DIR=/usr/share/openfoam
export FOAM_ETC=/usr/share/openfoam/etc
export WM_PROJECT=OpenFOAM

# Add Python path for this project
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/..:${PYTHONPATH}"

echo "OpenFOAM environment loaded"
