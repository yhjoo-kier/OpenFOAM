# OpenFOAM CFD Pipeline Docker Image
# Based on Ubuntu 24.04 with OpenFOAM, Gmsh, and Python visualization tools

FROM ubuntu:24.04

LABEL maintainer="OpenFOAM Pipeline"
LABEL description="CFD simulation pipeline with Gmsh, OpenFOAM, and visualization tools"

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Build tools
    build-essential \
    cmake \
    git \
    curl \
    wget \
    # OpenFOAM and CFD tools
    openfoam \
    gmsh \
    # Python
    python3 \
    python3-pip \
    python3-venv \
    # Visualization dependencies
    libgl1-mesa-glx \
    libglu1-mesa \
    libosmesa6-dev \
    libxrender1 \
    libxcursor1 \
    libxft2 \
    libxinerama1 \
    # Cleanup
    && rm -rf /var/lib/apt/lists/*

# Set OpenFOAM environment variables
ENV WM_PROJECT_DIR=/usr/share/openfoam
ENV FOAM_ETC=/usr/share/openfoam/etc
ENV WM_PROJECT=OpenFOAM
ENV PATH="/usr/bin:${PATH}"

# Create app directory
WORKDIR /app

# Install uv for Python package management
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Create virtual environment and install Python packages
RUN uv venv /app/.venv
ENV PATH="/app/.venv/bin:${PATH}"
ENV VIRTUAL_ENV="/app/.venv"

RUN uv pip install \
    gmsh \
    pyvista \
    matplotlib \
    numpy \
    vtk

# Copy project files
COPY src/ /app/src/
COPY scripts/ /app/scripts/
COPY cases/ /app/cases/

# Create output directory
RUN mkdir -p /app/results

# Set Python path
ENV PYTHONPATH="/app:${PYTHONPATH}"

# Default command
CMD ["python", "scripts/run_pipeline.py", "--help"]

# Usage examples:
# Build:
#   docker build -t openfoam-pipeline .
#
# Run interactive:
#   docker run -it --rm -v $(pwd)/results:/app/results openfoam-pipeline bash
#
# Run pipeline:
#   docker run -it --rm -v $(pwd)/results:/app/results openfoam-pipeline \
#       python scripts/run_pipeline.py --case-name heatsink
#
# Run with custom parameters:
#   docker run -it --rm -v $(pwd)/results:/app/results openfoam-pipeline \
#       python scripts/run_pipeline.py --heat-flux 10000 --end-time 500
