# OpenFOAM CFD Pipeline Docker Image
# Based on Ubuntu 24.04 with OpenFOAM, Gmsh, and Python visualization tools

FROM ubuntu:22.04

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
    sudo \
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

# Create non-root user 'ubuntu' with sudo privileges
ARG USERNAME=ubuntu
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set OpenFOAM environment variables
ENV WM_PROJECT_DIR=/usr/share/openfoam
ENV FOAM_ETC=/usr/share/openfoam/etc
ENV WM_PROJECT=OpenFOAM
ENV PATH="/usr/bin:${PATH}"

# Create app directory and set ownership
WORKDIR /app
RUN chown -R ubuntu:ubuntu /app

# Switch to ubuntu user for uv and claude code installation
USER ubuntu

# Install uv for Python package management (installs to ~/.local/bin)
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/home/ubuntu/.local/bin:${PATH}"

# Install Claude Code CLI (installs to ~/.local/bin)
RUN curl -fsSL https://claude.ai/install.sh | bash

# Create virtual environment and install Python packages
RUN uv venv /app/.venv
ENV PATH="/app/.venv/bin:${PATH}"
ENV VIRTUAL_ENV="/app/.venv"

RUN uv pip install \
    gmsh \
    pyvista \
    matplotlib \
    numpy \
    vtk \
    scipy

# Copy project files with correct ownership
COPY --chown=ubuntu:ubuntu src/ /app/src/
COPY --chown=ubuntu:ubuntu scripts/ /app/scripts/
COPY --chown=ubuntu:ubuntu cases/ /app/cases/

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
