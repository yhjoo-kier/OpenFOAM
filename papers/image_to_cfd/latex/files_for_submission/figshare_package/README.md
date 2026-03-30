# Benchmark Dataset and Code for "Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction"

## Description

This dataset accompanies the paper "Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction" submitted to Building and Environment. It contains the benchmark dataset, multi-view rendered input images, evaluation code, and pipeline source code needed to reproduce the results reported in the paper.

A rule-based benchmark of 20 indoor geometries spanning four complexity levels (A1: rectangular simple, A2: rectangular dense, A3: composite simple, A4: composite dense) is provided as ground-truth scene JSON files. Each geometry is rendered into five input view types (perspective, bird's-eye, floor plan, wireframe, section), yielding 100 evaluation images. The pipeline code converts a single 2D image into a steady-state CFD solution via VLM geometric abstraction and automated OpenFOAM meshing/solving.

## Contents

- `benchmark_scenes/` — 20 ground-truth scene JSON files defining room geometry, obstacles, and openings
- `benchmark_renderings/` — 100 multi-view input images (20 geometries × 5 view types: perspective, birdseye, floorplan, wireframe, section)
- `pipeline_code/` — Image-to-CFD pipeline entry point (requires OpenFOAM, Gmsh, and Gemini API access)
- `evaluation_code/` — Benchmark evaluation and aggregation scripts
- `vlm_prompt/` — VLM prompt template used for scene extraction

## Requirements

- OpenFOAM (v2312 or later)
- Gmsh (4.x)
- Python 3.10+
- Google Gemini API access

## License

CC-BY 4.0

## Related Publication

Joo, Y. (2026). Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction. Building and Environment. [DOI to be inserted upon acceptance]

## Funding

Korea Institute of Energy Research (C6-2419-63)
