#!/usr/bin/env python3
"""
MCP Server for Grid Independence Study.

This allows other AI agents to use the grid study framework as a tool
through the Model Context Protocol.

Usage:
    Run as MCP server:
        python -m grid_study.mcp_server

    Configure in Claude Desktop or other MCP clients:
        {
            "mcpServers": {
                "grid-study": {
                    "command": "python",
                    "args": ["-m", "grid_study.mcp_server"],
                    "cwd": "/workspaces/OpenFOAM/src"
                }
            }
        }
"""

import json
import sys
from pathlib import Path
from typing import Any

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from grid_study import GridStudyConfig, GridStudy, MeshLevel, run_grid_study


def handle_run_grid_study(params: dict) -> dict:
    """
    Run a grid independence study.

    Parameters:
        base_case (str): Path to base OpenFOAM case (required)
        adaptive (bool): Enable adaptive mode (default: True)
        threshold (float): Convergence threshold % (default: 1.0)
        max_cells (int): Maximum cells limit (default: 1000000)
        max_levels (int): Maximum levels (default: 8)
        metric_patch (str): Patch to monitor (default: "heat_source")
        metric_field (str): Field to extract (default: "T")
        metric_region (str): Region name (default: "solid")

    Returns:
        dict with: converged, stop_reason, levels, recommended_level, extrapolated_value
    """
    try:
        base_case = params.get("base_case")
        if not base_case:
            return {"error": "base_case is required"}

        # Validate path exists
        case_path = Path(base_case)
        if not case_path.exists():
            return {"error": f"Case path does not exist: {base_case}"}

        # Run study with adaptive mode
        analysis = run_grid_study(
            base_case=base_case,
            output_dir=params.get("output_dir", "grid_study_output"),
            study_name=params.get("study_name", "grid_study"),
            metric_patch=params.get("metric_patch", "heat_source"),
            metric_field=params.get("metric_field", "T"),
            metric_region=params.get("metric_region", "solid"),
            threshold=params.get("threshold", 1.0),
            adaptive=params.get("adaptive", True),
            max_cells=params.get("max_cells", 1_000_000),
            max_levels=params.get("max_levels", 8),
        )

        # Format results for AI consumption
        levels = []
        for i, lr in enumerate(analysis.level_results):
            level_info = {
                "name": lr.level_name,
                "cells": lr.num_cells,
                "metric_value": lr.metric_value,
            }
            # Add convergence info for levels after first
            if i > 0 and i - 1 < len(analysis.convergence_results):
                cr = analysis.convergence_results[i - 1]
                level_info["change_percent"] = cr.percent_change
                level_info["converged"] = cr.converged
            levels.append(level_info)

        return {
            "success": True,
            "converged": analysis.is_converged,
            "stop_reason": analysis.stop_reason,
            "levels": levels,
            "recommended_level": analysis.recommended_level,
            "extrapolated_value": analysis.extrapolated_value,
            "gci_fine": analysis.gci_fine,
        }

    except Exception as e:
        return {"error": str(e)}


def handle_check_case(params: dict) -> dict:
    """
    Check if a case is valid for grid study.

    Parameters:
        case_path (str): Path to OpenFOAM case

    Returns:
        dict with: valid, has_geometry_script, solver_type, regions
    """
    try:
        case_path = Path(params.get("case_path", ""))
        if not case_path.exists():
            return {"valid": False, "error": "Path does not exist"}

        # Check for required directories
        required = ["0", "constant", "system"]
        missing = [d for d in required if not (case_path / d).exists()]
        if missing:
            return {"valid": False, "error": f"Missing directories: {missing}"}

        # Check for geometry script
        geometry_script = None
        for script_name in ["generate_mesh.py", "create_geometry.py"]:
            script_path = case_path / "scripts" / script_name
            if script_path.exists():
                geometry_script = str(script_path)
                break

        # Check controlDict for solver
        control_dict = case_path / "system" / "controlDict"
        solver_type = "unknown"
        if control_dict.exists():
            content = control_dict.read_text()
            if "chtMultiRegionSimpleFoam" in content:
                solver_type = "chtMultiRegionSimpleFoam"
            elif "simpleFoam" in content:
                solver_type = "simpleFoam"

        # Check for regions
        regions = []
        constant_dir = case_path / "constant"
        for item in constant_dir.iterdir():
            if item.is_dir() and (item / "polyMesh").exists():
                regions.append(item.name)

        return {
            "valid": True,
            "has_geometry_script": geometry_script is not None,
            "geometry_script": geometry_script,
            "solver_type": solver_type,
            "regions": regions if regions else ["single_region"],
        }

    except Exception as e:
        return {"valid": False, "error": str(e)}


# Tool definitions for MCP
TOOLS = {
    "run_grid_study": {
        "description": """Run an automated grid independence study for OpenFOAM CFD cases.

This tool automatically:
1. Generates progressively finer meshes
2. Runs simulations at each mesh level
3. Extracts the monitored metric (e.g., temperature)
4. Checks for convergence (change < threshold)
5. Continues refining until converged or limits reached

Returns convergence status, recommended mesh level, and Richardson extrapolated value.""",
        "parameters": {
            "type": "object",
            "properties": {
                "base_case": {
                    "type": "string",
                    "description": "Path to the base OpenFOAM case directory",
                },
                "adaptive": {
                    "type": "boolean",
                    "description": "Enable adaptive refinement mode (recommended: true)",
                    "default": True,
                },
                "threshold": {
                    "type": "number",
                    "description": "Convergence threshold in percent (e.g., 1.0 for 1%)",
                    "default": 1.0,
                },
                "max_cells": {
                    "type": "integer",
                    "description": "Maximum number of cells before stopping",
                    "default": 1000000,
                },
                "max_levels": {
                    "type": "integer",
                    "description": "Maximum number of refinement levels",
                    "default": 8,
                },
                "metric_patch": {
                    "type": "string",
                    "description": "OpenFOAM patch name to monitor",
                    "default": "heat_source",
                },
                "metric_field": {
                    "type": "string",
                    "description": "Field to extract (T, p, U, etc.)",
                    "default": "T",
                },
                "metric_region": {
                    "type": "string",
                    "description": "Region name for multi-region cases",
                    "default": "solid",
                },
            },
            "required": ["base_case"],
        },
        "handler": handle_run_grid_study,
    },
    "check_grid_study_case": {
        "description": "Check if an OpenFOAM case is valid for grid study and get its configuration.",
        "parameters": {
            "type": "object",
            "properties": {
                "case_path": {
                    "type": "string",
                    "description": "Path to the OpenFOAM case directory",
                },
            },
            "required": ["case_path"],
        },
        "handler": handle_check_case,
    },
}


def main():
    """Run as MCP server (stdio transport)."""
    # Simple MCP-like protocol over stdin/stdout
    # For full MCP implementation, use the official MCP SDK

    print(json.dumps({
        "name": "grid-study",
        "version": "0.2.0",
        "description": "Grid Independence Study framework for OpenFOAM CFD simulations",
        "tools": {
            name: {
                "description": tool["description"],
                "parameters": tool["parameters"],
            }
            for name, tool in TOOLS.items()
        },
    }), flush=True)

    for line in sys.stdin:
        try:
            request = json.loads(line.strip())
            tool_name = request.get("tool")
            params = request.get("params", {})

            if tool_name in TOOLS:
                result = TOOLS[tool_name]["handler"](params)
            else:
                result = {"error": f"Unknown tool: {tool_name}"}

            print(json.dumps(result), flush=True)

        except json.JSONDecodeError:
            print(json.dumps({"error": "Invalid JSON"}), flush=True)
        except Exception as e:
            print(json.dumps({"error": str(e)}), flush=True)


if __name__ == "__main__":
    main()
