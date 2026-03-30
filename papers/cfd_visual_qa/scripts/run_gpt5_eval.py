#!/usr/bin/env python3
"""
GPT-5 (Codex CLI) blind evaluation of CFD flow field images.
Two modes: setup-conditioned and image-only, 30 items each.
"""

import json
import subprocess
import os
import random
import re
import time
from pathlib import Path

BASE = Path("/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa")
BENCHMARK = BASE / "benchmark"
IMAGES = BENCHMARK / "images"
CASES = BENCHMARK / "cases"
LABELS = BENCHMARK / "labels"

# Load the 30 items from Claude image-only eval
with open(LABELS / "vlm_eval_claude_image_only.json") as f:
    claude_data = json.load(f)

items = claude_data["items"]

# Temperature scenarios use V2, others use V6
TEMP_SCENARIOS = {"S1", "S6", "S9"}

def get_image_path(case_name):
    scenario = case_name.split("_")[0]
    if scenario in TEMP_SCENARIOS:
        return IMAGES / case_name / f"{case_name}_V2.png"
    return IMAGES / case_name / f"{case_name}_V6.png"

def get_setup_text(case_name):
    meta_path = CASES / case_name / "case_meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    # Build a setup description from metadata, without revealing case name
    parts = []
    if "scenario_name" in meta:
        name_map = {
            "lid_driven_cavity": "Lid-driven cavity flow",
            "natural_convection_tall": "Natural convection in a tall rectangular cavity (aspect ratio 2:1), heated left wall, cooled right wall",
            "channel_flow_laminar": "Laminar channel flow (10:1 aspect ratio)",
            "turbulent_channel_flow": "Turbulent channel flow (10:1 aspect ratio)",
            "backward_facing_step": "Backward-facing step flow",
            "differentially_heated_cavity": "Differentially heated square cavity (natural convection)",
            "mixed_convection_channel": "Mixed convection in a horizontal channel with bottom heating",
            "obstacle_channel_flow": "Channel flow past internal obstacles",
            "rayleigh_benard": "Rayleigh-Benard convection in a wide cavity (4:1 aspect ratio), heated from below",
            "ventilation_room": "Ventilation flow in a rectangular room (9x3 m) with inlet on upper-left wall",
        }
        sn = meta.get("scenario_name", "")
        parts.append(name_map.get(sn, sn))

    if "Re" in meta and meta["Re"] is not None:
        parts.append(f"Re = {meta['Re']}")
    if "Ra" in meta and meta.get("Ra") is not None:
        parts.append(f"Ra = {meta['Ra']}")
    if "Ri" in meta and meta.get("Ri") is not None:
        parts.append(f"Ri = {meta['Ri']}")
    if meta.get("turbulence"):
        parts.append("Turbulent flow (k-epsilon model)")
    else:
        parts.append("Laminar flow")

    # Add boundary conditions if available
    if "U_lid" in meta and meta.get("U_lid") is not None:
        parts.append(f"Lid velocity = {meta['U_lid']} m/s")
    if "U_inlet" in meta and meta.get("U_inlet") is not None:
        parts.append(f"Inlet velocity = {meta['U_inlet']} m/s")
    if "T_hot" in meta and meta.get("T_hot") is not None:
        parts.append(f"Hot wall T = {meta['T_hot']} K")
    if "T_cold" in meta and meta.get("T_cold") is not None:
        parts.append(f"Cold wall T = {meta['T_cold']} K")

    desc = meta.get("description", "")
    if desc:
        parts.append(f"Description: {desc}")

    return ". ".join(parts)

def generate_blind_codes(n):
    codes = set()
    while len(codes) < n:
        codes.add(f"CFD-{random.randint(0, 0xFFFF):04X}")
    return list(codes)

def run_codex(prompt, timeout=180):
    """Run codex exec and return the output."""
    try:
        result = subprocess.run(
            ["codex", "exec", prompt, "--full-auto"],
            capture_output=True, text=True, timeout=timeout,
            env={**os.environ, "TERM": "dumb"}
        )
        output = result.stdout.strip()
        if not output and result.stderr:
            output = result.stderr.strip()
        return output
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception as e:
        return f"ERROR: {e}"

def extract_codex_response(raw_output):
    """Extract the actual model response from codex CLI output.

    The codex CLI output has headers, then blocks prefixed with 'codex\\n'.
    The actual response is in the last codex block, before 'tokens used'.
    """
    if not raw_output or raw_output in ("TIMEOUT",) or raw_output.startswith("ERROR"):
        return raw_output

    # Split by 'tokens used' to get content before token count
    parts = raw_output.split("tokens used")
    content = parts[0]

    # Find all 'codex\n' blocks - the response is after the last one
    segments = content.split("\ncodex\n")
    if len(segments) > 1:
        response = segments[-1].strip()
    else:
        # Fallback: try to find content after 'viewed image' line
        lines = content.split("\n")
        response_lines = []
        found_image = False
        for line in lines:
            if "viewed image" in line:
                found_image = True
                continue
            if found_image and line.strip() and line.strip() != "codex":
                response_lines.append(line)
        response = "\n".join(response_lines).strip()

    return response if response else raw_output

def parse_verdict(raw_output):
    """Extract OK or Anomaly from codex output."""
    if not raw_output or raw_output == "TIMEOUT" or raw_output.startswith("ERROR"):
        return "error", raw_output

    response = extract_codex_response(raw_output)

    # Look for OK or Anomaly in the first few lines of the extracted response
    lines = response.strip().split("\n")
    for line in lines[:5]:
        line_clean = line.strip()
        line_upper = line_clean.upper()
        # Check for Anomaly first (more specific)
        if "ANOMALY" in line_upper:
            return "Anomaly", response
        # Check for OK (but not as part of another word like "looked")
        if re.match(r'^OK\b', line_upper):
            return "OK", response

    # Fallback: search entire response
    resp_upper = response.upper()
    if "ANOMALY" in resp_upper:
        return "Anomaly", response

    # Check for standalone OK
    m = re.search(r'\bOK\b', response)
    if m:
        return "OK", response

    return "unclear", response

def run_evaluation(mode="setup_conditioned"):
    """Run all 30 items in the given mode."""
    codes = generate_blind_codes(30)
    results = []

    print(f"\n{'='*60}")
    print(f"  GPT-5 Evaluation: {mode}")
    print(f"{'='*60}\n")

    for i, item in enumerate(items):
        case = item["case"]
        gt = item["gt"]
        code = codes[i]
        img_path = get_image_path(case)

        if not img_path.exists():
            print(f"  [{i+1}/30] {code}: IMAGE NOT FOUND {img_path}")
            results.append({
                "id": item["id"],
                "code": code,
                "case": case,
                "gt": gt,
                "verdict": "error",
                "match": False,
                "reasoning": f"Image not found: {img_path}"
            })
            continue

        if mode == "setup_conditioned":
            setup = get_setup_text(case)
            prompt = (
                f"Read the image at {img_path}. "
                f"Problem setup: {setup}. "
                f"Is this flow field physically plausible for the stated conditions? "
                f"Answer: OK (plausible) or Anomaly (non-physical pattern detected). "
                f"Provide your verdict on the first line, then brief reasoning."
            )
        else:
            prompt = (
                f"Read the image at {img_path}. "
                f"You have NO problem setup information. "
                f"Does this flow/temperature field look physically plausible based on visual appearance alone? "
                f"Answer: OK or Anomaly on the first line, then brief reasoning."
            )

        print(f"  [{i+1}/30] {code} ...", end=" ", flush=True)
        raw_output = run_codex(prompt)
        verdict, full_resp = parse_verdict(raw_output)

        # Determine match
        expected = "OK" if gt == "correct" else "Anomaly"
        match = (verdict == expected)

        print(f"verdict={verdict} (expected={expected}) {'MATCH' if match else 'MISS'}")

        results.append({
            "id": item["id"],
            "code": code,
            "case": case,
            "gt": gt,
            "verdict": verdict,
            "match": match,
            "reasoning": full_resp[:500] if full_resp else ""
        })

        # Small delay to avoid rate limiting
        time.sleep(1)

    return results

def compute_stats(results):
    total = len(results)
    correct = sum(1 for r in results if r["match"])

    # FP: predicted Anomaly but GT is correct
    fp = sum(1 for r in results if r["verdict"] == "Anomaly" and r["gt"] == "correct")
    fp_cases = [r["case"] for r in results if r["verdict"] == "Anomaly" and r["gt"] == "correct"]

    # FN: predicted OK but GT is error
    fn = sum(1 for r in results if r["verdict"] == "OK" and r["gt"] != "correct")
    fn_cases = [r["case"] for r in results if r["verdict"] == "OK" and r["gt"] != "correct"]

    # Errors (unclear/error verdicts)
    errors = sum(1 for r in results if r["verdict"] not in ("OK", "Anomaly"))
    error_cases = [r["case"] for r in results if r["verdict"] not in ("OK", "Anomaly")]

    return {
        "total": total,
        "correct": correct,
        "accuracy": correct / total if total > 0 else 0,
        "fp": fp,
        "fp_cases": fp_cases,
        "fn": fn,
        "fn_cases": fn_cases,
        "errors": errors,
        "error_cases": error_cases
    }

def save_results(results, stats, mode, output_path):
    data = {
        "evaluator": "gpt-5-codex-cli",
        "task": mode,
        "date": "2026-03-28",
        "total": stats["total"],
        "correct": stats["correct"],
        "accuracy": stats["accuracy"],
        "fp": stats["fp"],
        "fp_cases": stats["fp_cases"],
        "fn": stats["fn"],
        "fn_cases": stats["fn_cases"],
        "errors": stats["errors"],
        "error_cases": stats["error_cases"],
        "items": results
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"\nSaved to: {output_path}")

def print_summary(stats, mode):
    print(f"\n--- {mode} Summary ---")
    print(f"  Total: {stats['total']}")
    print(f"  Correct: {stats['correct']}")
    print(f"  Accuracy: {stats['accuracy']:.1%}")
    print(f"  False Positives: {stats['fp']} {stats['fp_cases']}")
    print(f"  False Negatives: {stats['fn']} {stats['fn_cases']}")
    print(f"  Errors/Unclear: {stats['errors']} {stats['error_cases']}")

if __name__ == "__main__":
    random.seed(42)

    # 1. Setup-conditioned evaluation
    sc_results = run_evaluation("setup_conditioned")
    sc_stats = compute_stats(sc_results)
    save_results(sc_results, sc_stats, "setup_conditioned",
                 LABELS / "vlm_eval_gpt5_setup_conditioned.json")
    print_summary(sc_stats, "setup_conditioned")

    # 2. Image-only evaluation
    io_results = run_evaluation("image_only")
    io_stats = compute_stats(io_results)
    save_results(io_results, io_stats, "image_only",
                 LABELS / "vlm_eval_gpt5_image_only.json")
    print_summary(io_stats, "image_only")

    # Final comparison
    print(f"\n{'='*60}")
    print(f"  FINAL COMPARISON")
    print(f"{'='*60}")
    print(f"  Setup-conditioned: {sc_stats['accuracy']:.1%} ({sc_stats['correct']}/{sc_stats['total']})")
    print(f"  Image-only:        {io_stats['accuracy']:.1%} ({io_stats['correct']}/{io_stats['total']})")
    print(f"{'='*60}")
