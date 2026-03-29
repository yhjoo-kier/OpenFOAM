#!/usr/bin/env python3
"""
Isolated Gemini API evaluation of CFD flow field images for CFD-VisQA benchmark.

Uses google-genai SDK to send base64-encoded images directly to Gemini.
The model receives ONLY image bytes + prompt text — NO filesystem access.

Two modes: setup-conditioned (30 items) and image-only (30 items) = 60 total calls.
"""

import json
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

from google import genai
from google.genai import types

# Paths
BASE = Path("/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa")
LABELS_DIR = BASE / "benchmark" / "labels"
IMAGES_DIR = BASE / "benchmark" / "images"
CASES_DIR = BASE / "benchmark" / "cases"

ITEMS_FILE = LABELS_DIR / "vlm_eval_claude_image_only.json"
OUTPUT_SETUP = LABELS_DIR / "vlm_eval_gemini_isolated_setup.json"
OUTPUT_IMAGEONLY = LABELS_DIR / "vlm_eval_gemini_isolated_imageonly.json"

# Temperature scenarios use V2 (temperature contour), others use V6 (velocity+streamlines)
TEMPERATURE_SCENARIOS = {"S1", "S6", "S9"}

# Model preference order
MODELS = ["gemini-2.5-pro", "gemini-2.0-flash-001"]


def get_image_path(case_name: str) -> Path:
    scenario = case_name.split("_")[0]
    suffix = "V2" if scenario in TEMPERATURE_SCENARIOS else "V6"
    return IMAGES_DIR / case_name / f"{case_name}_{suffix}.png"


def get_setup_text(case_name: str) -> str:
    meta_path = CASES_DIR / case_name / "case_meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    parts = [meta.get("description", "")]

    if "Re" in meta and meta["Re"] is not None:
        parts.append(f"Re = {meta['Re']}")
    if "Ra" in meta and meta["Ra"] is not None:
        parts.append(f"Ra = {meta['Ra']}")
    if "Ri" in meta and meta["Ri"] is not None:
        parts.append(f"Ri = {meta['Ri']}")
    if meta.get("turbulence"):
        turb_model = meta.get("turbulence_model", "unknown")
        parts.append(f"Turbulence model: {turb_model}")
    else:
        parts.append("Laminar flow")

    mesh_info = []
    if "mesh_nx" in meta:
        mesh_info.append(f"nx={meta['mesh_nx']}")
    if "mesh_ny" in meta:
        mesh_info.append(f"ny={meta['mesh_ny']}")
    if mesh_info:
        parts.append(f"Mesh: {', '.join(mesh_info)}")

    return ". ".join(p for p in parts if p)


def parse_verdict(output: str) -> str:
    """Parse OK or Anomaly from model output. Returns 'OK', 'Anomaly', or 'error'."""
    if not output or not output.strip():
        return "error"

    for line in output.strip().splitlines():
        line_clean = line.strip().lower()
        if not line_clean:
            continue
        if "anomaly" in line_clean:
            return "Anomaly"
        if line_clean.startswith("ok"):
            return "OK"
        if re.match(r"^(ok|anomaly)", line_clean):
            return "Anomaly" if "anomaly" in line_clean else "OK"

    # Fallback: search entire output
    text_lower = output.lower()
    anomaly_count = text_lower.count("anomaly")
    ok_count = len(re.findall(r'\bok\b', text_lower))

    if anomaly_count > 0 and anomaly_count >= ok_count:
        return "Anomaly"
    if ok_count > 0:
        return "OK"

    return "error"


def extract_reasoning(output: str) -> str:
    """Extract reasoning from output, trimming to reasonable length."""
    lines = output.strip().splitlines()
    reasoning_lines = []
    found_verdict = False
    for line in lines:
        stripped = line.strip()
        if not found_verdict:
            if stripped.lower().startswith(("ok", "anomaly")):
                found_verdict = True
                rest = re.sub(r'^(ok|anomaly)\s*', '', stripped, flags=re.IGNORECASE).strip()
                if rest:
                    reasoning_lines.append(rest)
                continue
            continue
        if stripped:
            reasoning_lines.append(stripped)

    result = " ".join(reasoning_lines)
    if len(result) > 500:
        result = result[:497] + "..."
    return result if result else output.strip()[:500]


def call_gemini(client: genai.Client, model: str, img_bytes: bytes, prompt: str) -> str:
    """Call Gemini API with image bytes and prompt. Returns response text."""
    response = client.models.generate_content(
        model=model,
        contents=[
            types.Part.from_bytes(data=img_bytes, mime_type="image/png"),
            types.Part.from_text(text=prompt),
        ],
    )
    return response.text


def call_with_fallback(client: genai.Client, img_bytes: bytes, prompt: str) -> str:
    """Try primary model, fall back to secondary on error."""
    for i, model in enumerate(MODELS):
        try:
            return call_gemini(client, model, img_bytes, prompt)
        except Exception as e:
            if i < len(MODELS) - 1:
                print(f"[{model} failed: {e}] falling back... ", end="", flush=True)
                time.sleep(2)
            else:
                raise


def run_evaluation(client: genai.Client, items: list, mode: str) -> dict:
    """Run one evaluation mode (setup or imageonly) on all items."""
    results = []
    correct = 0
    fp = 0
    fn = 0
    errors = 0
    fp_cases = []
    fn_cases = []
    error_cases = []

    total = len(items)

    for idx, item in enumerate(items):
        item_id = item["id"]
        case_name = item["case"]
        gt = item["gt"]
        image_path = get_image_path(case_name)

        if not image_path.exists():
            print(f"  [!] Image not found: {image_path}")
            error_cases.append(case_name)
            errors += 1
            results.append({
                "id": item_id,
                "case": case_name,
                "gt": gt,
                "verdict": "error",
                "match": False,
                "reasoning": f"Image not found: {image_path}",
            })
            continue

        # Read image as bytes
        img_bytes = image_path.read_bytes()

        # Build prompt
        if mode == "setup":
            setup_text = get_setup_text(case_name)
            prompt = (
                f"Problem setup: {setup_text}. "
                f"Is this flow field physically plausible for the stated conditions? "
                f"Answer on the first line: OK (plausible) or Anomaly (non-physical). "
                f"Then provide brief reasoning."
            )
        else:
            prompt = (
                "You have NO problem setup information. "
                "Based on visual appearance alone, does this flow/temperature field "
                "look physically plausible? "
                "Answer on the first line: OK or Anomaly. "
                "Then provide brief reasoning."
            )

        print(f"  [{mode}] Item {idx+1}/{total}: {case_name} (gt={gt}) ... ", end="", flush=True)

        try:
            raw_output = call_with_fallback(client, img_bytes, prompt)
            verdict = parse_verdict(raw_output)
            reasoning = extract_reasoning(raw_output) if verdict != "error" else raw_output[:500]
        except Exception as e:
            raw_output = f"ERROR: {e}"
            verdict = "error"
            reasoning = raw_output[:500]

        expected = "OK" if gt == "correct" else "Anomaly"
        match = (verdict == expected)

        if verdict == "error":
            errors += 1
            error_cases.append(case_name)
        elif match:
            correct += 1
        elif verdict == "Anomaly" and expected == "OK":
            fp += 1
            fp_cases.append(case_name)
        elif verdict == "OK" and expected == "Anomaly":
            fn += 1
            fn_cases.append(case_name)

        status = "MATCH" if match else ("MISS" if verdict != "error" else "ERROR")
        print(f"verdict={verdict} {status}")

        results.append({
            "id": item_id,
            "case": case_name,
            "gt": gt,
            "verdict": verdict,
            "match": match,
            "reasoning": reasoning,
        })

        # Rate limit: small delay between calls
        time.sleep(1)

    evaluated = total - errors
    accuracy = correct / evaluated if evaluated > 0 else 0.0

    return {
        "evaluator": f"gemini-isolated ({MODELS[0]})",
        "task": "setup_conditioned" if mode == "setup" else "image_only",
        "date": datetime.now().strftime("%Y-%m-%d"),
        "total": total,
        "correct": correct,
        "accuracy": accuracy,
        "fp": fp,
        "fp_cases": fp_cases,
        "fn": fn,
        "fn_cases": fn_cases,
        "errors": errors,
        "error_cases": error_cases,
        "items": results,
    }


def print_summary(setup_results: dict, imageonly_results: dict):
    """Print comparison table."""
    print("\n" + "=" * 70)
    print("SUMMARY: Gemini Isolated Evaluation Results")
    print("=" * 70)

    header = f"{'Metric':<25} {'Setup-Conditioned':>18} {'Image-Only':>18}"
    print(header)
    print("-" * 70)

    for key, label in [
        ("accuracy", "Accuracy"),
        ("correct", "Correct"),
        ("fp", "False Positives"),
        ("fn", "False Negatives"),
        ("errors", "Errors"),
    ]:
        v1 = setup_results[key]
        v2 = imageonly_results[key]
        if key == "accuracy":
            print(f"{label:<25} {v1:>17.1%} {v2:>17.1%}")
        else:
            s_total = setup_results["total"]
            i_total = imageonly_results["total"]
            print(f"{label:<25} {v1:>12}/{s_total} {v2:>12}/{i_total}")

    print("-" * 70)

    print(f"\nSetup FN cases ({setup_results['fn']}):")
    for c in setup_results["fn_cases"]:
        print(f"  - {c}")

    print(f"\nImage-only FN cases ({imageonly_results['fn']}):")
    for c in imageonly_results["fn_cases"]:
        print(f"  - {c}")

    if setup_results["fp_cases"]:
        print(f"\nSetup FP cases ({setup_results['fp']}):")
        for c in setup_results["fp_cases"]:
            print(f"  - {c}")

    if imageonly_results["fp_cases"]:
        print(f"\nImage-only FP cases ({imageonly_results['fp']}):")
        for c in imageonly_results["fp_cases"]:
            print(f"  - {c}")


def main():
    # Initialize Gemini client
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("ERROR: GEMINI_API_KEY environment variable not set")
        sys.exit(1)

    client = genai.Client(api_key=api_key)
    print(f"Gemini client initialized (models: {' -> '.join(MODELS)})")

    # Load items
    with open(ITEMS_FILE) as f:
        data = json.load(f)
    items = data["items"]
    print(f"Loaded {len(items)} evaluation items")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print()

    # Mode 1: Setup-conditioned
    print("=" * 60)
    print("MODE 1: Setup-Conditioned Evaluation")
    print("=" * 60)
    setup_results = run_evaluation(client, items, "setup")

    with open(OUTPUT_SETUP, "w") as f:
        json.dump(setup_results, f, indent=2)
    print(f"\nSaved: {OUTPUT_SETUP}")

    print()

    # Mode 2: Image-only
    print("=" * 60)
    print("MODE 2: Image-Only Evaluation")
    print("=" * 60)
    imageonly_results = run_evaluation(client, items, "imageonly")

    with open(OUTPUT_IMAGEONLY, "w") as f:
        json.dump(imageonly_results, f, indent=2)
    print(f"\nSaved: {OUTPUT_IMAGEONLY}")

    # Summary
    print_summary(setup_results, imageonly_results)


if __name__ == "__main__":
    main()
