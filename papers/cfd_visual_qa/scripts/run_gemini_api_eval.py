#!/usr/bin/env python3
"""
Gemini API evaluation of CFD flow field images for CFD-VisQA benchmark.

Uses google-genai SDK to send base64-encoded images directly to Gemini.
The model receives ONLY image bytes + prompt text — NO filesystem access.

Two modes: setup-conditioned (30 items) and image-only (30 items) = 60 total calls.

Output:
  vlm_eval_gemini_api_setup.json
  vlm_eval_gemini_api_imageonly.json
"""

import json
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

# Load API key from .env before importing genai
_env_path = Path("/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa/.env")
if _env_path.exists():
    for _line in _env_path.read_text().splitlines():
        _line = _line.strip()
        if _line and "=" in _line and not _line.startswith("#"):
            _k, _v = _line.split("=", 1)
            os.environ[_k.strip()] = _v.strip()

from google import genai
from google.genai import types

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa")
LABELS_DIR = BASE / "benchmark" / "labels"
IMAGES_DIR = BASE / "benchmark" / "images"
CASES_DIR = BASE / "benchmark" / "cases"

ITEMS_FILE = LABELS_DIR / "vlm_eval_claude_image_only.json"
OUTPUT_SETUP = LABELS_DIR / "vlm_eval_gemini_api_setup.json"
OUTPUT_IMAGEONLY = LABELS_DIR / "vlm_eval_gemini_api_imageonly.json"

# Scenarios that use V2 (temperature contour); others use V6 (velocity+streamlines)
TEMPERATURE_SCENARIOS = {"S1", "S6", "S9"}

PRIMARY_MODEL = "gemini-3.1-pro-preview"
FALLBACK_MODEL = "gemini-2.5-pro"

INTER_CALL_DELAY = 2   # seconds between calls
RATE_LIMIT_WAIT = 30   # seconds to wait on rate limit before retry


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_image_path(case_name: str) -> Path:
    scenario = case_name.split("_")[0]
    suffix = "V2" if scenario in TEMPERATURE_SCENARIOS else "V6"
    return IMAGES_DIR / case_name / f"{case_name}_{suffix}.png"


def get_setup_text(case_name: str) -> str:
    meta_path = CASES_DIR / case_name / "case_meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    parts = [meta.get("description", "")]
    if meta.get("Re") is not None:
        parts.append(f"Re = {meta['Re']}")
    if meta.get("Ra") is not None:
        parts.append(f"Ra = {meta['Ra']}")
    if meta.get("Ri") is not None:
        parts.append(f"Ri = {meta['Ri']}")
    if meta.get("turbulence"):
        turb_model = meta.get("turbulence_model", "unknown")
        parts.append(f"Turbulence model: {turb_model}")
    else:
        parts.append("Laminar flow")

    mesh_parts = []
    if "mesh_nx" in meta:
        mesh_parts.append(f"nx={meta['mesh_nx']}")
    if "mesh_ny" in meta:
        mesh_parts.append(f"ny={meta['mesh_ny']}")
    if mesh_parts:
        parts.append(f"Mesh: {', '.join(mesh_parts)}")

    return ". ".join(p for p in parts if p)


def parse_verdict(output: str) -> str:
    """Parse OK or Anomaly from first non-empty line. Returns 'OK', 'Anomaly', or 'error'."""
    if not output or not output.strip():
        return "error"

    for line in output.strip().splitlines():
        line_clean = line.strip().lower()
        if not line_clean:
            continue
        if "anomaly" in line_clean:
            return "Anomaly"
        if re.match(r"^ok\b", line_clean):
            return "OK"

    # Fallback: count occurrences in full text
    text_lower = output.lower()
    anomaly_count = text_lower.count("anomaly")
    ok_count = len(re.findall(r"\bok\b", text_lower))

    if anomaly_count > 0 and anomaly_count >= ok_count:
        return "Anomaly"
    if ok_count > 0:
        return "OK"

    return "error"


def extract_reasoning(output: str) -> str:
    """Return text after the first verdict line, trimmed to 500 chars."""
    lines = output.strip().splitlines()
    reasoning_lines = []
    found_verdict = False
    for line in lines:
        stripped = line.strip()
        if not found_verdict:
            if re.match(r"^(ok|anomaly)\b", stripped, re.IGNORECASE):
                found_verdict = True
                rest = re.sub(r"^(ok|anomaly)\s*[:\-–]?\s*", "", stripped, flags=re.IGNORECASE).strip()
                if rest:
                    reasoning_lines.append(rest)
            continue
        if stripped:
            reasoning_lines.append(stripped)

    result = " ".join(reasoning_lines)
    if len(result) > 500:
        result = result[:497] + "..."
    return result if result else output.strip()[:500]


def call_gemini(client: genai.Client, model: str, img_bytes: bytes, prompt: str) -> str:
    response = client.models.generate_content(
        model=model,
        contents=[
            types.Part.from_bytes(data=img_bytes, mime_type="image/png"),
            types.Part.from_text(text=prompt),
        ],
    )
    return response.text


def call_with_retry_and_fallback(
    client: genai.Client, img_bytes: bytes, prompt: str
) -> tuple[str, str]:
    """
    Try PRIMARY_MODEL; on rate-limit wait RATE_LIMIT_WAIT seconds and retry once.
    If still failing, fall back to FALLBACK_MODEL.
    Returns (response_text, model_used).
    """
    rate_limit_keywords = ("429", "quota", "rate", "resource_exhausted")

    for model in (PRIMARY_MODEL, FALLBACK_MODEL):
        for attempt in range(2):
            try:
                text = call_gemini(client, model, img_bytes, prompt)
                return text, model
            except Exception as e:
                err_str = str(e).lower()
                is_rate_limit = any(kw in err_str for kw in rate_limit_keywords)
                if is_rate_limit and attempt == 0:
                    print(f"\n  [rate-limit on {model}] waiting {RATE_LIMIT_WAIT}s ...", end="", flush=True)
                    time.sleep(RATE_LIMIT_WAIT)
                    continue
                elif model == PRIMARY_MODEL:
                    # Non-rate-limit error on primary, or second attempt — try fallback
                    print(f"\n  [{model} failed: {e}] trying {FALLBACK_MODEL} ...", end="", flush=True)
                    break
                else:
                    # Fallback also failed — raise to caller
                    raise

    raise RuntimeError(f"Both {PRIMARY_MODEL} and {FALLBACK_MODEL} failed")


# ---------------------------------------------------------------------------
# Evaluation runner
# ---------------------------------------------------------------------------

def run_evaluation(
    client: genai.Client, items: list, mode: str
) -> dict:
    """
    Run one evaluation mode over all items.
    mode: 'setup' | 'imageonly'
    Returns result dict ready to serialize as JSON.
    Raises RateLimitStop if we should stop early and save partial.
    """
    results = []
    correct = fp = fn = errors = 0
    fp_cases: list[str] = []
    fn_cases: list[str] = []
    error_cases: list[str] = []
    partial = False
    total = len(items)
    mode_label = "setup" if mode == "setup" else "image-only"

    for idx, item in enumerate(items):
        item_id = item["id"]
        case_name = item["case"]
        gt = item["gt"]
        image_path = get_image_path(case_name)

        if not image_path.exists():
            verdict = "error"
            reasoning = f"Image not found: {image_path}"
            match = False
            errors += 1
            error_cases.append(case_name)
            print(f"Item {idx+1}/{total} ({mode_label}) — error — {case_name}  [IMAGE MISSING]")
            results.append({"id": item_id, "case": case_name, "gt": gt,
                            "verdict": verdict, "match": match, "reasoning": reasoning})
            continue

        img_bytes = image_path.read_bytes()

        if mode == "setup":
            setup_text = get_setup_text(case_name)
            prompt = (
                f"Problem setup: {setup_text}. "
                "Is this flow field physically plausible for the stated conditions? "
                "Answer on the first line: OK (plausible) or Anomaly (non-physical). "
                "Then provide brief reasoning."
            )
        else:
            prompt = (
                "You have NO problem setup information. "
                "Based on visual appearance alone, does this flow/temperature field "
                "look physically plausible? "
                "Answer on the first line: OK or Anomaly. "
                "Then provide brief reasoning."
            )

        try:
            raw_output, model_used = call_with_retry_and_fallback(client, img_bytes, prompt)
            verdict = parse_verdict(raw_output)
            reasoning = extract_reasoning(raw_output) if verdict != "error" else raw_output[:500]
        except Exception as e:
            err_str = str(e).lower()
            is_rate_limit = any(kw in err_str for kw in ("429", "quota", "rate", "resource_exhausted"))
            if is_rate_limit:
                print(f"\n  [RATE LIMIT] saving partial results and stopping.")
                partial = True
                # Save what we have so far and return
                break
            verdict = "error"
            reasoning = f"ERROR: {e}"[:500]
            model_used = PRIMARY_MODEL
            error_cases.append(case_name)
            errors += 1

        if verdict != "error":
            expected = "OK" if gt == "correct" else "Anomaly"
            match = verdict == expected
            if match:
                correct += 1
            elif verdict == "Anomaly":
                fp += 1
                fp_cases.append(case_name)
            else:
                fn += 1
                fn_cases.append(case_name)
        else:
            match = False
            error_cases.append(case_name)
            errors += 1

        print(f"Item {idx+1}/{total} ({mode_label}) — {verdict} — {case_name}")

        results.append({
            "id": item_id,
            "case": case_name,
            "gt": gt,
            "verdict": verdict,
            "match": match,
            "reasoning": reasoning,
        })

        if idx < total - 1:
            time.sleep(INTER_CALL_DELAY)

    evaluated = len(results) - errors
    accuracy = correct / evaluated if evaluated > 0 else 0.0

    output = {
        "evaluator": PRIMARY_MODEL,
        "task": "setup_conditioned" if mode == "setup" else "image_only",
        "date": datetime.now().strftime("%Y-%m-%d"),
        "total": total,
        "evaluated": len(results),
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
    if partial:
        output["partial"] = True

    return output, partial


# ---------------------------------------------------------------------------
# Summary printer
# ---------------------------------------------------------------------------

def print_summary(setup: dict, imageonly: dict) -> None:
    print("\n" + "=" * 70)
    print("SUMMARY: Gemini API Evaluation Results")
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
        v1 = setup[key]
        v2 = imageonly[key]
        if key == "accuracy":
            print(f"{label:<25} {v1:>17.1%} {v2:>17.1%}")
        else:
            print(f"{label:<25} {v1:>12}/{setup['total']} {v2:>12}/{imageonly['total']}")
    print("-" * 70)

    if setup.get("fn_cases"):
        print(f"\nSetup FN cases ({setup['fn']}):")
        for c in setup["fn_cases"]:
            print(f"  - {c}")

    if imageonly.get("fn_cases"):
        print(f"\nImage-only FN cases ({imageonly['fn']}):")
        for c in imageonly["fn_cases"]:
            print(f"  - {c}")

    if setup.get("fp_cases"):
        print(f"\nSetup FP cases ({setup['fp']}):")
        for c in setup["fp_cases"]:
            print(f"  - {c}")

    if imageonly.get("fp_cases"):
        print(f"\nImage-only FP cases ({imageonly['fp']}):")
        for c in imageonly["fp_cases"]:
            print(f"  - {c}")

    for label, result in [("Setup", setup), ("Image-only", imageonly)]:
        if result.get("partial"):
            print(f"\n[!] {label} results are PARTIAL (rate-limited at item {result['evaluated']}/{result['total']})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        print("ERROR: GEMINI_API_KEY not set. Check .env file.")
        sys.exit(1)

    client = genai.Client(api_key=api_key)
    print(f"Gemini client initialized. Primary: {PRIMARY_MODEL}, Fallback: {FALLBACK_MODEL}")

    with open(ITEMS_FILE) as f:
        data = json.load(f)
    items = data["items"]
    print(f"Loaded {len(items)} items from {ITEMS_FILE.name}")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print()

    # --- Mode 1: Setup-conditioned ---
    print("=" * 60)
    print("MODE 1: Setup-Conditioned")
    print("=" * 60)
    setup_results, setup_partial = run_evaluation(client, items, "setup")
    with open(OUTPUT_SETUP, "w") as f:
        json.dump(setup_results, f, indent=2)
    print(f"\nSaved: {OUTPUT_SETUP}")

    if setup_partial:
        print("[!] Rate-limited during setup mode. Stopping.")
        sys.exit(0)

    print()

    # --- Mode 2: Image-only ---
    print("=" * 60)
    print("MODE 2: Image-Only")
    print("=" * 60)
    imageonly_results, imageonly_partial = run_evaluation(client, items, "imageonly")
    with open(OUTPUT_IMAGEONLY, "w") as f:
        json.dump(imageonly_results, f, indent=2)
    print(f"\nSaved: {OUTPUT_IMAGEONLY}")

    print_summary(setup_results, imageonly_results)


if __name__ == "__main__":
    main()
