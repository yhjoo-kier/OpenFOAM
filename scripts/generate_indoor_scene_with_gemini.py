#!/usr/bin/env python3
"""Generate indoor CFD scene JSON with Gemini and validate it.

Pipeline stage:
    prompt (+ optional images) -> Gemini CLI/API -> scene.json -> validator
"""

from __future__ import annotations

import argparse
import base64
import json
import mimetypes
import os
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

FALLBACK_MODELS = ["gemini-3-flash-preview", "gemini-2.5-flash-lite"]
API_BASE_URL = "https://generativelanguage.googleapis.com/v1beta/models"

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from validate_indoor_scene import format_report, validate_scene  # noqa: E402

DEFAULT_MODEL = "gemini-3.1-pro-preview"
DEFAULT_BACKEND = "api"
PRIMARY_MODEL_RETRIES = 5
RETRYABLE_ERROR_MARKERS = (
    "MODEL_CAPACITY_EXHAUSTED",
    "RESOURCE_EXHAUSTED",
    "status 429",
    "status 503",
    "UNAVAILABLE",
    "currently experiencing high demand",
)

PROMPT_TEMPLATE = """You are generating a structured indoor CFD scene specification.
Return ONLY valid JSON. Do not include markdown, comments, or explanation.

Interpret any provided image(s), photos, or renderings as references for a **simulation-ready abstraction**,
not as a requirement for literal photoreal reconstruction.
Preserve the important flow-relevant layout qualitatively while simplifying geometry into solver-friendly boxes/openings.

Use exactly this schema:
{{
  \"schema_version\": \"indoor_cfd_scene_v1\",
  \"units\": \"m\",
  \"coordinate_system\": {{
    \"origin\": \"room_min_corner\",
    \"x\": \"west_to_east\",
    \"y\": \"south_to_north\",
    \"z\": \"floor_to_ceiling\"
  }},
  \"room\": {{
    \"size\": {{
      \"Lx\": number,
      \"Ly\": number,
      \"Lz\": number
    }}
  }},
  \"obstacles\": [
    {{
      \"id\": string,
      \"type\": \"box\",
      \"min\": {{\"x\": number, \"y\": number, \"z\": number}},
      \"size\": {{\"dx\": number, \"dy\": number, \"dz\": number}}
    }}
  ],
  \"openings\": [
    {{
      \"id\": string,
      \"type\": \"inlet\" | \"outlet\",
      \"wall\": \"west\" | \"east\" | \"south\" | \"north\" | \"floor\" | \"ceiling\",
      \"center\": {{\"u\": number, \"v\": number}},
      \"size\": {{\"du\": number, \"dv\": number}}
    }}
  ],
  \"meta\": {{
    \"scene_type\": string,
    \"description\": string
  }}
}}

Constraints:
- Create one rectangular room.
- Create 3 to 5 non-overlapping box obstacles.
- Create exactly 1 inlet and 1 outlet.
- All obstacles must lie fully inside the room.
- Obstacles must not overlap each other.
- Openings must lie fully on their specified wall and remain within room bounds.
- Use units of meters.
- All lengths must be positive.
- Use concise IDs like obs_001, inlet_001, outlet_001.
- Favor solver-friendly abstraction over visual detail.
- If image cues are ambiguous, choose a conservative, simulation-stable layout.

Scenario:
{scenario}
"""


def ensure_gemini_available(backend: str) -> None:
    if backend == "cli":
        if shutil.which("gemini") is None:
            raise RuntimeError("gemini CLI not found in PATH")
        return
    if backend == "api":
        if not os.environ.get("GEMINI_API_KEY"):
            raise RuntimeError("GEMINI_API_KEY is not set for Gemini API backend")
        return
    raise RuntimeError(f"Unsupported backend: {backend}")


def build_prompt(scenario: str) -> str:
    return PROMPT_TEMPLATE.format(scenario=scenario.strip())


def run_gemini_cli(prompt: str, model: str) -> str:
    cmd = ["gemini", "--model", model, "--output-format", "json", prompt]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "Gemini CLI failed with exit code "
            f"{result.returncode}\nMODEL:{model}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return result.stdout.strip()


def _image_to_part(image_path: Path) -> dict:
    mime_type, _ = mimetypes.guess_type(str(image_path))
    if not mime_type:
        mime_type = "application/octet-stream"
    data = base64.b64encode(image_path.read_bytes()).decode("ascii")
    return {
        "inline_data": {
            "mime_type": mime_type,
            "data": data,
        }
    }


def run_gemini_api(prompt: str, model: str, image_paths: list[Path]) -> str:
    api_key = os.environ["GEMINI_API_KEY"]
    url = f"{API_BASE_URL}/{model}:generateContent?key={api_key}"

    parts: list[dict] = [{"text": prompt}]
    for image_path in image_paths:
        parts.append(_image_to_part(image_path))

    payload = {
        "contents": [
            {
                "role": "user",
                "parts": parts,
            }
        ],
        "generationConfig": {
            "temperature": 0.2,
            "responseMimeType": "application/json",
        },
    }

    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"Gemini API failed with status {exc.code}\nMODEL:{model}\nBODY:\n{body}"
        ) from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Gemini API request failed for model {model}: {exc}") from exc

    return raw.strip()


def _is_retryable_error(message: str) -> bool:
    return any(marker in message for marker in RETRYABLE_ERROR_MARKERS)


def _run_backend_request(prompt: str, model: str, backend: str, image_paths: list[Path]) -> str:
    if backend == "cli":
        if image_paths:
            raise RuntimeError("Gemini CLI backend does not yet support --image inputs in this script")
        return run_gemini_cli(prompt, model)
    if backend == "api":
        return run_gemini_api(prompt, model, image_paths)
    raise RuntimeError(f"Unsupported backend: {backend}")


def generate_with_fallback(
    prompt: str,
    primary_model: str,
    backend: str,
    image_paths: list[Path],
    allow_fallback: bool = True,
    primary_retries: int = PRIMARY_MODEL_RETRIES,
) -> tuple[str, str, list[dict]]:
    models = [primary_model]
    if allow_fallback:
        models.extend([m for m in FALLBACK_MODELS if m != primary_model])

    history: list[dict] = []
    last_error: RuntimeError | None = None

    for idx, model in enumerate(models):
        max_attempts = primary_retries if idx == 0 else 1
        for attempt in range(1, max_attempts + 1):
            try:
                raw = _run_backend_request(prompt, model, backend, image_paths)
                history.append({
                    "model": model,
                    "attempt": attempt,
                    "ok": True,
                })
                return raw, model, history
            except RuntimeError as exc:
                last_error = exc
                message = str(exc)
                retryable = _is_retryable_error(message)
                history.append({
                    "model": model,
                    "attempt": attempt,
                    "ok": False,
                    "retryable": retryable,
                    "error": message,
                })

                should_retry_same_model = idx == 0 and attempt < max_attempts and retryable
                if should_retry_same_model:
                    delay = min(2 ** (attempt - 1), 16)
                    print(
                        f"Primary model attempt {attempt}/{max_attempts} failed ({model}); retrying in {delay}s...",
                        file=sys.stderr,
                    )
                    time.sleep(delay)
                    continue

                if idx < len(models) - 1 and allow_fallback:
                    print(
                        f"Model {model} failed after {attempt} attempt(s); retrying with fallback model...",
                        file=sys.stderr,
                    )
                    break
                raise

    assert last_error is not None
    raise last_error


def _strip_code_fences(text: str) -> str:
    text = text.strip()
    if text.startswith("```"):
        lines = text.splitlines()
        if len(lines) >= 2 and lines[0].startswith("```") and lines[-1].strip() == "```":
            return "\n".join(lines[1:-1]).strip()
    return text


def _extract_text_from_candidate(payload: dict) -> str | None:
    candidates = payload.get("candidates")
    if not isinstance(candidates, list):
        return None
    for candidate in candidates:
        content = candidate.get("content")
        if not isinstance(content, dict):
            continue
        parts = content.get("parts")
        if not isinstance(parts, list):
            continue
        texts = [part.get("text") for part in parts if isinstance(part, dict) and isinstance(part.get("text"), str)]
        if texts:
            return "\n".join(texts).strip()
    return None


def parse_scene(raw_output: str) -> dict:
    try:
        payload = json.loads(raw_output)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Gemini output is not valid JSON: {exc}\nRaw output:\n{raw_output}") from exc

    if not isinstance(payload, dict):
        raise RuntimeError("Gemini output root must be a JSON object")

    if "response" in payload and isinstance(payload["response"], str):
        candidate = _strip_code_fences(payload["response"])
        try:
            scene = json.loads(candidate)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                "Gemini response field did not contain valid scene JSON: "
                f"{exc}\nResponse field:\n{payload['response']}"
            ) from exc
    else:
        candidate_text = _extract_text_from_candidate(payload)
        if candidate_text is not None:
            candidate = _strip_code_fences(candidate_text)
            try:
                scene = json.loads(candidate)
            except json.JSONDecodeError as exc:
                raise RuntimeError(
                    "Gemini candidate text did not contain valid scene JSON: "
                    f"{exc}\nCandidate text:\n{candidate_text}"
                ) from exc
        else:
            scene = payload

    if not isinstance(scene, dict):
        raise RuntimeError("Parsed scene root must be a JSON object")
    return scene


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate indoor_cfd_scene_v1 JSON with Gemini CLI or API and validate it"
    )
    parser.add_argument(
        "--scenario",
        type=str,
        default=(
            "Generate a simple mechanically ventilated room with internal obstacles "
            "representing furniture or equipment."
        ),
        help="Natural-language scenario instruction for Gemini",
    )
    parser.add_argument(
        "--backend",
        choices=["cli", "api"],
        default=DEFAULT_BACKEND,
        help=f"Gemini backend to use (default: {DEFAULT_BACKEND})",
    )
    parser.add_argument(
        "--image",
        dest="images",
        action="append",
        default=[],
        help="Optional image/photo/rendering path to use as multimodal input; repeatable",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=PROJECT_ROOT / "generated" / "indoor_scene.json",
        help="Where to save the generated JSON",
    )
    parser.add_argument(
        "--model",
        type=str,
        default=DEFAULT_MODEL,
        help=f"Gemini model name (default: {DEFAULT_MODEL})",
    )
    parser.add_argument(
        "--print-prompt",
        action="store_true",
        help="Print the full prompt before running Gemini",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with nonzero code if validation warnings exist",
    )
    parser.add_argument(
        "--no-fallback",
        action="store_true",
        help="Disable automatic fallback to alternate Gemini models",
    )
    args = parser.parse_args()

    image_paths = [Path(p).expanduser().resolve() for p in args.images]
    for image_path in image_paths:
        if not image_path.exists():
            raise RuntimeError(f"Image path does not exist: {image_path}")

    ensure_gemini_available(args.backend)
    prompt = build_prompt(args.scenario)

    if args.print_prompt:
        print("=" * 80)
        print(prompt)
        if image_paths:
            print("-" * 80)
            print("Images:")
            for image_path in image_paths:
                print(image_path)
        print("=" * 80)

    raw_output, used_model, generation_history = generate_with_fallback(
        prompt,
        args.model,
        backend=args.backend,
        image_paths=image_paths,
        allow_fallback=not args.no_fallback,
    )
    scene = parse_scene(raw_output)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        json.dump(scene, f, indent=2, ensure_ascii=False)
        f.write("\n")

    report = validate_scene(scene)
    print(format_report(report, args.output))

    fallback_used = used_model != args.model
    summary = {
        "ok": report.ok and (not args.strict or not report.warnings),
        "errors": report.errors,
        "warnings": report.warnings,
        "output": str(args.output),
        "backend": args.backend,
        "model": used_model,
        "requested_model": args.model,
        "fallback_used": fallback_used,
        "fallback_notice": None if not fallback_used else f"Requested model '{args.model}' was unavailable; used fallback model '{used_model}'.",
        "primary_model_retry_limit": PRIMARY_MODEL_RETRIES,
        "generation_history": generation_history,
        "scenario": args.scenario,
        "images": [str(p) for p in image_paths],
    }
    print("---")
    if fallback_used:
        print(
            f"NOTICE: Requested model '{args.model}' was unavailable after retries; "
            f"fallback model '{used_model}' was used.",
            file=sys.stderr,
        )
    print(json.dumps(summary, indent=2, ensure_ascii=False))

    if not report.ok:
        return 1
    if args.strict and report.warnings:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
