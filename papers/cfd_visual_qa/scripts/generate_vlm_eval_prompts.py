#!/usr/bin/env python3
"""Generate blind VLM evaluation prompts for CFD Visual QA benchmark.

Creates a directory of evaluation items, each with:
- Blind code (no case name exposure)
- Problem setup text
- Image path
- Question

Output can be fed to any VLM (Claude, Gemini, GPT-4o, etc.) via API or CLI.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
BENCH_DIR = PROJECT_ROOT / "papers" / "cfd_visual_qa" / "benchmark"
IMAGE_DIR = BENCH_DIR / "images"
LABELS_DIR = BENCH_DIR / "labels"

# Evaluation items — same as Day 1 expert labeling + extensions
# Each item: (case, viz, setup, question)
EVAL_ITEMS_DAY1 = [
    ("S5_correct_lam", "V6",
     "Lid-driven cavity (1m×1m). Re=100 laminar. Top wall moves right at U=1m/s. Other 3 walls stationary. Velocity contour + streamlines.",
     "Is the primary vortex location and rotation direction physically plausible for a lid-driven cavity flow?"),
    ("S5_E1_underconverged", "V6",
     "Lid-driven cavity (1m×1m). Re=10,000 turbulent. Top wall moves right at U=1m/s. k-ω SST. Velocity contour + streamlines.",
     "Has the flow pattern developed adequately for a lid-driven cavity?"),
    ("S5_E2_bc_swap", "V6",
     "Lid-driven cavity (1m×1m). Re=10,000 turbulent. Moving wall (lid) at U=1m/s. k-ω SST. Velocity contour + streamlines.",
     "Does the primary vortex rotation and high-velocity region match the lid movement direction?"),
    ("S5_E5_coarse_mesh", "V6",
     "Lid-driven cavity (1m×1m). Re=10,000 turbulent. Top lid moves right. k-ω SST. Velocity contour + streamlines.",
     "Is the vortex structure resolution adequate? Are corner secondary vortices resolved?"),
    ("S6_correct_Ra1e4", "V2",
     "Differentially heated cavity (1m×1m). Left wall hot (305K), right wall cold (295K), top/bottom adiabatic. Ra=10^4. Temperature contour.",
     "Is the temperature distribution and isotherm shape physically plausible for natural convection?"),
    ("S6_E8_gravity_flipped", "V2",
     "Differentially heated cavity (1m×1m). Left 305K, right 295K, top/bottom adiabatic. Ra=10^4. Temperature contour.",
     "Does the temperature distribution match the expected natural convection heat transfer pattern?"),
    ("S6_E2_bc_swap", "V2",
     "Differentially heated cavity (1m×1m). Wall temperature BCs applied for natural convection. Ra=10^4. Temperature contour.",
     "Are the hot/cold wall positions and convection circulation direction physically consistent?"),
    ("S9_correct_Ra5k", "V2",
     "Rayleigh-Bénard convection (4m×1m cavity). Bottom wall heated (310K), top wall cooled (300K), sidewalls adiabatic. Ra=5,000. Temperature contour.",
     "Are convection cells observed? Is the cell pattern plausible for Ra=5,000?"),
    ("S9_E8_gravity_flipped", "V2",
     "Rayleigh-Bénard (4m×1m). Bottom heated (310K), top cooled (300K). Ra=5,000. Temperature contour.",
     "Is this temperature distribution physically plausible for bottom-heated + top-cooled conditions?"),
    ("S9_E5_coarse_mesh", "V2",
     "Rayleigh-Bénard (4m×1m). Bottom 310K, top 300K. Ra=20,000. Temperature contour.",
     "Is the convection cell structure well-resolved? Any anomalies in the temperature distribution?"),
    ("S10_correct_turb", "V6",
     "Ventilated room (9m×3m). Slot inlet at top-left (0.168m, U=0.455m/s), slot outlet at bottom-right. Re≈5,000. k-ω SST. Velocity contour + streamlines.",
     "Is the inlet jet ceiling attachment (Coanda effect) and indoor circulation pattern physically plausible?"),
    ("S10_E1_underconverged", "V6",
     "Ventilated room (9m×3m). Top-left slot inlet, bottom-right outlet. Re≈5,000. k-ω SST. Velocity contour + streamlines.",
     "Has the flow pattern developed adequately for a ventilated room?"),
    ("S10_E4_wrong_turb_model", "V6",
     "Ventilated room (9m×3m). Top-left inlet, bottom-right outlet. Re≈5,000. Laminar solver. Velocity contour + streamlines.",
     "Is this flow pattern physically plausible at Re≈5,000?"),
    ("S2_correct_lam", "V6",
     "2D channel flow (10m×1m). Left inlet U=1m/s uniform, right outlet. Re=100 laminar. Top/bottom no-slip walls. Velocity contour + streamlines.",
     "Is the velocity profile development (uniform→parabolic) observed? Physically plausible?"),
    ("S2_E3_wrong_viscosity", "V6",
     "2D channel flow (10m×1m). Left inlet U=1m/s, right outlet. Laminar. Top/bottom no-slip. Velocity contour + streamlines.",
     "Is the velocity profile development rate and shape plausible for channel flow?"),
]


def generate_prompts(items: list[tuple], output_dir: Path) -> None:
    """Generate blind evaluation prompt files."""
    blind_map_path = LABELS_DIR / "blind_code_mapping.json"
    if blind_map_path.exists():
        blind = json.load(open(blind_map_path))["case_to_code"]
    else:
        import hashlib
        blind = {case: f"CFD-{hashlib.md5(case.encode()).hexdigest()[:4].upper()}"
                 for case, *_ in items}

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for i, (case, viz, setup, question) in enumerate(items):
        code = blind.get(case, f"CFD-{i:04d}")
        item_id = f"{i+1:02d}"

        # Copy image with blind name
        src = IMAGE_DIR / case / f"{case}_{viz}.png"
        dst = output_dir / f"{code}_{viz}.png"
        if src.exists():
            shutil.copy2(src, dst)

        # Generate prompt text
        prompt = (
            f"## Item {item_id}: {code} — {viz}\n\n"
            f"**Problem Setup:** {setup}\n\n"
            f"**Question:** {question}\n\n"
            f"**Image:** {code}_{viz}.png\n\n"
            f"Evaluate whether the flow/temperature field shown in the image is physically "
            f"plausible given the problem setup. Answer:\n"
            f"- 'OK' if plausible\n"
            f"- 'Anomaly: [explanation]' if something appears non-physical\n"
        )

        manifest.append({
            "id": int(item_id),
            "code": code,
            "case": case,  # stored in manifest but NOT in prompt
            "viz": viz,
            "setup": setup,
            "question": question,
            "image_file": f"{code}_{viz}.png",
        })

        (output_dir / f"{item_id}_{code}_prompt.txt").write_text(prompt)

    # Save manifest (contains case mapping — for scoring only)
    json.dump(manifest, open(output_dir / "manifest.json", "w"), indent=2)
    print(f"Generated {len(manifest)} evaluation prompts in {output_dir}")

    # Save combined prompt (for batch API calls)
    combined = "# CFD Visual QA — Blind VLM Evaluation\n\n"
    combined += "Evaluate each flow/temperature field image. For each item:\n"
    combined += "- Read the Problem Setup\n"
    combined += "- Examine the image\n"
    combined += "- Judge if physically plausible: 'OK' or 'Anomaly: [explanation]'\n\n"
    combined += "---\n\n"
    for item in manifest:
        prompt_file = output_dir / f"{item['id']:02d}_{item['code']}_prompt.txt"
        combined += prompt_file.read_text() + "\n---\n\n"
    (output_dir / "combined_prompt.txt").write_text(combined)
    print(f"Combined prompt: {output_dir / 'combined_prompt.txt'}")


def main():
    output_dir = BENCH_DIR / "vlm_eval_prompts" / "day01_15items"
    generate_prompts(EVAL_ITEMS_DAY1, output_dir)


if __name__ == "__main__":
    main()
