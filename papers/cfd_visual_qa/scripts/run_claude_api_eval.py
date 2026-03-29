#!/usr/bin/env python3
"""
Claude API evaluation for CFD Visual QA benchmark.
Two modes: setup-conditioned and image-only.
"""
import os
import json
import base64
import time
import datetime

# Load env
for line in open('/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa/.env'):
    line = line.strip()
    if '=' in line:
        k, v = line.split('=', 1)
        os.environ[k] = v

import anthropic
client = anthropic.Anthropic()

BASE = '/home/yhjoo/projects/OpenFOAM/papers/cfd_visual_qa/benchmark'
LABELS_PATH = f'{BASE}/labels/vlm_eval_claude_image_only.json'
IMAGES_DIR = f'{BASE}/images'
CASES_DIR = f'{BASE}/cases'
OUT_DIR = f'{BASE}/labels'

MODEL = 'claude-opus-4-6'

# Scenes where _V2.png is used instead of _V6.png
V2_SCENES = {'S1', 'S6', 'S9'}

def get_image_path(case_name: str) -> str:
    scene = case_name.split('_')[0]
    suffix = '_V2.png' if scene in V2_SCENES else '_V6.png'
    return f'{IMAGES_DIR}/{case_name}/{case_name}{suffix}'

def get_setup_text(case_name: str) -> str:
    meta_path = f'{CASES_DIR}/{case_name}/case_meta.json'
    with open(meta_path) as f:
        meta = json.load(f)
    return meta.get('description', json.dumps(meta))

def encode_image(img_path: str) -> str:
    with open(img_path, 'rb') as f:
        return base64.standard_b64encode(f.read()).decode('utf-8')

def parse_verdict(text: str) -> str:
    first_line = text.strip().split('\n')[0].strip()
    if first_line.lower().startswith('ok'):
        return 'OK'
    elif first_line.lower().startswith('anomaly'):
        return 'Anomaly'
    else:
        # Try to find in first line
        fl = first_line.lower()
        if 'anomaly' in fl:
            return 'Anomaly'
        elif 'ok' in fl:
            return 'OK'
        return 'unknown'

def call_api(img_data: str, prompt_text: str) -> str:
    response = client.messages.create(
        model=MODEL,
        max_tokens=300,
        messages=[{
            'role': 'user',
            'content': [
                {
                    'type': 'image',
                    'source': {
                        'type': 'base64',
                        'media_type': 'image/png',
                        'data': img_data,
                    }
                },
                {'type': 'text', 'text': prompt_text},
            ]
        }]
    )
    return response.content[0].text

def run_mode(items: list, mode: str, out_path: str):
    results = []
    rate_limited = False

    for i, item in enumerate(items, 1):
        case = item['case']
        gt = item['gt']
        item_id = item['id']

        img_path = get_image_path(case)
        if not os.path.exists(img_path):
            print(f"Item {i}/30 ({mode}) — ERROR: image not found — {case}")
            results.append({
                'id': item_id, 'case': case, 'gt': gt,
                'verdict': 'error', 'match': False,
                'reasoning': f'Image not found: {img_path}'
            })
            continue

        try:
            img_data = encode_image(img_path)

            if mode == 'setup':
                setup_text = get_setup_text(case)
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
                result_text = call_api(img_data, prompt)
            except (anthropic.RateLimitError, anthropic.APIStatusError) as e:
                if isinstance(e, anthropic.APIStatusError) and e.status_code != 429:
                    raise
                print(f"  Rate limit hit on item {i}/30, waiting 60s and retrying...")
                time.sleep(60)
                try:
                    result_text = call_api(img_data, prompt)
                except (anthropic.RateLimitError, anthropic.APIStatusError) as e2:
                    print(f"RATE LIMITED after item {i-1}/30. Partial results saved.")
                    rate_limited = True
                    break

            verdict = parse_verdict(result_text)
            reasoning = '\n'.join(result_text.strip().split('\n')[1:]).strip()

            # Determine match: gt 'correct' → OK is match; gt starts with 'error' → Anomaly is match
            if gt == 'correct':
                match = (verdict == 'OK')
            else:
                match = (verdict == 'Anomaly')

            print(f"Item {i}/30 ({mode}) — {verdict} — {case}")
            results.append({
                'id': item_id,
                'case': case,
                'gt': gt,
                'verdict': verdict,
                'match': match,
                'reasoning': reasoning,
            })

        except Exception as e:
            print(f"Item {i}/30 ({mode}) — error: {e} — {case}")
            results.append({
                'id': item_id, 'case': case, 'gt': gt,
                'verdict': 'error', 'match': False,
                'reasoning': str(e)
            })

        if i < len(items):
            time.sleep(2)

    # Compute stats
    valid = [r for r in results if r['verdict'] != 'error']
    correct_items = [r for r in valid if r['match']]
    n_correct = len(correct_items)
    n_total = len(valid)
    accuracy = n_correct / n_total if n_total > 0 else 0.0

    fp = len([r for r in valid if r['verdict'] == 'Anomaly' and r['gt'] == 'correct'])
    fn = len([r for r in valid if r['verdict'] == 'OK' and r['gt'] != 'correct'])
    fn_cases = [r['case'] for r in valid if r['verdict'] == 'OK' and r['gt'] != 'correct']

    output = {
        'evaluator': MODEL,
        'task': 'setup_conditioned' if mode == 'setup' else 'image_only',
        'date': datetime.date.today().isoformat(),
        'total': n_total,
        'correct': n_correct,
        'accuracy': accuracy,
        'fp': fp,
        'fn': fn,
        'fn_cases': fn_cases,
        'items': results,
    }
    if rate_limited:
        output['partial'] = True

    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved {len(results)} results to {out_path}")
    print(f"Accuracy: {n_correct}/{n_total} = {accuracy:.4f}")
    return rate_limited

def main():
    with open(LABELS_PATH) as f:
        label_data = json.load(f)
    items = label_data['items']

    print(f"Loaded {len(items)} items")
    print(f"Model: {MODEL}\n")

    # Mode 1: Setup-conditioned
    print("=" * 60)
    print("MODE: setup-conditioned")
    print("=" * 60)
    out_setup = f'{OUT_DIR}/vlm_eval_claude_api_setup.json'
    rate_limited = run_mode(items, 'setup', out_setup)
    if rate_limited:
        print("Stopping due to rate limit.")
        return

    print()

    # Mode 2: Image-only
    print("=" * 60)
    print("MODE: image-only")
    print("=" * 60)
    out_imageonly = f'{OUT_DIR}/vlm_eval_claude_api_imageonly.json'
    run_mode(items, 'imageonly', out_imageonly)

if __name__ == '__main__':
    main()
