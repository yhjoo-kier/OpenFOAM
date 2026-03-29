#!/usr/bin/env python3
"""Multi-trial VLM evaluation runner.
Usage: python run_multitrial_eval.py <model> <trial> [mode]
  model: claude | gemini | gpt5
  trial: 1 | 2 | 3
  mode: setup | imageonly | both (default: both)
"""
import os, json, base64, time, datetime, sys
from pathlib import Path

# Load env
for line in open(Path(__file__).parent.parent / '.env'):
    k, v = line.strip().split('=', 1)
    os.environ[k] = v

BASE = Path(__file__).resolve().parent.parent / 'benchmark'
ITEMS_PATH = BASE / 'labels' / 'vlm_eval_claude_image_only.json'
V2_SCENES = {'S1', 'S6', 'S9'}

def get_image_path(case):
    scene = case.split('_')[0]
    suffix = '_V2.png' if scene in V2_SCENES else '_V6.png'
    return str(BASE / 'images' / case / f'{case}{suffix}')

def get_setup_text(case):
    meta = json.load(open(BASE / 'cases' / case / 'case_meta.json'))
    return meta.get('description', json.dumps(meta))

def parse_verdict(text):
    fl = text.strip().split('\n')[0].lower()
    if 'anomaly' in fl: return 'Anomaly'
    if 'ok' in fl: return 'OK'
    return 'unknown'

def make_prompt(mode, case):
    if mode == 'setup':
        setup = get_setup_text(case)
        return (f"Problem setup: {setup}. Is this flow field physically plausible "
                "for the stated conditions? Answer on the first line: OK (plausible) "
                "or Anomaly (non-physical). Then provide brief reasoning.")
    else:
        return ("You have NO problem setup information. Based on visual appearance alone, "
                "does this flow/temperature field look physically plausible? "
                "Answer on the first line: OK or Anomaly. Then provide brief reasoning.")

# --- Model-specific API callers ---

def call_claude(img_path, prompt):
    import anthropic
    client = anthropic.Anthropic()
    img_data = base64.standard_b64encode(open(img_path, 'rb').read()).decode('utf-8')
    resp = client.messages.create(
        model='claude-opus-4-6', max_tokens=300,
        messages=[{'role':'user','content':[
            {'type':'image','source':{'type':'base64','media_type':'image/png','data':img_data}},
            {'type':'text','text':prompt}
        ]}])
    return resp.content[0].text

def call_gemini(img_path, prompt):
    from google import genai
    from google.genai import types
    client = genai.Client(api_key=os.environ['GEMINI_API_KEY'])
    img_bytes = open(img_path, 'rb').read()
    resp = client.models.generate_content(
        model='gemini-3.1-pro-preview',
        contents=[types.Part.from_bytes(data=img_bytes, mime_type='image/png'),
                  prompt])
    return resp.text

def call_gpt5(img_path, prompt):
    import subprocess
    full_prompt = f"Read the image at {img_path}. {prompt}"
    result = subprocess.run(
        ['codex', 'exec', full_prompt, '--full-auto'],
        capture_output=True, text=True, timeout=120)
    return result.stdout

def call_gpt5_api(img_path, prompt):
    from openai import OpenAI
    client = OpenAI(api_key=os.environ['OPENAI_API_KEY'])
    img_data = base64.standard_b64encode(open(img_path, 'rb').read()).decode('utf-8')
    resp = client.chat.completions.create(
        model='gpt-5.4',
        max_completion_tokens=300,
        messages=[{'role':'user','content':[
            {'type':'image_url','image_url':{'url':f'data:image/png;base64,{img_data}'}},
            {'type':'text','text':prompt}
        ]}])
    return resp.choices[0].message.content

CALLERS = {'claude': call_claude, 'gemini': call_gemini, 'gpt5': call_gpt5, 'gpt_api': call_gpt5_api}
MODEL_NAMES = {'claude': 'claude-opus-4-6', 'gemini': 'gemini-3.1-pro-preview', 'gpt5': 'gpt-5.4-codex', 'gpt_api': 'gpt-5.4'}
DELAY = {'claude': 3, 'gemini': 2, 'gpt5': 5, 'gpt_api': 2}

def run_eval(model_key, trial, mode):
    caller = CALLERS[model_key]
    model_name = MODEL_NAMES[model_key]
    delay = DELAY[model_key]

    with open(ITEMS_PATH) as f:
        items = json.load(f)['items']

    results = []
    for i, item in enumerate(items, 1):
        case, gt = item['case'], item['gt']
        img_path = get_image_path(case)
        prompt = make_prompt(mode, case)

        try:
            text = caller(img_path, prompt)
            verdict = parse_verdict(text)
            reasoning = text.strip().split('\n', 1)[-1].strip() if '\n' in text else ''
        except Exception as e:
            err_str = str(e)
            if '429' in err_str or 'rate' in err_str.lower() or 'quota' in err_str.lower():
                print(f"  RATE LIMITED at {i}/30 — saving partial results")
                break
            verdict, reasoning = 'error', err_str
            print(f"  {i}/30 ERROR: {err_str[:80]}")

        match = (verdict == 'OK') if gt == 'correct' else (verdict == 'Anomaly')
        results.append({'id':item['id'],'case':case,'gt':gt,'verdict':verdict,'match':match,'reasoning':reasoning})
        status = '✓' if match else '✗'
        print(f"  {i}/30 — {verdict} — {case} {status}")
        if i < len(items):
            time.sleep(delay)

    valid = [r for r in results if r['verdict'] not in ('error','unknown')]
    correct = sum(1 for r in valid if r['match'])
    n = len(valid)
    fp = sum(1 for r in valid if r['verdict']=='Anomaly' and r['gt']=='correct')
    fn = sum(1 for r in valid if r['verdict']=='OK' and r['gt']!='correct')

    output = {
        'evaluator': model_name, 'trial': trial,
        'task': 'setup_conditioned' if mode=='setup' else 'image_only',
        'date': datetime.date.today().isoformat(),
        'total': n, 'correct': correct,
        'accuracy': correct/n if n else 0,
        'fp': fp, 'fn': fn,
        'fn_cases': [r['case'] for r in valid if r['verdict']=='OK' and r['gt']!='correct'],
        'fp_cases': [r['case'] for r in valid if r['verdict']=='Anomaly' and r['gt']=='correct'],
        'items': results,
        'partial': len(results) < 30,
    }

    out_path = BASE / 'labels' / f'vlm_eval_{model_key}_api_t{trial}_{mode}.json'
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\n  Saved: {out_path.name}")
    acc_pct = correct/n*100 if n > 0 else 0.0
    print(f"  {model_name} Trial {trial} {mode}: {correct}/{n} = {acc_pct:.1f}% (FP={fp}, FN={fn})")
    return output

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    model_key = sys.argv[1]
    trial = int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else 'both'

    print(f"=== {MODEL_NAMES[model_key]} Trial {trial} ===")

    if mode in ('setup', 'both'):
        print(f"\n--- Setup-conditioned ---")
        run_eval(model_key, trial, 'setup')

    if mode in ('imageonly', 'both'):
        print(f"\n--- Image-only ---")
        run_eval(model_key, trial, 'imageonly')

    print("\nDone.")
