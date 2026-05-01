# Untracked artifact cleanup

> Date: 2026-05-02 KST

## Context

`~/projects/OpenFOAM` 작업 시작 전 `git status` 기준 untracked 항목이 605개 존재했다. 대부분 OpenFOAM/benchmark 실행 산출물, local agent/QC artifacts, LaTeX scratch, paper result candidates였다.

## Actions performed

### 1. Inventory

- Initial untracked total: 605
- Largest groups:
  - `cases/`: 401 entries, 약 8.8 GiB
  - `.omx/`: 77 entries
  - `benchmark/`: 71 entries
  - `papers/`: 33 entries, 약 609 MiB
  - root scratch images / LaTeX outputs

### 2. `.gitignore` strengthened

Added ignore rules for clearly reproducible/generated artifacts:

- generated OpenFOAM/evaluation cases
  - `cases/eval_*/`
  - `cases/phase2_ref_*/`
  - `cases/phase2_pred_*/`
  - `cases/grid_test_*/`
  - `cases/bench_*_ms*/`
  - `cases/demo_case*/`
  - `cases/cfd_metrics_*.json`
- benchmark generated evaluation roots
  - `benchmark/evaluations_*/`
  - `benchmark/manifests/grid_test_logs/`
- local/QC artifacts
  - `.omx/`
  - root-level QC scratch images
  - root-level `main.{aux,log,out,spl,pdf}`
- paper generated artifacts
  - `papers/**/.env`
  - `papers/**/benchmark/cases/`
  - `papers/**/results/paperbanana_*/`
  - `papers/**/results/paperbanana*.pdf`
  - `papers/**/latex/files_for_submission/figshare_package/`
  - `papers/**/latex/files_for_submission/*.zip`
  - `papers/**/latex/files_for_submission/texput.log`
  - `*:Zone.Identifier`

No files were deleted. No commit was made.

## Result

After the `.gitignore` update:

- Remaining untracked entries: 21
- Tracked change: `.gitignore`

Remaining untracked items are intentionally left visible because they look like possible research records, lightweight summaries, scripts, labels, or submission artifacts that need a human decision before ignoring/deleting.

## Remaining review candidates

### Likely commit/review candidates

- `docs/26-03-20_posthoc_scaled_full_aggregate.md`
- `docs/26-03-20_posthoc_scaled_full_comparison.md`
- `docs/26-03-20_scale_hint_post_view_guarded_gate_note.md`
- `scripts/launch_all_grid_tests.sh`
- `scripts/run_grid_test.sh`
- `benchmark/manifests/*.json` posthoc/grid summary files
- `papers/cfd_visual_qa/benchmark/labels/*.json`
- `papers/cfd_visual_qa/docs/handoff_2026-03-29.md`
- `papers/cfd_visual_qa/docs/supervisor_review_log.md`
- `papers/cfd_visual_qa/docs/superpowers/`

### Submission files confirmed by user

- `papers/cfd_visual_qa/latex/files_for_submission/declarationStatement.docx`
- `papers/image_to_cfd/latex/files_for_submission/declarationStatement.docx`

User confirmed these two `.docx` files are Elsevier submission files. They should be preserved and treated as commit candidates, not ignored/deleted.

## Second-pass review (2026-05-02)

Remaining untracked entries were reviewed at a lightweight level.

- Secret-like pattern scan on remaining untracked text files: **0 findings** for common OpenAI/Anthropic/Google key patterns or generic key/token assignments.
- `papers/cfd_visual_qa/docs/handoff_2026-03-29.md` mentions that API keys exist in `.env`, but does not contain the key values. `.env` itself is ignored.
- Remaining files are small enough to be versioned except the two `.docx` files, which are submission artifacts intentionally retained.

## Recommended commit bundles (do not commit without explicit user request)

1. **Repository hygiene**
   - `.gitignore`
   - `docs/26-05-02_untracked_artifact_cleanup.md`

2. **Image-to-CFD benchmark/posthoc records**
   - `docs/26-03-20_posthoc_scaled_full_aggregate.md`
   - `docs/26-03-20_posthoc_scaled_full_comparison.md`
   - `docs/26-03-20_scale_hint_post_view_guarded_gate_note.md`
   - `benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
   - `benchmark/manifests/evaluation_batch_summary_posthoc_scaled_full.json`
   - `benchmark/manifests/phase1_grid_test_results.json`
   - `benchmark/manifests/posthoc_scaled_batch_A1.json`
   - `benchmark/manifests/posthoc_scaled_batch_A2.json`
   - `benchmark/manifests/posthoc_scaled_batch_A3.json`
   - `benchmark/manifests/posthoc_scaled_batch_A4.json`
   - `benchmark/manifests/posthoc_scaled_vs_baseline_comparison.json`

3. **Grid independence rerun scripts**
   - `scripts/launch_all_grid_tests.sh`
   - `scripts/run_grid_test.sh`

4. **CFD Visual QA labels and paper-management records**
   - `papers/cfd_visual_qa/benchmark/labels/vlm_eval_gemini_isolated_imageonly.json`
   - `papers/cfd_visual_qa/benchmark/labels/vlm_eval_gemini_isolated_setup.json`
   - `papers/cfd_visual_qa/benchmark/labels/vlm_eval_gpt54_raw.json`
   - `papers/cfd_visual_qa/docs/handoff_2026-03-29.md`
   - `papers/cfd_visual_qa/docs/superpowers/`
   - `papers/cfd_visual_qa/docs/supervisor_review_log.md`

5. **Elsevier submission declarations**
   - `papers/cfd_visual_qa/latex/files_for_submission/declarationStatement.docx`
   - `papers/image_to_cfd/latex/files_for_submission/declarationStatement.docx`
