# CFD Visual QA — Supervisor Review Log

> Daily supervisory reviews for the autonomous research workflow.
> Default cadence: once per day (target 20:30 KST).
> Purpose: track progress, identify research/paper risks, and suggest next priorities.

---

## 2026-04-03 20:30 KST — Supervisor Review

### Progress
- The project remains in the correct late-stage phase: manuscript hardening and submission packaging, not further benchmark expansion. Recent timestamps show active updates concentrated in manuscript-facing assets on 2026-04-01, including `latex/main.tex`, blinded/unblinded PDFs, submission metadata, cover-letter materials, highlights, graphical abstract assets, and synchronized figure copies under `latex/figures/`.
- There is a real milestone here: the work is no longer just a promising benchmark with internal notes. It has been shaped into a submission-grade bundle with manuscript PDF, blinded version, submission metadata, highlights, cover letter, and figure package.
- The strongest paper-facing narrative appears coherent in the files that matter most for submission. `docs/cover_letter.md`, `docs/highlights.md`, and the late-stage LaTeX bundle align around the API-isolated evaluation story, the setup-conditioned vs image-only rank reversal, and the contamination warning as a methodological contribution.
- However, the repository still contains clearly stale manuscript-adjacent text that directly conflicts with the submission narrative. In particular, `docs/paper_draft_v1.md` still presents the superseded high-accuracy Claude story (99.6% / majority vote 100%) and broader benchmark-era framing, while `docs/handoff_2026-03-29.md` explicitly states that those earlier subagent results were contaminated and that the API-isolated rerun is the valid paper story.

### Strengths
- **The paper thesis is now genuinely sharper than the original leaderboard-style framing.** The interesting result is not that one model got a high score, but that evaluator rankings reverse when setup text is removed, which says something meaningful about what these systems are actually doing.
- **Submission-facing assets are mostly aligned.** The cover letter and highlights are now consistent with the contamination-aware, setup-conditioned interpretation rather than the earlier near-perfect-Claude narrative.
- **The benchmark still looks paper-worthy.** Ten canonical CFD scenarios, controlled error types, expert comparison, and dual-condition evaluation together remain a credible benchmark/protocol package.
- **The project has already passed the most dangerous scientific integrity test.** A flattering but compromised result was replaced with a weaker but more defensible API-isolated result, and the paper still retained a publishable story.

### Concerns/Risks
- **Major paper-risk: stale manuscript text is still present in a dangerously reusable location.** `docs/paper_draft_v1.md` is not just an old log; it is a full paper draft with outdated claims, old contribution wording, and superseded results. This creates a real risk of accidental reuse in revisions, rebuttals, abstracts, captions, or response letters.
- **Canonicalization is still not fully solved.** The project clearly has at least two competing narratives inside the repo: the older benchmark-era “Claude nearly perfect” story and the newer API-isolated rank-reversal story. Even if LaTeX is intended as canonical, the repo does not yet enforce that strongly enough.
- **Title/framing risk remains.** Submission-facing text still uses the broad title “Can Vision-Language Models Evaluate CFD?” That remains catchier than the evidence is broad. The stronger, safer evidence is about setup-conditioned physical plausibility assessment of CFD visualizations on this benchmark.
- **State/provenance hygiene is still weak.** `research_state.json` remains historically inconsistent with the actual project state and still reflects the earlier fast-moving benchmark phase structure rather than a clean submission-era record.
- **Hidden reviewer risk: the stale draft could leak old numbers back into future edits.** This is exactly the kind of self-inflicted inconsistency that produces avoidable reviewer suspicion.

### Next Priorities
1. **Quarantine or archive `docs/paper_draft_v1.md` immediately.** Move it out of the active manuscript path, or add an unmistakable “OBSOLETE — superseded by API-isolated LaTeX manuscript” banner at the top.
2. **Make one submission-canonical evidence path explicit.** The repo should clearly indicate that the LaTeX manuscript plus its current referenced result files are the only valid paper-facing sources.
3. **Do a repo-wide stale-claim sweep.** Search for outdated percentages, majority-vote claims, old contamination-paradox wording, and pre-rerun result summaries across `docs/`, `results/`, and benchmark summaries.
4. **Tighten title precision if possible before submission.** A more specific title would reduce reviewer attack surface and better match the actual scope of the evidence.
5. **Repair archival state files.** Update `research_state.json` so the phase history and evaluation status no longer contradict the project’s true completion state.

### Paper Verdict
- **Overall verdict:** the project is in the appropriate final consolidation phase and has reached a meaningful manuscript-packaging milestone.
- **Supervisor judgment:** there is a **major paper-risk worth surfacing now**: stale manuscript text with superseded near-perfect results is still sitting in the active docs area, which undermines canonical evidence hygiene.
- **Reviewer-style bottom line:** the main remaining threat is not novelty, experimental scope, or lack of a paper story. It is internal inconsistency. If stale draft isolation and canonicalization are cleaned promptly, this remains a credible submission-grade benchmark/protocol paper.

---

## 2026-04-01 20:30 KST — Supervisor Review

### Progress
- The project is still in the right macro-phase: late manuscript hardening rather than benchmark expansion. Recent artifacts indicate active submission-style polishing, including maintained LaTeX outputs, updated result figures, graphical-abstract exploration, and a repo that now clearly contains both manuscript-facing assets and historical research traces.
- A meaningful milestone has already been reached at the evidence level: the paper is no longer just a benchmark build with impressive but shaky model scores. It now has a stronger story around setup-conditioned validation, rank reversal under image-only ablation, and the methodological necessity of API isolation.
- However, the archival state remains visibly split across at least two result eras. The LaTeX manuscript presents the newer API-isolated story, while multiple docs under `docs/` and some handoff/planning artifacts still preserve older 99.6--100% Claude narratives, contaminated-evaluation interpretations, or stale pilot framing.
- `research_state.json` remains structurally inconsistent with the actual project history: it marks `current_phase` as `P5_paper_draft` while `P4_full_eval` is `not_started`, and it still carries older totals/daily snapshots rather than a fully reconciled submission-era state.

### Strengths
- **The paper thesis is now genuinely interesting.** The strongest contribution is not a leaderboard result but a methodological/behavioral finding: setup-conditioned scientific validation behaves differently from image-only visual judgment, and evaluator ranking reverses when textual setup is removed.
- **The benchmark contribution remains paper-worthy.** Ten canonical CFD scenarios, controlled error types, expert comparison, and dual-condition evaluation together make a coherent benchmark/protocol package.
- **The contamination issue was discovered early enough to strengthen rather than sink the paper.** Catching and replacing contaminated evaluation with API-isolated evaluation materially improves credibility.
- **The current evidence supports a narrower, more defensible claim.** The manuscript appears strongest when framed as setup-conditioned CFD visual validation on this benchmark, with model-specific capability differences rather than broad claims about all VLMs.

### Concerns/Risks
- **Major paper-risk: stale manuscript-adjacent artifacts still conflict with the current paper story.** Files such as `docs/paper_draft_v1.md`, `docs/paper_draft_results_discussion.md`, `docs/paper_draft_intro_related.md`, `docs/26-03-26_paper_outline.md`, and `docs/handoff_2026-03-29.md` still contain superseded claims or numbers (e.g., Claude 99.6--100%, GPT-5.4 at 80.0%, older contamination framing). Even if these are "historical," they are close enough to the manuscript workflow to create real risk of accidental reuse in revision, rebuttal, captions, or cover-letter updates.
- **Major design/framing risk: the title still overclaims relative to the evidence.** `Can Vision-Language Models Evaluate CFD?` remains catchy but still invites a reviewer to test the broadest interpretation. The actual evidence supports a more precise claim about physical plausibility assessment of CFD visualizations under setup-conditioned evaluation, with strong model dependence.
- **Canonical-data risk: `canonical_results_summary.json` is itself stale relative to the newer API-isolated narrative.** It still reports the older near-perfect Claude setup-conditioned result and different image-only values, which means the file named as canonical is not actually canonical anymore. That is a hidden but serious provenance problem.
- **Project-state/provenance risk remains unresolved.** `research_state.json` and daily logs still encode the accelerated exploratory history rather than a clean submission-era state. This is not fatal for the paper text, but it weakens auditability and future revision/rebuttal hygiene.
- **Reviewer-positioning risk: the current conclusion can still be misread as "VLMs lack visual physics understanding" in a too-general sense.** The safer interpretation is that, on this benchmark, the strongest observed VLM capability is setup--image cross-referencing; pure image-only physical plausibility remains substantially weaker and model-dependent.

### Next Priorities
1. **Fix canonicalization immediately.** Either replace `benchmark/labels/canonical_results_summary.json` with the API-isolated final numbers or create a new clearly named submission-canonical summary and point every manuscript-facing asset to it.
2. **Quarantine stale drafts.** Move superseded Markdown drafts/handoffs/outlines into an `archive/` or explicitly mark them as obsolete so they cannot be mistaken for current paper sources.
3. **Tighten the title and headline language.** Consider wording that foregrounds setup-conditioned physical plausibility assessment rather than the broader `evaluate CFD` phrasing.
4. **Repair provenance files.** Update `research_state.json` so full evaluation is actually marked complete and the current phase/history matches what the repo now contains.
5. **Do one final manuscript-consistency audit.** Check abstract, cover letter, highlights, graphical abstract text, captions, and supplementary descriptions against the same final result snapshot.
6. **Preempt reviewer pushback more explicitly.** In Discussion/Limitations, directly acknowledge model-count limits, single-expert baseline, benchmark templating concerns, and the intended real-world nature of setup-conditioned validation.

### Paper Verdict
- **Overall verdict:** the project remains in an appropriate final consolidation phase, and the benchmark/paper package is strong enough to justify continued manuscript polishing rather than further benchmark growth.
- **Supervisor judgment:** there is a **major paper-risk that should be surfaced**: the repo still lacks true canonicalization because stale drafts and even the nominal `canonical_results_summary.json` conflict with the current API-isolated manuscript story.
- **Reviewer-style bottom line:** this is close to a credible submission-grade benchmark/protocol paper, but the remaining danger is not lack of novelty. It is evidence hygiene. Until canonical results, title precision, and stale-draft isolation are cleaned up, the paper remains unnecessarily exposed to self-inflicted inconsistency.

---

## 2026-03-31 20:30 KST — Supervisor Review

### Progress
- The project appears to have reached a genuine manuscript-consolidation stage rather than open-ended exploration. The LaTeX manuscript (`latex/main.tex`) and compiled PDF were updated on 2026-03-30, indicating the paper has been carried into a submission-style artifact rather than remaining only in Markdown drafts.
- The current manuscript story is substantially sharper than the early benchmark-only framing: it centers on API-isolated evaluation, contamination risk, and the rank reversal between setup-conditioned and image-only CFD visual QA.
- Recent activity on 2026-03-31 shows continued figure/protocol polishing work under `results/paperbanana_eval_protocol_*`, which suggests the team is improving paper communication quality rather than merely adding more raw experiments.
- The benchmark asset base remains strong and stable: 60 CFD cases, 258 images, 279 QA items, and 80 expert labels are already sufficient for a focused benchmark/protocol paper if the claims remain disciplined.

### Strengths
- **The paper now has a clear, reviewer-relevant thesis.** The strongest contribution is no longer “a model got a high score,” but rather: setup-conditioned scientific visual validation behaves differently from image-only judgment, and evaluation isolation materially changes the reported capability.
- **The LaTeX manuscript seems internally aligned with the updated contamination-aware narrative.** The abstract and contribution statements consistently emphasize API isolation, model-specific setup dependence, and rank reversal, which is much stronger than the earlier over-simple high-accuracy headline.
- **The current experimental framing is appropriate for a benchmark/protocol paper.** Ten canonical scenarios, controlled error taxonomy, expert comparison, and dual-condition evaluation together form a coherent study design.
- **The contamination finding remains a standout methodological contribution.** Many papers would have hidden this; here it is becoming part of the scientific contribution, which increases credibility if presented carefully.

### Concerns/Risks
- **Major paper-risk: stale parallel artifacts still likely exist in the repo.** Older Markdown drafts and logs remain alongside the LaTeX paper. If any of those are reused during revision, caption extraction, or response writing, the project could accidentally reintroduce outdated numbers or pre-isolation claims. The canonical source must be explicit.
- **Major design/interpretation risk: the title still overgeneralizes.** `Can Vision-Language Models Evaluate CFD?` is attractive, but the evidence more precisely supports: frontier VLMs show model-dependent capability for **setup-conditioned CFD visual validation** on this benchmark. As written, the title risks inviting reviewers to challenge the broad “evaluate CFD” framing.
- **Evidence breadth risk remains.** The core API-isolated comparison is compelling, but still modest in scale for sweeping generalization across all VLMs. The manuscript should continue to avoid category-wide claims and foreground benchmark/protocol contribution over model leaderboard rhetoric.
- **Bookkeeping inconsistency remains a hidden archival risk.** `research_state.json` still reflects an awkward phase history (e.g., P5 in progress while P4 is not started), which may not hurt the submitted paper directly but does weaken project auditability and could confuse later revisions or rebuttal prep.
- **Recent PaperBanana-style figure-generation outputs suggest a communication sprint, but also a possible distraction risk.** At this stage, figure aesthetics should not outrun evidence hardening. Better-looking protocol figures help, but they are secondary to claim precision, canonicalization, and denominator clarity.

### Next Priorities
1. **Declare one canonical paper source immediately.** Treat `latex/main.tex` + its referenced result files as the only submission-valid artifact, and mark older Markdown draft files as historical/outdated.
2. **Tighten the title and headline claim if possible.** Consider wording closer to setup-conditioned validation / physical plausibility assessment, rather than the broader “evaluate CFD” phrasing.
3. **Do a final repo-wide stale-number audit.** Search for old claims, old sample counts, contaminated-result percentages, and outdated captions across `docs/`, `results/`, and any release-facing text.
4. **Repair state provenance files.** Update `research_state.json` and, if needed, the daily/supervisor logs so the archival record reflects that full evaluation and manuscript consolidation actually occurred.
5. **Keep new work narrowly scoped to paper-hardening.** Prioritize abstract/title/discussion claim discipline, figure-caption consistency, supplement planning, and reproducibility notes over any new exploratory experiments.
6. **Preempt reviewer pushback in Discussion.** Explicitly acknowledge model-count limits, single-expert baseline, and the fact that setup-conditioned validation is the intended task rather than pure image-only physical intelligence.

### Paper Verdict
- **Overall verdict:** the project is in an appropriate late-stage paper-consolidation phase and has crossed the threshold of a meaningful manuscript milestone.
- **Supervisor judgment:** no direction change is needed, but there is a **major paper-risk worth surfacing**: the manuscript is now strong enough that the main remaining danger is not novelty or lack of results, but accidental overclaiming / stale-artifact leakage from older drafts.
- **Reviewer-style bottom line:** this looks like a credible benchmark-and-methodology paper with a sharp empirical message. I would support moving forward, provided the team now treats canonicalization, title/claim precision, and archival cleanup as top-priority near-final tasks.

---

## 2026-03-25 20:30 KST — Supervisor Review

### Progress
- Phase-level execution speed is far ahead of the original directive schedule. All 10 scenarios appear generated with 60 CFD cases, 258 images, and 279 QA items already prepared.
- Pilot-to-preliminary transition produced meaningful early evidence rather than only infrastructure: 35 expert labels collected, blind-code evaluation pipeline built, Google Drive DOCX labeling workflow established, and at least two VLM blind evaluation files produced.
- A genuinely useful methodological correction was made early: contaminated evaluation was recognized and replaced with blind-coded evaluation. This is a strong sign that the workflow is capable of self-correction.
- Current state bookkeeping is internally inconsistent: `current_phase` is `P2_preliminary_eval`, `P3_full_expansion` is marked `complete`, and `current_day` remains `1` despite work spanning beyond the original day plan. This is not fatal, but it weakens auditability and later paper-method reproducibility if not cleaned now.

### Strengths
- The project has already passed the “toy only” stage. There is a nontrivial benchmark structure across multiple CFD phenomena, multiple error types, and multiple visualization modes.
- The most promising result so far is not the raw 100% Claude score itself, but the observation that blind evaluation plus explicit problem setup can expose subtle inconsistencies that even an expert missed. That is a publishable angle if framed carefully.
- Error taxonomy quality is decent. The current set mixes obvious failures (gravity flip, under-convergence) with subtle but physically important ones (BC swap, wrong viscosity, wrong turbulence model, coarse mesh).
- Expert feedback is already informative enough to improve benchmark design: some question prompts lack critical setup information (e.g. lid direction, Reynolds number), and some visualizations are too dense to support fair judgment. This is exactly the kind of iteration a serious benchmark paper needs.

### Concerns/Risks
- **Major paper-risk: evidence is still too thin for a strong claim.** The headline result (Claude blind 15/15) is based on a very small pilot subset. It is encouraging, but not yet robust evidence for a paper-level conclusion.
- **Claim inflation risk:** if the paper starts sounding like “VLM surpasses experts in CFD QA,” reviewers will immediately attack sample size, single-expert dependence, and task construction. The safer claim is currently: “blind, setup-aware VLM evaluation shows promise on a pilot benchmark and can catch subtle inconsistencies missed by human inspection.”
- **Ground-truth / protocol risk:** several “subtle” detections depend strongly on the textual setup being explicit. That is acceptable, but the paper must clearly define whether the task is (a) image-only visual QA, or (b) setup-conditioned visual QA. Mixing the two will create reviewer pushback.
- **Benchmark validity risk:** some error classes are not uniformly detectable across scenarios (e.g. coarse mesh, wrong viscosity, wrong turbulence model). That is scientifically fine, but it means the benchmark difficulty is highly scenario-dependent. This must be analyzed explicitly rather than averaged away.
- **Labeling design risk:** only 35 expert labels are complete against a target of ~1,000–1,200 QA. If the paper relies on expert-vs-VLM comparison, label coverage and sampling policy need to be specified now, not later.
- **State-management risk:** `research_state.json` and `daily_log.md` are not fully synchronized. If autonomous progress continues at this pace without bookkeeping cleanup, later reconstruction of what was evaluated when may become messy.

### Next Priorities
1. **Freeze the task definition now.** Decide explicitly whether the main benchmark is “setup-conditioned visual QA” or “image-only visual QA.” My recommendation: make setup-conditioned QA the main task, and treat image-only judgment as a secondary harder variant.
2. **Expand blind evaluation before expanding claims.** Run a larger blind subset across both stronger and weaker models with the same protocol, and report abstentions/no-response separately rather than hiding them under valid-only accuracy.
3. **Standardize prompt completeness.** Every item should include the minimum physics context needed for fair judgment (e.g. Re, BC direction, heated wall identity, lid direction, expected regime). Missing setup information currently confounds some results.
4. **Define paper-safe metrics.** Beyond raw accuracy, track per-error-type recall, per-scenario accuracy, no-response rate, and calibration/uncertainty if available. Reviewers will want to know where the models fail, not just the average score.
5. **Repair bookkeeping immediately.** Synchronize `research_state.json`, daily log, and phase progression so the autonomous record matches the actual workflow.
6. **Use expert labels strategically.** Given the labeling budget, prioritize a stratified subset covering all scenarios, all main error families, and both easy/subtle cases, rather than spreading labels too thinly.

### Paper Verdict
- **Overall verdict:** promising and worth continuing, but **not yet ready for a strong paper claim**.
- **Phase appropriateness:** yes, the work is in a reasonable preliminary-evaluation stage even though data generation has already raced ahead. The current bottleneck is no longer case generation; it is evaluation rigor and claim discipline.
- **Supervisor judgment:** continue in the current direction **without changing the core topic**, but tighten protocol definition before presenting milestone claims externally.
- **Reviewer-style bottom line:** the project currently has a believable benchmark-and-analysis paper path, but the main hidden danger is overclaiming from a small pilot with mixed task definitions. If that is controlled early, the work has real paper potential.

---

## 2026-03-26 20:30 KST — Supervisor Review

### Progress
- The project has advanced from pilot-stage benchmarking to a materially paper-shaped state in just two days: 60 CFD cases, 258 images, 279 QA items, 80 expert labels, expanded Claude blind evaluation, Gemini retry handling, and a full 7-section paper draft are now on record.
- Relative to yesterday’s review, the most important milestone is that the benchmark is no longer supported only by a 15–30 item pilot. The daily log now reports **80 expert-labeled items completed** and **75 Claude blind-evaluated items**, which is a meaningful step toward paper-grade evidence.
- The main task definition has improved. The directive now clearly positions **setup-conditioned visual QA as the primary task** and image-only QA as a secondary harder variant. This resolves one of yesterday’s most important protocol ambiguities.
- However, the project record is currently fragmented: `research_state.json`, `daily_log.md`, and the paper draft do not describe the same evaluation snapshot. The draft still presents a 30-item result table, while the daily log claims 75-item Claude evaluation and 80-item expert coverage.

### Strengths
- Fast execution has not been purely operational; there is real methodological maturation. The blind protocol, contamination-aware framing, setup-conditioned task definition, and explicit no-response handling are all paper-relevant improvements.
- The benchmark has a credible structure: multiple canonical CFD scenarios, multiple visualization types, multiple controlled error classes, and a task formulation tied to an actual automation use case (result validation).
- The emerging result pattern is interesting enough for a paper: VLM strength seems to come from **systematic setup-image cross-checking**, whereas the human expert is more vulnerable on plausible-looking but setup-inconsistent errors.
- The project already contains strong paper assets beyond raw scores: a benchmark taxonomy, a reusable evaluation protocol, example figures, and a concrete lesson about evaluation contamination.

### Concerns/Risks
- **Major paper-risk: internal result inconsistency.** Right now there are at least three “truths” in the repo: (1) 30-item pilot-style paper draft, (2) 75-item Claude / 80-item expert daily-log result, and (3) `research_state.json` still indicating `current_phase: P2_preliminary_eval` while P3 is complete and P5 is in progress. A reviewer will not see the logs, but this inconsistency is dangerous because it can leak into tables, captions, and methods.
- **Major interpretation risk: the strongest narrative is still based mainly on Claude.** Gemini has reliability problems and GPT-5.4 image evaluation failed. That means the paper currently supports a claim about **at least one frontier VLM** more strongly than a broad claim about “VLMs” in general.
- **Protocol risk: 3-trial independence is now part of the directive, but the manuscript text still reads mostly like single-pass evaluation.** I can see trial2/trial3 files in the benchmark folder, but the methodology/results sections do not yet integrate this into the headline metrics. If three-trial majority-vote is the official protocol, the paper must report it clearly; if not, the directive and evaluation assets are drifting apart.
- **Benchmark novelty risk:** as currently written, the paper leans heavily on benchmark creation + one strong model result. That can be publishable, but only if the benchmark itself is framed as the main contribution. If the paper instead over-centers “Claude beats expert,” reviewers may dismiss it as underpowered and model-specific.
- **Ground-truth subtlety risk:** some error classes (wrong viscosity, BC swap, reversed lid/gravity) are not “visual errors” in the pure image sense; they are **setup-inconsistency errors**. This is okay, but the paper must say so explicitly and avoid implying that the image alone encodes the anomaly.
- **Human baseline risk:** there is still only one expert baseline, and the expert labeling route mixes convenience sampling and staged collection. This is acceptable for an initial benchmark paper, but the sampling and scope must be described honestly.

### Next Priorities
1. **Freeze a single canonical evaluation snapshot for the paper.** Decide whether the headline result is 30 items, 75 items, or another finalized subset—and make `paper_draft_*`, tables, and logs all match that snapshot.
2. **Integrate the multi-trial protocol explicitly.** If trial2/trial3 are real and usable, report majority-vote accuracy, consistency buckets (3/3, 2/3, 1/3), and single-trial variance. This would materially strengthen the paper.
3. **Narrow the claim language.** Prefer “frontier VLMs show promise for setup-conditioned CFD validation” over “VLMs can evaluate CFD” unless multiple robust models support the broader claim.
4. **Reframe the contribution hierarchy.** Put the benchmark/protocol first, model comparison second. The benchmark is the stable contribution; any one model’s score is less stable.
5. **Update the Results and Discussion sections before further expansion.** Right now the manuscript is the weakest-synchronized artifact in the repo. Fix tables, sample counts, and wording before adding more experiments.
6. **Audit stratification of the 75–80 evaluated items.** Confirm coverage across scenarios, error types, and severity levels; otherwise a strong overall accuracy could still be hiding a biased subset.
7. **Handle Gemini and GPT carefully in the narrative.** Gemini’s no-response behavior is itself a result; GPT-5.4 image failure is a tooling limitation, not necessarily a model weakness. Keep those distinctions clean.

### Paper Verdict
- **Overall verdict:** real milestone reached, but the paper is entering a dangerous phase where writing can outrun protocol consolidation.
- **Phase appropriateness:** yes—moving into paper drafting is appropriate **only because** data generation and initial evaluation are already unusually far along. But from this point, rigor and synchronization matter more than speed.
- **Supervisor judgment:** continue in the same direction; no core topic pivot is needed. However, I would treat **result synchronization + trial-protocol integration** as near-blocking paper tasks before any aggressive external claim.
- **Reviewer-style bottom line:** this now looks like a potentially publishable benchmark/protocol paper with an interesting early finding, not just a speculative pilot. The hidden danger is not lack of novelty; it is letting an attractive headline outrun a fully consistent evidence package.

---

## 2026-03-28 20:30 KST — Supervisor Review

### Progress
- The project has now crossed from a fast autonomous prototype into a credible paper package: the benchmark asset set is substantial (60 CFD cases, 258 images, 279 QA), the paper draft is materially synchronized with a canonical result summary, and the three-trial Claude evaluation has been consolidated into a single reproducible JSON summary.
- A meaningful maturation since the previous review is that the manuscript is no longer anchored to the earlier 30-item pilot framing. The current draft and `canonical_results_summary.json` both describe the 75-item / 3-trial Claude result, the 30-item / 3-trial Gemini result, and the 80-item expert baseline.
- The benchmark contribution is also clearer now. The work is not merely “Claude got a high score”; it includes an explicit setup-conditioned task definition, a blind protocol, an image-only ablation, per-error-type recall, and a contamination-paradox finding that is methodologically interesting in its own right.
- One workflow inconsistency remains: `research_state.json` marks `current_phase` as `P5_paper_draft` while `P4_full_eval` remains `not_started`. Given the repo contents, P4 has effectively already happened. This is a bookkeeping issue, but it is the sort of inconsistency that can later leak into the Methods section or a response-to-reviewers.

### Strengths
- **The paper now has a stable main contribution.** The strongest publishable unit is the benchmark plus evaluation protocol for setup-conditioned CFD visual validation, not any single model score.
- **The setup-conditioned vs image-only gap is a strong result.** It clearly establishes what kind of reasoning is being tested: not generic image description, but consistency checking between physical setup and simulation output.
- **The multi-trial Claude result is substantially stronger than a single-shot score.** Seventy-five items across three independent trials with near-perfect consistency is much more defensible than the earlier pilot snapshot, especially because the only observed Claude error is a stochastic false positive rather than a repeated miss on an error case.
- **The expert comparison is informative rather than decorative.** The key signal is not “human bad, AI good,” but that experts and VLMs fail differently—especially on BC swap and wrong-viscosity cases that remain visually plausible.
- **The benchmark design looks extensible.** Ten scenarios and a controlled error taxonomy create a scaffold that can later support additional models, more experts, or 3D/transient extensions.

### Concerns/Risks
- **Major paper-risk: claim scope can still outrun model coverage.** The evidence is strong for **Claude on this benchmark** and suggestive for Gemini, but not yet strong enough for broad wording like “VLMs can evaluate CFD” without qualification. What the current evidence robustly supports is: **a frontier VLM can perform setup-conditioned CFD visual validation surprisingly well on this benchmark**.
- **Reviewer skepticism risk around benchmark structure.** Even with the blind protocol, reviewers may worry that the error taxonomy is too templated or that the benchmark rewards text-image consistency checking more than deeper physical understanding. This does not invalidate the work, but the manuscript should preempt it by explicitly positioning setup-conditioned validation as the intended real-world task.
- **Human baseline interpretation risk.** The single-expert result is useful, but any framing that sounds like a general human-vs-AI superiority claim will likely trigger pushback. The comparison should be presented as a structured baseline with complementary failure modes.
- **Evaluation asymmetry risk.** Claude has 75 items and 3 trials, Gemini has 30 items and 3 trials, and the expert has 80 items and 1 trial. This is acceptable for an initial paper, but the asymmetry should be acknowledged explicitly as a limitation rather than hidden behind one aggregated table.
- **State/provenance risk remains nonzero.** `research_state.json`, daily log chronology, and manuscript timeline still look like an accelerated autonomous project rather than a fully cleaned archival record. Before submission, the repo needs one evidence-freeze pass so that counts, dates, and protocol language are perfectly aligned.
- **Novelty framing risk.** If the introduction leans too hard on “first benchmark” alone, reviewers may ask whether this is mainly a benchmark paper or mainly a model-capability paper. The safer and stronger framing is a benchmark/protocol paper with a compelling early empirical result.

### Next Priorities
1. **Freeze the paper-safe claim language.** Replace any broad wording with benchmark-first, model-specific wording.
2. **Do one final evidence-consistency audit.** Align `research_state.json`, manuscript tables/captions, figure labels, and any release text to the same canonical snapshot (`canonical_results_summary.json`).
3. **Strengthen reviewer preemption in the Discussion.** Add a short paragraph explicitly addressing likely criticisms: templated errors, single-expert baseline, model-count asymmetry, and setup-conditioned nature of the task.
4. **Clarify contribution hierarchy in Abstract / Introduction / Conclusion.** Lead with benchmark + protocol + analysis; present the high Claude score as an empirical finding enabled by the benchmark, not the sole reason the paper matters.
5. **Prepare a minimal appendix or supplement plan.** Even if not written today, identify which artifacts should be supplementary: prompt templates, blind-code procedure, representative error examples, and trial-consistency details.
6. **Update project state formally.** Mark full evaluation as completed in the state record, or revise the phase labels so the project history reflects what actually happened.

### Paper Verdict
- **Overall verdict:** this is now a **meaningful milestone**. The project looks paper-credible as a benchmark/protocol study with a striking initial result, not just an autonomous exploration artifact.
- **Phase appropriateness:** yes. The right phase now is paper consolidation and evidence hardening, not further uncontrolled benchmark growth.
- **Supervisor judgment:** no core direction change is needed. Continue on the present track, with emphasis on claim discipline and archival cleanup rather than new expansion.
- **Reviewer-style bottom line:** I would not block this from moving toward manuscript polishing. But I would strongly insist that the paper be written as **“benchmark + setup-conditioned validation protocol + model-specific findings”** rather than **“VLMs broadly solve CFD validation.”** That distinction will strongly affect reviewer reception.

---

## 2026-03-29 20:30 KST — Supervisor Review

### Progress
- A major methodological correction occurred since the previous review: the earlier near-perfect subagent VLM results were found to be contaminated by filesystem-access leakage, and the benchmark was re-evaluated via direct API calls with base64 image transfer and zero ground-truth file access. This is an important credibility upgrade, even though it weakened the headline score.
- The project now has a much stronger and more defensible central story: **rank reversal** between setup-conditioned and image-only evaluation. In the updated isolated evaluation, Claude Opus 4.6 performs well with setup text (88.9%) but drops sharply without setup (33.3%), GPT-5.4 shows moderate setup-conditioned performance (80.0%) and weak image-only performance (43.3%), Gemini remains less stable, and the expert baseline sits between these behaviors (73.8% setup-conditioned, 66.7% image-only).
- Figure production and manuscript packaging appear substantially more mature. The session handoff reports a rewritten 18-page LaTeX manuscript, 10 figures, API-level raw evaluation files, and a committed branch state. This is a meaningful milestone: the work is no longer just a benchmark build plus draft notes; it is now a coherent paper candidate with a sharper thesis.
- However, the project documentation is now split across generations of the story. `docs/paper_draft_results_discussion.md` still contains the older 30-item / 100%-Claude narrative, while the 2026-03-29 handoff and raw API files describe the new isolated-evaluation results. Unless the repo clearly marks LaTeX as canonical and deprecates stale Markdown sections, internal inconsistency may still leak into the final paper.

### Strengths
- **The contamination was caught before submission.** That is one of the most important positive signals in the entire project. A weaker but valid result is vastly better than a spectacular but compromised one.
- **The revised paper thesis is stronger scientifically.** “Setup-image cross-referencing is a model-specific capability, and image-only plausibility judgment behaves differently” is more interesting and more publishable than a simple “Claude beats expert” claim.
- **The new results support a real benchmark/protocol contribution.** The benchmark now demonstrates not only performance ranking, but also how evaluation framing changes model behavior. That is a meaningful methodological insight for scientific-vision benchmarking.
- **The expert baseline became more informative after adding image-only evaluation.** The expert drops less from setup-conditioned to image-only than Claude does, which supports the interpretation that the expert relies more on gestalt plausibility while Claude benefits disproportionately from explicit setup-image matching.
- **Raw provenance is improving.** The presence of isolated API trial files, a unified evaluation runner, and a handoff summarizing the exact rerun basis all help the project move toward paper-auditability.

### Concerns/Risks
- **Major paper-risk: the paper narrative has changed, but stale artifacts remain.** Any residual text, caption, or table that still reports the older contaminated 99.6–100% Claude story is now dangerous. This is the highest immediate manuscript risk.
- **Major interpretation risk: sample size is still modest for a strong cross-model generalization claim.** The updated core comparison appears to be 30 items × 2 conditions × 3 trials for the API-isolated study. That is enough for a strong pilot-to-initial-paper finding, but still too small for broad claims about “VLMs” as a category.
- **Design risk: evaluator asymmetry remains.** The expert uses 80 setup-conditioned items and 30 image-only items; Claude/GPT/Gemini appear compared on 30-item setup and image-only subsets. This is acceptable if clearly stated, but it complicates direct headline comparisons and must be disclosed carefully.
- **Novelty/framing risk:** if the introduction still frames the work mainly as “first benchmark for VLM-based CFD visual QA,” reviewers may ask whether the stronger contribution is actually the benchmark, the setup-conditioned protocol, or the rank-reversal finding. The hierarchy should be explicit: benchmark + protocol first, empirical rank-reversal second.
- **Evidence-freeze risk:** the repo now contains multiple result eras (contaminated blind eval, expanded subagent eval, API-isolated rerun). Before submission, the project needs one clearly labeled archival cut that says which results are valid, which are deprecated, and which file is canonical.
- **Hidden reviewer risk:** some reviewers may interpret the strong setup-conditioned performance as mostly text-guided verification rather than deep physical understanding. This is not fatal, but the paper must proactively frame the task as *simulation validation against stated setup*, not pure image intelligence.

### Next Priorities
1. **Purge or clearly deprecate stale result text.** Any Markdown draft, figure caption, or log-derived manuscript snippet still carrying the contaminated high-accuracy story should be marked outdated or removed from the paper path.
2. **Freeze one canonical evidence package.** I recommend treating the API-isolated rerun as the only paper-valid evaluation set and explicitly naming it as such in the manuscript and repo.
3. **Clarify comparison denominators everywhere.** For each evaluator and condition, state item count, trials, and whether accuracy is majority-vote or per-trial mean. Do not allow readers to compare 80-item expert numbers directly against 30-item model numbers without context.
4. **Strengthen the Discussion around task nature.** Add a concise paragraph arguing that setup-conditioned CFD validation is the intended real-world use case, and that the image-only drop is therefore not an embarrassment but a revealing ablation.
5. **Tighten the claim language further.** The safe paper claim is now: frontier multimodal models show nontrivial capability for **setup-conditioned** CFD validation, but this capability is strongly model-dependent and degrades sharply without setup context.
6. **Add a short provenance note or appendix table.** Explicitly summarize the contaminated evaluation discovery, why it was invalid, and why the API-isolated rerun resolves that concern. Reviewers will appreciate the honesty.
7. **Clean the project state files.** `research_state.json`, daily logs, and paper-draft docs should point to the same phase history and result snapshot, or at minimum distinguish “historical log” from “submission-canonical result.”

### Paper Verdict
- **Overall verdict:** this is a **meaningful milestone worth surfacing**. The project just survived the kind of methodological stress test that kills weak papers: a flattering result was invalidated, and the work still retained a publishable, arguably stronger scientific story afterward.
- **Phase appropriateness:** yes. The correct phase now is not further benchmark expansion; it is manuscript hardening, provenance cleanup, and reviewer-preemption.
- **Supervisor judgment:** no core direction change is needed, but I do consider **canonical-result cleanup and stale-claim removal** to be near-blocking tasks before any submission or public sharing.
- **Reviewer-style bottom line:** after the API-isolated rerun, this looks less like “VLMs solved CFD QA” and more like a sharper, more credible paper about **setup-conditioned scientific visual validation**, model-specific multimodal reasoning, and the danger of evaluation contamination. That is a better paper than the old headline, provided the evidence package is cleaned thoroughly.

---

## 2026-04-02 20:30 KST — Supervisor Review

### Progress
- The project is clearly in the correct late-stage phase: submission packaging and evidence hardening, not benchmark expansion. Since the previous review, manuscript-facing artifacts were updated again on 2026-04-01, including `latex/main.tex`, blinded/unblinded PDFs, submission metadata, highlights, cover-letter materials, and refreshed result figures.
- A meaningful milestone has been reached: the paper now exists as a coherent submission bundle rather than a loose research draft. The main manuscript, blinded manuscript, compiled PDFs, and submission-side metadata are all present and recently updated.
- The main paper narrative is now internally aligned across the submission-facing files that matter most. `latex/main.tex`, `docs/highlights.md`, `docs/cover_letter.md`, and `latex/files_for_submission/submission_metadata.md` all consistently describe the API-isolated evaluation story with the setup-conditioned vs image-only rank reversal and the contamination warning.
- However, the project’s archival/control files remain inconsistent with the manuscript-facing story. `research_state.json` still reports `current_phase = P5_paper_draft` while `P4_full_eval` is `not_started`, and the file still encodes the older near-perfect Claude-era interpretation. The nominally canonical file `benchmark/labels/canonical_results_summary.json` is also stale and conflicts directly with the manuscript now under submission preparation.

### Strengths
- **The manuscript package is now real and submission-oriented.** This is no longer just a benchmark repo with promising notes; the project now has the artifacts expected for an actual journal submission workflow.
- **The core scientific story is sharper than the older leaderboard framing.** The strongest contribution is the discovery that setup-conditioned CFD validation and image-only plausibility judgment are different tasks with different evaluator rankings.
- **Methodological credibility improved rather than weakened after the contamination discovery.** The paper’s insistence on API isolation is now a genuine contribution, not just a cleanup detail.
- **The current experimental scope is appropriate for the claimed paper type** if the claims stay benchmark/protocol-first: 10 scenarios, controlled error taxonomy, expert comparison, dual-condition evaluation, and isolation-aware methodology together form a coherent EAAI-style contribution.

### Concerns/Risks
- **Major paper-risk: the repo still lacks a true canonical result source.** `benchmark/labels/canonical_results_summary.json` still reports the superseded 99.6% Claude narrative and outdated image-only values, while the manuscript now claims the API-isolated 88.9% / 33.3% story. A file explicitly named “canonical” disagreeing with the paper is a serious provenance flaw.
- **Major title/framing risk remains unresolved.** The manuscript title still asks, `Can Vision-Language Models Evaluate Computational Fluid Dynamics?` This is stronger and broader than the evidence actually supports. The paper more defensibly answers whether VLMs can perform **setup-conditioned physical plausibility assessment of CFD visualizations on this benchmark**.
- **State/provenance hygiene is still weak.** `research_state.json`, daily history, and canonical summary files remain out of sync with the final paper package. This may not block submission mechanically, but it creates avoidable rebuttal/revision risk and weakens auditability.
- **Residual reviewer risk: the conclusion can still be overread as a statement about “visual physics understanding.”** The evidence instead supports that the best-performing model is strong at setup--image cross-referencing, while pure image-only judgment remains weak to moderate and model-dependent.
- **Experimental asymmetry still needs careful handling.** The manuscript compares 30-item common subsets across models while the expert also has an 80-item broader setup-conditioned evaluation. This is acceptable, but must remain crystal-clear everywhere to avoid denominator confusion.

### Next Priorities
1. **Repair canonicalization immediately.** Replace or deprecate `benchmark/labels/canonical_results_summary.json` so the repo has one explicit submission-canonical result source matching the LaTeX manuscript.
2. **Tighten the title and headline phrasing.** A more precise title would reduce reviewer attack surface and better match the benchmark’s actual evidentiary scope.
3. **Synchronize provenance files.** Update `research_state.json` so phase completion and current status reflect the real project history.
4. **Do one final manuscript-adjacent audit.** Recheck abstract, title, highlights, cover letter, submission metadata, and figure captions against the same final numeric snapshot and denominator policy.
5. **Protect the benchmark-first framing.** Keep the paper positioned as a benchmark/protocol paper with model-specific findings, not a broad claim that VLMs generally “evaluate CFD.”

### Paper Verdict
- **Overall verdict:** the project has reached a meaningful manuscript milestone and is in the correct final consolidation phase.
- **Supervisor judgment:** there is still a **major paper-risk worth surfacing**—not lack of novelty, but lack of canonical evidence hygiene. Until the repo’s nominally canonical result/state files are reconciled with the submission manuscript, the paper remains unnecessarily exposed to internal inconsistency.
- **Reviewer-style bottom line:** this is close to submission-ready as a benchmark/protocol paper with a strong methodological story. The remaining high-value work is not more experimentation; it is canonicalization, title precision, and provenance cleanup.

---


## 2026-04-05 11:05 KST — Supervisor Review

### Progress
- The project is still in the appropriate final phase: manuscript hardening and submission-package maintenance, not further benchmark expansion. Recent artifact timestamps remain concentrated around 2026-04-01 manuscript assets (`latex/main.tex`, compiled PDFs, submission metadata, cover letter, highlights, graphical abstract, and synchronized figure copies), with the last supervisory update on 2026-04-03.
- A meaningful milestone remains intact: the work has a real submission-grade package, and the LaTeX manuscript appears internally synchronized around the API-isolated story rather than the older near-perfect subagent-evaluation narrative.
- The current manuscript-facing evidence is sharper than the early benchmark-era claim set. `latex/main.tex` now explicitly centers the setup-conditioned vs image-only rank reversal, uses the API-isolated 88.9% / 33.3% Claude numbers, and frames contamination as a methodological lesson rather than a side note.
- However, the repository still contains unresolved evidence-hygiene conflicts. `benchmark/labels/canonical_results_summary.json` continues to present the superseded 99.6% / 100% Claude story and is directly inconsistent with the current manuscript title-page narrative and results tables. `research_state.json` also remains historically inconsistent, still encoding an earlier project-phase interpretation.

### Strengths
- **The paper thesis is now scientifically sharper than the original leaderboard story.** The strongest contribution is the task/mechanism distinction: setup-conditioned CFD validation and image-only plausibility judgment are not the same capability, and evaluator ranking reverses between them.
- **Submission-facing manuscript assets look coherent.** The LaTeX manuscript appears aligned around the API-isolated result set, the contamination warning, and the narrower setup-conditioned validation claim.
- **The benchmark itself still looks paper-worthy.** Ten CFD scenarios, controlled error taxonomy, expert comparison, and dual-condition evaluation remain a credible benchmark/protocol contribution.
- **The project survived the integrity stress test.** A flattering contaminated result was replaced with a weaker but defensible one, and the paper still retained a publishable story. That improves trustworthiness rather than weakening it.

### Concerns/Risks
- **Major paper-risk: the repo still does not have a true canonical result source.** A file explicitly named `canonical_results_summary.json` still advertises the superseded near-perfect Claude narrative. This is now the clearest hidden paper risk because it can leak stale numbers back into future manuscript revisions, rebuttals, metadata forms, or figure captions.
- **State/provenance hygiene is still weak.** `research_state.json` still reports a phase history that no longer matches the actual completion state of the project. That is less dangerous than stale result numbers, but it still weakens auditability.
- **Title/framing risk remains nonzero.** The manuscript title is improved by adding “A Benchmark for Physical Plausibility Assessment of Flow Field Visualizations,” but `Can Vision-Language Models Evaluate Computational Fluid Dynamics?` still invites reviewers to read the scope more broadly than the evidence comfortably supports.
- **Hidden reviewer risk: denominator and protocol confusion can still reappear if stale files are reused.** The manuscript seems careful, but the repository still contains multiple result eras, which raises avoidable consistency risk during revision.

### Next Priorities
1. **Repair canonicalization immediately.** Replace or deprecate `benchmark/labels/canonical_results_summary.json` so the repo has one explicit submission-canonical result source matching `latex/main.tex`.
2. **Update archival state files.** Reconcile `research_state.json` with the actual project completion/evaluation status so the internal record matches the submission-stage reality.
3. **Do one more repo-wide stale-claim sweep.** Search for old 99.6% / 100% Claude figures, contaminated-evaluation wording, and pre-rerun summaries across `docs/`, `benchmark/`, and any submission-support files.
4. **Consider one last title-precision pass.** If feasible before submission, a slightly narrower title/subtitle would reduce reviewer attack surface and better match the evidence scope.
5. **Keep work focused on evidence hygiene, not new experiments.** The main remaining threat is inconsistency, not missing novelty.

### Paper Verdict
- **Overall verdict:** the project remains in the correct final consolidation phase and already has a meaningful manuscript-package milestone.
- **Supervisor judgment:** there is still a **major paper-risk worth surfacing**: canonical evidence hygiene is not fully repaired because the nominal canonical summary file still conflicts with the manuscript’s valid API-isolated story.
- **Reviewer-style bottom line:** this remains a credible submission-grade benchmark/protocol paper, and the main remaining danger is self-inflicted inconsistency rather than scientific weakness. I would not ask for a direction change, but I would treat canonical result cleanup as near-blocking before any further submission or revision activity.


## 2026-04-05 20:30 KST — Supervisor Review

### Progress
- The project remains in the correct final-stage phase: manuscript consolidation, submission packaging, and evidence hygiene rather than further benchmark growth. File activity is still concentrated around the 2026-04-01 submission bundle (`latex/main.tex`, PDFs, submission metadata, cover letter, highlights, and synchronized figures), with no sign that new experiments are currently the bottleneck.
- The manuscript-facing story still appears coherent. `latex/main.tex` is aligned around the API-isolated evaluation narrative, the setup-conditioned vs image-only rank reversal, and the contamination warning as a methodological contribution rather than an embarrassment to hide.
- The benchmark and paper package remain substantively strong enough for a submission-grade benchmark/protocol paper: 10 scenarios, controlled CFD error taxonomy, expert baseline, dual-condition evaluation, and a credible methodological lesson about evaluation isolation.
- However, the repository is still carrying multiple live result eras in easy-to-reuse locations. In addition to `benchmark/labels/canonical_results_summary.json` still reporting the superseded 99.6% / majority-vote-100% Claude story, several manuscript-adjacent draft files under `docs/` still contain outdated “100%” era phrasing. A newer file, `benchmark/labels/paper_canonical_apiisolated.json`, appears to hold the correct API-isolated narrative, but the repo still does not clearly elevate that as the single canonical source.

### Strengths
- **The main scientific contribution is still sharp and defensible.** The strongest result is the capability/mechanism distinction: setup-conditioned CFD validation and image-only plausibility judgment are different tasks, and evaluator rankings reverse between them.
- **The submission-facing manuscript seems internally consistent.** The LaTeX paper, cover-letter/support files, and figure set appear anchored to the API-isolated story rather than the older contaminated near-perfect result.
- **The benchmark contribution remains paper-worthy independent of any one model score.** The scenario coverage, error taxonomy, expert comparison, and protocol design support a benchmark/protocol contribution even if reviewer enthusiasm for specific model scores varies.
- **The project has already done the hard scientific-hygiene correction.** Detecting contamination and preserving the paper with a weaker but valid result is a major credibility positive.

### Concerns/Risks
- **Major paper-risk: canonicalization is still visibly unresolved.** The repo currently exposes both an outdated `canonical_results_summary.json` and a newer `paper_canonical_apiisolated.json`. That naming collision itself is a provenance hazard. If someone asks “what is the canonical result file?”, the wrong answer is still too easy.
- **Major paper-risk: stale manuscript-adjacent drafts are still reusable.** Files such as `docs/paper_draft_v1.md`, `docs/paper_draft_methodology.md`, `docs/paper_draft_intro_related.md`, `docs/paper_draft_results_discussion.md`, and `docs/26-03-26_paper_outline.md` still contain older 100%-era language. Even if intended as historical traces, they remain close enough to the paper workflow to leak stale claims back into revisions or response writing.
- **State/provenance hygiene remains weak.** `research_state.json` still reports an outdated phase history (`P4_full_eval` not started while the manuscript clearly reflects completed evaluation and submission-stage work). This is not the main scientific risk, but it is still archival debt.
- **Title/framing risk remains secondary but real.** `Can Vision-Language Models Evaluate Computational Fluid Dynamics?` is still broader than the evidence sweet spot. The actual strongest evidence is about setup-conditioned physical plausibility assessment on this benchmark.

### Next Priorities
1. **Make one file unambiguously canonical.** Either rename/promote `paper_canonical_apiisolated.json` as the sole canonical result file or explicitly deprecate `canonical_results_summary.json` with a banner and replacement pointer.
2. **Quarantine stale draft files.** Move outdated manuscript-adjacent Markdown drafts into an `archive/` area or add unmistakable `OBSOLETE / superseded by API-isolated LaTeX manuscript` banners.
3. **Repair `research_state.json`.** Update the project state so the phase history and evaluation status reflect the actual completed benchmark/evaluation/manuscript state.
4. **Do a repo-wide stale-claim audit one more time.** Search for 99.6%, “100%,” “majority vote 100%,” and pre-rerun language across `docs/`, `benchmark/`, and submission-support files.
5. **If there is time before submission, narrow the headline/title slightly.** This is less urgent than canonicalization, but it would reduce reviewer attack surface.

### Paper Verdict
- **Overall verdict:** the project is still in the correct late paper-hardening phase and remains a credible submission-grade benchmark/protocol paper.
- **Supervisor judgment:** there is still a **major paper-risk worth surfacing now**: unresolved canonicalization plus stale manuscript-adjacent drafts continue to expose the project to self-inflicted inconsistency.
- **Reviewer-style bottom line:** novelty and experimental design are no longer the main threat. The main threat is provenance hygiene. If the canonical result path is made explicit and stale drafts are quarantined, the paper should be on much safer ground for submission and later revision.

---

## 2026-04-06 20:30 KST — Supervisor Review

### Progress
- The project is still in the correct final-stage phase: manuscript hardening and submission-evidence maintenance, not further benchmark expansion. Recent timestamps remain concentrated on the 2026-04-01 LaTeX/submission bundle, and no new experimental assets suggest that further data growth is the current bottleneck.
- The manuscript-facing core remains coherent. `latex/main.tex`, `latex/main_blinded.tex`, submission metadata, and cover-letter/highlights assets all consistently use the API-isolated story: Claude 88.9% setup-conditioned vs 33.3% image-only, explicit rank reversal, and contamination as a methodological lesson.
- A meaningful paper milestone remains intact: the submission bundle exists, the paper argument is scientifically sharper than the older leaderboard-like framing, and the benchmark/protocol package is already strong enough for a submission-grade paper if evidence hygiene is cleaned.
- However, today’s repo sweep confirms that the same provenance problem is still live rather than merely historical. The file explicitly named `benchmark/labels/canonical_results_summary.json` still contains the superseded 99.6% / majority-vote-100% Claude narrative, while `benchmark/labels/paper_canonical_apiisolated.json` contains the paper-valid API-isolated numbers. Several Markdown drafts under `docs/` also still contain reusable outdated claims.

### Strengths
- **The scientific story is now sharper and more defensible than the original one.** The strongest contribution is not that one VLM scored highly, but that setup-conditioned CFD validation and image-only plausibility judgment are different tasks with different evaluator rankings.
- **Submission-facing assets are internally aligned.** The active LaTeX manuscript and submission-support documents consistently reflect the API-isolated evaluation rather than the contaminated near-perfect story.
- **The benchmark remains paper-worthy independent of one model’s score.** Ten scenarios, controlled error taxonomy, expert comparison, and dual-condition evaluation still support a credible benchmark/protocol contribution.
- **The project already passed the hard integrity test.** The attractive but compromised result was replaced with a weaker valid one, and the paper still retained a publishable message.

### Concerns/Risks
- **Major paper-risk: the repo still presents two competing “truths,” and the one named canonical is wrong for the paper.** This is now the clearest hidden reviewer/revision risk because stale numbers can easily leak back into rebuttal text, abstract edits, captions, metadata, or future revisions.
- **Major paper-risk: stale manuscript-adjacent drafts remain too reusable.** `docs/paper_draft_v1.md`, `docs/paper_draft_intro_related.md`, `docs/paper_draft_methodology.md`, `docs/paper_draft_results_discussion.md`, and `docs/26-03-26_paper_outline.md` still contain superseded 100%-era language and old contamination-era framing.
- **State/provenance hygiene is still unresolved.** `research_state.json` still encodes an earlier accelerated benchmark phase history and does not reflect the real completion state of evaluation/manuscript hardening.
- **Title precision remains a secondary risk.** The current title is stronger than the evidence sweet spot; the safest scope remains setup-conditioned physical plausibility assessment on this benchmark rather than broad CFD evaluation.

### Next Priorities
1. **Resolve the canonical file conflict immediately.** Either replace/deprecate `benchmark/labels/canonical_results_summary.json` and point to `paper_canonical_apiisolated.json`, or rename the API-isolated file so the repo has one unmistakable canonical source.
2. **Quarantine obsolete drafts.** Archive or banner stale Markdown paper-draft files so they cannot be mistaken for active manuscript sources.
3. **Repair `research_state.json`.** Update project phase/evaluation status so the archival record matches the actual submission-stage reality.
4. **Run one final stale-claim sweep before any revision/submission touch.** Specifically search 99.6%, 100%, majority-vote-100%, old GPT-5.4 80.0%, and pre-rerun contamination wording across `docs/`, `benchmark/`, and submission-support assets.
5. **If there is still time, tighten title scope slightly.** This remains lower priority than canonicalization, but it would reduce reviewer attack surface.

### Paper Verdict
- **Overall verdict:** the project remains in the right final consolidation phase and still looks submission-grade from a benchmark/protocol perspective.
- **Supervisor judgment:** there is still a **major paper-risk that should be surfaced**: canonical evidence hygiene is not repaired yet, and the repo still exposes stale manuscript claims in easy-to-reuse locations.
- **Reviewer-style bottom line:** novelty and experimental design are already good enough; the remaining threat is self-inflicted inconsistency. If canonicalization and stale-draft isolation are cleaned, the paper should be on much firmer ground for submission and revision.
