# Gemini strategy

## Current position

This project already uses Gemini successfully in a CLI-driven prototype loop.
That was the right choice for rapid validation.

## Recommended next step

Migrate the **pipeline-default path** toward Gemini API integration while keeping CLI support for ad-hoc experiments.

## Why API fits this framework better

Use API for:
- image + text multimodal input
- strict JSON output requirements
- structured retries after invalid scene generation
- provenance logging of prompt, model, and response
- batch or background execution
- tighter integration with project scripts and result summaries

Keep CLI for:
- quick prompt experiments
- manual debugging
- one-off comparisons across Gemini models

## Suggested interface direction

A future scene-generation module should support both modes:
- `--backend cli`
- `--backend api`

Potential later default:
- API for normal project runs
- CLI as explicit fallback/debug mode

## Migration principle

Do not rewrite the whole framework at once.
Start by replacing only the Gemini scene-generation stage with an API-backed implementation that preserves the downstream JSON contract.
