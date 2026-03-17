#!/usr/bin/env bash
set -euo pipefail

# Wrapper: run a long command and proactively notify via OpenClaw Telegram when done.
#
# Usage:
#   ./scripts/run_with_notify.sh \
#     --target 8292931055 \
#     --label "OpenFOAM steady" \
#     --log /home/yhjoo/projects/OpenFOAM/logs/steady_$(date +%F_%H%M).log \
#     -- bash -lc 'export PATH=/app/.venv/bin:$PATH; bash ./run_water_case_steady.sh'

OPENCLAW_BIN_DEFAULT="$HOME/.local/share/fnm/node-versions/v22.16.0/installation/bin/openclaw"
OPENCLAW_BIN="${OPENCLAW_BIN:-$OPENCLAW_BIN_DEFAULT}"
CHANNEL="telegram"
TARGET=""
LABEL="Long run"
LOG_PATH=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target) TARGET="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --log) LOG_PATH="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Unknown arg: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$TARGET" ]]; then
  echo "--target is required" >&2
  exit 2
fi

if [[ $# -eq 0 ]]; then
  echo "No command provided after --" >&2
  exit 2
fi

mkdir -p "$(dirname "${LOG_PATH:-/tmp/dummy}")" 2>/dev/null || true
START_TS=$(date '+%F %T')
START_EPOCH=$(date +%s)
HOST=$(hostname)

notify() {
  local text="$1"
  if [[ -x "$OPENCLAW_BIN" ]]; then
    "$OPENCLAW_BIN" message send --channel "$CHANNEL" --target "$TARGET" --message "$text" >/dev/null 2>&1 || true
  fi
}

notify "🚀 [$LABEL] 시작됨 (${START_TS}, host=${HOST})"

set +e
if [[ -n "$LOG_PATH" ]]; then
  "$@" 2>&1 | tee "$LOG_PATH"
  RUN_EXIT=${PIPESTATUS[0]}
else
  "$@"
  RUN_EXIT=$?
fi
set -e

END_TS=$(date '+%F %T')
END_EPOCH=$(date +%s)
ELAPSED=$((END_EPOCH-START_EPOCH))

if [[ $RUN_EXIT -eq 0 ]]; then
  MSG="✅ [$LABEL] 완료 (${END_TS}) | elapsed=${ELAPSED}s"
  [[ -n "$LOG_PATH" ]] && MSG+=" | log=${LOG_PATH}"
  notify "$MSG"
  exit 0
else
  MSG="❌ [$LABEL] 실패 (${END_TS}) | exit=${RUN_EXIT} | elapsed=${ELAPSED}s"
  [[ -n "$LOG_PATH" ]] && MSG+=" | log=${LOG_PATH}"
  notify "$MSG"
  exit $RUN_EXIT
fi
