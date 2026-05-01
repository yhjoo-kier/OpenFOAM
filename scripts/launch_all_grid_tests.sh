#!/bin/bash
# Launch all 24 grid independence runs with controlled parallelism
# 3 cases x 4 mesh sizes x 2 types (ref + pred) = 24 runs
# Batch 1 (coarse, fast): 0.35 + 0.25  -> 12 runs, ~6 parallel
# Batch 2 (medium):       0.18         ->  6 runs, ~6 parallel
# Batch 3 (fine, slow):   0.10         ->  6 runs, ~6 parallel
cd /home/yhjoo/projects/OpenFOAM

LOGDIR="benchmark/manifests/grid_test_logs"
mkdir -p "$LOGDIR"
echo "=== Grid test started at $(date) ===" > "$LOGDIR/master.log"

run_cfd() {
    local name=$1
    local scene=$2
    local mesh=$3
    local log="$LOGDIR/${name}.log"
    local start_ts=$(date +%s)

    echo "[$(date +%H:%M:%S)] START $name mesh=$mesh" >> "$LOGDIR/master.log"
    python3 scripts/run_indoor_stabilized.py \
        --scenario "$scene" \
        --name "$name" \
        --mesh-size "$mesh" \
        --skip-mesh-ladder \
        --solver-timeout 3600 \
        --disable-repair \
        > "$log" 2>&1
    local rc=$?
    local end_ts=$(date +%s)
    local elapsed=$((end_ts - start_ts))
    echo "[$(date +%H:%M:%S)] DONE  $name rc=$rc elapsed=${elapsed}s" >> "$LOGDIR/master.log"
    echo "$rc" > "$LOGDIR/${name}.rc"
    echo "$elapsed" > "$LOGDIR/${name}.elapsed"
}

CASES=("a1_01" "a3_03" "a4_03")

launch_batch() {
    local mesh=$1
    local ml=$2
    local PIDS=()
    local NAMES=()

    for case_id in "${CASES[@]}"; do
        # Reference
        local ref_name="grid_test_${case_id}_ref_${ml}"
        local ref_scene="benchmark/scenes/${case_id}.json"
        run_cfd "$ref_name" "$ref_scene" "$mesh" &
        PIDS+=($!)
        NAMES+=("$ref_name")

        # Predicted
        local pred_name="grid_test_${case_id}_pred_${ml}"
        local pred_scene="benchmark/evaluations_posthoc_scaled_longest_span/bench_${case_id}/floorplan/predicted_scene.json"
        run_cfd "$pred_name" "$pred_scene" "$mesh" &
        PIDS+=($!)
        NAMES+=("$pred_name")
    done

    echo "  Batch mesh=$mesh: waiting for ${#PIDS[@]} jobs..."
    for i in "${!PIDS[@]}"; do
        wait ${PIDS[$i]}
        rc=$?
        if [ $rc -ne 0 ]; then
            echo "  FAILED: ${NAMES[$i]} (rc=$rc)"
        else
            echo "  OK: ${NAMES[$i]}"
        fi
    done
}

echo "=== Batch 1: mesh=0.35 (6 runs) ==="
launch_batch 0.35 035

echo "=== Batch 2: mesh=0.25 (6 runs) ==="
launch_batch 0.25 025

echo "=== Batch 3: mesh=0.18 (6 runs) ==="
launch_batch 0.18 018

echo "=== Batch 4: mesh=0.10 (6 runs) ==="
launch_batch 0.10 010

echo ""
echo "=== ALL BATCHES COMPLETE at $(date) ==="
cat "$LOGDIR/master.log"
