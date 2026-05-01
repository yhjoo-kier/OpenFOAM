#!/bin/bash
# Grid independence test runner
# Runs reference + predicted CFD for 3 cases x 4 mesh sizes = 24 total runs

set -e
cd /home/yhjoo/projects/OpenFOAM

CASES=("a1_01" "a3_03" "a4_03")
MESHES=("0.35" "0.25" "0.18" "0.10")
MESH_LABELS=("035" "025" "018" "010")

LOGDIR="/home/yhjoo/projects/OpenFOAM/benchmark/manifests/grid_test_logs"
mkdir -p "$LOGDIR"

RUN_SCRIPT="python3 scripts/run_indoor_stabilized.py"

run_one() {
    local case_id=$1
    local mesh_size=$2
    local mesh_label=$3
    local run_type=$4  # ref or pred
    local scene_json=$5
    local run_name="grid_test_${case_id}_${run_type}_${mesh_label}"
    local logfile="$LOGDIR/${run_name}.log"

    echo "[$(date +%H:%M:%S)] START $run_name (mesh=$mesh_size)"
    $RUN_SCRIPT \
        --scenario "$scene_json" \
        --name "$run_name" \
        --mesh-size "$mesh_size" \
        --skip-mesh-ladder \
        --solver-timeout 3600 \
        --disable-repair \
        > "$logfile" 2>&1
    local rc=$?
    echo "[$(date +%H:%M:%S)] DONE  $run_name rc=$rc"
    return $rc
}

# Launch specified run type
CASE=$1
MESH_IDX=$2
TYPE=$3  # ref or pred

case_id="${CASES[$CASE]}"
mesh_size="${MESHES[$MESH_IDX]}"
mesh_label="${MESH_LABELS[$MESH_IDX]}"

if [ "$TYPE" = "ref" ]; then
    scene_json="benchmark/scenes/${case_id}.json"
else
    scene_json="benchmark/evaluations_posthoc_scaled_longest_span/bench_${case_id}/floorplan/predicted_scene.json"
fi

run_one "$case_id" "$mesh_size" "$mesh_label" "$TYPE" "$scene_json"
