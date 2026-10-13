#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
mkdir -p logs
setsid nohup bash "${SCRIPT_DIR}/comp_top10_worker.sh" 0 pos comp01 comp02 comp03 comp04 comp05 comp06 comp07 comp08 comp09 comp10 > logs/workerA.log 2>&1 < /dev/null &
setsid nohup bash "${SCRIPT_DIR}/comp_top10_worker.sh" 1 neg comp01 comp04 comp07 comp10 known_peptide > logs/workerB.log 2>&1 < /dev/null &
setsid nohup bash "${SCRIPT_DIR}/comp_top10_worker.sh" 2 neg comp02 comp05 comp08 > logs/workerC.log 2>&1 < /dev/null &
echo launched
