#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
mkdir -p logs
setsid nohup bash "${SCRIPT_DIR}/comp_top10_worker.sh" 3 neg comp03 comp06 comp09 > logs/workerD.log 2>&1 < /dev/null &
echo launched
