#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"
if [[ -f .project_config ]]; then
  source .project_config
fi
if [[ "$(command -v qcc)" != "$PWD/basilisk/src/qcc" || ! -f basilisk/.comphy-lock ]]; then
  echo "A pinned project-local Basilisk compiler is required." >&2
  exit 1
fi

qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove --params "${1:-parameters/single-drop.cfg}"
