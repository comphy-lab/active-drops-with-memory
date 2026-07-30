#!/bin/sh
set -eu

qcc -O2 -Wall -disable-dimensions dropMove.c -o dropMove -lm
./dropMove --params "${1:-parameters/single-drop.cfg}"
