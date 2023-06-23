#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

pushd frontend
trunk build --release
popd

cargo run --bin pgr-server --no-default-features --release --  --addr 0.0.0.0 --port 3000 --static-dir ./dist --frg-file --data-path-prefix $HOME/Sandbox/pgr-tk-data/HGRP-y1-evaluation-set_fragdb 
