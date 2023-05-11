#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

pushd frontend
trunk build --release
popd

cargo run --bin pgr-server --release --  --addr 0.0.0.0 --port 3000 --static-dir ./dist --data-path-prefix /wd/pgr-tk-demo-data/data/pgr-tk-HGRP-y1-evaluation-set-v0
