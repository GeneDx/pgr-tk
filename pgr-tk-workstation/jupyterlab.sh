#!/bin/bash
mkdir -p /wd/results/
mkdir -p /wd/code/
mkdir -p /wd/data//
ln -sf /wd/* /
. /opt/conda/etc/profile.d/conda.sh
jupyter-lab --ip="*" --allow-root --no-browser --port 8899 --NotebookApp.disable_check_xsrf=True /wd
