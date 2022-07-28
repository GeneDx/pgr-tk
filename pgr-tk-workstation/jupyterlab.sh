#!/bin/bash
ln -s /wd/* /
. /opt/conda/etc/profile.d/conda.sh
jupyter-lab --ip="*" --allow-root --no-browser --port 8899 --NotebookApp.disable_check_xsrf=True /wd
