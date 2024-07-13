#!/bin/bash
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if /sbin/ldconfig -p | grep -q libnvidia-ml; then
    exec "${script_dir}/LightwaveExplorer"
else
    . /app/opt/intel/oneapi/setvars.sh
    exec "${script_dir}/LightwaveExplorerNoCuda"
fi