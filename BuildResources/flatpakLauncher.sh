#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/opt/cuda/lib64
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NVML="libnvidia-ml.so.1"
OLD_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"


for dir in $(echo $OLD_LD_LIBRARY_PATH | tr ":" "\n"); do
    if [ -f "$dir/$LIBRARY" ]; then
        echo "Sorry, this was just a test."
        exit 0
    fi
done
exec "${script_dir}/LightwaveExplorer"