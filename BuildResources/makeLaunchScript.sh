rm -f LightwaveExplorerLauncher.sh
touch LightwaveExplorerLauncher.sh
echo "#!/bin/bash -l" >>  LightwaveExplorerLauncher.sh
echo ". ${ONEAPI_ROOT}/setvars.sh" >> LightwaveExplorerLauncher.sh
echo "export SYCL_CACHE_PERSISTENT=1" >>  LightwaveExplorerLauncher.sh
echo '"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/app/opt/cuda/lib64' >> LightwaveExplorerLauncher.sh
echo 'script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"' >> LightwaveExplorerLauncher.sh
echo 'exec "${script_dir}/LightwaveExplorer"' >>  LightwaveExplorerLauncher.sh