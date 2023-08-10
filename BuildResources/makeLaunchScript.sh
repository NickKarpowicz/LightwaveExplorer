rm -f LightwaveExplorerLauncher.sh
touch LightwaveExplorerLauncher.sh
echo "#!/bin/bash -l" >>  LightwaveExplorerLauncher.sh
echo ". ${ONEAPI_ROOT}/setvars.sh" >> LightwaveExplorerLauncher.sh
echo "export SYCL_CACHE_PERSISTENT=1" >>  LightwaveExplorerLauncher.sh
echo "exec /usr/local/bin/LightwaveExplorer" >>  LightwaveExplorerLauncher.sh