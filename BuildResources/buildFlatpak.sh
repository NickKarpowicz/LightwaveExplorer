rm -rf build
mkdir build
cd build
. /opt/intel/oneapi/setvars.sh
cmake -DONEAPI_ROOT=${ONEAPI_ROOT} -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 ..
make
mkdir -p appdir/appdir/lib
cmake --install . --prefix=appdir/appdir
cp ../BuildResources/io.github.NickKarpowicz.LightwaveExplorer.yml ./

touch appdir/appdir/bin/LWELauncher.sh
echo "#!/bin/bash" >> appdir/appdir/bin/LWELauncher.sh
echo "env OCL_ICD_FILENAMES=libintelocl_emu.so:libalteracl.so:libintelocl.so LightwaveExplorer" >> appdir/appdir/bin/LWELauncher.sh
chmod +x appdir/appdir/bin/LWELauncher.sh

cp -a ${ONEAPI_ROOT}/compiler/latest/linux/compiler/lib/intel64/* appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/x64/* appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/libpi_level_zero.so appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/libpi_opencl.so appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/sycl.conf appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/libOpenCL.so* appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/compiler/latest/linux/lib/libsycl.so* appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/tbb/latest/lib/intel64/gcc4.8/libtbb.so* appdir/appdir/lib/
cp -a ${ONEAPI_ROOT}/tbb/latest/lib/intel64/gcc4.8/libtbbmalloc.so* appdir/appdir/lib/
cp -a /usr/lib/libcuda.so* appdir/appdir/lib/
cp -a /usr/lib/libnvidia-ml.so* appdir/appdir/lib/
cp -a /opt/cuda/targets/x86_64-linux/lib/libcufft* appdir/appdir/lib/

flatpak-builder build-dir io.github.NickKarpowicz.LightwaveExplorer.yml --force-clean
flatpak-builder --user --install --force-clean build-dir  io.github.NickKarpowicz.LightwaveExplorer.yml