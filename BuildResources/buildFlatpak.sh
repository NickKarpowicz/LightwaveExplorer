rm -rf build
mkdir build
cd build
. ~/intel/oneapi/setvars.sh
cmake -DONEAPI_ROOT=${ONEAPI_ROOT} -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 ..
make
mkdir -p appdir/appdir/lib
cmake --install . --prefix=appdir/appdir
cp ../BuildResources/io.github.NickKarpowicz.LightwaveExplorer.yml ./

cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libintelocl.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/libpi_level_zero.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/libpi_opencl.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/sycl.conf appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/../clbltfnshared.rtl appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/__ocl_svml_l9.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/cl.cfg appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/cl.fpga_emu.cfg appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/clbltfnl9.rtl appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/cllibrary.rtl appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/cllibraryl9.o appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libOclCpuBackEnd.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libOclCpuBackEnd.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libclang_compiler.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libclang_compiler.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcpu_device.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcpu_device.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcpu_device_emu.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcpu_device_emu.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libtask_executor.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/libimf.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/libintlc.so*
cp -a /home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/libiomp5.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/libirng.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64/libsvml.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/libOpenCL.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/libsycl.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcl_logger.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcl_logger_emu.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libcommon_clang.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libintelocl.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libintelocl_emu.so appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libtask_executor.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/compiler/latest/linux/lib/x64/libtask_executor_emu.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/tbb/latest/lib/intel64/gcc4.8/libtbb.so* appdir/appdir/lib/
cp -a /home/nick/intel/oneapi/tbb/latest/lib/intel64/gcc4.8/libtbbmalloc.so* appdir/appdir/lib/
cp -a /usr/lib/libcuda.so* appdir/appdir/lib/
cp -a /usr//lib/libnvidia-ml.so* appdir/appdir/lib/
cp -a /opt/cuda/targets/x86_64-linux/lib/libcufft* appdir/appdir/lib/

flatpak-builder build-dir io.github.NickKarpowicz.LightwaveExplorer.yml --force-clean
flatpak-builder --user --install --force-clean build-dir  io.github.NickKarpowicz.LightwaveExplorer.yml