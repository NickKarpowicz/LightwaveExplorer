#!/bin/bash -l
cd ..
make
cd AppImage
rm -rf LightwaveExplorer.AppDir
mkdir -p LightwaveExplorer.AppDir
mkdir -p LightwaveExplorer.AppDir/usr
mkdir -p LightwaveExplorer.AppDir/usr/lib
mkdir -p LightwaveExplorer.AppDir/usr/bin
mkdir -p LightwaveExplorer.AppDir/usr/share
mkdir -p LightwaveExplorer.AppDir/usr/share/glib-2.0
mkdir -p LightwaveExplorer.AppDir/usr/share/LightwaveExplorer
cp /usr/local/cuda-12.0/targets/x86_64-linux/lib/libcufft.so.11 LightwaveExplorer.AppDir/usr/lib/libcufft.so.11
cp /usr/lib/x86_64-linux-gnu/libnvidia-ml.so.1 LightwaveExplorer.AppDir/usr/lib/libnvidia-ml.so.1
cp /usr/local/cuda-12.0/targets/x86_64-linux/lib/libcudart.so.12 LightwaveExplorer.AppDir/usr/lib/libcudart.so.12
cp /usr/lib/x86_64-linux-gnu/libgomp.so.1 LightwaveExplorer.AppDir/usr/lib/libgomp.so.1
cp /usr/lib/x86_64-linux-gnu/libgtk-4.so.1 LightwaveExplorer.AppDir/usr/lib/libgtk-4.so.1
cp /usr/lib/x86_64-linux-gnu/libpangocairo-1.0.so.0 LightwaveExplorer.AppDir/usr/lib/libpangocairo-1.0.so.0
cp /usr/lib/x86_64-linux-gnu/libpango-1.0.so.0 LightwaveExplorer.AppDir/usr/lib/libpango-1.0.so.0
cp /usr/lib/x86_64-linux-gnu/libharfbuzz.so.0 LightwaveExplorer.AppDir/usr/lib/libharfbuzz.so.0
cp /usr/lib/x86_64-linux-gnu/libgdk_pixbuf-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libgdk_pixbuf-2.0.so.0
cp /usr/lib/x86_64-linux-gnu/libcairo-gobject.so.2 LightwaveExplorer.AppDir/usr/lib/libcairo-gobject.so.2
cp /usr/lib/x86_64-linux-gnu/libcairo.so.2 LightwaveExplorer.AppDir/usr/lib/libcairo.so.2
cp /usr/lib/x86_64-linux-gnu/libgobject-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libgobject-2.0.so.0
cp /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libglib-2.0.so.0
cp /home/nick/intel/oneapi/tbb/2021.8.0/lib/intel64/gcc4.8/*.so* LightwaveExplorer.AppDir/usr/lib/
cp -r /home/nick/intel/oneapi/compiler/2023.0.0/linux/lib/* LightwaveExplorer.AppDir/usr/lib/
cp -r /home/nick/intel/oneapi/compiler/2023.0.0/linux/lib/x64/* LightwaveExplorer.AppDir/usr/lib/
rm -rf LightwaveExplorer.AppDir/usr/lib/x64
rm -rf LightwaveExplorer.AppDir/usr/lib/clc
rm -rf LightwaveExplorer.AppDir/usr/lib/oclfpga
rm -rf LightwaveExplorer.AppDir/usr/lib/clang


rm LightwaveExplorer.AppDir/usr/lib/libOclCpuBackEnd_emu.so.2022.15.12.0
rm LightwaveExplorer.AppDir/usr/lib/icx-lto.so
cp -r /home/nick/intel/oneapi/compiler/2023.0.0/linux/compiler/lib/intel64_lin/*.so* LightwaveExplorer.AppDir/usr/lib/
cp /usr/lib/x86_64-linux-gnu/libm.so.6 LightwaveExplorer.AppDir/usr/lib/libm.so.6
cp /usr/lib/x86_64-linux-gnu/libstdc++.so.6 LightwaveExplorer.AppDir/usr/lib/libstdc++.so.6
cp /usr/lib/x86_64-linux-gnu/libc.so.6 LightwaveExplorer.AppDir/usr/lib/libc.so.6
cp /usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2 LightwaveExplorer.AppDir/usr/lib/ld-linux-x86-64.so.2
cp -r /usr/share/glib-2.0 LightwaveExplorer.AppDir/usr/share
cp ../LightwaveExplorer LightwaveExplorer.AppDir/usr/bin/LightwaveExplorer
cp ../CrystalDatabase.txt LightwaveExplorer.AppDir/usr/bin/CrystalDatabase.txt
cp ../DefaultValues.ini LightwaveExplorer.AppDir/usr/bin/DefaultValues.ini
cp ico512.png LightwaveExplorer.AppDir/LightwaveExplorer.png
cp AppRun LightwaveExplorer.AppDir/AppRun

touch LightwaveExplorer.AppDir/usr/bin/LaunchLWE.sh
echo "#!/bin/bash -l" >>  LightwaveExplorer.AppDir/usr/bin/LaunchLWE.sh
echo "export OCL_ICD_FILENAMES=libintelocl_emu.so:libalteracl.so:libintelocl.so" >>  LightwaveExplorer.AppDir/usr/bin/LaunchLWE.sh
echo "LightwaveExplorer" >>  LightwaveExplorer.AppDir/usr/bin/LaunchLWE.sh
chmod +x LightwaveExplorer.AppDir/usr/bin/LaunchLWE.sh

touch LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "[Desktop Entry]" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "Exec=LaunchLWE.sh" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "Icon=LightwaveExplorer" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "Type=Application" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop
echo "Categories=Utility" >> LightwaveExplorer.AppDir/LightwaveExplorer.desktop


ARCH=x86_64 ./appimagetool-x86_64.AppImage LightwaveExplorer.AppDir

rm -rf LightwaveExplorer.AppDir/
cd ..
make clean