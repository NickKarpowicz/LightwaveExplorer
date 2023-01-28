#!/bin/bash -l

#Target AppDir
APP=AppDir

#prepare the AppDir structure
rm -rf appimage-build
rm -rf $APP

cd ..
cd ..
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=AppDir/usr -DONEAPI_ROOT=${ONEAPI_ROOT} -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.0/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 ..
make install
mkdir -p $APP/usr/share/icons

cp ../BuildResources/AppImage/ico512.png $APP/usr/share/icons/LightwaveExplorer.png
cp ../BuildResources/AppImage/AppImageBuilder2.yml AppImageBuilder.yml

#optional launcher file, only used with appimage-builder --generate it so that the environment is set
#probably not needed anymore
# touch $APP/usr/bin/LaunchLWE.sh
# echo ". "
# echo "#!/bin/bash -l" >>  $APP/usr/bin/LaunchLWE.sh
# echo "export OCL_ICD_FILENAMES=libintelocl_emu.so:libalteracl.so:libintelocl.so" >>  $APP/usr/bin/LaunchLWE.sh
# echo 'HERE="$(dirname "$(readlink -f "${0}")")"' >> $APP/usr/bin/LaunchLWE.sh
# echo 'APPBASE="$HERE/../.."' >> $APP/usr/bin/LaunchLWE.sh
# echo "export LD_LIBRARY_PATH=$NEWPATH"'/usr/lib:/home/nick/intel/oneapi/compiler/latest/linux/lib:/home/nick/intel/oneapi/compiler/latest/linux/compiler/lib/intel64:/home/nick/intel/oneapi/tbb/latest/lib/intel64/gcc4.8:/home/nick/intel/oneapi/compiler/latest/linux/lib/x64:$LD_LIBRARY_PATH' >> $APP/usr/bin/LaunchLWE.sh
# echo 'echo "lp is : $PAFF"' >> $APP/usr/bin/LaunchLWE.sh
# echo "exec ./LightwaveExplorer" >>  $APP/usr/bin/LaunchLWE.sh
# chmod +x $APP/usr/bin/LaunchLWE.sh

# #buid the .desktop file
# touch $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "[Desktop Entry]" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "Name=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# #note: use this instead to trick appimage-builder --generate into doing something useful
# #echo "Exec=LaunchLWE.sh" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "Exec=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "Icon=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "Type=Application" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
# echo "Categories=Utility" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

# #build the image
appimage-builder --skip-tests






