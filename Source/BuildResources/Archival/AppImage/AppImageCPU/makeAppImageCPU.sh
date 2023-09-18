#!/bin/bash -l

#Target AppDir
APP=AppDir

#prepare the AppDir structure
rm -rf appimage-build
rm -rf $APP

cd ..
cd ..
cd ..
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=AppDir/usr ..
make install
mkdir -p $APP/usr/share/icons
cp ../BuildResources/AppImage/ico512.png $APP/usr/share/icons/LightwaveExplorer.png
cp ../BuildResources/AppImage/AppImageCPU/AppImageBuilder2.yml AppImageBuilder.yml
appimage-builder --skip-tests






