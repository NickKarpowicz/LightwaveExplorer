#!/bin/bash -l

#Target AppDir
APP=AppDir

#Target executable
BIN=../../../build/LightwaveExplorer

cd ..
cd ..
cd ..
./BuildResources/buildLinuxCPU.sh
cd BuildResources/AppImage/AppImageCPU

#prepare the AppDir structure
rm -rf appimage-build
rm -rf $APP
mkdir -p $APP
mkdir -p $APP/usr
mkdir -p $APP/usr/lib
mkdir -p $APP/usr/bin
mkdir -p $APP/usr/share
mkdir -p $APP/usr/share/glib-2.0
mkdir -p $APP/usr/share/LightwaveExplorer
mkdir -p $APP/usr/share/icons
cp $BIN $APP/usr/bin
cp ../CrystalDatabase.txt $APP/usr/bin
cp ../DefaultValues.ini $APP/usr/bin
cp ../CrystalDatabase.txt $APP/
cp ../DefaultValues.ini $APP/
cp ico512.png $APP/usr/share/icons/LightwaveExplorer.png
cp -r /usr/share/glib-2.0 $APP/usr/share/
wget https://github.com/AppImage/AppImageKit/releases/download/continuous/AppRun-x86_64
cp AppRun-x86_64 $APP/AppRun
rm AppRun-x86_64
cp AppRun.env $APP/AppRun.env


#buid the .desktop file
touch $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "[Desktop Entry]" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Icon=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Type=Application" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Categories=Utility" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

#build the image
appimage-builder --skip-tests

#clean
cd .. 
cd ..
cd ..
cp GPLsource.tar BuildResources/AppImage/AppImageCPU/
rm -rf build
cd BuildResources/AppImage/AppImageCPU
rm -rf AppDir
rm -rf appimage-build






