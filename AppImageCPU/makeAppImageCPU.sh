#!/bin/bash -l

#Target AppDir
APP=AppDir

#Target executable
BIN=../LightwaveExplorer

cd ..
make cpuonly
cd AppImageCPU

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

wget https://github.com/AppImage/AppImageKit/releases/download/continuous/AppRun-x86_64
cp AppRun-x86_64 $APP/AppRun
rm AppRun-x86_64
cp AppRun.env $APP/AppRun.env

NEWLPATH=":/usr/lib:"

#optional launcher file, only used with appimage-builder --generate it so that the environment is set
#probably not needed anymore
touch $APP/usr/bin/LaunchLWE.sh
echo "#!/bin/bash -l" >>  $APP/usr/bin/LaunchLWE.sh
echo 'HERE="$(dirname "$(readlink -f "${0}")")"' >> $APP/usr/bin/LaunchLWE.sh
echo 'APPBASE="$HERE/../.."' >> $APP/usr/bin/LaunchLWE.sh
echo "export LD_LIBRARY_PATH=$NEWPATH"'$LD_LIBRARY_PATH' >> $APP/usr/bin/LaunchLWE.sh
echo 'echo "lp is : $PAFF"' >> $APP/usr/bin/LaunchLWE.sh
echo "exec LightwaveExplorer" >>  $APP/usr/bin/LaunchLWE.sh
chmod +x $APP/usr/bin/LaunchLWE.sh

#buid the .desktop file
touch $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "[Desktop Entry]" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

#note: use this instead to trick appimage-builder --generate into doing something useful
echo "Exec=LaunchLWE.sh" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

echo "Icon=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Type=Application" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Categories=Utility" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

#build the image
appimage-builder --skip-tests

#clean
cd .. 
make clean
cd AppImage
rm -rf AppDir
rm -rf appimage-build






