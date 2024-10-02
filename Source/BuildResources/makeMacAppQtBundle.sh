#!/bin/bash -l
if [ ! -e "$1" ]; then
    echo "please give the path to the Qt bin folder, e.g. /Users/nick/Qt/6.7.2/macos/bin"
    exit 1
fi
QTPATH=$1
BIN=LightwaveExplorer
APP=build/${BIN}.app
BINPATH=${APP}/Contents/MacOS/${BIN}

#Shared libraries location
LIBS="/Users"

#build executable
rm -rf build
mkdir build
cd build
${QTPATH}/qt-cmake .. -G Ninja
cmake --build . --config Release
cd ..
mkdir ${APP}/Contents/Resources
cp CrystalDatabase.txt $APP/Contents/Resources
cp Source/BuildResources/DefaultValues.ini $APP/Contents/Resources
cp Source/BuildResources/Licenses.txt $APP/Contents/Resources
cp Source/BuildResources/AppIcon.icns $APP/Contents/Resources
${QTPATH}/macdeployqt ${APP} -dmg
mv build/${BIN}.dmg build/${BIN}MacOS.dmg
