brew install cmake make pkgconfig qt cairo wget
mkdir LightwaveExplorerBuild
cd LightwaveExplorerBuild
if [ -n "$1" ]; then
    git clone --depth 1 --branch "$1" https://github.com/NickKarpowicz/LightwaveExplorer
else
    git clone --depth 1 https://github.com/NickKarpowicz/LightwaveExplorer
fi
cd LightwaveExplorer

BIN=LightwaveExplorer
APP=build/${BIN}.app
BINPATH=${APP}/Contents/MacOS/${BIN}

#build executable
mkdir build
cd build
cmake -DCMAKE_OSX_DEPLOYMENT_TARGET=$(sw_vers -productVersion) -DCMAKE_CXX_FLAGS="-O3 -march=native" ..
make
cd ..

#copy in the databases and icons
mkdir $APP/Contents/Resources/
cp CrystalDatabase.txt $APP/Contents/Resources
cp Source/BuildResources/DefaultValues.ini $APP/Contents/Resources
cp Source/BuildResources/Licenses.txt $APP/Contents/Resources
cp Source/BuildResources/AppIcon.icns $APP/Contents/Resources

cp -r $APP /Applications/
cd ../..
rm -rf LightwaveExplorerBuild
