brew install cmake make pkgconfig qt cairo wget
mkdir LightwaveExplorerBuild
cd LightwaveExplorerBuild
git clone --depth 1 https://github.com/NickKarpowicz/LightwaveExplorer
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

#complete the directory structure of the .app
mkdir $APP/Contents/Resources/
mkdir $APP/Contents/Resources/lib
mkdir $APP/Contents/Resources/bin
mkdir $APP/Contents/Resources/etc
mkdir $APP/Contents/Resources/share

#copy in the databases and icons
cp CrystalDatabase.txt $APP/Contents/Resources
cp Source/BuildResources/DefaultValues.ini $APP/Contents/Resources
cp Source/BuildResources/Licenses.txt $APP/Contents/Resources
cp Source/BuildResources/AppIcon.icns $APP/Contents/Resources

cp -r $APP /Applications/
