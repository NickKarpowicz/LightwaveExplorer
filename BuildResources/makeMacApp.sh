#!/bin/bash -l
BIN=LightwaveExplorer
APP=build/${BIN}.app
BINPATH=${APP}/Contents/MacOS/${BIN}

#install the required packages from homebrew
brew install --quiet cmake make llvm fftw gtk4 fmt pkgconfig libomp

#Homebrew libraries location
LIBS="$(brew --prefix)"

#find-and-replace to use fftw, replace std::format .etc with fmt::format since clang doesn't have it yet
sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LightwaveExplorerUtilities.h
sed -i'.bak' 's/fftw3_mkl.h/fftw3.h/g' LWEActiveDeviceCPU.h 
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
sed -i'.bak' 's!<format>!<fmt/format.h>!g ; s/std::format/fmt::format/g ; s/std::vformat/fmt::vformat/g ; s/std::make_format_args/fmt::make_format_args/g' LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp

#use .mm extension for the interface to use objective-c++; allows direct calls to Cocoa for file dialogs
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm

#build executable
cp BuildResources/CMakeListsMac.txt CMakeLists.txt
rm -rf build
mkdir build
cd build
cmake ..
make
cd ..

#restore the original source and clean up
cp BuildResources/AppImage/AppImageCPU/COPYING COPYING
tar cf GPLsource.tar COPYING *.cpp *.cu *.h LightwaveExplorerGTK/* DlibLibraryComponents/* BuildResources/macplistbase.plist BuildResources/CMakeListsMac.txt BuildResources/AppIcon.icns
rm COPYING
rm LightwaveExplorerUtilities.h
rm LWEActiveDeviceCPU.h
rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h
rm LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.mm
mv LightwaveExplorerUtilities.h.bak LightwaveExplorerUtilities.h
mv LWEActiveDeviceCPU.h.bak LWEActiveDeviceCPU.h
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.cpp
mv LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h.bak LightwaveExplorerGTK/LightwaveExplorerFrontendGTK.h

#complete the directory structure of the .app
mkdir $APP/Contents/Resources/
mkdir $APP/Contents/Resources/lib
mkdir $APP/Contents/Resources/bin
mkdir $APP/Contents/Resources/etc
mkdir $APP/Contents/Resources/share

#copy in the databases and icons
cp CrystalDatabase.txt $APP/Contents/Resources
cp DefaultValues.ini $APP/Contents/Resources
cp BuildResources/AppIcon.icns $APP/Contents/Resources

#Functions:

#for a given binary, copy its dynamic link dependencies to the $APP/Contents/Resources/lib folder
copySharedLibraries(){
    NLIBS=$(otool -L $1 | grep "$LIBS" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(otool -L $1 | grep "$LIBS" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        cp -n $CURRENT $APP/Contents/Resources/lib
    done
}

#for a given binary, redirect its search path for its dependencies to the $APP/Contents/Resources/lib folder
rehomeSharedLibraries(){
    OTOUT=$(otool -L $1)
    NLIBS=$(echo "$OTOUT" | grep "$LIBS" | wc -l)
    CURRENTBASE=$(basename $1)
    install_name_tool -id "@executable_path/../Resources/lib/$CURRENTBASE" $1
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$OTOUT" | grep "$LIBS" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        CURRENTBASE=$(basename $CURRENT)
        install_name_tool -change "$CURRENT" "@executable_path/../Resources/lib/$CURRENTBASE" $1
    done
}

#go through all of the files in the $APP/Contents/Resources/lib folder and copy dependencies there
copyLibraryDependencies(){
    NLIBSD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSD"
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copySharedLibraries "$APP/Contents/Resources/lib/$CURRENTD"
    done
    return $NLIBSD
}

#apply rehomeSharedLibraries to everything in the bundle's /lib
redirectLibraryDependencies(){
    NLIBSR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    for((k=1; k<=$NLIBSR; k++))
    do
        CURRENTR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$k 'FNR==i')
        rehomeSharedLibraries "$APP/Contents/Resources/lib/$CURRENTR"
    done
}

#if the binary has fixed library paths, send them to the bundle /lib
fixBinPaths(){
    OTOUT=$(otool -l $BINPATH | grep path | grep -v @exec | awk '{print $2}')
    NLIBS=$(echo "$OTOUT" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$OTOUT" | awk -v i=$i 'FNR==i')
        install_name_tool -rpath "$CURRENT" "@executable_path/../Resources/lib" $BINPATH
    done
}

#first dependency pass, direct dependencies of the executable
copySharedLibraries $BINPATH

#apply to the copied libraries repeatedly until the number of libraries converges
echo "Checking  local dependencies"
LASTCOPY="$(copyLibraryDependencies)"
echo "$LASTCOPY libraries"
CURRENTCOPY="$(copyLibraryDependencies)"
echo "$CURRENTCOPY libraries"
while [ $LASTCOPY -lt $CURRENTCOPY ]
do
    LASTCOPY=$CURRENTCOPY
    CURRENTCOPY="$(copyLibraryDependencies)"
    echo "$CURRENTCOPY libraries"
done
echo "I think I have them all."

#correct paths of all the files
echo "Redirecting dependencies to .app bundle"
rehomeSharedLibraries $BINPATH
echo "Fixing intra-library dependencies"
redirectLibraryDependencies
echo "Fixing rpath"
fixBinPaths