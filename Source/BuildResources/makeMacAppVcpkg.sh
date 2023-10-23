#!/bin/bash -l
BIN=LightwaveExplorer
APP=build/${BIN}.app
BINPATH=${APP}/Contents/MacOS/${BIN}

#Shared libraries location
LIBS="/Users"

#build executable
rm -rf build
mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../../vcpkg/scripts/buildsystems/vcpkg.cmake .. -G Ninja
ninja
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
cp Source/BuildResources/AppIcon.icns $APP/Contents/Resources

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

#for a given binary, redirect its search path for its dependencies to the $APP/Contents/Resources/lib folder
rehomeSharedLibrariesFromOpt(){
    OTOUT=$(otool -L $1)
    NLIBS=$(echo "$OTOUT" | grep "$LIBS" | wc -l)
    CURRENTBASE=$(basename $1)
    install_name_tool -id "@executable_path/../Resources/lib/$CURRENTBASE" $1
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$OTOUT" | grep "opt" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
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
        rehomeSharedLibrariesFromOpt "$APP/Contents/Resources/lib/$CURRENTR"

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