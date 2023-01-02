#!/bin/bash -l
APP=LightwaveExplorer.app
BIN=LightwaveExplorer
BIN86=LightwaveExplorer_x86
BINARM=LightwaveExplorer_Arm64

#detect local architecture
LOCALARCH=$(arch)

#Homebrew libraries location, obviously people other than me will have to change this
ARMLIBS="/Users/nick/arm-target"
LIBS="/usr/local"
TARGETLIBS="$APP/Contents/Resources"

#compile
make macARMonIntel

#set up the directory structure of the .app
rm -rf $APP
mkdir $APP
mkdir $APP/Contents/
mkdir $APP/Contents/MacOS/
mkdir $APP/Contents/Resources/
mkdir $APP/Contents/Resources/lib
mkdir $APP/Contents/Resources/lib_arm64
mkdir $APP/Contents/Resources/lib_x86
mkdir $APP/Contents/Resources/bin
mkdir $APP/Contents/Resources/etc
mkdir $APP/Contents/Resources/share

#copy in the executable, databases and icons
cp CrystalDatabase.txt $APP/Contents/Resources
cp DefaultValues.ini $APP/Contents/Resources
cp MacResources/AppIcon.icns $APP/Contents/Resources
cp MacResources/macplistbase.plist $APP/Contents/info.plist

#for a given binary, copy its dynamic link dependencies to the $APP/Contents/Resources/lib folder
copySharedLibraries(){
    NLIBS=$(otool -L $1 | grep "$LIBS" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(otool -L $1 | grep "$LIBS" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        cp -n $CURRENT $APP/Contents/Resources/lib_x86
    done
}

#for a given binary, copy its dynamic link dependencies to the $APP/Contents/Resources/lib folder
copyArmSharedLibraries(){
    NLIBS=$(otool -L $1 | grep "$ARMLIBS" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(otool -L $1 | grep "$ARMLIBS" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        cp -n $CURRENT $APP/Contents/Resources/lib_arm64
    done
}

universalizeSharedLibraries(){
    NLIBS=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | wc -l | tr -d '[:blank:]')
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        CURRENTTARGET=$(echo "$CURRENT" | xargs basename)
        lipo -create -output $TARGETLIBS/lib/$CURRENTTARGET $APP/Contents/Resources/lib_x86/$CURRENTTARGET $APP/Contents/Resources/lib_arm64/$CURRENTTARGET
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
rehomeArmSharedLibraries(){
    OTOUT=$(otool -L $1)
    NLIBS=$(echo "$OTOUT" | grep "$ARMLIBS" | wc -l)
    CURRENTBASE=$(basename $1)
    install_name_tool -id "@executable_path/../Resources/lib/$CURRENTBASE" $1
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$OTOUT" | grep "$ARMLIBS" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
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
        copyArmSharedLibraries "$APP/Contents/Resources/lib/$CURRENTD"
    done
    return $NLIBSD
}

#go through all of the files in the $APP/Contents/Resources/lib folder and copy dependencies there
copyArmLibraryDependencies(){
    NLIBSD=$(ls -1a $APP/Contents/Resources/lib_arm64 | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSD"
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/Contents/Resources/lib_arm64 | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copyArmSharedLibraries "$APP/Contents/Resources/lib_arm64/$CURRENTD"
    done
    return $NLIBSD
}

#go through all of the files in the $APP/Contents/Resources/lib folder and copy dependencies there
copyx86LibraryDependencies(){
    NLIBSD=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSD"
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copySharedLibraries "$APP/Contents/Resources/lib_x86/$CURRENTD"
    done
    return $NLIBSD
}

redirectArmLibraryDependencies(){
    NLIBSR=$(ls -1a $APP/Contents/Resources/lib_arm64 | grep "lib" | wc -l | tr -d '[:blank:]')
    for((k=1; k<=$NLIBSR; k++))
    do
        CURRENTR=$(ls -1a $APP/Contents/Resources/lib_arm64 | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$k 'FNR==i')
        rehomeArmSharedLibraries "$APP/Contents/Resources/lib_arm64/$CURRENTR"
    done
}

redirectx86LibraryDependencies(){
    NLIBSR=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | wc -l | tr -d '[:blank:]')
    for((k=1; k<=$NLIBSR; k++))
    do
        CURRENTR=$(ls -1a $APP/Contents/Resources/lib_x86 | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$k 'FNR==i')
        rehomeSharedLibraries "$APP/Contents/Resources/lib_x86/$CURRENTR"
    done
}

#first dependency pass, direct dependencies of the executable
copySharedLibraries $BIN86
copyArmSharedLibraries $BINARM

#make two passes through the dependencies of the copied libraries, if the list is growing
#do a while loop until it stops growing
echo "Checking arm local dependencies"
LASTCOPY="$(copyArmLibraryDependencies)"
echo "$LASTCOPY libraries"
CURRENTCOPY="$(copyArmLibraryDependencies)"
echo "$CURRENTCOPY libraries"
while [ $LASTCOPY -lt $CURRENTCOPY ]
do
    LASTCOPY=$CURRENTCOPY
    CURRENTCOPY="$(copyArmLibraryDependencies)"
    echo "$CURRENTCOPY libraries"
done
echo "I think I have them all."

echo "Checking x86 local dependencies"
LASTCOPY="$(copyx86LibraryDependencies)"
echo "$LASTCOPY libraries"
CURRENTCOPY="$(copyx86LibraryDependencies)"
echo "$CURRENTCOPY libraries"
while [ $LASTCOPY -lt $CURRENTCOPY ]
do
    LASTCOPY=$CURRENTCOPY
    CURRENTCOPY="$(copyx86LibraryDependencies)"
    echo "$CURRENTCOPY libraries"
done
echo "I think I have them all."

#correct paths of all the copied files
echo "Redirecting dependencies to App folder"
rehomeSharedLibraries $BIN86
rehomeArmSharedLibraries $BINARM
redirectArmLibraryDependencies
redirectx86LibraryDependencies
install_name_tool -rpath /usr/local/Cellar/llvm/15.0.6/lib @executable_path/../Resources/lib $BIN86
install_name_tool -rpath /usr/local/Cellar/llvm/15.0.6/lib @executable_path/../Resources/lib $BINARM

#make universal binaries
universalizeSharedLibraries
lipo -create -output "$APP/Contents/MacOS/$BIN" $BIN86 $BINARM


#remove objects and raw binaries
rm -rf $APP/Contents/Resources/lib_arm64
rm -rf $APP/Contents/Resources/lib_x86
rm $BIN86
rm $BINARM
make clean



