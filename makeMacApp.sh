#!/bin/bash -l
APP=LightwaveExplorer.app
BIN=LightwaveExplorer

#brew install make gcc fftw gtk4

#detect local architecture
LOCALARCH=$(arch)

#Homebrew libraries location, 
# on intel it will be /usr/local
# on Arm64 it will be /opt/homebrew
LIBS="/usr/local"

make mac

#set up the directory structure of the .app
rm -rf $APP
mkdir $APP
mkdir $APP/Contents/
mkdir $APP/Contents/MacOS/
mkdir $APP/Contents/Resources/
mkdir $APP/Contents/Resources/lib
mkdir $APP/Contents/Resources/bin
mkdir $APP/Contents/Resources/etc
mkdir $APP/Contents/Resources/share

#copy in the executable, databases and icons
cp $BIN $APP/Contents/MacOS
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
copyx86LibraryDependencies(){
    NLIBSD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSD"
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copySharedLibraries "$APP/Contents/Resources/lib/$CURRENTD"
    done
    return $NLIBSD
}

redirectx86LibraryDependencies(){
    NLIBSR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    for((k=1; k<=$NLIBSR; k++))
    do
        CURRENTR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$k 'FNR==i')
        rehomeSharedLibraries "$APP/Contents/Resources/lib/$CURRENTR"
    done
}

#first dependency pass, direct dependencies of the executable
copySharedLibraries $BIN

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
rehomeSharedLibraries $BIN
redirectx86LibraryDependencies
install_name_tool -rpath /usr/local/Cellar/llvm/15.0.6/lib @executable_path/../Resources/lib $BIN
cp $BIN $APP/Contents/MacOS/
make clean
