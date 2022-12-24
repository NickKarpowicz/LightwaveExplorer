#!/bin/bash -l
APP=LightwaveExplorer.app
BIN=../LightwaveExplorer
rm -rf $APP
mkdir $APP
mkdir $APP/Contents/
mkdir $APP/Contents/MacOS/
mkdir $APP/Contents/Resources/
mkdir $APP/Contents/Resources/lib
mkdir $APP/Contents/Resources/bin
mkdir $APP/Contents/Resources/etc
mkdir $APP/Contents/Resources/share

cp $BIN $APP/Contents/MacOS
cp ../CrystalDatabase.txt $APP/Contents/Resources
cp ../DefaultValues.ini $APP/Contents/Resources
cp AppIcon.icns $APP/Contents/Resources
cp macplistbase.plist $APP/Contents/info.plist

copySharedLibraries(){
    NLIBS=$(otool -L $1 | grep "/usr/local" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(otool -L $1 | grep "/usr/local" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        cp -n $CURRENT $APP/Contents/Resources/lib
    done
}
#install_name_tool -change "/usr/local/opt/glib/lib/libgobject-2.0.0.dylib" "@executable_path/../Resources/lib/libgobject-2.0.0.dylib" Resources/lib/libatk-1.0.0.dylib

rehomeSharedLibraries(){
    OTOUT=$(otool -L $1)
    NLIBS=$(echo "$OTOUT" | grep "/usr/local" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$OTOUT" | grep "/usr/local" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$i 'FNR==i')
        CURRENTBASE=$(basename $CURRENT)
        install_name_tool -change "$CURRENT" "@executable_path/../Resources/lib/$CURRENTBASE" $1
        echo "b: $CURRENTBASE $1"
    done
}

copyLibraryDependencies(){
    NLIBSD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSD"
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copySharedLibraries $APP/Contents/Resources/lib/$CURRENTD
    done
    return $NLIBSD
}

redirectLibraryDependencies(){
    NLIBSR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | wc -l | tr -d '[:blank:]')
    echo "$NLIBSR"
    for((k=1; k<=$NLIBSR; k++))
    do
        CURRENTR=$(ls -1a $APP/Contents/Resources/lib | grep "lib" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$k 'FNR==i')
        rehomeSharedLibraries $APP/Contents/Resources/lib/$CURRENTR
    done

}

copySharedLibraries $BIN

echo "Checking local dependencies"
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

echo "Redirecting dependencies to App folder"
redirectLibraryDependencies





