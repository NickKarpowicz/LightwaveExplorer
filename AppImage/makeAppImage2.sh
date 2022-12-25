#!/bin/bash -l

#Target AppDir
APP=LightwaveExplorer.AppDir

#Target executable
BIN=../LightwaveExplorer

#Search paths
SEARCHPATHS="/opt/intel /usr"

# cd ..
# make cpuonly
# cd AppImage

#prepare the AppDir structure
rm -rf $APP
mkdir -p $APP
mkdir -p $APP/usr
mkdir -p $APP/usr/lib
mkdir -p $APP/usr/bin
mkdir -p $APP/usr/share
mkdir -p $APP/usr/share/glib-2.0
mkdir -p $APP/usr/share/LightwaveExplorer

cp $BIN $APP/usr/bin
cp ../CrystalDatabase.txt $APP/usr/bin
cp ../DefaultValues.ini $APP/usr/bin
cp ico512.png $APP/LightwaveExplorer.png
cp AppRun $APP/AppRun
cp -r /usr/share/glib-2.0 LightwaveExplorer.AppDir/usr/share
cp /usr/local/cuda-12.0/targets/x86_64-linux/lib/libcufft.so.11 LightwaveExplorer.AppDir/usr/lib/libcufft.so.11
cp /usr/local/cuda-12.0/targets/x86_64-linux/lib/libcudart.so.12 LightwaveExplorer.AppDir/usr/lib/libcudart.so.12
cp /usr/lib/x86_64-linux-gnu/libgomp.so.1 LightwaveExplorer.AppDir/usr/lib/libgomp.so.1
cp /usr/lib/x86_64-linux-gnu/libgtk-4.so.1 LightwaveExplorer.AppDir/usr/lib/libgtk-4.so.1
cp /usr/lib/x86_64-linux-gnu/libpangocairo-1.0.so.0 LightwaveExplorer.AppDir/usr/lib/libpangocairo-1.0.so.0
cp /usr/lib/x86_64-linux-gnu/libpango-1.0.so.0 LightwaveExplorer.AppDir/usr/lib/libpango-1.0.so.0
cp /usr/lib/x86_64-linux-gnu/libharfbuzz.so.0 LightwaveExplorer.AppDir/usr/lib/libharfbuzz.so.0
cp /usr/lib/x86_64-linux-gnu/libgdk_pixbuf-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libgdk_pixbuf-2.0.so.0
cp /usr/lib/x86_64-linux-gnu/libcairo-gobject.so.2 LightwaveExplorer.AppDir/usr/lib/libcairo-gobject.so.2
cp /usr/lib/x86_64-linux-gnu/libcairo.so.2 LightwaveExplorer.AppDir/usr/lib/libcairo.so.2
cp /usr/lib/x86_64-linux-gnu/libgobject-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libgobject-2.0.so.0
cp /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0 LightwaveExplorer.AppDir/usr/lib/libglib-2.0.so.0

# cp $ONEAPI_ROOT/tbb/2021.8.0/lib/intel64/gcc4.8/*.so* LightwaveExplorer.AppDir/usr/lib/
# cp -r $ONEAPI_ROOT/compiler/2023.0.0/linux/lib/* LightwaveExplorer.AppDir/usr/lib/
# cp -r $ONEAPI_ROOT/compiler/2023.0.0/linux/lib/x64/* LightwaveExplorer.AppDir/usr/lib/
# rm -rf LightwaveExplorer.AppDir/usr/lib/x64
# rm -rf LightwaveExplorer.AppDir/usr/lib/clc
# rm -rf LightwaveExplorer.AppDir/usr/lib/oclfpga
# rm -rf LightwaveExplorer.AppDir/usr/lib/clang
# rm LightwaveExplorer.AppDir/usr/lib/libOclCpuBackEnd_emu.so.2022.15.12.0
# rm LightwaveExplorer.AppDir/usr/lib/icx-lto.so
# cp -r $ONEAPI_ROOT/compiler/2023.0.0/linux/compiler/lib/intel64_lin/*.so* LightwaveExplorer.AppDir/usr/lib/

copySharedLibraries(){
    LIBLIST1=$(ldd $1 | grep "=>" | awk ' {print $3 }')
    NLIBS=$(echo "$LIBLIST1" | wc -l)
    for((i=1; i<=$NLIBS; i++))
    do
        CURRENT=$(echo "$LIBLIST1" | awk -v i=$i 'FNR==i')
        cp -n $CURRENT $APP/usr/lib
    done
}

copyLibraryDependencies(){
    LSLIST=$(ls -1a $APP/usr/lib | grep ".so")
    NLIBSD=$(ls -1a $APP/usr/lib | grep ".so" | wc -l )
    for((j=1; j<=$NLIBSD; j++))
    do
        CURRENTD=$(ls -1a $APP/usr/lib | grep ".so" | sed 's/([^)]*)//g' | tr -d '[:blank:]' | awk -v i=$j 'FNR==i')
        copySharedLibraries "$APP/usr/lib/$CURRENTD"
    done
    return $NLIBSD
}

copySharedLibraries $BIN

echo "Checking local dependencies"
LASTCOPY="$(copyLibraryDependencies)"
CURRENTCOPY="$(copyLibraryDependencies)"
while [[ $LASTCOPY -lt $CURRENTCOPY ]]
do
    LASTCOPY="$CURRENTCOPY"
    CURRENTCOPY="$(copyLibraryDependencies)"
done
echo "I think I have them all."

touch $APP/usr/bin/LaunchLWE.sh
echo "#!/bin/bash -l" >>  $APP/usr/bin/LaunchLWE.sh
echo "export OCL_ICD_FILENAMES=libintelocl_emu.so:libalteracl.so:libintelocl.so" >>  $APP/usr/bin/LaunchLWE.sh
echo "LightwaveExplorer" >>  $APP/usr/bin/LaunchLWE.sh
chmod +x $APP/usr/bin/LaunchLWE.sh

touch $APP/LightwaveExplorer.desktop
echo "[Desktop Entry]" >> $APP/LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> $APP/LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> $APP/LightwaveExplorer.desktop
echo "Icon=LightwaveExplorer" >> $APP/LightwaveExplorer.desktop
echo "Type=Application" >> $APP/LightwaveExplorer.desktop
echo "Categories=Utility" >> $APP/LightwaveExplorer.desktop


ARCH=x86_64 ./appimagetool-x86_64.AppImage $APP




