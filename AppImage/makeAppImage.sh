#!/bin/bash -l

#Target AppDir
APP=AppDir

#Target executable
BIN=../LightwaveExplorer

cd ..
make nocuda
cd AppImage

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
cp ico512.png $APP/usr/share/icons/LightwaveExplorer.png

wget https://github.com/AppImage/AppImageKit/releases/download/continuous/AppRun-x86_64
cp AppRun-x86_64 $APP/AppRun
rm AppRun-x86_64
cp AppRun.env $APP/AppRun.env

NEWLPATH=":/usr/lib:"

#copy nvidia files manually from the output of tracefile.pl
USEDFILES=$(cat traceFileOutput.txt | grep -Fv cache | grep -Fv proc | grep -Fv run | grep nvidia)
NUSEDFILES=$(echo "$USEDFILES" | wc -l)
for((n=1; n<=NUSEDFILES; n++))
do
    CURRENTF=$(echo "$USEDFILES" | awk -v i=$n 'FNR==i')
    cp -p $CURRENTF "$APP/usr/lib"
done

#do the same for OneAPI, adding the paths to the library path
USEDFILES=$(cat traceFileOutput.txt | grep -Fv cache | grep -Fv proc | grep -Fv run | grep oneapi)
NUSEDFILES=$(echo "$USEDFILES" | wc -l)
for((n=1; n<=NUSEDFILES; n++))
do
    CURRENTF=$(echo "$USEDFILES" | awk -v i=$n 'FNR==i')
    CURRENTTARGETF=$(echo "$CURRENTF" | sed 's:/home/nick/intel/oneapi:/oneapi:g') 
    CURRENTTARGETD=$(dirname "$CURRENTTARGETF")
    mkdir -p "$APP$CURRENTTARGETD"
    cp -p $CURRENTF "$APP$CURRENTTARGETF"
    NEWPATH+='$APPBASE'"$CURRENTTARGETD:"
done
echo "$NEWPATH"

#remove path duplicates
 if [ -n "$NEWPATH" ]; then
  old_PATH=$NEWPATH:; NEWPATH=
  while [ -n "$old_PATH" ]; do
    x=${old_PATH%%:*}      
    case $NEWPATH: in
      *:"$x":*) ;;          
      *) NEWPATH=$NEWPATH:$x;;   
    esac
    old_PATH=${old_PATH#*:}
  done
  NEWPATH=${NEWPATH#:}
  unset old_PATH x
fi

#optional launcher file, only used with appimage-builder --generate it so that the environment is set
#probably not needed anymore
touch $APP/usr/bin/LaunchLWE.sh
echo "#!/bin/bash -l" >>  $APP/usr/bin/LaunchLWE.sh
echo "export OCL_ICD_FILENAMES=libintelocl_emu.so:libalteracl.so:libintelocl.so" >>  $APP/usr/bin/LaunchLWE.sh
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
#echo "Exec=LaunchLWE.sh" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

echo "Icon=LightwaveExplorer" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Type=Application" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop
echo "Categories=Utility" >> $APP/io.github.nickkarpowicz.LightwaveExplorer.desktop

#build the image
appimage-builder --skip-tests
cd .. make clean






