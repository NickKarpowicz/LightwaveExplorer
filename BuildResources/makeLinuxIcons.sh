#!/bin/bash

# Path to the 512x512 pixel PNG icon
icon_path="../BuildResources/AppImage/ico512.png"
sizes=("256x256" "128x128" "64x64" "48x48" "32x32" "24x24" "16x16")
mkdir icons/hicolor/512x512/apps
cp ${icon_path} icons/hicolor/512x512/apps/
for size in "${sizes[@]}"; do
  output_path="icons/hicolor/${size}/apps/io.github.NickKarpowicz.LightwaveExplorer.png"
  convert "${icon_path}" -resize "${size}" "${output_path}"
done

#buid the .desktop file
touch io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "[Desktop Entry]" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Icon=io.github.NickKarpowicz.LightwaveExplorer" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Type=Application" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Categories=Utility" >> $APP/io.github.NickKarpowicz.LightwaveExplorer.desktop

# # Update icon cache
# sudo gtk-update-icon-cache /usr/share/icons/hicolor
