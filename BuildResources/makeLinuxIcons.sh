#!/bin/bash

# Path to the 512x512 pixel PNG icon
icon_path="../BuildResources/AppImage/ico512.png"
sizes=("256x256" "128x128" "64x64" "48x48" "32x32" "24x24" "16x16")
mkdir -p icons/hicolor/512x512/apps
cp ${icon_path} icons/hicolor/512x512/apps/io.github.NickKarpowicz.LightwaveExplorer.png
for size in "${sizes[@]}"; do
  mkdir -p icons/hicolor/${size}/apps
  output_path="icons/hicolor/${size}/apps/io.github.NickKarpowicz.LightwaveExplorer.png"
  convert "${icon_path}" -resize "${size}" "${output_path}"
done

#buid the .desktop file
touch io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "[Desktop Entry]" >> io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Name=LightwaveExplorer" >> io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Exec=LightwaveExplorer" >> io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Icon=io.github.NickKarpowicz.LightwaveExplorer" >> io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Type=Application" >> io.github.NickKarpowicz.LightwaveExplorer.desktop
echo "Categories=Utility" >> io.github.NickKarpowicz.LightwaveExplorer.desktop

# # Update icon cache
# sudo gtk-update-icon-cache /usr/share/icons/hicolor
