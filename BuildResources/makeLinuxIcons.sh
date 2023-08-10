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

# # Update icon cache
# sudo gtk-update-icon-cache /usr/share/icons/hicolor