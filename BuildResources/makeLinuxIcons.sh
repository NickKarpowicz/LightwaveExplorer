#!/bin/bash
# Path to the 512x512 pixel PNG icon
icon_path="../BuildResources/AppImage/ico512.png"
mkdir -p icons/hicolor/512x512/apps
cp ${icon_path} icons/hicolor/512x512/apps/io.github.NickKarpowicz.LightwaveExplorer.png
mkdir -p icons/hicolor/256x256/apps
convert "${icon_path}" -resize "256x256" "icons/hicolor/256x256/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/128x128/apps
convert "${icon_path}" -resize "128x128" "icons/hicolor/128x128/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/64x64/apps
convert "${icon_path}" -resize "64x64" "icons/hicolor/64x64/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/48x48/apps
convert "${icon_path}" -resize "48x48" "icons/hicolor/48x48/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/32x32/apps
convert "${icon_path}" -resize "32x32" "icons/hicolor/32x32/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/24x24/apps
convert "${icon_path}" -resize "24x24" "icons/hicolor/24x24/apps/io.github.NickKarpowicz.LightwaveExplorer.png"

mkdir -p icons/hicolor/16x16/apps
convert "${icon_path}" -resize "16x16" "icons/hicolor/16x16/apps/io.github.NickKarpowicz.LightwaveExplorer.png"
