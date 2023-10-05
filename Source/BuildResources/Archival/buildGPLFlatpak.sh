rm -rf build
mkdir build
cd build
cp ../BuildResources/io.github.NickKarpowicz.LightwaveExplorerGPL.yml ./
cmake ..
make
mkdir -p appdir/appdir
cmake --install . --prefix=appdir/appdir
flatpak-builder build-dir io.github.NickKarpowicz.LightwaveExplorerGPL.yml --force-clean
flatpak-builder --user --install build-dir --force-clean io.github.NickKarpowicz.LightwaveExplorerGPL.yml