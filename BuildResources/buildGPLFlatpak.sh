rm -rf build
mkdir build
cd build
git clone https://github.com/davisking/dlib --branch v19.24
mkdir dlibtmp
cp -rf dlib/dlib dlibtmp/dlib
rm -rf dlib
mv dlibtmp dlib
cp ../BuildResources/io.NickKarpowicz.LightwaveExplorerGPL.yml ./
flatpak-builder build-dir io.NickKarpowicz.LightwaveExplorerGPL.yml --force-clean
#flatpak-builder --user --install --force-clean build-dir  io.NickKarpowicz.LightwaveExplorerGPL.yml