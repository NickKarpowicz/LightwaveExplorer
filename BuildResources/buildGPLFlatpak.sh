rm -rf build
mkdir build
cd build
git clone https://github.com/davisking/dlib
mkdir dlibtmp
cp -rf dlib/dlib dlibtmp/dlib
rm -rf dlib
mv dlibtmp dlib
cp ../BuildResources/io.NickKarpowicz.LightwaveExplorer.yml ./
flatpak-builder build-dir io.NickKarpowicz.LightwaveExplorer.yml --force-clean
flatpak-builder --user --install --force-clean build-dir  io.NickKarpowicz.LightwaveExplorer.yml