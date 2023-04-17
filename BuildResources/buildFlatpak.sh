rm -rf build
mkdir build
cd build
git clone https://github.com/davisking/dlib --branch v19.24
mkdir dlibtmp
cp -rf dlib/dlib dlibtmp/dlib
rm -rf dlib
mv dlibtmp dlib
cp ../../l_BaseKit_p_2023.1.0.46401_offline.sh l_BaseKit_p_2023.1.0.46401_offline.sh
cp ../BuildResources/io.NickKarpowicz.LightwaveExplorer.yml ./
flatpak-builder build-dir io.NickKarpowicz.LightwaveExplorer.yml --force-clean
flatpak-builder --user --install --force-clean build-dir  io.NickKarpowicz.LightwaveExplorer.yml