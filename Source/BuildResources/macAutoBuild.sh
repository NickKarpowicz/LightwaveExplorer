brew install cmake make pkgconfig qt cairo wget
mkdir LightwaveExplorerBuild
cd LightwaveExplorerBuild
git clone --depth 1 https://github.com/NickKarpowicz/LightwaveExplorer >& /dev/null
git clone --depth 1 --branch v19.24.2 https://github.com/davisking/dlib >& /dev/null
cd LightwaveExplorer
./Source/BuildResources/makeMacAppQt.sh
cp -r build/LightwaveExplorer.app /Applications/
