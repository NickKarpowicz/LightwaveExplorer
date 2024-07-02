brew install cmake ninja make pkgconfig libomp qt cairo fmt wget
mkdir LightwaveExplorerBuild
cd LightwaveExplorerBuild
BASE_DIR=$(pwd)
echo "Building fftw..."
curl -O http://fftw.org/fftw-3.3.10.tar.gz >& /dev/null
tar -xvf fftw-3.3.10.tar.gz >& /dev/null
cd fftw-3.3.10 
curl -O https://raw.githubusercontent.com/andrej5elin/howto_fftw_apple_silicon/main/fftw-3-3-10-configure-diff.txt >& /dev/null
patch configure fftw-3-3-10-configure-diff.txt >& /dev/null

./configure --enable-threads --enable-neon --enable-armv8-cntvct-el0 --enable-silent-rules >& /dev/null
make >& /dev/null
make DESTDIR=$BASE_DIR/fftw install >& /dev/null

./configure --enable-threads --enable-neon --enable-armv8-cntvct-el0 --enable-float --enable-silent-rules >& /dev/null
make clean >& /dev/null
make >& /dev/null
make DESTDIR=$BASE_DIR/fftw install >& /dev/null

cd ..
rm -rf fftw-3.3.10

git clone -b QT_interface --single-branch --depth 1 https://github.com/NickKarpowicz/LightwaveExplorer >& /dev/null
git clone --depth 1 --branch v19.24.2 https://github.com/davisking/dlib >& /dev/null


cd LightwaveExplorer
./Source/BuildResources/makeMacAppQt.sh
cp -r build/LightwaveExplorer.app /Applications/
