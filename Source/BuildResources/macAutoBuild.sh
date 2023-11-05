brew install cmake ninja make pkgconfig libomp libtool autoconf automake
BASE_DIR=$(pwd)
echo $BASE_DIR
curl -O http://fftw.org/fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar.gz
cd fftw-3.3.10 
curl -O https://raw.githubusercontent.com/andrej5elin/howto_fftw_apple_silicon/main/fftw-3-3-10-configure-diff.txt
patch configure fftw-3-3-10-configure-diff.txt

./configure --enable-threads --enable-neon --enable-armv8-cntvct-el0 --enable-silent-rules
make
make DESTDIR=$BASE_DIR/fftw install

./configure --enable-threads --enable-neon --enable-armv8-cntvct-el0 --enable-float --enable-silent-rules
make clean
make
make DESTDIR=$BASE_DIR/fftw install

cd ..
rm -rf fftw-3.3.10

git clone --depth 1 https://github.com/NickKarpowicz/LightwaveExplorer
git clone --depth 1 --branch v19.24.2 https://github.com/davisking/dlib
git clone https://github.com/microsoft/vcpkg
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install gtk gcem fmt

cd LightwaveExplorer
./Source/BuildResources/makeMacAppVcpkg.sh
cp -r build/LightwaveExplorer.app /Applications/
