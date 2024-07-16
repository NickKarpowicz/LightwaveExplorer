$RootDir = (Get-Item -Path ".\").FullName
$BuildDir = "$RootDir\build"

if (-Not (Test-Path -Path $BuildDir)) {
    New-Item -ItemType Directory -Path $BuildDir
}

Set-Location -Path $BuildDir

cmake -DMAKESYCL=1 .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -T "Intel(R) oneAPI DPC++ Compiler 2024"
cmake --build . --config Release

cmake --fresh .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -DCMAKE_CUDA_ARCHITECTURES="75;86"
cmake --build . --config Release
