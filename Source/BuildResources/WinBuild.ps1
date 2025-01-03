$RootDir = (Get-Item -Path ".\").FullName
$BuildDir = "$RootDir\build"

if (-Not (Test-Path -Path $BuildDir)) {
    New-Item -ItemType Directory -Path $BuildDir
}

Set-Location -Path $BuildDir

cmake -DMAKESYCL=1 .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -T "Intel(R) oneAPI DPC++ Compiler 2025"
cmake --build . --config Release

cmake --fresh .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE="C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" -DCMAKE_CUDA_ARCHITECTURES="75;86"
cmake --build . --config Release

Copy-Item -Path Release -Destination LightwaveExplorerWin64 -Recurse 
Remove-Item .\LightwaveExplorerWin64\LightwaveExplorerSYCL.exp -Force
Remove-Item .\LightwaveExplorerWin64\LightwaveExplorerSYCL.lib

if ($args.Count -gt 0){
    $sevenZipPath = $args[0]
    & $sevenZipPath a -tzip LightwaveExplorerWin64.zip LightwaveExplorerWin64\*
}

