<p align="center"><img src="Source/BuildResources/icons/hicolor/512x512/apps/io.github.NickKarpowicz.LightwaveExplorer.png" width="256" height="256"></p>

## <p align="center">Lightwave Explorer</p>
<p align="center">Nick Karpowicz<br>
 Max Planck Institute of Quantum Optics</p>

---
New! 

Publication!
- N. Karpowicz, Open source, heterogenous, nonlinear-optics simulation. [*Optics Continuum* **2**, 2244-2254 (2023).](https://opg.optica.org/optcon/fulltext.cfm?uri=optcon-2-11-2244&id=540999)
  
Tutorials on YouTube!
- Tutorial 1: <a href="https://youtu.be/J1-dh6V5flQ">Introduction and quick start, with walkthrough of simulating SHG and OPA</a></li>
- Tutorial 2: <a href="https://youtu.be/7osRWaI91nk">Understanding and ensuring convergence</a>

---
### Latest release: 2023.09
<p>This release adds the following fixes and improvements:</p>
  <ul>
  <li>Adds a new CPU-based propagation mode on Linux/Windows using C++ std::execution::par instead of OpenMP, which on some computers is nearly twice as fast</li>
  <li>MacOS version now uses Grand Central Dispatch for multithreading, with a slight performance improvement.
  <li>Applies PEMDAS rules to the input arguments of sequence functions. In short: previously an input of 1+2*3 would give 9 because the operators were resolved from left to right, now it has the correct behavior (i.e. 7).</li>
  <li>Multiple small optimizations to propagation kernels with a small (percent-level) performance improvement</li>
  <li>Various bug fixes and stability enhancements (mostly regarding mutex locks).</li>
  <li>A button to collapse the data input panel so that the plots and images fill the whole window</li>
  </ul>
  
---


### What and why

Lightwave explorer is an open source nonlinear optics simulator, intended to be fast, visual, and flexible for students and researchers to play with ultrashort laser pulses and nonlinear optics without having to buy a laser first.

<p style="text-align: center;"><img src="Documentation/Images/Interface_screenshot.png"></p>

The simulation was written CUDA in order to run quickly on modern graphics cards. I've subsequently generalized it so that it can be run in two other ways: SYCL on CPUs and Intel GPUs, and using OpenMP to run on CPUs. Accordingly, I hope that the results are fast enough that even complicated systems can be simulated within a human attention span.

---

#### Main goals:
 - _Easily extensible database of materials:_ Eveything the program knows about nonlinear materials comes from a human-readable text file giving the appropriate coefficients and tensors. If you want to use a new material, or you've done a measurement in a new range where typical extrapolations from older data isn't relevant, it's easy to add and correct. There are places for references for the key parameters, and these references are stored in the saved simulation results for future reference. Especially if you have simulations that you checked against experiments, I'd be very happy for you to add your crystal definitions to the central database in the project Github.
 - _Accurate modeling of nonlinear optics_ using multiple, user-selectable physical models, including the unidirectional nonlinear wave equation and finite-difference time-domain approaches. This allows calculations that accommodate large systems where forward-propagation is an appropriate assumption, but also of etalon effects in thin crystals where reflections cannot be neglected.
 - _Efficient code so that complicated systems can be simulated in 3D:_ Real laser pulses can be messy, and if they weren't so before a nonlinear crystal, there's a good chance they are after (but not always). If things are slow, it's hard to go beyond one dimension on tolerable time scales, and then you miss out on the whole weird world of spatiotemporal couplings. Here you have options for rather fast simulations when there's a symmetry to apply (e.g. cylindrical or along one Cartesian dimension), alongside fully 3D propagation. Runs natively on both GPU and CPU to make use of whatever you have to work with.
 - _A graphical interface that lets you see what you're doing:_ A lot of us think in visual terms. Being able to adjust and scan parameters and immediately see what happens can really make it easier to understand what you're looking at. 
 - _A flexible sequence mode:_ By stringing together elements, not just nonlinear crystals but also spherical or parabolic mirrors, apertures, filters, free space propagation and other elements, simulate how  one interaction affects another. Sequences of events can be scripted and even programmed with loop functions to see how things change over the course of repeated interactions.
 - _Fitting modes:_ Sometimes the data that we measure depends in an interesting way on a parameter, and we'd actually like to go back and figure out what that parameter was from the data. Solving this kind of inverse problem can be tough when the parameter lives inside a partial differential equation, but by simulating the whole thing and doing a fit, you have a chance to do it! The fitting algorithm can be used to narrow down a huge space of variables to come at your best estimation of what was happening in an experiment, or to adjust your experimental system to maximize output at a given frequency.
 - _A Python module for easy postprocessing of the results:_ I hope that you get something interesting out that you want to plot and maybe publish. One of the nicest platforms for making nice plots is Python in my opinion (that's why the documentation is in a Jupyter notebook), so purely out of self interest I tried to make it easy to load the results in Python. The module also has some functions related to typical operations you'd like to do on the data to make it easy for all of us. The program also gives you a Matlab loading script for those who want to use that.
 - _Multiplatform:_ Works on Windows, Linux, and Mac.
 - _Command line interface for running on Linux/clusters:_ the simulation core can be compiled as a command line application to be controlled via the SLURM system. The GUI app can automatically configure the SLURM script, as well. I use this to run it on the clusters of the Max Planck Society, and other institutes and universities likely have similar systems. This lets you do a lot more if your personal resources are limited but you want to run simulations on a large grid or cover a lot of different parameters!

---

  ### Publications
  Lightwave Explorer has been used to perform the nonlinear optics simulations in the following papers!
  - Maciej Kowalczyk, *et al.*, Ultra-CEP-stable single-cycle pulses at 2.2 µm. [*Optica* **10**, 801-811 (2023).](https://opg.optica.org/optica/fulltext.cfm?uri=optica-10-6-801)
  - Najd Altwaijry, *et al.*, Broadband Photoconductive Sampling in Gallium Phosphide. [*Advanced Optical Materials* **11**, 2202994 (2023).](https://onlinelibrary.wiley.com/doi/full/10.1002/adom.202202994)
  - Hadil Kassab, *et al.*, In-line synthesis of multi-octave phase-stable infrared light, [*Optics Express* **31**, 24862 (2023).](https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-15-24862)

---
  ### Installation on a Windows PC
  You can get the [latest version directly from the Github releases](https://github.com/NickKarpowicz/LightwaveExplorer/releases/tag/2023.09.01), where there's a LightwaveExplorerWin64.zip file that you can just download, extract, and run.

  If you want to use SYCL for propagation, you need to install the [Intel® oneAPI DPC++/C++ Compiler Runtime for Windows](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html).

  The Python module for working with the results is [here](https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/Documentation/LightwaveExplorer.py) in this repo; I'd recommend putting it somewhere in your Python path if you're going to work with it a lot, otherwise just copy it into your working folder.

---
  ### Installation on Linux
  The easiest way to install on Linux is using the Flatpak, which is [available on Flathub!](https://flathub.org/apps/io.github.NickKarpowicz.LightwaveExplorer)

  If your system's graphical installer integrates Flathub, you might already be able to find it there (it's quite new, as of August 2023), or you can install it from the command line with:
  ```
  flatpak install flathub io.github.NickKarpowicz.LightwaveExplorer
  ```

  Once it's installed it should show up in your system menu in the Science category. Otherwise you can also run it with 
  
  ```
  flatpak run io.github.NickKarpowicz.LightwaveExplorer
  ```

  Of course you can also build it yourself using the instructions below. That's currently what you'll need to do to use the SYCL propagator, since that's kinda still in beta.

---

### Installation on Mac

The Mac version is also available [directly from the Github relases](https://github.com/NickKarpowicz/LightwaveExplorer/releases/tag/2023.09.01). In there, there's a LightwaveExplorerMacOS.zip file that you can just download and run. You might have to right-click on it and run it from there, and accept a prompt that you actually want to run it.

This version makes use of the FFTW library for Fourier transforms and is therefore released under the GNU Public License v3.

The application bundle contains all the required files. If you want to edit the crystal database or default settings, open the app as a folder (right click or control-click on the app and select "show package contents") - You will find them in the Resources folder.

I can't give you a binary that is native to the new Apple M1 and M2 chips, due to how Apple locks down their hardware, and how I refuse to pay them $100/year. However, you can compile it yourself using the instructions below (it's actually pretty easy).

---
  ### How do I know which configuration to run?
  At the bottom of the window, you'll see two pulldown menus marked "Config" - these let you choose whether the simulation runs in CUDA, SYCL, or OpenMP. It should start set to the fastest option for your system, but if you don't have the right drivers/runtimes, it might step down to something else that is present. OpenMP is typically the slowest (except on Linux), but will run on basically any system.

  - If you have an Intel CPU, chances are it will have an integrated GPU, which can be used by SYCL. In my experience, the ones that show up as "Iris Graphics" are actually much faster than running the code on the actual CPU, while the ones named "HD Graphics" are sometimes slower. Just try them both.

  - The SYCL code compiles itself for your specific hardware the first time it runs. So, the first run will be slower than the rest - subsequent runs will be faster!

  - The second pulldown is for offloading work onto other parts of your system. For example, if you are running a big batch of simulations on a CUDA-capable GPU, you can send some of them to the CPU to work on. This is what the second menu and following number do: chose the offload target, and the number to send to it. You don't have to use this if you don't want to.

  - Basically, the order in which you should choose the backends is: CUDA if you have something which supports it. If not, SYCL if you're on windows and/or have a supported Intel integrated GPU. You're running on CPU in Linux, the GPL3 version might be faster.

  ---
  ### Compilation in Visual Studio
  LWE was developed in Visual Studio, so you'll need that. You will also need [vcpkg](https://vcpkg.io), and use that to install dlib (make sure you get the 64-bit version, not the 32-bit one).
  
   Next, install the [CUDA development kit](https://developer.nvidia.com/cuda-downloads) from NVIDIA, and [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) (including the Math Kernel Library and the DPC++ compiler).
   
   Next, you'll need a compiled version of [GTK4](http://gtk.org). The resulting compiled thing should be kept in a folder next to the LightwaveExplorer folder (e.g. they're both in the same parent folder).


  ---
### Compiling the GUI app on Linux (Easy CPU-only version)
The easiest version to compile on Linux is the GPL3 version, which doesn't include the CUDA or OneAPI propagators. This means it will _only_ run on CPU, but if you don't have a compatible GPU anyway, it makes use of FFTW for the FFTs, which may be faster on your hardware in any case.

The prerequisite packages are: gcc, cmake, GTK4, and FFTW (plus git to download the repo). Their exact names in your package manager may vary... 

If you are on an Ubuntu-based distro, you can use this to grab everything:

```
sudo apt install gcc git cmake libgtk-4-1 libgtk-4-dev libfftw3-3 libfftw3-dev
```

On OpenSUSE Tumbleweed, I needed:
```
sudo zypper install git gcc-c++ cmake gtk4-devel fftw-devel
```

Once you have that, type the following into the terminal:

```
git clone https://github.com/NickKarpowicz/LightwaveExplorer
mkdir LightwaveExplorer/build
cd LightwaveExplorer/build
cmake ..
make
```

It should then spend a bit of time building and finally produce a LightwaveExplorer executable in the build directory.

You can install the application in your default location (probably /usr/local) with the command:

```
sudo cmake --install .
```

If you want to install it somewhere else, append --prefix "/where/you/want/it/to/go"

Installing will also place the CrystalDatabase.txt and DefaultValues.ini text files in the /share/LightwaveExplorer folder alongside the /bin folder where the binary ends up. You can edit these freely to add crystals or changes the values that populate the program's interface when it starts.

### Compiling the GUI app on Linux (CUDA and SYCL version)

  You'll need everything required to build the GPL3 version above, except for FFTW, which isn't used in this version. I'd recommend building the one above first to make sure it works. Next, install the prerequisites for this version:

  - [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html) - I recommend using the online installer so you can choose what to install (the full thing is quite large). You'll at least need the DPC++ compiler, the Math Kernel Library, and Thread Building Blocks (and their prerequisites).
  - [NVIDIA CUDA](https://developer.nvidia.com/cuda-downloads) - this might already be in your package manager, but I'd recommend at least version 11.6.

  Now that you have everything, in order to build the full version, first you have to set the OneAPI environment variables, typically with:
  ```
  . ~/intel/oneapi/setvars.sh
  ```
  if you installed OneAPI as a normal user or
  ```
  . /opt/intel/oneapi/setvars.sh
  ```
  if you installed as root.

  Then, build the executable with:
  ```
  git clone https://github.com/NickKarpowicz/LightwaveExplorer
  mkdir LightwaveExplorer/build
  cd LightwaveExplorer/build
  cmake -DMAKEFULL=TRUE -DCMAKE_CXX_COMPILER=icpx -DCMAKE_CUDA_HOST_COMPILER=clang++ -DCMAKE_CUDA_COMPILER=nvcc -DCMAKE_CUDA_ARCHITECTURES=75 ..
  make
  ```
  
  If it doesn't fail, you should now have an executable file named LightwaveExplorer in the build folder. You can install using the same process as the CPU-only version above.

---

  ### Compiling on Mac

  The first thing you'll need is [Homebrew](https://brew.sh/). If you go there, you'll see a command that you have to run in the terminal. Just paste it and follow the instructions.

  I also made a build script that you can run in the same way; just copy and paste the command below that matches your system and it will compile everything it needs and put the application in your Applications folder. It will take a while, so go get a coffee!

  Apple Silicon (M1, M2, .etc) version:

  ```
  curl -s https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/Source/BuildResources/macAutoBuild.sh | zsh -s
  ```

  Intel version:
  ```
  curl -s https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/Source/BuildResources/macAutoBuildIntel.sh | zsh -s
  ```
  
---
  ### Compilation on clusters
  
  A script is provided to compile the CUDA command line version on Linux. This is made specifically to work on the clusters of the MPCDF but will likely work with small modifications on other distributions depending on the local environment. The CUDA development kit and Intel OneAPI should be available in advance. With these prerequisites, the following command should work:
  ```
curl -s https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/Source/BuildResources/compileCommandLineLWEfromRepos.sh | tcsh -s
 ```
 On other clusters you might have to instead dowload the script (e.g. with wget) and change it to suit that system before you run it.

 If you have the GUI version installed locally, you can set up your calculation and then generate a SLURM script to run on the cluster (it will tell you what to do).

 ---

  ### Libraries used
Thanks to the original authors for making their work available! They are all freely available, but of course have their own licenses .etc.
  - [NVIDIA CUDA](https://developer.nvidia.com/cuda-toolkit): This provides the basic CUDA runtime, compiler, and cuFFT, for running the simulations on NVIDIA GPUs, and is the basis of the fastest version of this code.
  - [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html), specifically the [Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.cw3ci4): This is used for performing fast fourier transforms when running in CPU mode. The DPC++ compiler allows the program to run on both CPUs and a wider range of GPUs, including the integrated ones on Intel chips. I found that on my rather old laptop, SYCL on the GPU is several times faster than running on CPU, so it's useful even for systems without dedicated GPUs.
  - [Dlib](http://dlib.net/): This library is the basis of the optimization routines. I make use of the global optimization functions for the fitting/optimization modes. The library is [available on Github](https://github.com/davisking/dlib), and their excellent documentation and further information is on the [main project website](http://dlib.net/).
  - [GTK](https://www.gtk.org): The new version of the user interface uses GTK 4; this is why it looks pretty much the same on Windows, Linux, and Mac. It was pretty easy to get working cross-platform, which again is nice for the goal that everybody should be able to reproduce calculations in LWE.
  - [FFTW](https://www.fftw.org/): This is used for Fast Fourier Transforms in the GPL 3.0 version (i.e. the CPU-only Linux and Mac versions). On a given CPU this is on average the fastest FFT you can find.
  
  ---

  ### Programming note

  The code is written in a "trilingual" way - a single core code file is compiled (after some includes and preprocessor definitions) by the three different compilers, Nvidia nvcc, a c++ compiler (either Microsoft's, g++, or clang++ have all worked), and Intel dpc++. 

  Although CUDA was the initial platform and what I use (and test) most extensively, I've added two additional languages for those who don't have an Nvidia graphics card. 
  
  One is in c++, with multithreading done with OpenMP. 
  
  The other language is SYCL. This also allows the simulation to run on the CPU and should allow it to run on Intel's graphics cards, as well as the integrated graphics of many Intel CPUs. The same language should be able to run on AMD cards, but support for the DPC++ toolchain with the HipSYCL backend is quite new, and I don't have an AMD card to test it on. 
    
  The different architectures are using the same algorithm, aside from small differences in their floating point math and intrinsic functions. So when I make changes or additions, there will never be any platform gaining over the other (again, reproducibility by anyone is part of the goals here).
