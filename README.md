## <p style="text-align: center;">Lightwave Explorer</p>

<p style="text-align: center;">Nick Karpowicz</p>
<p style="text-align: center;">Max Planck Institute of Quantum optics</p>

<p style="text-align: center;"><img src="/Documentation/Images/LWEicon.png" width="200" height="200"></p>

<p style="text-align: center;">(icon made by Stable Diffusion)</p>

---


### What and why

Lightwave explorer is an open source nonlinear optics simulator, intended to be fast, visual, and flexible for students and researchers to play with ultrashort laser pulses and nonlinear optics without having to buy a laser first.

<p style="text-align: center;"><img src="/Documentation/Images/LinuxScreenshot.png"></p>

The simulation was written CUDA in order to run quickly on modern graphics cards. I've subsequently generalized it so that it can be run in two other ways: SYCL on CPUs and Intel GPUs, and using OpenMP to run on CPUs. Accordingly, I hope that the results are fast enough that even complicated system can be simulated within a human attention span.

---

#### Main goals:
 - _Easily extensible database of materials:_ Eveything the program knows about nonlinear materials comes from a human-readable text file giving the appropriate coefficients and tensors. If you want to use a new material, or you've done a measurement in a new range where typical extrapolations from older data isn't relevant, it's easy to add and correct. There are places for references for the key parameters, and these references are stored in the saved simulation results for future reference. Especially if you have simulations that you checked against experiments, I'd be very happy for you to add your crystal definitions to the central database in the project Github.
 - _Efficient code so that complicated systems can be simulated in 3D:_ Real laser pulses can be messy, and if they weren't so before a nonlinear crystal, there's a good chance they are after (but not always). If things are slow, it's hard to go beyond one dimension on tolerable time scales, and then you miss out on the whole weird world of spatiotemporal couplings. Here you have options for rather fast simulations when there's a symmetry to apply (e.g. cylindrical or along one Cartesian dimension), alongside fully 3D propagation. Runs natively on both GPU and CPU to make use of whatever you have to work with.
 - _A graphical interface that lets you see what you're doing:_ A lot of us think in visual terms. Being able to adjust and scan parameters and immediately see what happens can really make it easier to understand what you're looking at. 
 - _A flexible sequence mode:_ By stringing together elements, not just nonlinear crystals but also spherical or parabolic mirrors, apertures, filters, free space propagation and other elements, simulate how  one interaction affects another. Sequences of events can be scripted and even programmed with loop functions to see how things change over the course of repeated interactions.
 - _Fitting modes:_ Sometimes the data that we measure depends on an interesting way on a parameter, and we'd actually like to go back and figure out what that parameter was from the data. Solving this kind of inverse problem can be tough when the parameter lives inside a partial differential equation, but by simulating the whole thing and doing a fit, you have a chance to do it! The fitting algorithm can be used to narrow down a huge space of variables to come at your best estimation of what was happening in an experiment, or to adjust your experimental system to maximize output at a given frequency.
 - _A Python module for easy postprocessing of the results:_ I hope that you get something interesting out that you want to plot and maybe publish. One of the nicest platforms for making nice plots is Python in my opinion (that's why the documentation is in a Jupyter notebook), so purely out of self interest I tried to make it easy to load the results in Python. The module also has some functions related to typical operations you'd like to do on the data to make it easy for all of us. The program also gives you a Matlab loading script for those who want to use that.
 - _Command line interface for running on Linux/clusters:_ The main application runs on Windows, but the simulation core can be compiled on Linux. I use this to run it on the clusters of the Max Planck Society, and other institutes and universities likely have similar systems. This lets you do a lot more if your personal resources are limited but you want to run simulations on a large grid or cover a lot of different parameters!

---

  ### Installation on a Windows PC
  In order to install and run Lightwave Explorer on Windows, just download the file LightwaveExplorerGTK.zip or LightwaveExplorerGTK.7z (the two files have the same contents, I just upload both because not everyone can open 7z) from this [shared volume on the Max Planck Computing and Data Facility DataShare](https://datashare.mpcdf.mpg.de/s/oJj9eFYDBFmViFP).

  The Python module for working with the results is also in that folder for convenience; I'd recommend putting it somewhere in your Python path if you're going to work with it a lot, otherwise just copy it into your working folder. It's also in this repo if you think of any improvements.

  You can also download the version using the Win32 interface rather than the newer interface, in case there are some bugs in the latter which I haven't found yet.
  
  To use the SYCL version, you might also need to install the [Intel® oneAPI DPC++/C++ Compiler Runtime for Windows](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html).

---

  ### Installation on Linux
  The process to get it running on Linux is a bit more involved at the moment, so this is just for the adventurous. 

  First, you'll need to install some stuff.

  - [GTK 4](https://www.gtk.org/docs/installations/linux/)
  - [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
  - [NVIDIA CUDA (optional)](https://developer.nvidia.com/cuda-downloads)

  You'll also need your basic compilation stuff, which I assume you already installed at some point (e.g. sudo apt install build-essential git g++).

  Once those are installed, cross your fingers and enter this in the terminal with your remaining hand:

  ```
  git clone https://github.com/NickKarpowicz/LightwaveExplorer
  git clone https://github.com/davisking/dlib
  cd LightwaveExplorer
  . /opt/intel/oneapi/setvars.sh
  make
  sudo make install
  ```

  That should have done it. If you don't want to install CUDA (i.e. you don't have an NVIDIA board so why bother), you can replace "make" with "make nocuda" and still use SYCL. If you don't want CUDA or SYCL, you can use "make cpuonly".

  This will copy the application binary to /usr/bin/LightwaveExplorer and the text files that the program uses in /usr/shared/LightwaveExplorer. If you want them somewhere else, edit the makefile before you run "make install". In the end, you should be able to call it from anywhere just typing LightwaveExplorer.

  If anyone knows how to load all of this stuff into a package or Flatpak to obviate the need to download 20 gigabytes of development tools, let me know; Linux isn't really my specialty.

---

  ### How do I know which configuration to run?
  At the bottom of the window, you'll see two pulldown menus marked "Config" - these let you choose whether the simulation runs in CUDA, SYCL, or OpenMP. It should start set to the fastest option for your system, but if you don't have the right drivers/runtimes, it might step down to something else that is present. OpenMP is typically the slowest, but will run on basically any system.

  - If you have an Intel CPU, chances are it will have an integrated GPU, which can be used by SYCL. In my experience, the ones that show up as "Iris Graphics" are actually much faster than running the code on the actual CPU, while the ones named "HD Graphics" are sometimes slower. Just try them both.

  - The SYCL code compiles itself for your specific hardware the first time it runs. So, the first run will be slower than the rest - don't reject it because it that; subsequent runs will be faster!

  - The second pulldown is for offloading work onto other parts of your system. For example, if you are running a big batch of simulations on a CUDA-capable GPU, you can send some of them to the CPU to work on. This is what the second menu and following number do: chose the offload target, and the number to send to it. You don't have to use this if you don't want to.

  - Basically, the order in which you should choose the backends is: CUDA if you have something which supports it, SYCL otherwise, and OpenMP if those options aren't available.

  ---
  
  ### Compilation in Visual Studio
  LWE was developed in Visual Studio. If you clone this repo, if you have installed the latest CUDA development kit from NVIDIA, oneAPI from Intel (including the Math Kernel Library and the DPC++ compiler), clone also the dlib repo side-by-side with LWE and it should compile directly.

  ---

  ### Compilation on clusters
  
  A script is provided to compile the CUDA command line version on Linux. This is made specifically to work on the clusters of the MPCDF but will likely work with small modifications on other distributions depending on the local environment. The CUDA development kit and Intel OneAPI should be available in advance. With these prerequisites, entering the following should work:
  ```
wget https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/compileCommandLineLWEfromRepos.sh

chmod +x compileCommandLineLWEfromRepos.sh

./compileCommandLineLWEfromRepos.sh 
 ```
 
 If you have the GUI version installed locally, you can set up your calculation and then generate a SLURM script to run on the cluster (it will tell you what to do).

 ---

  ### Libraries used
 For some of the mathematical operations, I use the following libraries. Thanks to the original authors for making their work available! In order to compile Lightwave Explorer, you'll need them. They are all freely available, but of course have their own licenses .etc.
  - [NVIDIA CUDA](https://developer.nvidia.com/cuda-toolkit): This provides the basic CUDA runtime, compiler, and cuFFT (which does the Fourier transforms here).
  - [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html), specifically the [Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.cw3ci4): This is used for performing fast fourier transforms when running in CPU mode. The DPC++ compiler allows the program to run on both CPUs and a wider range of GPUs, including the integrated ones on Intel chips. I found that on my rather old laptop, SYCL on the GPU is several times faster than running on CPU, so it's useful even for systems without dedicated GPUs.
  - [Dlib](http://dlib.net/): This library is the basis of the optimization routines. I make use of the global optimization functions for the fitting/optimization modes. The library is [available on Github](https://github.com/davisking/dlib), and their excellent documentation and further information is on the [main project website](http://dlib.net/).
  
  ---

  ### Programming note
  Although CUDA was the initial platform and what I use (and test) most extensively, I've added two additional languages for those who don't have an Nvidia graphics card. One is in c++, with multithreading done with OpenMP. The other (just added) language is SYCL. This also allows the simulation to run on the CPU and should allow it to run on Intel's graphics cards, although I have not tested this yet. The same language should be able to run on AMD cards, but their support is not there yet. The code is written in a "trilingual" way - a single core code file is compiled (after some includes and preprocessor definitions) by the three different compilers, Nvidia nvcc, Microsoft Visual Studio, and Intel dpc++. 
  
  Because of this, the different architectures are using the same algorithm, aside from small differences in their floating point math and intrinsic functions. So when I make changes or additions, there will never be any platform gaining over the other (again, reproducibility by anyone is part of the goals here).
