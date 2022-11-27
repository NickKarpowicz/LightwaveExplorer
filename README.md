## <p style="text-align: center;">Lightwave Explorer</p>

<p style="text-align: center;">Nick Karpowicz</p>
<p style="text-align: center;">Max Planck Institute of Quantum optics</p>

<p style="text-align: center;"><img src="/Documentation/Images/LWEicon.png" width="200" height="200"></p>

<p style="text-align: center;">(icon made by Stable Diffusion)</p>


### What and why

Lightwave explorer is an open source nonlinear optics simulator, intended to be fast, visual, and flexible for students and researchers to play with ultrashort laser pulses and nonlinear optics without having to buy a laser first.

<p style="text-align: center;"><img src="/Documentation/Images/Screenshot1122.png" width="941" height="495"></p>

The simulation was written CUDA in order to run quickly on modern graphics cards. I've subsequently generalized it so that it can be run in two other ways: SYCL on CPUs and Intel GPUs, and using OpenMP to run on CPUs. Accordingly, I hope that the results are fast enough that even complicated system can be simulated within a human attention span.

#### Main goals:
 - _Easily extensible database of materials:_ Eveything the program knows about nonlinear materials comes from a human-readable text file giving the appropriate coefficients and tensors. If you want to use a new material, or you've done a measurement in a new range where typical extrapolations from older data isn't relevant, it's easy to add and correct. There are places for references for the key parameters, and these references are stored in the saved simulation results for future reference. Especially if you have simulations that you checked against experiments, I'd be very happy for you to add your crystal definitions to the central database in the project Github.
 - _Efficient code so that complicated systems can be simulated in 3D:_ Real laser pulses can be messy, and if they weren't so before a nonlinear crystal, there's a good chance they are after (but not always). If things are slow, it's hard to go beyond one dimension on tolerable time scales, and then you miss out on the whole weird world of spatiotemporal couplings. Here you have options for rather fast simulations when there's a symmetry to apply (e.g. cylindrical or along one Cartesian dimension), alongside fully 3D propagation. Runs natively on both GPU and CPU to make use of whatever you have to work with.
 - _A graphical interface that lets you see what you're doing:_ A lot of us think in visual terms. Being able to adjust and scan parameters and immediately see what happens can really make it easier to understand what you're looking at. 
 - _A flexible sequence mode:_ By stringing together elements, not just nonlinear crystals but also spherical or parabolic mirrors, apertures, filters, free space propagation and other elements, simulate how  one interaction affects another. Sequences of events can be scripted and even programmed with loop functions to see how things change over the course of repeated interactions.
 - _Fitting modes:_ Sometimes the data that we measure depends on an interesting way on a parameter, and we'd actually like to go back and figure out what that parameter was from the data. Solving this kind of inverse problem can be tough when the parameter lives inside a partial differential equation, but by simulating the whole thing and doing a fit, you have a chance to do it! The fitting algorithm can be used to narrow down a huge space of variables to come at your best estimation of what was happening in an experiment, or to adjust your experimental system to maximize output at a given frequency.
 - _A Python module for easy postprocessing of the results:_ I hope that you get something interesting out that you want to plot and maybe publish. One of the nicest platforms for making nice plots is Python in my opinion (that's why the documentation is in a Jupyter notebook), so purely out of self interest I tried to make it easy to load the results in Python. The module also has some functions related to typical operations you'd like to do on the data to make it easy for all of us. The program also gives you a Matlab loading script for those who want to use that.
 - _Command line interface for running on Linux/clusters:_ The main application runs on Windows, but the simulation core can be compiled on Linux. I use this to run it on the clusters of the Max Planck Society, and other institutes and universities likely have similar systems. This lets you do a lot more if your personal resources are limited but you want to run simulations on a large grid or cover a lot of different parameters!

  ### Installation (Windows)
  In order to install and run Lightwave Explorer on Windows, just download the file LightwaveExplorerWin64.zip from this [shared volume on the Max Planck Computing and Data Facility DataShare](https://datashare.mpcdf.mpg.de/s/oJj9eFYDBFmViFP).
  
  To use the SYCL version, you will also need to install the [Intel® oneAPI DPC++/C++ Compiler Runtime for Windows](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html).
  
  ### Compilation
  LWE was developed in Visual Studio. If you clone this repo, if you have installed the latest CUDA development kit from NVIDIA, oneAPI from Intel (including the Math Kernel Library and the DPC++ compiler), clone also the dlib repo side-by-side with LWE and it should compile directly.
  
  A script is provided to compile the CUDA command line version on Linux. This is made specifically to work on the clusters of the MPCDF but will likely work with small modifications on other distributions depending on the local environment. The CUDA development kit and Intel OneAPI should be available in advance. With these prerequisites, entering the following should work:
  ```
          wget https://raw.githubusercontent.com/NickKarpowicz/LightwaveExplorer/master/compileCommandLineLWEfromRepos.sh

          chmod +x compileCommandLineLWEfromRepos.sh

          ./compileCommandLineLWEfromRepos.sh 
 ```
 
 I have also compiled the command line version on an Intel Mac using a similar procedure. If you wish to do the same, I can send you the configuration files for VS Code, which is likely the easiest way.
 
 The Windows binary does work (minus aesthetics) using Wine, as tested under both Linux and on an Intel MacBook, but as far as I can tell only using the c++ (CPU) propgation.
 
  ### Libraries used
 For some of the mathematical operations, I use the following libraries. Thanks to the original authors for making their work available! In order to compile Lightwave Explorer, you'll need them. They are all freely available, but of course have their own licenses .etc.
  - [NVIDIA CUDA](https://developer.nvidia.com/cuda-toolkit): This provides the basic CUDA runtime, compiler, and cuFFT (which does the Fourier transforms here).
  - [Intel OneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html), specifically the [Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.cw3ci4): This is used for performing fast fourier transforms when running in CPU mode. The DPC++ compiler allows the program to run on both CPUs and a wider range of GPUs, including the integrated ones on Intel chips. I found that on my rather old laptop, SYCL on the GPU is several times faster than running on CPU, so it's useful even for systems without dedicated GPUs.
  - [Dlib](http://dlib.net/): This library is the basis of the optimization routines. I make use of the global optimization functions for the fitting/optimization modes. The library is [available on Github](https://github.com/davisking/dlib), and their excellent documentation and further information is on the [main project website](http://dlib.net/).
  
  ### Programming note
  Although CUDA was the initial platform and what I use (and test) most extensively, I've added two additional languages for those who don't have an Nvidia graphics card. One is in c++, with multithreading done with OpenMP. The other (just added) language is SYCL. This also allows the simulation to run on the CPU and should allow it to run on Intel's graphics cards, although I have not tested this yet. The same language should be able to run on AMD cards, but their support is not there yet. The code is written in a "trilingual" way - a single core code file is compiled (after some includes and preprocessor definitions) by the three different compilers, Nvidia nvcc, Microsoft Visual Studio, and Intel dpc++. 
  
  Because of this, the different architectures are using the same algorithm, aside from small differences in their floating point math and intrinsic functions. So when I make changes or additions, there will never be any platform gaining over the other (again, reproducibility by anyone is part of the goals here).
