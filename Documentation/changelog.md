### Changes in 2026.1
New features:
- Support for complex input beams, on either a Hermite-Gauss or Laguerre-Gauss basism with frequency-dependent beamwaist, rotation, 3D offset and angle.
- offsetAperture() sequence function for including a non-centered iris
Bug fixes:
- Will load the last crystal in the database even if there isn't a newline at the end
- Circular polarization settings are explained in the tooltip, and values of outside of -1..1 are no longer passed to the simulation with odd results
Programming things:
- The CPU FFT library was switched to PocketFFT. This leads to much better performance on the AMD systems I tested; slightly worse on Intel. However, it also makes maintenance and compiling much easier.

### Changes in 2025.6:
Bug fixes:
- Fix a bug where setting a non-colinear angle would cause pulse injection to fail in FDTD simulations
- Fix a bug where a failed load of a custom FDTD grid wouldn't cause the simulation to throw an exception
- The Windows version is now built by a Github action for better reproducibility
- The Intel Mac version is deprecated since Apple stopped updates for them. If someone still wants a build, tell me. That means that to use it on Mac, you should follow the compilation directions.

### Changes in 2025.5:
Bug fixes + some crystals:
- Fixes bug where files could not be overwritten after loading
- Fixes behavior of default values within for-loop sequences
- Add Lorentzian model of GaSe, Gaussian model of YAG

### Changes in 2025.4:
Bug fix release:
- Corrects the carrier-envelope phase unit on the GUI
- Fixes batch size of 1 refusing to run
- Fixes a crash when rendering large 2D beam views
- Fixes a crash in FDTD mode for large grids and short time lengths
- Fixes GUI lock-up when rendering large beam views

### Changes in 2025.3:
- Add polarizer() sequence function
- Use OneMath 0.7 in SYCL version

### Changes in 2025.2:
- Added new beam visualization mode, where you can make color images of the light field that results from the simulation
- Add new controls for the appearance of plots, which you can access by clicking the ↔️ button on the bottom of the window.
- Fixed a bug where a FROG file loaded for pulse 2 wouldn't be used.
- Fixed a bug where very large grids whose dimensions were divisible by 4 but not by 8 would look weird.
- The LightwaveExplorer python module is now on pip! i.e. you can run "pip install LightwaveExplorer" to get the latest version.

### Changes in 2025.1:
- Fixed a bug introduced in 2025.0 where on CUDA the plasma would not affect fields in the y-direction
- Fixed a bug where cancelling a batch job caused the subsequent batch to cancel itself
- Axes and scales are now shown on the space and momentum images, and they're also saved when an SVG is produced.
- Plot font size scaling works better with high-DPI displays
- Refactored device structures to more efficiently and consistently allocate memory
- NOTE: to read the files produced with 2025.0 and later in python, you need to update the LightwaveExplorer module (now available on pip!), many thanks to Florian Lindinger for the improved version that works with new and old files!

### Changes in 2025.0:

 - FDTD mode is more accessible: besides the options for crystals in the Propagation pulldown menu, you can also load a custom grid containing different nonlinear materials at different positions using the fdtd() function
 - FDTD reflection mode for nonlinear reflection from a surface, using the fdtdReflection() sequence command.
 - Improved convergence for the nonlinear absorption/plasma response. Note that the changes I made here will re-scale the values of the nonlinear absorption parameter (now they are consistent with what was in the paper). The conversion factor for older values is $\left(\frac{1}{2\pi\epsilon_0}\right)^{\frac{1}{N_p-1}}$ where $N_p$ is the number of photons in the multiphoton transition (divide the old value by this number).
