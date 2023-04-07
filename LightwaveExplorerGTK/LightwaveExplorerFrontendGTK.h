#pragma once
#include <sstream>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>
#include "../LightwaveExplorerDevices/LightwaveExplorerUtilities.h"
#include "LightwaveExplorerGraphicalClasses.h"
#undef __noinline__
#include <gtk/gtk.h>
#if defined __linux__
#include<fmt/format.h>
#define Sformat fmt::format
#define Svformat fmt::vformat
#define Smake_format_args fmt::make_format_args
#else
#include <format>
#define Sformat std::format
#define Svformat std::vformat
#define Smake_format_args std::make_format_args
#endif

#include <thread>
#include <chrono>
#include <locale>
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCPU.h"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCounter.h"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreFP32.cuh"
#include "../LightwaveExplorerDevices/LightwaveExplorerCoreCPUFP32.h"

//conditional includes and definitions
//Apple libraries
#ifdef __APPLE__
#include <mach-o/dyld.h>
#import<Cocoa/Cocoa.h>
#endif

//if not running on a cpu-only build load CUDA and SYCL code
#ifndef CPUONLY
//if using CUDA include CUDA libraries and the 
//CUDA version of the simulation
#ifndef NOCUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <nvml.h>
#include "../LightwaveExplorerCore.cuh"
#endif
//if on windows, include the header for the Windows version
//of the SYCL code and the Windows headers needed to check
//if the Intel DPC++ runtime is present
#ifdef _WIN32
#include "../LightwaveExplorerDevices/LightwaveExplorerSYCL.h"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
//on Linux, load the Linux versions of the SYCL code
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerDPCPPlibFP32.h"
#endif
#endif


void openFileDialogCallback(GtkWidget* widget, gpointer pathTarget);
void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget);
void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawField1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawProgress(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void checkLibraryAvailability();
void setInterfaceValuesToActiveValues();
int insertAfterCharacter(std::string& s, char target, std::string appended);
int insertAfterCharacterExcept(std::string& s, char target, std::string appended, std::string exclude);
int formatSequence(std::string& s);
void mainSimThread(int pulldownSelection, int pulldownSelection2, bool use64bitFloatingPoint);
void launchRunThread();
void independentPlotQueue();
void loadCallback(GtkWidget* widget, gpointer pathTarget);
void launchFitThread();
void fittingThread(int pulldownSelection, bool use64bitFloatingPoint);
void stopButtonCallback();
void svgCallback();
void createRunFile();
static void buttonAddSameCrystal();
static void buttonAddDefault();
static void buttonAddRotation();
static void buttonAddPulse();
static void buttonAddMirror();
static void buttonAddFilter();
static void buttonAddLinear();
static void buttonAddAperture();
static void buttonAddFarFieldAperture();
static void buttonAddForLoop();
bool updateDisplay();
void destroyMainWindowCallback();

