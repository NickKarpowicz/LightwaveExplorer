#pragma once
#include <sstream>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <chrono>
#include <locale>
#include <functional>
#undef __noinline__
#include <QApplication>
#include <QWidget>
#include <QObject>
#include <QTextEdit>
#include <QPushButton>
#include <QMainWindow>
#include <QLineEdit>
#include <QComboBox>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QProgressBar>
#include <QSlider>
#include <QLabel>
#include <QImage>
#include <QPainter>
#include <QThread>
#include <QFileDialog>
#include <cairo.h>
#include "../LightwaveExplorerUtilities.h"
std::mutex GTKmutex;
#include "LightwaveExplorerPlots.h"
#include "../Devices/LightwaveExplorerCoreCPU.h"
#include "../Devices/LightwaveExplorerCoreCounter.h"
#include "../Devices/LightwaveExplorerCoreCPUFP32.h"
class LWEGui;
using CairoFunction = std::function<void(cairo_t*,int,int,LWEGui&)>;

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
#include "../Devices/LightwaveExplorerCoreFP32.cuh"
#include "../LightwaveExplorerCore.cuh"
#endif
//if on windows, include the header for the Windows version
//of the SYCL code and the Windows headers needed to check
//if the Intel DPC++ runtime is present
#ifdef _WIN32
#include "../Devices/LightwaveExplorerSYCL.h"
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
//on Linux, load the Linux versions of the SYCL code
#ifndef NOSYCL
#include "../Devices/LightwaveExplorerSYCLLinux.h"
#include "../Devices/LightwaveExplorerSYCLLinuxFP32.h"
#endif
#endif
#endif

#if defined __linux__
#include<fmt/format.h>
#define Sformat fmt::format
#define Svformat fmt::vformat
#define Smake_format_args fmt::make_format_args
#elif defined __APPLE__
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

// void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawField1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
// void drawProgress(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
bool sliderResponseToArrows(void* widget, unsigned int keyValue);
std::string checkLibraryAvailability(simulationBatch& theSim);
void setInterfaceValuesToActiveValues();
int insertAfterCharacter(std::string& s, char target, std::string appended);
int insertAfterCharacterExcept(std::string& s, char target, std::string appended, std::string exclude);
int formatSequence(std::string& s);
void readDefaultValues(simulationBatch& sim, crystalDatabase& db);
void launchRunThread();
void independentPlotQueue();
void loadCallback();
void savePathCallback();
void waveform1PathCallback();
void waveform2PathCallback();
void fittingPathCallback();
void loadDatabaseCallback();
void saveRunFileCallback();
void loadFromPath(const std::string path);
void saveThread();
void launchFitThread();
void fittingThread(int pulldownSelection, bool use64bitFloatingPoint);
void stopButtonCallback();
void svgCallback();
void dataPanelCollapseCallback();
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

void drawTimeImage1(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawField1Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawField2Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawSpectrum1Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawSpectrum2Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawTimeImage2(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawFourierImage1(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawFourierImage2(cairo_t* cr, int width, int height, LWEGui& theGui);

class simulationRun {
public:
    std::function<unsigned long(simulationParameterSet*)> normalFunction;
    std::function<unsigned long(simulationParameterSet*)> sequenceFunction;
    int assignedGPU{};
    bool forceCPU = false;
    bool useOpenMP = false;
    simulationRun(int pulldownSelection, bool use64bitFloatingPoint, simulationBatch& theSim){
        sequenceFunction = use64bitFloatingPoint ? 
            &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
        normalFunction = use64bitFloatingPoint ? 
            &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
        int assignedGPU = 0;
        bool forceCPU = false;
        bool useOpenMP = false;
    #ifdef CPUONLY
        useOpenMP = true;
    #endif
        [[maybe_unused]]int SYCLitems = 0;
        #if !defined(CPUONLY)
        if (theSim.base().syclGPUCount == 0) {
            SYCLitems = (int)theSim.base().SYCLavailable;
        }
        else {
            SYCLitems = 3;
        }
        #endif
        #if !defined(CPUONLY) && !defined(NOCUDA)
        if (pulldownSelection < theSim.base().cudaGPUCount) {
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequence;
                normalFunction = &solveNonlinearWaveEquation;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceFP32;
                normalFunction = &solveNonlinearWaveEquationFP32;
            }

            assignedGPU = pulldownSelection;
        }
        #endif
        #ifndef NOSYCL
        if (pulldownSelection == theSim.base().cudaGPUCount && theSim.base().SYCLavailable) {
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
            }

        }
        else if (pulldownSelection == theSim.base().cudaGPUCount + 1 && SYCLitems > 1) {
            forceCPU = 1;
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
            }
        }
        else if (pulldownSelection == theSim.base().cudaGPUCount + 2 && SYCLitems > 1) {
            assignedGPU = 1;
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
            }
        }
        else if (pulldownSelection == (theSim.base().cudaGPUCount + SYCLitems + 1)){
            useOpenMP = true;
        }
        #endif
    };
};

void mainSimThread(LWEGui& theGui, simulationRun theRun, simulationRun theOffloadRun);

