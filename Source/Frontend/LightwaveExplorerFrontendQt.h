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
#include <QSettings>
#include <QPalette>
#include <QColor>
#include <QStyleFactory>
#include <QStyle>
#include <QStyleHints>
#include <QStack>
#include <QSyntaxHighlighter>
#include <QRegularExpression>
#include <QTextCharFormat>
#include <cairo.h>
#include "../LightwaveExplorerUtilities.h"
std::mutex GTKmutex;
#include "LightwaveExplorerPlots.h"
#include "../Devices/LightwaveExplorerCoreCPU.h"
#include "../Devices/LightwaveExplorerCoreCounter.h"
#include "../Devices/LightwaveExplorerCoreCPUFP32.h"
#include "LWEVisualizationsCPU.hpp"
class LWEGui;
using CairoFunction = std::function<void(cairo_t*,int,int,LWEGui&)>;

//conditional includes and definitions
//Apple libraries
#ifdef __APPLE__
#include <mach-o/dyld.h>
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

std::string checkLibraryAvailability(simulationBatch& theSim);
void setInterfaceValuesToActiveValues();
int insertAfterCharacter(std::string& s, char target, std::string appended);
int insertAfterCharacterExcept(std::string& s, char target, std::string appended, std::string exclude);
int formatSequence(std::string& s);
void readDefaultValues(simulationBatch& sim, crystalDatabase& db);
void loadFromPath(const std::string path);
void saveThread();
bool updateDisplay();
void drawBeamImage(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawTimeImage1(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawField1Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawField2Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawSpectrum1Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawSpectrum2Plot(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawTimeImage2(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawFourierImage1(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawFourierImage2(cairo_t* cr, int width, int height, LWEGui& theGui);
void drawBeamThread(LWEGui& theGui, VisualizationConfig config);
class simulationRun {
public:
    std::function<unsigned long(simulationParameterSet*)> normalFunction;
    std::function<unsigned long(simulationParameterSet*)> sequenceFunction;
    std::function<unsigned long(simulationParameterSet*)> fittingFunction;
    int assignedGPU{};
    bool forceCPU = false;
    bool useOpenMP = false;
    simulationRun(int pulldownSelection, bool use64bitFloatingPoint, simulationBatch& theSim){
        sequenceFunction = use64bitFloatingPoint ? 
            &solveNonlinearWaveEquationSequenceCPU : &solveNonlinearWaveEquationSequenceCPUFP32;
        normalFunction = use64bitFloatingPoint ? 
            &solveNonlinearWaveEquationCPU : &solveNonlinearWaveEquationCPUFP32;
        fittingFunction = use64bitFloatingPoint ?
            &runDlibFittingCPU : &runDlibFittingCPUFP32;
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
                fittingFunction = &runDlibFitting;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceFP32;
                normalFunction = &solveNonlinearWaveEquationFP32;
                fittingFunction = &runDlibFittingFP32;
            }

            assignedGPU = pulldownSelection;
        }
        #endif
        #ifndef NOSYCL
        if (pulldownSelection == theSim.base().cudaGPUCount && theSim.base().SYCLavailable) {
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
                fittingFunction = &runDlibFittingSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
                fittingFunction = &runDlibFittingSYCLFP32;
            }

        }
        else if (pulldownSelection == theSim.base().cudaGPUCount + 1 && SYCLitems > 1) {
            forceCPU = 1;
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
                fittingFunction = &runDlibFittingSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
                fittingFunction = &runDlibFittingSYCLFP32;
            }
        }
        else if (pulldownSelection == theSim.base().cudaGPUCount + 2 && SYCLitems > 1) {
            assignedGPU = 1;
            if (use64bitFloatingPoint) {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCL;
                normalFunction = &solveNonlinearWaveEquationSYCL;
                fittingFunction = &runDlibFittingSYCL;
            }
            else {
                sequenceFunction = &solveNonlinearWaveEquationSequenceSYCLFP32;
                normalFunction = &solveNonlinearWaveEquationSYCLFP32;
                fittingFunction = &runDlibFittingSYCL;
            }
        }
        else if (pulldownSelection == (theSim.base().cudaGPUCount + SYCLitems + 1)){
            useOpenMP = true;
        }
        #endif
    };
};

QFont getEmojiFont() {
    QFont emojiFont;
#if defined(Q_OS_WIN)
    emojiFont.setFamily("Segoe UI Emoji");
#elif defined(Q_OS_MAC)
    emojiFont.setFamily("Apple Color Emoji");
#elif defined(Q_OS_LINUX)
    emojiFont.setFamily("Noto Color Emoji");
#else
    // Fallback option
    emojiFont.setFamily("Segoe UI Emoji");
#endif
    return emojiFont;
}

void mainSimThread(LWEGui& theGui, simulationRun theRun, simulationRun theOffloadRun);
void fittingThread(LWEGui& theGui,  simulationRun theRun);
void createRunFile(LWEGui& theGui);



