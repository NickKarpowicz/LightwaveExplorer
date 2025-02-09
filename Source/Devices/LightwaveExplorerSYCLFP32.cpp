// dllmain.cpp : Defines the entry point for the DLL application.
#define RUNONSYCL
#define LWEFLOATINGPOINT 32
#include "../LightwaveExplorerUtilities.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"
#include "LightwaveExplorerSYCL.h"
#include "LightwaveExplorerCore.cu"
