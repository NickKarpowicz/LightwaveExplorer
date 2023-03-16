#define RUNONSYCL
#define LWEFLOATINGPOINT 32
#define isnan(x) std::isnan(x)
#include "LightwaveExplorerDevices/LightwaveExplorerUtilities.h"
#include "LightwaveExplorerDPCPPlib.h"
#include "LightwaveExplorerCore.cu"
