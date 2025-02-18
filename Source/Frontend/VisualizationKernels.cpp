#include "../LightwaveExplorerTrilingual.h"
#include "../DeviceFunctions.hpp"


using namespace deviceFunctions;
namespace kernelNamespace{

};

namespace {
    void renderBeamPower(ActiveDevice& d, const VisualizationConfig& config){
        
    }
}
unsigned long renderVisualizationX(ActiveDevice& d, const VisualizationConfig& config){
    std::lock_guard<std::mutex> lock(d.visualization->memoryMutex);
    switch(config.type){
        case VisualizationType::beamPower:
            break;
    }
}