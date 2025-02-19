#include "../LightwaveExplorerTrilingual.h"
#include "../DeviceFunctions.hpp"
#include "../LightwaveExplorerInterfaceClasses.hpp"

using namespace deviceFunctions;
namespace kernelNamespace{
struct beamPowerRenderKernel{
    const float* Et;
    float* workspace;
    const VisualizationConfig config;
    deviceFunction void operator()(const int64_t i){
        //i has range Ngrid * 2
        //each thread deals with one point in field
        //get the two field polarization components, square and add
        //atomic add this to correct pixel.
        const int64_t t_field = i % config.Nt;
        const int64_t x_field = i / config.Nt;
        const int64_t y_field = 0;
        const int64_t x_image = x_field;
        //const int64_t y_image = i / (config.Nx * config.Nt);
        atomicAdd(workspace[x_image], Et[i]*Et[i]);
    }
};
};

namespace {
    void renderBeamPower(ActiveDevice& d, const VisualizationConfig& config){
        d.visualization.
    }
}
unsigned long renderVisualizationX(ActiveDevice& d, const VisualizationConfig& config){
    std::lock_guard<std::mutex> lock(d.visualization->memoryMutex);
    switch(config.type){
        case VisualizationType::beamPower:
            renderBeamPower(d,config);
            break;
    }
}