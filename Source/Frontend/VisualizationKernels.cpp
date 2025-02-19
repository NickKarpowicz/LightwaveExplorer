#ifdef LWEFLOATINGPOINT
#undef LWEFLOATINGPOINT
#endif
#define LWEFLOATINGPOINT 32
#include "../LightwaveExplorerTrilingual.h"
#include "../DeviceFunctions.hpp"
#include "../LightwaveExplorerInterfaceClasses.hpp"

using namespace deviceFunctions;
namespace kernelNamespace{
struct beamPowerRenderKernel{
    const float* Et;
    float* workspace;
    const VisualizationConfig config;
    deviceFunction void operator()(const int64_t i) const {
        //i has range Ngrid * 2
        //each thread deals with one point in field
        //get the two field polarization components, square and add
        //atomic add this to correct pixel.
        //const int64_t t_field = i % config.Nt;
        const int64_t x_field = i / config.Nt;
        //const int64_t y_field = 0;
        const int64_t x_image = x_field;
        //const int64_t y_image = i / (config.Nx * config.Nt);
        atomicAdd(&workspace[x_image], Et[i]*Et[i]);
    }
};

struct floatLinetoImageKernel{
    const float* line;
    uint8_t* image;
    const VisualizationConfig config;
    deviceFunction void operator()(const int64_t i) const {
        const int64_t x = i % config.Nx;
        const int64_t y = i / config.Nx;
        float mid = static_cast<float>(config.Nx)/2.0f;
        const int64_t r = static_cast<int64_t>(mid + deviceFPLib::hypot(static_cast<float>(x) - mid, static_cast<float>(y) - mid)); 
        //stride in image is 4, r g b a
        uint8_t pt = (r < config.Nx) ? static_cast<uint8_t>(std::clamp(1e-18f * line[r], 0.f, 255.f)) : 0U;
        image[4*i] = pt;
        image[4*i + 1] = 0U;
        image[4*i + 2] = pt;
    }
};
};
using namespace kernelNamespace;

namespace {
    void renderBeamPower(ActiveDevice& d, const VisualizationConfig config){
        d.visualization->gridFrequency.initialize_to_zero();
        d.deviceLaunch(config.Ngrid/64, 64, 
        beamPowerRenderKernel{
                d.visualization->gridTimeSpace.device_ptr() + 2*config.Ngrid*config.simIndex, 
                reinterpret_cast<float*>(d.visualization->gridFrequency.device_ptr()),
                 config});
        d.deviceLaunch((config.width * config.height) / 64, 64, 
        floatLinetoImageKernel{
            reinterpret_cast<const float*>(d.visualization->gridFrequency.device_ptr()), 
            d.visualization->imageGPU.device_ptr(),
            config});
        d.visualization->syncImages();
        std::cout << "out: " << d.visualization->gridTimeSpace.device_ptr() + 2 * config.Ngrid * config.simIndex << std::endl;
    }
}
unsigned long renderVisualizationX(VisualizationConfig config){
    ActiveDevice d(config.width, config.height, config.sCPU);
    std::lock_guard<std::mutex> lock(d.visualization->memoryMutex);
    switch(config.type){
        case VisualizationType::beamPower:
            renderBeamPower(d,config);
            break;
        case VisualizationType::beamFalseColor:
            break;
        default:
            return 1;
    }
    if(config.result_pixels != nullptr){
        config.result_pixels->assign(d.visualization->image.begin(), d.visualization->image.end());
    }
    return 0;
}