#ifdef LWEFLOATINGPOINT
#undef LWEFLOATINGPOINT
#endif
#define LWEFLOATINGPOINT 32
#include "../LightwaveExplorerTrilingual.h"
#include "../DeviceFunctions.hpp"
#include "../LightwaveExplorerInterfaceClasses.hpp"

namespace deviceFunctions{

	deviceFunction static inline int64_t expandCylindricalBeamAtRadius(
		const float r, 
		const int64_t Nspace) {

        //is remainder of r closer to 1/4 or 3/4
        bool negativeSide = fmodf(r,1.0f)-0.5 > 0.0f;
        int64_t index = static_cast<int64_t>((negativeSide ? 
            static_cast<float>(Nspace)/2.0f - r + 0.25f:
            static_cast<float>(Nspace)/2.0f + r + 0.25f));
		return index;
	}

    deviceFunction static inline float findCenter(const int64_t Nspace, const int64_t stride, const float* data){
        float beamTotal{};
        float beamCenter{};
            for (int64_t j{}; j < Nspace; ++j) {
                float a = (data[j * stride]);
                beamTotal += a;
                beamCenter += static_cast<float>(j) * a;
            }
            if (beamTotal > 0.0f) {
                beamCenter /= beamTotal;
            }
            return beamCenter;
    }

    deviceFunction static inline float findCenter(const int64_t Nspace, const int64_t stride, const deviceComplex* data){
        float beamTotal{};
        float beamCenter{};
            for (int64_t j{}; j < Nspace; ++j) {
                float a = cModulusSquared(data[j * stride]);
                beamTotal += a;
                beamCenter += static_cast<float>(j) * a;
            }
            if (beamTotal > 0.0f) {
                beamCenter /= beamTotal;
            }
            return beamCenter;
    }
}
using namespace deviceFunctions;

namespace kernelNamespace{

// struct findBeamCentersKernel{

// };
    struct beamPowerRenderKernel2D{
        const float* Et;
        float* red;
        float* green;
        float* blue;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const int64_t x_field = (i / config.Nt) * (config.width/config.Nx);
            atomicAdd(&red[x_field], Et[i]*Et[i]);
            atomicAdd(&green[x_field], Et[i + config.Ngrid] * Et[i + config.Ngrid]);
            atomicAdd(&blue[x_field], Et[i]*Et[i] + Et[i + config.Ngrid] * Et[i + config.Ngrid]);
        }
    };

    struct falseColorKernel2D{
        const deviceComplex* Ef;
        float* red;
        float* green;
        float* blue;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            float f = (i+1) * config.df;
            
        }
    }

    struct floatLinetoImageKernel{
        const float* red;
        const float* green;
        const float* blue;
        uint8_t* image;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float mid = static_cast<float>(config.Nx)/2.0f;
            const float x = static_cast<float>(i % config.Nx)-mid;
            const float y = static_cast<float>(i / config.Nx)-mid;
            const float r = abs(hypotf(x,y));
            const int64_t idx = std::clamp(expandCylindricalBeamAtRadius(r, config.Nx),int64_t{},config.Nx);
            //stride in image is 4, r g b a

            image[4*i] = static_cast<uint8_t>(
                std::clamp(
                    1e-18f * blue[idx],
                    0.f, 255.f));;
            image[4*i + 1] = static_cast<uint8_t>(
                std::clamp(
                    1e-18f * green[idx],
                    0.f, 255.f));;
            image[4*i + 2] = static_cast<uint8_t>(
                std::clamp(
                    1e-18f * red[idx],
                    0.f, 255.f));;
        }
    };

    struct floatLinetoImageKernelCartesian{
        const float* red;
        const float* green;
        const float* blue;
        uint8_t* image;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float mid = static_cast<float>(config.Nx)/2.0f;
            const float center = findCenter(config.Nx, 1, blue);
            const float x = static_cast<float>(i % config.Nx)-center;
            const float y = static_cast<float>(i / config.Nx)-mid;
            const float r = abs(hypotf(x,y));
            const int64_t idx1 = static_cast<int64_t>(roundf(r+center));
            const int64_t idx2 = static_cast<int64_t>(roundf(r-center));
            bool idx1Valid = idx1 < config.Nx && idx1>=0;
            bool idx2Valid = idx2 < config.Nx && idx2>=0;
            if(idx1Valid && idx2Valid){
                image[4*i] = static_cast<uint8_t>(
                    std::clamp(
                        0.5e-18f * (blue[idx1]+ blue[idx2]),
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        0.5e-18f * (green[idx1]+green[idx2]),
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        0.5e-18f * (red[idx1] + red[idx2]),
                        0.f, 255.f));;
            }
            else if(idx1Valid){
                image[4*i] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (blue[idx1]),
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (green[idx1]),
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (red[idx1]),
                        0.f, 255.f));;
            }
            else if(idx2Valid){
                image[4*i] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (blue[idx2]),
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (green[idx2]),
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        1e-18f * (red[idx2]),
                        0.f, 255.f));;
            }
            else{
                image[4*i] = {};
                image[4*i+1] = {};
                image[4*i+2] = {};
            }

        }
    };
};
using namespace kernelNamespace;

namespace {
    static void renderFalseColor(ActiveDevice& d, const VisualizationConfig config){
        d.fft(d.visualization->gridTimeSpace.device_ptr(),
            d.visualization->gridFrequency.device_ptr(),
            deviceFFT::D2Z_1D);
    }
    static void renderBeamPower(ActiveDevice& d, const VisualizationConfig config){
        d.visualization->gridFrequency.initialize_to_zero();

        switch(config.mode){
            case CoordinateMode::cartesian2D:
                d.deviceLaunch(config.Ngrid/64, 64, 
                    beamPowerRenderKernel2D{
                            d.visualization->gridTimeSpace.device_ptr(), 
                            d.visualization->red.device_ptr(),
                            d.visualization->green.device_ptr(),
                            d.visualization->blue.device_ptr(),
                            config});
                    d.deviceLaunch((config.width * config.height) / 64, 64, 
                        floatLinetoImageKernelCartesian{
                            d.visualization->red.device_ptr(),
                            d.visualization->green.device_ptr(),
                            d.visualization->blue.device_ptr(),
                            d.visualization->imageGPU.device_ptr(),
                            config});
                break;
            case CoordinateMode::radialSymmetry:
                d.deviceLaunch(config.Ngrid/64, 64, 
                    beamPowerRenderKernel2D{
                            d.visualization->gridTimeSpace.device_ptr(), 
                            d.visualization->red.device_ptr(),
                            d.visualization->green.device_ptr(),
                            d.visualization->blue.device_ptr(),
                            config});
                    d.deviceLaunch((config.width * config.height) / 64, 64, 
                        floatLinetoImageKernel{
                            d.visualization->red.device_ptr(),
                            d.visualization->green.device_ptr(),
                            d.visualization->blue.device_ptr(),
                            d.visualization->imageGPU.device_ptr(),
                            config});
                break;
            default:
                break;
        }


        d.visualization->syncImages();
    }
}
unsigned long renderVisualizationX(VisualizationConfig config){
    ActiveDevice d(config.width, config.height, config.sCPU);
    std::lock_guard<std::mutex> lock(d.visualization->memoryMutex);
    d.visualization->fetchSim(config.sCPU, config.simIndex);
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