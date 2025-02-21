#ifdef LWEFLOATINGPOINT
#undef LWEFLOATINGPOINT
#endif
#define LWEFLOATINGPOINT 32
#include "../LightwaveExplorerTrilingual.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"

namespace deviceFunctions{
    deviceFunction static inline float deviceNorm(const deviceComplex& x){
        return x.real() * x.real() + x.imag() * x.imag();
    }
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
                float a = deviceNorm(data[j * stride]);
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
            const float A2 = Et[i]*Et[i];
            const float B2 = Et[i + config.Ngrid] * Et[i + config.Ngrid];
            atomicAdd(&red[x_field], A2);
            atomicAdd(&green[x_field], B2);
            atomicAdd(&blue[x_field], A2 + B2);
        }
    };

    struct falseColorPrerenderKernel2D{
        const deviceComplex* Ef;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float f = (i+1) * config.df;

            const float scale = 1e-18 * config.df / (config.Nt);
            const float filter_red = scale * config.red_amplitude * expf(-(f-config.red_f0)*(f-config.red_f0)/(2*config.red_sigma*config.red_sigma));
            const float filter_green = scale * config.green_amplitude * expf(-(f-config.green_f0)*(f-config.green_f0)/(2*config.green_sigma*config.green_sigma));
            const float filter_blue = scale * config.blue_amplitude * expf(-(f-config.blue_f0)*(f-config.blue_f0)/(2*config.blue_sigma*config.blue_sigma));
            
            int64_t dataSize = (config.Nt/2)*(config.Nx);
            float* red = workspace + i;
            float* green = red + dataSize;
            float* blue = green + dataSize;
            float* red_centers = blue + dataSize;
            float* green_centers = red_centers + config.Nt/2;
            float* blue_centers = green_centers + config.Nt/2;

            for(int j{}; j<config.Nx; ++j){
                float Ef2 = deviceNorm(Ef[(i+1) + j*config.Nf]) 
                + deviceNorm(Ef[(i+1) + j*config.Nf + config.Ngrid_complex]);
                int offset = j * config.Nt/2;
                red[offset] = filter_red * Ef2;
                green[offset]= filter_green * Ef2;
                blue[offset] = filter_blue * Ef2;
            }
            red_centers[0] = findCenter(config.Nx, config.Nt/2, red);
            green_centers[0] = findCenter(config.Nx, config.Nt/2, green);
            blue_centers[0] = findCenter(config.Nx, config.Nt/2, blue);
        }
    };

    struct falseColorPrerenderKernelCylindrical{
        const deviceComplex* Ef;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float f = (i+1) * config.df;

            const float scale = 1e-18 * config.df / (config.Nt);
            const float filter_red = scale * config.red_amplitude * expf(-(f-config.red_f0)*(f-config.red_f0)/(2*config.red_sigma*config.red_sigma));
            const float filter_green = scale * config.green_amplitude * expf(-(f-config.green_f0)*(f-config.green_f0)/(2*config.green_sigma*config.green_sigma));
            const float filter_blue = scale * config.blue_amplitude * expf(-(f-config.blue_f0)*(f-config.blue_f0)/(2*config.blue_sigma*config.blue_sigma));
            
            float* red_line = workspace;
            float* green_line = red_line + config.Nx;
            float* blue_line = green_line + config.Nx; 

            for(int j{}; j<config.Nx; ++j){
                float Ef2 = deviceNorm(Ef[(i+1) + j*config.Nf])
                 + deviceNorm(Ef[(i+1) + j*config.Nf + config.Ngrid_complex]);
                atomicAdd(&red_line[j], filter_red * Ef2);
                atomicAdd(&green_line[j], filter_green * Ef2);
                atomicAdd(&blue_line[j], filter_blue * Ef2);
            }
        }
    };

    struct falseColorRender2DCartesianKernel{
        float* red;
        float* green;
        float* blue;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            //i range 0..Nx*Nx*Nt/2
            const float mid = static_cast<float>(config.Nx)/2.0f;
            const int64_t Nf = config.Nt/2; //not +1 because I skip 0
            const int64_t f_ind = i % Nf;
            
            int64_t dataSize = Nf*(config.Nx);
            float* red_data = workspace;
            float* green_data = red_data + dataSize;
            float* blue_data = green_data + dataSize;
            float* red_centers = blue_data + dataSize;
            float* green_centers = red_centers + config.Nt/2;
            float* blue_centers = green_centers + config.Nt/2;
            auto addPoint = [&](float* channel, const float* data, const float* centers, const int x_ind, const float y, const int64_t img_ind){
                    const float center = centers[f_ind];
                    const float x = static_cast<float>(x_ind)-center;
                    const float r = hypotf(x,y);
                    const int64_t idx1 = static_cast<int64_t>((r+center));
                    const int64_t idx2 = static_cast<int64_t>((r-center));
                    const bool idx1Valid = (idx1 < config.Nx) && (idx1>=0);
                    const bool idx2Valid = (idx2 < config.Nx) && (idx2>=0);

                    if(idx1Valid && idx2Valid){
                        atomicAdd(&channel[img_ind], data[f_ind + idx1*(config.Nt/2)] + data[f_ind + idx2*(config.Nt/2)]); 
                    }
                    else if(idx1Valid){
                        atomicAdd(&channel[img_ind], data[f_ind + idx1*(config.Nt/2)]); 
                    }
                    else if(idx2Valid){
                        atomicAdd(&channel[img_ind], data[f_ind + idx2*(config.Nt/2)]); 
                    }
            };
            
            for(int x_ind{}; x_ind<config.Nx; ++x_ind){
                for(int y_ind{}; y_ind<config.Nx; ++y_ind){
                    const float y = static_cast<float>(y_ind)-mid;
                    const int64_t img_ind = x_ind + y_ind * config.Nx;
                    addPoint(red, red_data, red_centers, x_ind, y, img_ind);
                    addPoint(green, green_data, green_centers, x_ind, y, img_ind);
                    addPoint(blue, blue_data, blue_centers, x_ind, y, img_ind);
                }
            }
        }
    };

    struct floatImagesToPixelsKernel{
        const float* red;
        const float* green;
        const float* blue;
        uint8_t* pixels;
        float multiplier = 1.0f;
        deviceFunction void operator()(const int64_t i) const {
            pixels[4*i] = static_cast<uint8_t>(std::clamp(blue[i], 0.0f, 255.f));
            pixels[4*i + 1] = static_cast<uint8_t>(std::clamp(green[i], 0.0f, 255.f));
            pixels[4*i + 2] = static_cast<uint8_t>(std::clamp(red[i], 0.0f, 255.f));
        }
    };

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
            const float center = config.Nx/2;//findCenter(config.Nx, 1, blue);
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
    static void floatImageToPixels(ActiveDevice& d, const VisualizationConfig config, float multiplier=1.0f){
        d.deviceLaunch((config.width * config.height)/64, 64, 
            floatImagesToPixelsKernel{
                    d.visualization->red.device_ptr(), 
                    d.visualization->green.device_ptr(),
                    d.visualization->blue.device_ptr(),
                    d.visualization->imageGPU.device_ptr(),});
    }
    static void renderFalseColor(ActiveDevice& d, const VisualizationConfig config){

        d.fft(d.visualization->gridTimeSpace.device_ptr(),
            d.visualization->gridFrequency.device_ptr(),
            deviceFFT::D2Z_1D);
        d.visualization->gridTimeSpace.initialize_to_zero();
        switch(config.mode){
            case CoordinateMode::cartesian2D:
                d.deviceLaunch(config.Nt/8, 4, 
                    falseColorPrerenderKernel2D{
                            d.visualization->gridFrequency.device_ptr(), 
                            d.visualization->gridTimeSpace.device_ptr(),
                            config});
                d.deviceLaunch(config.Nt/8, 4, 
                    falseColorRender2DCartesianKernel{
                            d.visualization->red.device_ptr(),
                            d.visualization->green.device_ptr(),
                            d.visualization->blue.device_ptr(),
                            d.visualization->gridTimeSpace.device_ptr(),
                            config});
                floatImageToPixels(d, config, 1.0f);
                break;
            case CoordinateMode::radialSymmetry:
                d.deviceLaunch(config.Nt/8, 4, 
                    falseColorPrerenderKernelCylindrical{
                            d.visualization->gridFrequency.device_ptr(), 
                            d.visualization->gridTimeSpace.device_ptr(),
                            config});
                d.deviceLaunch((config.width * config.height) / 64, 64, 
                floatLinetoImageKernelCartesian{
                    d.visualization->gridTimeSpace.device_ptr(),
                    d.visualization->gridTimeSpace.device_ptr() + config.Nx,
                    d.visualization->gridTimeSpace.device_ptr() + 2*config.Nx,
                    d.visualization->imageGPU.device_ptr(),
                config});
                break;
            default:
                break;

            
        }
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
    }
}
unsigned long renderVisualizationX(VisualizationConfig config){
    ActiveDevice d(config.width, config.height, config.sCPU);
    std::lock_guard<std::mutex> lock(d.visualization->memoryMutex);
    d.visualization->fetchSim(config.sCPU, config.simIndex);
    switch(config.type){
        case VisualizationType::beamPower:
            renderBeamPower(d, config);
            break;
        case VisualizationType::beamFalseColor:
            renderFalseColor(d, config);
            break;
        default:
            return 1;
    }
    d.visualization->syncImages(config.result_pixels);
    return 0;
}