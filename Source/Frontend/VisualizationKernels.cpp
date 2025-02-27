#ifdef LWEFLOATINGPOINT
#undef LWEFLOATINGPOINT
#endif
#define LWEFLOATINGPOINT 32
#include "../LightwaveExplorerTrilingual.h"
#include "../LightwaveExplorerInterfaceClasses.hpp"
#include "../LightwaveExplorerHelpers.h"
ActiveDevice* d_vis = nullptr;

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

    deviceFunction static inline int64_t indexShiftFFT(const int64_t i, const int64_t Nx){
        return (i < Nx/2) ?
            Nx/2 + i:
            i - Nx/2;
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
            const float A2 = 1e-18f * Et[i]*Et[i];
            const float B2 = 1e-18f * Et[i + config.Ngrid] * Et[i + config.Ngrid];
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

            const float scale = 1.0f / ((1e-9 * config.df) * config.Nt * config.Nt);
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

    struct falseColorFarfieldPrerenderKernel2D{
        const deviceComplex* Ef;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float f = (i+1) * config.df;
            const float k0 = twoPi<float>()*f/lightC<float>()/config.dk;
            const float scale = 1.0f/ ((1e-20 * config.dk * config.dk) * config.df * config.Ngrid * config.Ngrid);
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
                const float theta = j*config.dTheta - config.maxAngle;
                const int64_t k_x = config.Nx/2 + 
                    static_cast<int64_t>(round(k0*deviceFPLib::sin(theta)));
                if(k_x > 0 && k_x < config.Nx){
                    int64_t idx_k = indexShiftFFT(k_x, config.Nx);
                    float Ef2 = deviceNorm(Ef[(i+1) + idx_k*config.Nf]) 
                    + deviceNorm(Ef[(i+1) + idx_k*config.Nf + config.Ngrid_complex]);
                    int offset = j * config.Nt/2;
                    red[offset] = filter_red * Ef2;
                    green[offset]= filter_green * Ef2;
                    blue[offset] = filter_blue * Ef2;
                }
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

            const float scale = 1.0 / ((1e-9 * config.df) * config.Nt * config.Nt);
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

    struct falseColorFarfieldPrerenderKernelCylindrical{
        const deviceComplex* Ef;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const float f = (i+1) * config.df;
            const float k0 = twoPi<float>()*f/lightC<float>()/config.dk;
            const float scale = 1.0f / ((1e-20 * config.dk * config.dk) * config.Ngrid * config.Ngrid * config.df);
            const float filter_red = scale * config.red_amplitude * expf(-(f-config.red_f0)*(f-config.red_f0)/(2*config.red_sigma*config.red_sigma));
            const float filter_green = scale * config.green_amplitude * expf(-(f-config.green_f0)*(f-config.green_f0)/(2*config.green_sigma*config.green_sigma));
            const float filter_blue = scale * config.blue_amplitude * expf(-(f-config.blue_f0)*(f-config.blue_f0)/(2*config.blue_sigma*config.blue_sigma));
            
            float* red_line = workspace;
            float* green_line = red_line + config.Nx;
            float* blue_line = green_line + config.Nx; 

            for(int j{}; j<config.Nx; ++j){
                const float theta = j*config.dTheta - config.maxAngle;
                const int64_t k_x = config.Nx/2 + 
                    static_cast<int64_t>(round(k0*deviceFPLib::sin(theta)));
                if(k_x > 0 && k_x < config.Nx){
                    int64_t idx_k = indexShiftFFT(k_x, config.Nx);
                    float Ef2 = deviceNorm(Ef[(i+1) + idx_k*config.Nf])
                    + deviceNorm(Ef[(i+1) + idx_k*config.Nf + config.Ngrid_complex]);
                    atomicAdd(&red_line[j], filter_red * Ef2);
                    atomicAdd(&green_line[j], filter_green * Ef2);
                    atomicAdd(&blue_line[j], filter_blue * Ef2);
                }  
            }
        }
    };

    struct falseColor3DKernel{
        const deviceComplex* Ef;
        float* red;
        float* green;
        float* blue;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const int64_t f_index = (i % (config.Nt/2)) + 1;
            const int64_t x_index = (i / (config.Nt/2)) % config.Nx;
            const int64_t y_index = (i / (config.Nt/2)) / config.Nx;
            const int64_t field_index = f_index + (x_index + y_index * config.Nx) * config.Nf;
            const int64_t image_index = x_index + config.Nx * y_index;
            const float f = f_index * config.df;

            const float scale = 1.0f /(1e-9 * config.df* config.Nt * config.Nt);
            const float filter_red = scale * config.red_amplitude * expf(-(f-config.red_f0)*(f-config.red_f0)/(2*config.red_sigma*config.red_sigma));
            const float filter_green = scale * config.green_amplitude * expf(-(f-config.green_f0)*(f-config.green_f0)/(2*config.green_sigma*config.green_sigma));
            const float filter_blue = scale * config.blue_amplitude * expf(-(f-config.blue_f0)*(f-config.blue_f0)/(2*config.blue_sigma*config.blue_sigma));

            float Ef2 = deviceNorm(Ef[field_index]) + deviceNorm(Ef[field_index + config.Ngrid_complex]);
            atomicAdd(&red[image_index],filter_red * Ef2);
            atomicAdd(&green[image_index],filter_green * Ef2);
            atomicAdd(&blue[image_index],filter_blue * Ef2);
        }
    };

    struct falseColorFarfield3DKernel{
        const deviceComplex* Ef;
        float* red;
        float* green;
        float* blue;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            //relative to fixed max angle
            //find angle of current point
            //find corresponding kx value omega/c * sin theta
            //round(kx/dk)+Nx/2 is field index
            const int64_t f_index = (i % (config.Nt/2)) + 1;
            const float f = f_index * config.df;
            const float w = twoPi<float>() * f;
            const int64_t x_index = (i / (config.Nt/2)) % config.Nx;
            const int64_t y_index = (i / (config.Nt/2)) / config.Nx;
            const float angle_x = x_index * config.dTheta - config.maxAngle;
            const float angle_y = y_index * config.dTheta - config.maxAngle;
            const float kx = deviceFPLib::sin(angle_x)*w/lightC<float>();
            const float ky = deviceFPLib::sin(angle_y)*w/lightC<float>();
            const int64_t kx_index = round(kx/config.dk) + config.Nx/2;
            const int64_t ky_index = round(ky/config.dk) + config.Nx/2;
            if(kx_index >= config.Nx 
                || ky_index >= config.Nx
                || kx_index<0
                || ky_index<0) return;
            const int64_t field_index = f_index + (indexShiftFFT(kx_index, config.Nx) + indexShiftFFT(ky_index, config.Ny) * config.Nx) * config.Nf;
            const int64_t image_index = x_index + config.Nx * y_index;
            
            const float scale = 1.0f / (deviceFPLib::pow(config.dk * 3e-5, 4)* config.Ngrid * config.Ngrid);
            const float filter_red = scale * config.red_amplitude * expf(-(f-config.red_f0)*(f-config.red_f0)/(2*config.red_sigma*config.red_sigma));
            const float filter_green = scale * config.green_amplitude * expf(-(f-config.green_f0)*(f-config.green_f0)/(2*config.green_sigma*config.green_sigma));
            const float filter_blue = scale * config.blue_amplitude * expf(-(f-config.blue_f0)*(f-config.blue_f0)/(2*config.blue_sigma*config.blue_sigma));

            float Ef2 = deviceNorm(Ef[field_index]) + deviceNorm(Ef[field_index + config.Ngrid_complex]);
            atomicAdd(&red[image_index],filter_red * Ef2);
            atomicAdd(&green[image_index],filter_green * Ef2);
            atomicAdd(&blue[image_index],filter_blue * Ef2);
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
            const int64_t Nf = config.Nt/2; 
            const int64_t f_ind = i % Nf; //not +1 because I skip 0
            
            int64_t dataSize = Nf*(config.Nx);
            float* red_data = workspace;
            float* green_data = red_data + dataSize;
            float* blue_data = green_data + dataSize;
            float* red_centers = blue_data + dataSize;
            float* green_centers = red_centers + config.Nt/2;
            float* blue_centers = green_centers + config.Nt/2;
            auto addPoint = [&](
                float* channel, 
                const float* data, 
                const float* centers, 
                const int x_ind, 
                const float y, 
                const int64_t img_ind){
                    const float center = centers[f_ind];
                    const float x = static_cast<float>(x_ind)-center;
                    const float r = hypotf(x,y);
                    const int64_t idx1 = static_cast<int64_t>((r+center));
                    const bool idx1Valid = (idx1 < config.Nx) && (idx1>=0);
                    if(idx1Valid){
                        atomicAdd(&channel[img_ind], data[f_ind + idx1*(config.Nt/2)]); 
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

    struct falseColorRender2DCartesianNoRotationKernel{
        float* red;
        float* green;
        float* blue;
        float* workspace;
        const VisualizationConfig config;
        deviceFunction void operator()(const int64_t i) const {
            const int64_t Nf = config.Nt/2;
            const int64_t f_ind = i % Nf;
            
            int64_t dataSize = Nf*(config.Nx);
            float* red_data = workspace;
            float* green_data = red_data + dataSize;
            float* blue_data = green_data + dataSize;
        
            for(int x_ind{}; x_ind<config.Nx; ++x_ind){
                for(int y_ind{}; y_ind<config.Nx; ++y_ind){
                    const int64_t img_ind = x_ind + y_ind * config.Nx;
                    atomicAdd(&red[img_ind], red_data[f_ind + x_ind*(config.Nt/2)]);
                    atomicAdd(&green[img_ind], green_data[f_ind + x_ind*(config.Nt/2)]);
                    atomicAdd(&blue[img_ind], blue_data[f_ind + x_ind*(config.Nt/2)]);
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
                    blue[idx],
                    0.f, 255.f));;
            image[4*i + 1] = static_cast<uint8_t>(
                std::clamp(
                    green[idx],
                    0.f, 255.f));;
            image[4*i + 2] = static_cast<uint8_t>(
                std::clamp(
                    red[idx],
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
                        0.5f*(blue[idx1]+ blue[idx2]),
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        0.5f*(green[idx1]+green[idx2]),
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        0.5f*(red[idx1] + red[idx2]),
                        0.f, 255.f));;
            }
            else if(idx1Valid){
                image[4*i] = static_cast<uint8_t>(
                    std::clamp(
                        blue[idx1],
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        green[idx1],
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        red[idx1],
                        0.f, 255.f));;
            }
            else if(idx2Valid){
                image[4*i] = static_cast<uint8_t>(
                    std::clamp(
                        blue[idx2],
                        0.f, 255.f));;
                image[4*i + 1] = static_cast<uint8_t>(
                    std::clamp(
                        green[idx2],
                        0.f, 255.f));;
                image[4*i + 2] = static_cast<uint8_t>(
                    std::clamp(
                        red[idx2],
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
    static void floatImageToPixels(const VisualizationConfig config, float multiplier=1.0f){
        d_vis->deviceLaunch((config.width * config.height)/64, 64, 
            floatImagesToPixelsKernel{
                    d_vis->visualization->red.device_ptr(), 
                    d_vis->visualization->green.device_ptr(),
                    d_vis->visualization->blue.device_ptr(),
                    d_vis->visualization->imageGPU.device_ptr(),});
    }
    static void renderFalseColor(const VisualizationConfig config){
        d_vis->fft(d_vis->visualization->gridTimeSpace.device_ptr(),
            d_vis->visualization->gridFrequency.device_ptr(),
            deviceFFT::D2Z_1D);
        d_vis->visualization->gridTimeSpace.initialize_to_zero();
        switch(config.mode){
            case CoordinateMode::cartesian2D: [[fallthrough]];
            case CoordinateMode::FDTD2D:
                d_vis->deviceLaunch(config.Nt/8, 4, 
                    falseColorPrerenderKernel2D{
                            d_vis->visualization->gridFrequency.device_ptr(), 
                            d_vis->visualization->gridTimeSpace.device_ptr(),
                            config});
                if(config.rotate2D){
                    d_vis->deviceLaunch(config.Nt/8, 4, 
                        falseColorRender2DCartesianKernel{
                                d_vis->visualization->red.device_ptr(),
                                d_vis->visualization->green.device_ptr(),
                                d_vis->visualization->blue.device_ptr(),
                                d_vis->visualization->gridTimeSpace.device_ptr(),
                                config});
                }
                else{
                    d_vis->deviceLaunch(config.Nt/8, 4, 
                        falseColorRender2DCartesianNoRotationKernel{
                                d_vis->visualization->red.device_ptr(),
                                d_vis->visualization->green.device_ptr(),
                                d_vis->visualization->blue.device_ptr(),
                                d_vis->visualization->gridTimeSpace.device_ptr(),
                                config});
                }
                
                floatImageToPixels(config, 1.0f);
                d_vis->visualization->red.initialize_to_zero();
                d_vis->visualization->green.initialize_to_zero();
                d_vis->visualization->blue.initialize_to_zero();
                break;
            case CoordinateMode::radialSymmetry:
                d_vis->deviceLaunch(config.Nt/8, 4, 
                    falseColorPrerenderKernelCylindrical{
                            d_vis->visualization->gridFrequency.device_ptr(), 
                            d_vis->visualization->gridTimeSpace.device_ptr(),
                            config});
                d_vis->deviceLaunch((config.width * config.height) / 64, 64, 
                floatLinetoImageKernel{
                    d_vis->visualization->gridTimeSpace.device_ptr(),
                    d_vis->visualization->gridTimeSpace.device_ptr() + config.Nx,
                    d_vis->visualization->gridTimeSpace.device_ptr() + 2*config.Nx,
                    d_vis->visualization->imageGPU.device_ptr(),
                config});
                break;
            case CoordinateMode::cartesian3D: [[fallthrough]];
            case CoordinateMode::FDTD3D:
                d_vis->deviceLaunch(config.Ngrid/64, 32, falseColor3DKernel{
                    d_vis->visualization->gridFrequency.device_ptr(),
                    d_vis->visualization->red.device_ptr(),
                    d_vis->visualization->green.device_ptr(),
                    d_vis->visualization->blue.device_ptr(),
                    config});
                floatImageToPixels(config, 1.0f);
                    d_vis->visualization->red.initialize_to_zero();
                    d_vis->visualization->green.initialize_to_zero();
                    d_vis->visualization->blue.initialize_to_zero();
            default:
                break;
        }
    }

    static void renderFalseColorFarfield(const VisualizationConfig config){
        d_vis->visualization->gridTimeSpace.initialize_to_zero();
        switch(config.mode){
            case CoordinateMode::cartesian2D: [[fallthrough]];
            case CoordinateMode::FDTD2D:
                d_vis->deviceLaunch(config.Nt/8, 4, 
                    falseColorFarfieldPrerenderKernel2D{
                            d_vis->visualization->gridFrequency.device_ptr(), 
                            d_vis->visualization->gridTimeSpace.device_ptr(),
                            config});
                if(config.rotate2D){
                    d_vis->deviceLaunch(config.Nt/8, 4, 
                        falseColorRender2DCartesianKernel{
                                d_vis->visualization->red.device_ptr(),
                                d_vis->visualization->green.device_ptr(),
                                d_vis->visualization->blue.device_ptr(),
                                d_vis->visualization->gridTimeSpace.device_ptr(),
                                config});
                }
                else{
                    d_vis->deviceLaunch(config.Nt/8, 4, 
                        falseColorRender2DCartesianNoRotationKernel{
                                d_vis->visualization->red.device_ptr(),
                                d_vis->visualization->green.device_ptr(),
                                d_vis->visualization->blue.device_ptr(),
                                d_vis->visualization->gridTimeSpace.device_ptr(),
                                config});
                }
                
                floatImageToPixels(config, 1.0f);
                d_vis->visualization->red.initialize_to_zero();
                d_vis->visualization->green.initialize_to_zero();
                d_vis->visualization->blue.initialize_to_zero();
                break;
            case CoordinateMode::radialSymmetry:
                d_vis->deviceLaunch(config.Nt/8, 4, 
                    falseColorFarfieldPrerenderKernelCylindrical{
                            d_vis->visualization->gridFrequency.device_ptr(), 
                            d_vis->visualization->gridTimeSpace.device_ptr(),
                            config});
                d_vis->deviceLaunch((config.width * config.height) / 64, 64, 
                floatLinetoImageKernelCartesian{
                    d_vis->visualization->gridTimeSpace.device_ptr(),
                    d_vis->visualization->gridTimeSpace.device_ptr() + config.Nx,
                    d_vis->visualization->gridTimeSpace.device_ptr() + 2*config.Nx,
                    d_vis->visualization->imageGPU.device_ptr(),
                config});
                break;
            case CoordinateMode::cartesian3D: [[fallthrough]];
            case CoordinateMode::FDTD3D:
                d_vis->deviceLaunch(config.Ngrid/64, 32, falseColorFarfield3DKernel{
                    d_vis->visualization->gridFrequency.device_ptr(),
                    d_vis->visualization->red.device_ptr(),
                    d_vis->visualization->green.device_ptr(),
                    d_vis->visualization->blue.device_ptr(),
                    config});
                floatImageToPixels(config, 1.0f);
                    d_vis->visualization->red.initialize_to_zero();
                    d_vis->visualization->green.initialize_to_zero();
                    d_vis->visualization->blue.initialize_to_zero();
            default:
                break;
        }
    }

    static void renderBeamPower(const VisualizationConfig config){
        d_vis->visualization->gridFrequency.initialize_to_zero();

        switch(config.mode){
            case CoordinateMode::cartesian2D:
                d_vis->deviceLaunch(config.Ngrid/64, 64, 
                    beamPowerRenderKernel2D{
                            d_vis->visualization->gridTimeSpace.device_ptr(), 
                            d_vis->visualization->red.device_ptr(),
                            d_vis->visualization->green.device_ptr(),
                            d_vis->visualization->blue.device_ptr(),
                            config});
                d_vis->deviceLaunch((config.width * config.height) / 64, 64, 
                    floatLinetoImageKernelCartesian{
                        d_vis->visualization->red.device_ptr(),
                        d_vis->visualization->green.device_ptr(),
                        d_vis->visualization->blue.device_ptr(),
                        d_vis->visualization->imageGPU.device_ptr(),
                        config});
                break;
            case CoordinateMode::radialSymmetry:
                d_vis->deviceLaunch(config.Ngrid/64, 64, 
                    beamPowerRenderKernel2D{
                            d_vis->visualization->gridTimeSpace.device_ptr(), 
                            d_vis->visualization->red.device_ptr(),
                            d_vis->visualization->green.device_ptr(),
                            d_vis->visualization->blue.device_ptr(),
                            config});
                d_vis->deviceLaunch((config.width * config.height) / 64, 64, 
                        floatLinetoImageKernel{
                            d_vis->visualization->red.device_ptr(),
                            d_vis->visualization->green.device_ptr(),
                            d_vis->visualization->blue.device_ptr(),
                            d_vis->visualization->imageGPU.device_ptr(),
                            config});
                break;
            default:
                break;
        }
        d_vis->visualization->red.initialize_to_zero();
        d_vis->visualization->green.initialize_to_zero();
        d_vis->visualization->blue.initialize_to_zero();
    }
}
unsigned long renderVisualizationX(VisualizationConfig config){
    if(d_vis == nullptr){
        d_vis = new ActiveDevice(config.width, config.height, config.sCPU);
    }
    else{
        d_vis->reset(config.sCPU);
        d_vis->visualization->setImageDimensions(config.width, config.height);
    }
    d_vis->visualization->fetchSim(config.sCPU, config.simIndex);
    switch(config.type){
        case VisualizationType::beamPower:
            renderBeamPower(config);
            break;
        case VisualizationType::beamFalseColor:
            renderFalseColor(config);
            break;
        case VisualizationType::farFieldFalseColor:
            renderFalseColorFarfield(config);
            break;
        default:
            return 1;
    }
    d_vis->visualization->syncImages(config.result_pixels);
    d_vis->visualization->imageGPU.initialize_to_zero();
    return 0;
}