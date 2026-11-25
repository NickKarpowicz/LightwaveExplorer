#include "vkFFT.h" // IWYU pragma: keep
#include "pocketfft_hdronly.h"
#include <cstddef>
#include <vpp/vpp.hpp>
#include <iostream>
#include <concepts>
enum class FourierKind{R2C,C2C};

template <typename Backend, std::floating_point F>
struct Fourier {
    void forward(void* in, void* out){
        static_cast<Backend*>(this)->forward_impl(in, out);
    }
    void backward(void* in, void* out){
        static_cast<Backend*>(this)->backward_impl(in, out);
    }
};

template <std::floating_point F>
struct CpuFourier : public Fourier<CpuFourier<F>, F>{
    const ptrdiff_t fp_size;
    const ptrdiff_t complex_size;
    FourierKind fft_kind;
   	pocketfft::shape_t shape;
	pocketfft::shape_t axes;
	pocketfft::stride_t stride_time;
	pocketfft::stride_t stride_freq;
	CpuFourier(
        size_t size_x = 0,
        size_t size_y = 1,
        size_t size_z = 1,
        size_t batch = 1,
        FourierKind kind = FourierKind::C2C
	) : fp_size(static_cast<ptrdiff_t>(sizeof(F))),
	complex_size(2*fp_size),
	fft_kind(kind){
        uint64_t dim = 1;
        ptrdiff_t size_kx = kind == FourierKind::R2C ? size_x / 2u + 1u : size_x;
        if(size_y>1) dim = 2;
        if(size_z>1) dim = 3; //I'm assuming size_y != 1
        ptrdiff_t x_element_size = kind == FourierKind::R2C ? fp_size : complex_size;
        switch (dim){
            case 1:{
                shape = pocketfft::shape_t{size_x};
                axes=pocketfft::shape_t{0};
                stride_time = pocketfft::stride_t{x_element_size};
                stride_freq = pocketfft::stride_t{complex_size};
                }
                break;
            case 2:{
                shape = pocketfft::shape_t{size_y, size_x};
                axes=pocketfft::shape_t{0, 1};
                stride_time = pocketfft::stride_t{
                    x_element_size * static_cast<ptrdiff_t>(size_x),
                    x_element_size};
                stride_freq = pocketfft::stride_t{
                    complex_size * size_kx,
                    complex_size};
                }
                break;
            case 3:{
                shape = pocketfft::shape_t{size_z, size_y, size_x};
                axes=pocketfft::shape_t{0, 1, 2};
                stride_time = pocketfft::stride_t{
                    x_element_size * static_cast<ptrdiff_t>(size_x) * static_cast<ptrdiff_t>(size_y),
                    x_element_size * static_cast<ptrdiff_t>(size_x),
                    x_element_size};
                stride_freq = pocketfft::stride_t{
                    complex_size * size_kx * static_cast<ptrdiff_t>(size_y),
                    complex_size * size_kx,
                    complex_size};
                }
                break;
        }
	}

	void forward_impl(void* in, void* out){
	    switch(fft_kind){
			case FourierKind::C2C:
                pocketfft::c2c(
                    shape,
                    stride_time,
                    stride_freq,
                    axes,
                    pocketfft::FORWARD,
   					reinterpret_cast<std::complex<F>*>(in),
   					reinterpret_cast<std::complex<F>*>(out),
   					1.0f);
                return;
			case FourierKind::R2C:
                pocketfft::r2c(
                    shape,
                    stride_time,
                    stride_freq,
                    axes,
                    pocketfft::FORWARD,
   					reinterpret_cast<F*>(in),
   					reinterpret_cast<std::complex<F>*>(out),
   					1.0f);
                return;
		}
	}

	void backward_impl(void* in, void* out){
	    switch(fft_kind){
			case FourierKind::C2C:
                pocketfft::c2c(
                    shape,
                    stride_time,
                    stride_freq,
                    axes,
                    pocketfft::BACKWARD,
   					reinterpret_cast<std::complex<F>*>(in),
   					reinterpret_cast<std::complex<F>*>(out),
   					1.0f);
                return;
			case FourierKind::R2C:
                pocketfft::r2c(
                    shape,
                    stride_time,
                    stride_freq,
                    axes,
                    pocketfft::BACKWARD,
   					reinterpret_cast<F*>(in),
   					reinterpret_cast<std::complex<F>*>(out),
   					1.0f);
                return;
		}
	}
};

template <std::floating_point F>
struct VulkanFourier : public Fourier<VulkanFourier<F>, F>{
    uint64_t buffer_size;
    uint64_t result_buffer_size;
    vpp::VulkanDevice& device;
    VkFence vk_fence;
    VkPhysicalDevice vk_physical_device;
    VkDevice vk_device;
    VkQueue vk_queue;
    VkCommandPool vk_command_pool;
    VkFFTConfiguration configuration;
    VkFFTApplication app;
    ~VulkanFourier(){
        deleteVkFFT(&app);
    }
    VulkanFourier(
        vpp::VulkanDevice& in_device,
        size_t size_x = 0,
        size_t size_y = 1,
        size_t size_z = 1,
        size_t batch = 1,
        FourierKind kind = FourierKind::C2C) :
            buffer_size(size_x*size_y*size_z*batch*(kind==FourierKind::R2C ? 1u : 2u)*sizeof(F)),
            result_buffer_size(kind==FourierKind::R2C ?
               (size_x/2u + 1u)*size_y*size_z*batch*2*sizeof(F)
               : size_x*size_y*size_z*batch*2*sizeof(F)),
            device(in_device),
            vk_fence(static_cast<VkFence>(device.context->device.createFence({}))),
            vk_physical_device(static_cast<VkPhysicalDevice>(device.context->physical_device)),
            vk_device(static_cast<VkDevice>(device.context->device)),
            vk_queue(static_cast<VkQueue>(device.context->queue)),
            vk_command_pool(static_cast<VkCommandPool>(device.context->command_pool)),
            configuration({}),
            app({}) {
                uint64_t dim = 1;
                if(size_y>1) dim = 2;
                if(size_z>1) dim = 3; //I'm assuming size_y != 1
                std::cout << "dim: " << dim << '\n';
                std::cout << "buffer_size: " << buffer_size << '\n';
                std::cout << "buffer_size_out: " << result_buffer_size << '\n';
                configuration.physicalDevice = &vk_physical_device;
                configuration.device = &vk_device;
                configuration.queue = &vk_queue;
                configuration.commandPool = &vk_command_pool;
                configuration.fence = &vk_fence;
                configuration.FFTdim = dim;
                configuration.size[0] = size_x;
                configuration.size[1] = dim > 1 ? size_y : 0;
                configuration.size[2] = dim > 2 ? size_z : 0;
                configuration.numberBatches = batch;
                configuration.bufferSize = &result_buffer_size;
                configuration.inputBufferSize = &buffer_size;
                configuration.performR2C = kind == FourierKind::R2C ? 1 : 0;
                configuration.doublePrecision = sizeof(F) == 8 ? 1 : 0;
                configuration.isInputFormatted=1;
                VkFFTResult result = initializeVkFFT(&app, configuration);
                std::cout << "create result: " << getVkFFTErrorString(result) << '\n';
            }
    void forward_impl(void* in, void* out){
        VkFFTLaunchParams launch{};
        launch.inputBuffer = reinterpret_cast<VkBuffer*>(in);
        launch.buffer = reinterpret_cast<VkBuffer*>(out);
        launch.commandBuffer = device.get_raw_active_command_buffer();
        VkFFTResult result = VkFFTAppend(&app,-1,&launch);
        std::cout << "append result: " << getVkFFTErrorString(result) << '\n';
    }
    void backward_impl(void* in, void* out){
        VkFFTLaunchParams launch{};
        launch.inputBuffer = reinterpret_cast<VkBuffer*>(in);
        launch.buffer = reinterpret_cast<VkBuffer*>(out);
        launch.commandBuffer = device.get_raw_active_command_buffer();
        VkFFTResult result = VkFFTAppend(&app,1,&launch);
        std::cout << "append result: " << getVkFFTErrorString(result) << '\n';
    }
};
