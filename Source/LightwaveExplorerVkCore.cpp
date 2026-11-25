#include "Fourier.hpp"
#include "vkFFT/vkFFT_AppManagement/vkFFT_DeleteApp.h"
#include "vkFFT/vkFFT_AppManagement/vkFFT_InitializeApp.h"
#include "vkFFT/vkFFT_Structs/vkFFT_Structs.h"
#include "vpp/vpp.hpp"
#include <iostream>
#include <vector>
#include <complex>
#include <print>
#include <vpp/device_vulkan.hpp>
#include <vulkan/vulkan_core.h>

template <typename T>
void print_matrix(const std::vector<std::complex<T>>& data, const size_t x_length, const size_t y_length, const size_t pages, const float factor){

    for(size_t p = 0; p<pages;p++){
        size_t offset = p * y_length * x_length;
        for(size_t y = 0; y<y_length; y++){
            for(size_t x = 0; x<x_length; x++){
                float value = factor * std::abs(data[x + y*x_length + offset]);
                if(value>1) std::cout << std::format("\x1b[91m{:2.0f}\x1b[0m", value);
                else std::cout << std::format("\x1b[90m{:2.0f}\x1b[0m", value);
                //std::cout << std::format("{:2.1f} ", factor * data[x + y*x_length].real());
            }
            std::cout << "\n";
        }
        std::cout << "==========\n";
    }
}

int main(){
    vpp::VulkanDevice device;
    auto fence = static_cast<VkFence>(device.context->device.createFence({}));
    VkPhysicalDevice vkfft_physical_device = static_cast<VkPhysicalDevice>(device.context->physical_device);
    VkDevice vkfft_device = static_cast<VkDevice>(device.context->device);
    VkQueue vkfft_queue = static_cast<VkQueue>(device.context->queue);
    VkCommandPool vkfft_pool = static_cast<VkCommandPool>(device.context->command_pool);

    uint32_t x_dim = 16;
    uint32_t y_dim = 4;
    uint64_t buffer_size = x_dim*y_dim*sizeof(float);

    uint32_t kx_dim = x_dim/2 + 1;
    uint32_t ky_dim = y_dim;
    uint64_t result_buffer_size = kx_dim * ky_dim * 2 * sizeof(float);

    std::vector<float> data_local(x_dim*y_dim,0.0f);
    std::vector<std::complex<float>> result_cpu(result_buffer_size,0.0f);
    data_local[10] = 1.0f;
    data_local[27] = 1.0f;
    auto data_vk = device.create_buffer_from_vec(data_local);
    auto result_vk = device.create_buffer_from_vec(result_cpu);

    std::cout << "testing Fourier\n";
    VkFFTConfiguration configuration = {};
    configuration.physicalDevice=&vkfft_physical_device;
    configuration.device = &vkfft_device;
    configuration.queue = &vkfft_queue;
    configuration.commandPool = &vkfft_pool;
    configuration.fence = &fence;
    configuration.FFTdim=2;
    configuration.size[0] = x_dim;
    configuration.size[1] = y_dim;
    configuration.bufferSize = &result_buffer_size;
    configuration.inputBufferSize = &buffer_size;
    configuration.performR2C=1;
    configuration.isInputFormatted=1;
    VkFFTApplication app = {};
    VkFFTResult res = initializeVkFFT(&app, configuration);
    std::cout << "create result: " << getVkFFTErrorString(res) << '\n';

    VkFFTLaunchParams launch = {};
    launch.inputBuffer = reinterpret_cast<VkBuffer*>(&(data_vk.buffer));
    launch.buffer = reinterpret_cast<VkBuffer*>(&(result_vk.buffer));
    launch.commandBuffer = device.get_raw_active_command_buffer();
    auto result_fft = VkFFTAppend(&app,-1,&launch);
    device.go();
    device.wait();
    std::cout << "fft result: " << getVkFFTErrorString(result_fft) << '\n';
    auto result = result_vk.to_vec<std::complex<float>>();
    print_matrix(result, kx_dim, ky_dim, 1, 1);
    deleteVkFFT(&app);
    return res;
}
