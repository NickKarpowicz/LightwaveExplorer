#include "Fourier.hpp"
#include "vkFFT/vkFFT_AppManagement/vkFFT_DeleteApp.h"
#include "vkFFT/vkFFT_AppManagement/vkFFT_InitializeApp.h"
#include "vkFFT/vkFFT_Structs/vkFFT_Structs.h"
#include "vpp/vpp.hpp"
#include <iostream>
#include <vpp/device_vulkan.hpp>
#include <vulkan/vulkan_core.h>
int main(){
    vpp::VulkanDevice device;
    auto fence = static_cast<VkFence>(device.context->device.createFence({}));
    VkPhysicalDevice vkfft_physical_device = static_cast<VkPhysicalDevice>(device.context->physical_device);
    VkDevice vkfft_device = static_cast<VkDevice>(device.context->device);
    VkQueue vkfft_queue = static_cast<VkQueue>(device.context->queue);
    VkCommandPool vkfft_pool = static_cast<VkCommandPool>(device.context->command_pool);

    uint32_t x_dim = 256;
    uint32_t y_dim = 256;
    uint64_t buffer_size = x_dim*y_dim*sizeof(float);
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
    configuration.bufferSize = &buffer_size;
    VkFFTApplication app = {};
    VkFFTResult res = initializeVkFFT(&app, configuration);
    std::cout << "result: " << getVkFFTErrorString(res) << '\n';
    deleteVkFFT(&app);
    return res;
}
