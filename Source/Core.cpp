#include "Fourier.hpp"
#include "LightwaveExplorerInterfaceClasses.hpp"
#include "vpp/vpp.hpp"
#include <iostream>
#include <vector>
#include <complex>
#include <concepts>
#include <vpp/kernel.hpp>

template <std::floating_point F>
uint32_t run_fitting(const SimulationParameterSet& launch_parameters){
    return 0u;
}

template <std::floating_point T>
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

template <std::floating_point F>
struct TestKernel : public vpp::VppKernel<TestKernel<F>>{
    F* buffer;
    void main(const uint32_t i) const {
        buffer[i] *= 2.0f;
    }
};
std::vector<std::complex<float>> forward_pocketfft(std::vector<float>& data_local, uint32_t x_dim, uint32_t y_dim){
    CpuFourier<float> fft(x_dim,y_dim,1,1,FourierKind::R2C);
    std::vector<std::complex<float>> result(y_dim * (x_dim/2u + 1u));
    fft.forward(data_local.data(), result.data());
	return result;
}

std::vector<std::complex<float>> test_vkfft(std::vector<float>& data_local, uint32_t x_dim, uint32_t y_dim){
    vpp::VulkanDevice device;
    VulkanFourier<float> fft(device,x_dim,y_dim,1,1,FourierKind::R2C);
    auto data_vk = device.create_buffer_from_vec(data_local);
    auto result_vk = device.create_buffer(fft.result_buffer_size);
    fft.forward(&(data_vk.buffer), &(result_vk.buffer));
    device.launch_kernel<24>(TestKernel<float>{.buffer=data_vk.as<float>()},10);
    device.go();
    device.wait();
    return result_vk.to_vec<std::complex<float>>();
}

int main(){
    uint32_t x_dim = 24;
    uint32_t y_dim = 10;
    uint32_t kx_dim = x_dim/2 + 1;
    uint32_t ky_dim = y_dim;

    std::vector<float> data_local(x_dim*y_dim,0.0f);
    data_local[10] = 1.0f;
    data_local[192] = 1.0f;

    auto result = test_vkfft(data_local, x_dim, y_dim);
    auto result_pocket = forward_pocketfft(data_local, x_dim, y_dim);
    print_matrix(result_pocket, kx_dim, ky_dim, 1, 1);
    print_matrix(result, kx_dim, ky_dim, 1, 1);
}
