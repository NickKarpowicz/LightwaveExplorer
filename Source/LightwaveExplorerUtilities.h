#pragma once
#include <fstream>
#include <iostream>

#include <mutex>
#include <shared_mutex>
#include <algorithm>
#include <utility>
#include <filesystem>
#include <ranges>
#include <string_view>
#include <charconv>
#include "LightwaveExplorerHelpers.h"
#include "DataStructures.hpp"
#include <miniz/miniz.h>

static const unsigned int threadsPerBlock = 64;
static const unsigned int minGridDimension = 8;


std::string     getBasename(const std::string& fullPath);
double          cModulusSquared(const std::complex<double>& x);
std::size_t     findParenthesesClosure(std::string& a);
double          parameterStringToDouble(const std::string& ss, const double* iBlock, const double* vBlock);
void            stripWhiteSpace(std::string& s);
void            stripLineBreaks(std::string& s);
int             interpretParameters(const std::string& cc, const int n, const double *iBlock, const double *vBlock, double *parameters, bool* defaultMask);

template <FPType T, int number_of_modes, int max_expansion_order>
BeamSpecification<T, number_of_modes, max_expansion_order> beam_specification_from_string(const std::string& str){
    BeamSpecification<T, number_of_modes, max_expansion_order> spec;
    for(int mode = 0; mode < number_of_modes; mode++){

    }

}

template<typename T>
static void zipIntoMemory(std::string zipPath, std::string filename, T* data, std::size_t dataSize){
	mz_zip_archive zip = {};
	mz_zip_reader_init_file(&zip, zipPath.c_str(), 0);
	int fileIndex = mz_zip_reader_locate_file(&zip, filename.c_str(), nullptr, 0);
	if(fileIndex < 0) return;
	mz_zip_archive_file_stat fileStat = {};
    mz_zip_reader_file_stat(&zip, fileIndex, &fileStat);
	mz_zip_reader_extract_to_mem(&zip, fileIndex, data, dataSize, 0);
	mz_zip_reader_end(&zip);
}

template<typename T>
static void zipIntoMemory(std::string zipPath, std::string filename, std::vector<T>& data){
	mz_zip_archive zip = {};
	mz_zip_reader_init_file(&zip, zipPath.c_str(), 0);
	int fileIndex = mz_zip_reader_locate_file(&zip, filename.c_str(), nullptr, 0);
	if(fileIndex < 0) return;
	mz_zip_archive_file_stat fileStat = {};
    mz_zip_reader_file_stat(&zip, fileIndex, &fileStat);
	data = std::vector<T>(fileStat.m_uncomp_size/sizeof(T));
	mz_zip_reader_extract_to_mem(&zip, fileIndex, data.data(), data.size(), 0);
	mz_zip_reader_end(&zip);
}

static bool zipContainsFile(std::string zipPath, std::string filename){
	mz_zip_archive zip = {};
	mz_zip_reader_init_file(&zip, zipPath.c_str(), 0);
	int fileIndex = mz_zip_reader_locate_file(&zip, filename.c_str(), nullptr, 0);
	mz_zip_reader_end(&zip);
	return fileIndex > 0;
}

[[maybe_unused]] static std::string zipGetBasename(const std::string& zipPath){
	mz_zip_archive zip = {};
	mz_zip_reader_init_file(&zip, zipPath.c_str(), 0);
	auto nfiles = mz_zip_reader_get_num_files(&zip);
	if(nfiles > 0){
		char name[1024]={};
		mz_zip_reader_get_filename(&zip, 0, name, 1024);
		std::string interiorName(name);
		mz_zip_reader_end(&zip);
		return std::filesystem::path(interiorName).stem().string();
	}
	else return "ERROR";
}

template <typename T>
std::vector<std::vector<double>> parse_string_to_vecs(const std::string& txt)
{
    auto outer = txt | std::ranges::views::split(';')
                 | std::ranges::views::transform([](auto&& row_range) {
                       std::string_view row_sv{
                           &*row_range.begin(),
                           static_cast<std::size_t>(row_range.end() - row_range.begin())
                       };
                       auto inner_view = row_sv
                                       | std::ranges::views::split(' ')
                                       | std::ranges::views::filter([](auto&& sub) { return !sub.empty(); })
                                       | std::ranges::views::transform(
                                           [](auto&& token_range) {
                                                T val;
                                                std::string_view token_sv{
                                                    &*token_range.begin(),
                                                    static_cast<std::size_t>(token_range.end() -
                                                                            token_range.begin())
                                                };
                                                std::from_chars(token_sv.data(), token_sv.data() + token_sv.size(), val);
                                                return val;
                                            });

                       return std::vector<T>{ std::ranges::begin(inner_view),
                                                    std::ranges::end(inner_view) };
                   });

    auto result = std::vector<std::vector<T>>{
        std::ranges::begin(outer), std::ranges::end(outer)
    };
    if(result.size() == 0){
        throw std::runtime_error("Unparsable or empty string.");
    }
    if(result.size() > 1){
        result.erase(std::remove_if(result.begin(), result.end(), [](auto&& c){ return c.empty(); }));
    }
    return result;
}
