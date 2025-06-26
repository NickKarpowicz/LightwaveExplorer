#pragma once
#include <fstream>
#include <iostream>

#include <mutex>
#include <shared_mutex>
#include <algorithm>
#include <utility>
#include <filesystem>
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
