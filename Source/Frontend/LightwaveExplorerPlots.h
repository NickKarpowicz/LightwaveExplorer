#include <sstream>
#include <fstream>
#include <complex>
#include <vector>
#include <array>
#include <string>
#include <optional>
#include <algorithm>
#include <thread>
#include <mutex>
#include <gcem.hpp>
#include <cairo.h>
#include "../LightwaveExplorerHelpers.h"
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

//temporary set of macros until std::format is on all platforms
#if defined __linux__
#include<fmt/format.h>
#define Sformat fmt::format
#define Svformat fmt::vformat
#define Smake_format_args fmt::make_format_args
#elif defined __APPLE__
#include<fmt/format.h>
#define Sformat fmt::format
#define Svformat fmt::vformat
#define Smake_format_args fmt::make_format_args
#else
#include <format>
#define Sformat std::format
#define Svformat std::vformat
#define Smake_format_args std::make_format_args
#endif

//Limit the number of threads used to draw the interface if the processor supports a lot
const int interfaceThreads = 
maxN(2, minN(4, static_cast<int>(std::thread::hardware_concurrency() / 2)));

constexpr inline uint8_t clampcast(const double d){
    return static_cast<uint8_t>(std::clamp(d, 0.0, 255.0));
}
enum ColorMap {
    greyscale = 0,
    jet = 1,
    jet_white = 2,
    purple = 3,
    cyan_magenta = 4
};

consteval std::array<std::array<uint8_t, 3>, 256> createColormap(const ColorMap cm) {
    std::array<std::array<uint8_t, 3>, 256> colorMap{};
    constexpr double oneOver255 = 1.0 / 255.0;
    for (int j = 0; j < 256; ++j) {
        const double nval = static_cast<double>(j) * oneOver255;
        switch (cm) {
        case ColorMap::greyscale:
            colorMap[j][0] = static_cast<uint8_t>(j);
            colorMap[j][1] = static_cast<uint8_t>(j);
            colorMap[j][2] = static_cast<uint8_t>(j);
            break;
        case ColorMap::jet:
            colorMap[j][0] = clampcast(255. * gcem::cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = clampcast(255. * gcem::cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = clampcast(255. * gcem::sin(vPi<double>() * nval / 2.));
            break;
        case ColorMap::jet_white:
            colorMap[j][0] = clampcast(255. * gcem::cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = clampcast(255. * gcem::cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = clampcast(255. * gcem::sin(vPi<double>() * nval / 2.));
            if (nval < 0.02) {
                colorMap[j][0] = 255u;
                colorMap[j][1] = 128u;
                colorMap[j][2] = 128u;
            }
            if (nval < 0.01) {
                colorMap[j][0] = 255u;
                colorMap[j][1] = 255u;
                colorMap[j][2] = 255u;
            }
            break;
        case ColorMap::purple:
            {
            double dj = static_cast<double>(j);
            colorMap[j][0] = clampcast(255. *
                (0.998 * gcem::exp(-gcem::pow(7.7469e-03 * (dj - 160.), 6))
                    + 0.22 * gcem::exp(-gcem::pow(0.016818 * (dj - 305.), 4))));
            colorMap[j][1] = clampcast(255. *
                (0.022 * gcem::exp(-gcem::pow(0.042045 * (dj - 25.), 4))
                    + 0.11 * gcem::exp(-gcem::pow(0.015289 * (dj - 120.), 4))
                    + 1 * gcem::exp(-gcem::pow(4.6889e-03 * (dj - 400.), 6))));
            colorMap[j][2] = clampcast(255. *
                (gcem::exp(-gcem::pow(3.1101e-03 * (dj - 415.), 10))));
            }
            break;
        case ColorMap::cyan_magenta:
            colorMap[j][0] = clampcast(255. * (1.0 * gcem::exp(-gcem::pow(4.5 * (nval - 0.05), 2))
                + 1.00 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][1] = clampcast(255. * (0.95 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][2] = clampcast(255. * (0.9 * gcem::exp(-gcem::pow(4.5 * (nval - 0.05), 2))
                + 0.2 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
        }
    }
    return colorMap;
}

static constexpr std::array<std::array<std::array<uint8_t, 3>, 256>, 5> LweColorMaps{
        createColormap(ColorMap::greyscale), 
        createColormap(ColorMap::jet), 
        createColormap(ColorMap::jet_white), 
        createColormap(ColorMap::purple), 
        createColormap(ColorMap::cyan_magenta)};



class LweColor {
public:
    double r;
    double g;
    double b;
    double a;
    constexpr LweColor(double rIn, double gIn, double bIn, double aIn) 
        : r(rIn), g(gIn), b(bIn), a(aIn) {}
    [[nodiscard]] constexpr int rHex() const noexcept { return static_cast<int>(15 * r); }
    [[nodiscard]] constexpr int gHex() const noexcept { return static_cast<int>(15 * g); }
    [[nodiscard]] constexpr int bHex() const noexcept { return static_cast<int>(15 * b); }
    [[nodiscard]] constexpr int aHex() const noexcept { return static_cast<int>(15 * a); }
    void setCairo(cairo_t* cr) const { cairo_set_source_rgb(cr, r, g, b); }
    void setCairoA(cairo_t* cr) const { cairo_set_source_rgba(cr, r, g, b, a); }
};

//function to stream png data into a vector
static cairo_status_t cairoWritePNGtoVector(
    void *closure, 
    const unsigned char *data, 
    unsigned int length) {
    std::vector<uint8_t> *vector = static_cast<std::vector<uint8_t>*>(closure);
    vector->insert(vector->end(), data, data + length);
    return CAIRO_STATUS_SUCCESS;
}

template <typename T, std::size_t dimension>
struct DataSlice{
    T* data = nullptr;
    std::array<size_t,dimension> shape = {};
};

template <typename T>
class LweImage {
    public:
        int width = 0;
        int height = 0;
        std::optional<DataSlice<T,2>> data;
        std::optional<DataSlice<std::complex<T>,2>> complexData;
        ColorMap colorMap = ColorMap::purple;
        bool logScale = false;
        bool hasFullSizeRenderedImage = false;
        T logMin = 0;
        std::vector<uint8_t> pixels;
        size_t yDim(){
            if (data.has_value()) return data.value().shape[1];
            else if (complexData.has_value()) return complexData.value().shape[1];
            else return height;
        }
        size_t xDim(){
            if (data.has_value()) return data.value().shape[0];
            else if (complexData.has_value()) return complexData.value().shape[0];
            else return width;
        }
        void render() {
            if(hasFullSizeRenderedImage) return;
            int64_t plotSize = static_cast<int64_t>(width) * static_cast<int64_t>(height);
            std::vector<float> plotarr2(plotSize);
            if(data.has_value()){
                linearRemapDoubleToFloat(
                    data.value().data, 
                    static_cast<int>(data.value().shape[1]), 
                    static_cast<int>(data.value().shape[0]),
                    plotarr2.data(),
                    height,
                    width);
                createBitmapFromArray(width, height, plotarr2.data(), colorMap);
            }
            else if(complexData.has_value()){
                linearRemapZToLogFloatShift(
                    complexData.value().data, 
                    static_cast<int>(complexData.value().shape[1]),
                    static_cast<int>(complexData.value().shape[0]),
                    plotarr2.data(),
                    height,
                    width,
                    logMin);
                createBitmapFromArray(width, height, plotarr2.data(), colorMap);
            }
        }
    
        [[nodiscard]] constexpr T cModulusSquared(const std::complex<T>& x) const {
            return x.real() * x.real() + x.imag() * x.imag();
        }
    
        void linearRemapZToLogFloatShift(
            const std::complex<T>* A, 
            const int nax, 
            const int nay, 
            float* B, 
            const int nbx, 
            const int nby, 
            const T logMin) {
            const int div2 = nax / 2;
    #pragma omp parallel for num_threads(interfaceThreads)
            for (int i = 0; i < nbx; ++i) {
                float f = i * (nax / (float)nbx);
                int nx0 = minN((int)f, nax);
                nx0 -= div2 * ((nx0 >= div2) - (nx0 < div2));
                nx0 *= nay;
                for (int j = 0; j < nby; ++j) {
                    f = (j * (nay / (float)nby));
                    int ny0 = minN(nay, (int)f);
                    B[i * nby + j] = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
                }
            }
        }
    
        void linearRemapZToLogFloat(
            const std::complex<T>* A, 
            const int nax, 
            const int nay, 
            float* B, 
            const int nbx, 
            const int nby, 
            const T logMin) {
    #pragma omp parallel for num_threads(interfaceThreads)
            for (int i = 0; i < nbx; ++i) {
                float f = i * (nax / (float)nbx);
                int Ni = (int)f;
                int nx0 = nay * minN(Ni, nax);
                for (int j = 0; j < nby; ++j) {
                    f = (j * (nay / (float)nby));
                    int Nj = (int)f;
                    int ny0 = minN(nay, Nj);
                    B[i * nby + j] = (float)log10(cModulusSquared(A[ny0 + nx0]) + logMin);
                }
            }
        }
    
        void linearRemapDoubleToFloat(
            const T* A, 
            const int nax, 
            const int nay, 
            float* B, 
            const int nbx, 
            const int nby) {
    #pragma omp parallel for num_threads(interfaceThreads)
            for (int i = 0; i < nbx; ++i) {
                int Ni = (int)(i * (nax / (float)nbx));
                int nx0 = nay * minN(Ni, nax);
                for (int j = 0; j < nby; ++j) {
                    int Nj = (int)((j * (nay / (float)nby)));
                    int ny0 = minN(nay, Nj);
                    B[i * nby + j] = (float)A[ny0 + nx0];
                }
            }
        }
    
        void createBitmapFromArray(
            const int Nx,
            const int Ny,
            const float* data,
            const int cm){
            // creating input
            const int64_t Ntot = static_cast<int64_t>(Nx) * static_cast<int64_t>(Ny);
            pixels.resize(4 * Ntot);
    
            const std::array<std::array<uint8_t, 3>, 256> colorMap = LweColorMaps[cm];
            constexpr int stride = 4;
            //Find the image maximum and minimum
            float imin = data[0];
            float imax = data[0];
            for (int64_t i = 1; i < Ntot; ++i) {
                if (data[i] > imax) imax = data[i];
                if (data[i] < imin) imin = data[i];
            }
            if (cm == 4) {
                imax = maxN(imax, -imin);
                imin = minN(imin, -imax);
            }
    
            if (imin != imax) {
    #pragma omp parallel for num_threads(interfaceThreads)
                for (int p = 0; p < Ntot; p++) {
                    uint8_t currentValue = 
                        clampcast(255.0 * (data[p] - imin) / (imax - imin));
                    pixels[stride * p] = colorMap[currentValue][0];
                    pixels[stride * p + 1] = colorMap[currentValue][1];
                    pixels[stride * p + 2] = colorMap[currentValue][2];
                }
            }
        }
        void drawArrayAsBitmap(
            cairo_t* cr, 
            const int Nx, 
            const int Ny, 
            const float* data, 
            const int cm,
            const double x_offset = 0.0,
            const double y_offset = 0.0) {
            if (Nx * Ny == 0) return;
            createBitmapFromArray(Nx, Ny, data, cm);
            std::unique_lock GTKlock(GTKmutex);
            const int caiStride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, Nx);
            cairo_surface_t* cSurface = cairo_image_surface_create_for_data(
                pixels.data(),
                CAIRO_FORMAT_RGB24, 
                Nx, Ny, 
                caiStride);
            cairo_set_source_surface(cr, cSurface, x_offset, y_offset);
            cairo_paint(cr);
            cairo_surface_finish(cSurface);
            cairo_surface_destroy(cSurface);
        }

        void drawRenderedPixels(
            cairo_t* cr, 
            const double x_offset = 0.0,
            const double y_offset = 0.0,
            int drawingwidth = 0,
            int drawingheight = 0) {
            if (height * width == 0) return;
            bool reshapeImage = drawingwidth != 0 || drawingheight != 0;
            const int caiStride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, width);
            cairo_surface_t* cSurface = cairo_image_surface_create_for_data(
                pixels.data(),
                CAIRO_FORMAT_RGB24, 
                width, height, 
                caiStride);
            if(reshapeImage){
                cairo_save(cr);
                double x_scale = static_cast<double>(drawingwidth)/width;
                double y_scale = static_cast<double>(drawingheight)/height;
                cairo_scale(
                    cr, 
                    x_scale, 
                    y_scale);
                    cairo_set_source_surface(cr, cSurface, x_offset/x_scale, y_offset/y_scale);
                    cairo_paint(cr);
                    cairo_restore(cr);
                }
                else{
                    cairo_set_source_surface(cr, cSurface, x_offset, y_offset);
                    cairo_paint(cr);
                }

            cairo_surface_finish(cSurface);
            cairo_surface_destroy(cSurface);
        }

        std::string encodeBase64(const std::vector<uint8_t>& data) {
            static const char* base64_chars = 
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
            std::string encoded;
            size_t i = 0;
            while (i < data.size()) {
                uint32_t value = 0;
                size_t bytes = std::min(size_t(3), data.size() - i);
                for (size_t j = 0; j < bytes; ++j) {
                    value |= static_cast<uint32_t>(data[i + j]) << (16 - j * 8);
                }
                for (size_t j = 0; j < 4; ++j) {
                    if (j <= bytes) {
                        encoded += base64_chars[(value >> (18 - j * 6)) & 0x3F];
                    } else {
                        encoded += '=';
                    }
                }
                i += bytes;
            }
            return encoded;
        }


        std::string pngFromRenderedPixels(
            cairo_t* cr,
            const double x_offset = 0.0,
            const double y_offset = 0.0,
            const int width_in_svg = 320,
            const int height_in_svg = 240){
                std::vector<uint8_t> png;
                const int caiStride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, width);
                cairo_surface_t* cSurface = cairo_image_surface_create_for_data(
                    pixels.data(),
                    CAIRO_FORMAT_RGB24, 
                    width, 
                    height, 
                    caiStride);
                cairo_set_source_surface(cr, cSurface, x_offset, y_offset);
                cairo_paint(cr);
                cairo_surface_write_to_png_stream(cSurface, cairoWritePNGtoVector, &png);
                cairo_surface_finish(cSurface);
                cairo_surface_destroy(cSurface);
                auto pngString = encodeBase64(png);
                std::string tag = Sformat(
                    "<image x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" "
                    "preserveAspectRatio=\"none\" xlink:href=\"data:image/png;base64,{}\"/>", 
                    x_offset, 
                    y_offset, 
                    width_in_svg, 
                    height_in_svg, 
                    pngString);
                return tag;
        }
};

template<typename T>
class LwePlot {
public:
    bool makeSVG = false;
    bool markers = true;
    double width = 0;
    double height = 0;
    double fontSize = 14.0;
    std::vector<T*> data;
    std::vector<T*> dataX;
    std::optional<LweImage<T>*> image;
    std::string xLabel;
    std::string yLabel;
    bool logScale = false;
    T logMin = 0;
    int dataType = 0;
    T dx = 1.0;
    T x0 = 0.0;
    int64_t Npts = 0;
    T unitY = 1.0;
    std::optional<T> forcedYmin;
    std::optional<T> forcedYmax;
    std::optional<T> forcedXmin;
    std::optional<T> forcedXmax;
    LweColor axisColor = LweColor(0.8, 0.8, 0.8, 0.8);
    LweColor textColor = LweColor(0.8, 0.8, 0.8, 0.8);
    LweColor backgroundColor = LweColor(0.0, 0.0, 0.0, 0.0);
    std::vector<LweColor> lineColors;
    std::string SVGString;

    int plot(cairo_t* cr) {
        if (image.has_value()) Npts = 5;
        if (Npts == 0) return 1;
        std::unique_lock GTKlock(GTKmutex);
        int64_t iMin = 0;
        int64_t iMax = Npts;
        cairo_font_extents_t fe{};
        cairo_text_extents_t te{};
        cairo_set_font_size(cr, fontSize);
        cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_font_extents(cr, &fe);
        T x1, y1, x2, y2;

        double layoutTop = 0.0;
        double layoutBottom = 0.0;
        double layoutLeft = 0.0;
        double layoutRight = 0.0;

        //Fixed parameters affecting aesthetics
        double radius = 2;
        double lineWidth = 1.5;
        double axisSpaceX = 5.36 * fontSize;
        double axisSpaceY = 2.5 * fontSize;
        double axisLabelSpaceX = 1.5 * fontSize;

        //get limits and make the plotting arrays
        T maxY = -1.0e300;
        T minY = 0.0;
        T maxX = 0.0;
        T minX = 0.0;
        T currentY;
        T currentX;
        std::vector<double> xValues(Npts + 2, 0.0);
        if(!image.has_value()){
            for (int i = 0; i < Npts; ++i) {
                if (dataX.size()) currentX = (double)dataX[0][i];
                else { currentX = (double)(i * dx + x0); }
                if (i == 0) {
                    minX = currentX;
                    maxX = currentX;
                }
                xValues[i] = currentX;
                if (forcedXmin.has_value() && (currentX < forcedXmin.value())) {
                    iMin = i + 1;
                }
                if (forcedXmax.has_value() && (currentX > forcedXmax.value())) {
                    iMax = i;
                    break;
                }
                maxX = maxN(currentX, maxX);
                minX = minN(currentX, minX);
            }
            if (iMin >= iMax || iMin >= Npts) return -1;
            for(size_t j=0; j<data.size(); ++j){
                for (int64_t i = iMin; i < iMax; ++i) {
                    if (logScale) { currentY = (double)log10(data[j][i]); }
                    else { currentY = (double)data[j][i]; }
                    maxY = maxN(currentY, maxY);
                    minY = minN(currentY, minY);
                }
            }
    
            if (minY == maxY) {
                minY = -1;
                maxY = 1;
            }

            if (forcedYmin.has_value()) {
                minY = forcedYmin.value() * unitY;
                if (logScale) minY = log10(minY);
            }
            if (forcedYmax.has_value()) {
                maxY = forcedYmax.value() * unitY;
                if (logScale) maxY = log10(maxY);
            }
            if (forcedXmin.has_value()) {
                minX = forcedXmin.value();
            }
            if (forcedXmax.has_value()) {
                maxX = forcedXmax.value();
            }
        }
        else{
            //drawing an image: require input x and y limits
            if(forcedXmin.has_value() 
            && forcedXmin.has_value()
            && forcedYmin.has_value()
            && forcedYmax.has_value()){
                minX = forcedXmin.value();
                maxX = forcedXmax.value();
                maxY = forcedYmax.value() * unitY;
                minY = forcedYmin.value() * unitY;
            }
            else{
                minX = -1.0;
                maxX = 1.0;
                minY = -1.0;
                maxY = 1.0;
            }
        }

        //Tickmark labels
        constexpr int NyTicks = 3;
        std::string messageBuffer;
        const double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
        const double xTicks1[3] = { 
            minX + 0.25 * (maxX - minX), 
            minX + 0.5 * (maxX - minX), 
            minX + 0.75 * (maxX - minX) };

        //Start SVG file if building one
        auto constexpr SVGh = [&](const double x) {
            return static_cast<int>(15 * x);
        };
        if (makeSVG) {
            SVGString.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
                "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
                "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
            SVGString.append(Sformat("<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}"
                "\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink="
                "\"http://www.w3.org/1999/xlink\">\n",
                width, height, width, height));
            SVGString.append(Sformat("<rect fill=\"#{:x}{:x}{:x}\" stroke="
                "\"#000\" x=\"0\" y=\"0\" width=\"{}\" height=\"{}\"/>\n",
                SVGh(backgroundColor.r), SVGh(backgroundColor.g), SVGh(backgroundColor.b), width, height));
        }
        
        cairo_rectangle(cr, 0, 0, width, height);
        backgroundColor.setCairo(cr);
        cairo_fill(cr);
        width -= axisSpaceX;
        height -= axisSpaceY;
        if(image.has_value()){
            if(makeSVG){
                image.value()->width = image.value()->xDim();
                image.value()->height = image.value()->yDim();
                image.value()->render();
                std::string imagePNG = image.value()->pngFromRenderedPixels(cr,axisSpaceX,0.0,width,height);
                SVGString.append(imagePNG);
            }
            if(image.value()->hasFullSizeRenderedImage){
                image.value()->drawRenderedPixels(cr,axisSpaceX,0.0,width,height);
            } 
            else{
                image.value()->width = width;
                image.value()->height = height;
                image.value()->render();
                image.value()->drawRenderedPixels(cr,axisSpaceX,0.0);
            }

        }
        double scaleX = width / (maxX - minX);
        double scaleY = height / (maxY - minY);
        LweColor currentColor = textColor;
        currentColor = textColor;

        //lambdas for writing components of SVG file
        auto SVGstdline = [&]() {
            if (makeSVG)SVGString.append(Sformat("<line x1=\"{}\" y1=\"{}\" "
                "x2=\"{}\" y2=\"{}\" stroke=\"#{:x}{:x}{:x}\" stroke-width=\"{}\"/>\n", 
                x1, y1, x2, y2, 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth));
        };

        auto SVGstartPolyLine = [&]() {
            SVGString.append(Sformat("<polyline points=\""));
        };

        auto SVGendPolyLine = [&]() {
            SVGString.append(Sformat("\" stroke=\"#{:x}{:x}{:x}\" stroke-width="
                "\"{}\" fill=\"none\"/>\n", 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth));
        };

        auto SVGaddXYtoPolyLine = [&](const double a, const double b) {
            SVGString.append(Sformat("{},{} ", a, b));
        };

        auto SVGstdcircle = [&]() {
            if (makeSVG)SVGString.append(Sformat("<circle cx=\"{}\" cy=\"{}\" r=\"{}\""
                " stroke=\"none\" fill=\"#{:x}{:x}{:x}\" />\n", 
                x1, y1, radius, 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex()));
        };

        auto SVGstartgroup = [&]() {
            if (makeSVG)SVGString.append("<g>\n");
        };

        auto SVGendgroup = [&]() {
            if (makeSVG)SVGString.append("</g>\n");
        };

        auto SVGcentertext = [&]() {
            if (makeSVG)SVGString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}"
                "\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" "
                "text-anchor=\"middle\">\n{}\n</text>\n", 
                fontSize - 1, 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 
                0.5 * (layoutLeft + layoutRight), 
                0.5 * (layoutBottom + layoutTop - te.height), 
                messageBuffer));
        };

        auto SVGrighttext = [&]() {
            if (makeSVG){
                cairo_text_extents(cr, messageBuffer.c_str(), &te);
                SVGString.append(Sformat("<text font-family=\"Arial\" "
                    "font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" "
                    "y=\"{}\">\n{}\n</text>\n", 
                    fontSize - 1, 
                    currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 
                    layoutRight - te.x_advance - 0.25*fontSize, 0.5 * (layoutBottom + layoutTop - te.height), 
                    messageBuffer));
            }
        };

        cairo_set_line_width(cr, lineWidth);
        auto cairoLine = [&]() {
            currentColor.setCairo(cr);
            cairo_move_to(cr, x1, y1);
            cairo_line_to(cr, x2, y2);
            cairo_stroke(cr);
        };
        auto cairoRightText = [&]() {
            currentColor.setCairo(cr);
            cairo_text_extents(cr, messageBuffer.c_str(), &te);
            cairo_move_to(cr, layoutRight - te.x_advance - 0.25 * fontSize, 
                0.5 * (layoutBottom + layoutTop - te.height));
            cairo_show_text(cr, messageBuffer.c_str());
        };
        auto cairoCenterText = [&]() {
            currentColor.setCairo(cr);
            cairo_text_extents(cr, messageBuffer.c_str(), &te);
            cairo_move_to(cr, 
                0.5 * (layoutLeft + layoutRight - te.x_advance), 
                0.5 * (layoutBottom + layoutTop - te.height));
            cairo_show_text(cr, messageBuffer.c_str());
        };

        auto cairoVerticalText = [&]() {
            currentColor.setCairo(cr);
            cairo_text_extents(cr, messageBuffer.c_str(), &te);
            cairo_move_to(cr, 0.0, height);
            cairo_rotate(cr, -vPi<double>() / 2);
            cairo_rel_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), fontSize);
            cairo_show_text(cr, messageBuffer.c_str());
            cairo_rotate(cr, vPi<double>() / 2);
        };

        currentColor = textColor;
        //y-tick text labels
        for (int i = 0; i < NyTicks; ++i) {
            double ytVal = yTicks1[i] / unitY;
            if (logScale) {
                switch (i) {
                case 0:
                    ytVal = pow(10.0, maxY) / unitY;
                    break;
                case 1:
                    ytVal = pow(10.0, minY + 0.5 * (maxY - minY)) / unitY;
                    break;
                case 2:
                    ytVal = pow(10.0, minY) / unitY;
                }
            }
            if (abs(ytVal) > 10.0 || abs(ytVal) < 0.01) {
                messageBuffer = Sformat("{:.2g}", ytVal);
            }
            else {
                messageBuffer = Sformat("{:4.4g}", ytVal);
            }
            if(messageBuffer=="0.0e+00") messageBuffer = "0";
            layoutLeft = axisLabelSpaceX;
            layoutTop = (i * (0.5 * (height)));
            if (i == 2) layoutTop -= 8.0f;
            if (i == 1) layoutTop -= 6.0f;
            layoutBottom = layoutTop + axisSpaceY;
            layoutRight = axisSpaceX;
            cairoRightText();
            SVGrighttext();
        }
        //y-axis name
        if (yLabel.size()) {
            messageBuffer.assign(yLabel);
            layoutLeft = 0;
            layoutTop = height;
            layoutBottom = height + axisSpaceY;
            layoutRight = height;

            cairoVerticalText();
            if (makeSVG)SVGString.append(Sformat("<text font-family=\"Arial\" font-size"
                "=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" text-anchor="
                "\"middle\" transform=\"translate({}, {}) rotate(-90)\">\n{}\n</text>\n", 
                fontSize, 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 
                0.5 * (layoutLeft + layoutRight), 
                layoutTop + fontSize, 
                -(layoutLeft + layoutRight), 
                height, messageBuffer));
        }

        //x-axis name
        if (xLabel.size()) {
            layoutLeft = axisSpaceX;
            layoutTop = height + 2.8 * fontSize;
            layoutBottom = height + axisSpaceY;
            layoutRight = axisSpaceX + width;
            messageBuffer.assign(xLabel);
            cairoCenterText();
            SVGcentertext();
        }

        //x-axis tick labels
        for (int i = 0; i < 3; ++i) {
            messageBuffer.assign(Sformat("{}", (int)round(xTicks1[i])));
            if(messageBuffer=="0.0e+00") messageBuffer = "0";
            layoutLeft = axisSpaceX + 0.25 * width * ((int64_t)(i)+1) - 0.5 * axisSpaceX;
            layoutTop = height + 3;
            layoutBottom = height + axisSpaceY;
            layoutRight = layoutLeft + axisSpaceX;

            cairoCenterText();
            SVGcentertext();
        }

        //Draw axes and tickmarks
        SVGstartgroup();
        currentColor = axisColor;
        x1 = axisSpaceX;
        y1 = height;
        x2 = scaleX * (maxX - minX) + axisSpaceX;
        y2 = height;
        cairoLine();
        SVGstdline();
        x1 = axisSpaceX;
        y1 = height;
        x2 = x1;
        y2 = 0.0;
        cairoLine();
        SVGstdline();
        for (int i = 0; i < 2; ++i) {
            y1 = height - scaleY * (yTicks1[i] - minY);
            y2 = y1;
            x1 = axisSpaceX;
            x2 = x1 + 10.;
            cairoLine();
            SVGstdline();
        }
        for (int i = 0; i < 3; ++i) {
            y1 = height;
            y2 = y1 - 10.;
            x1 = axisSpaceX + scaleX * (xTicks1[i] - minX);
            x2 = x1;
            cairoLine();
            SVGstdline();
        }
        SVGendgroup();

        //Lambdas for plotting a line
        std::vector<double> scaledX(Npts);
        std::vector<double> scaledY(Npts);
        auto getNewScaledXY = [&](const std::vector<double>& xValues, const T* y) {
            if (logScale) {
#pragma omp parallel for num_threads(interfaceThreads)
                for (int i = 0; i < Npts; i++) {
                    scaledX[i] = scaleX * (xValues[i] - minX) + axisSpaceX;
                    scaledY[i] = height - scaleY * (static_cast<double>(log10(y[i])) - static_cast<double>(minY));
                    if (std::isnan(scaledY[i]) || std::isinf(scaledY[i])) scaledY[i] = maxX * 2.;
                }
            }
            else {
#pragma omp parallel for num_threads(interfaceThreads)
                for (int i = 0; i < Npts; i++) {
                    scaledX[i] = scaleX * (xValues[i] - minX) + axisSpaceX;
                    scaledY[i] = height - scaleY * (static_cast<double>(y[i]) - static_cast<double>(minY));
                    if (std::isnan(scaledY[i]) || std::isinf(scaledY[i])) scaledY[i] = maxX * 2.;
                }
            }
        };

        auto plotCairoDots = [&]() {
            currentColor.setCairo(cr);
            for (int64_t i = iMin; i < iMax; ++i) {
                if (scaledY[i] <= height) {
                    cairo_arc(cr, scaledX[i], scaledY[i], radius, 0, twoPi<double>());
                    cairo_fill(cr);
                }
            }
        };
        //I would think it would be faster to only call cairo_fill() 
        //at the end, but this requires calling cairo_cloase_path()
        //in the loop, which seems to be even slower....
        auto plotCairoPolyline = [&]() {
            currentColor.setCairo(cr);
            cairo_move_to(cr, scaledX[iMin], scaledY[iMin]);
            for (int64_t i = iMin + 1; i < iMax; ++i) {
                if (scaledY[i - 1] <= height) {
                    if (scaledY[i] <= height) {
                        cairo_line_to(cr, scaledX[i], scaledY[i]);
                    }
                    else {
                        cairo_line_to(cr, 
                            scaledX[i - 1] + 
                            (height - scaledY[i - 1]) 
                            / ((scaledY[i] - scaledY[i - 1]) / (scaledX[i] - scaledX[i - 1])), 
                            height);
                    }
                }
                else if (scaledY[i] <= height) {
                    cairo_move_to(cr, 
                        scaledX[i - 1] + 
                        (height - scaledY[i - 1]) 
                        / ((scaledY[i] - scaledY[i - 1]) 
                            / (scaledX[i] - scaledX[i - 1])), 
                        height);
                    cairo_line_to(cr, scaledX[i], scaledY[i]);
                }
            }
            cairo_stroke(cr);
        };

        auto plotSVGPolyline = [&]() {
            bool lineOn = false;
            x2 = scaledX[iMin];
            y2 = scaledY[iMin];
            if (y2 <= height) {
                SVGstartPolyLine();
                SVGaddXYtoPolyLine(x2, y2);
                lineOn = true;
            }
            for (int64_t i = iMin + 1; i < iMax; ++i) {
                x1 = x2;
                x2 = scaledX[i];
                y1 = y2;
                y2 = scaledY[i];
                if (y1 <= height) {
                    if (y2 <= height) {
                        if (!lineOn) {
                            SVGstartPolyLine();
                            lineOn = true;
                        }
                        SVGaddXYtoPolyLine(x2, y2);
                    }
                    else {
                        x1 += (height - y1) / ((y2 - y1) / (x2 - x1));
                        SVGaddXYtoPolyLine(x1, height);
                        SVGendPolyLine();
                        lineOn = false;
                    }
                }
                else if (y2 <= height) {
                    if (!lineOn) {
                        x1 += (height - y1) / ((y2 - y1) / (x2 - x1));
                        SVGstartPolyLine();
                        SVGaddXYtoPolyLine(x1, height);
                        lineOn = true;
                    }
                    SVGaddXYtoPolyLine(x2, y2);
                }
                else if (lineOn) {
                    SVGendPolyLine();
                    lineOn = false;
                }
            }
            if (lineOn) {
                SVGendPolyLine();
            }
        };

        auto plotSVGDots = [&]() {
            SVGstartgroup();
            for (int64_t i = iMin; i < iMax - 1; ++i) {
                x1 = scaledX[i];
                y1 = scaledY[i];
                if (y1 <= height && !std::isnan(y1) && !std::isinf(y1)) {
                    SVGstdcircle();
                }
            }
            SVGendgroup();
        };

        //Plot the lines
        if(!image.has_value()){
            for(size_t i = 0; i<data.size(); ++i){
                getNewScaledXY(xValues, data[i]);
                currentColor = lineColors[i % lineColors.size()];
                plotCairoPolyline();
                if (markers)plotCairoDots();
                if (makeSVG) {
                    plotSVGPolyline();
                    if (markers)plotSVGDots();
                }
            }
        }

        if (makeSVG) {
            SVGString.append("</svg>");
        }
        return 0;
    }
};

