#pragma once
#include <gtk/gtk.h>
#include <sstream>
#include <fstream>
#include <complex>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <gcem.hpp>
#include "../LightwaveExplorerDevices/LightwaveExplorerHelpers.h"

//GLOBAL VARIABLE: GTK MUTEX
std::mutex GTKmutex;
//temporary set of macros until std::format is on all platforms
#if defined __linux__
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

constexpr std::array<std::array<uint8_t, 3>, 256> createColormap(const int cm) {
    std::array<std::array<uint8_t, 3>, 256> colorMap{};
    const double oneOver255 = 1.0f / 255.0f;
    for (int j = 0; j < 256; ++j) {
        double nval;
        switch (cm) {
        case 0:
            colorMap[j][0] = static_cast<uint8_t>(j);
            colorMap[j][1] = static_cast<uint8_t>(j);
            colorMap[j][2] = static_cast<uint8_t>(j);
            break;
        case 1:
            nval = j * oneOver255;
            colorMap[j][0] = static_cast<uint8_t>(255. * gcem::cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = static_cast<uint8_t>(255. * gcem::cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = static_cast<uint8_t>(255. * gcem::sin(vPi<double>() * nval / 2.));
            break;
        case 2:
            nval = j * oneOver255;
            colorMap[j][0] = static_cast<uint8_t>(255. * gcem::cos(vPi<double>() * nval / 2.));
            colorMap[j][1] = static_cast<uint8_t>(255. * gcem::cos(vPi<double>() * (nval - 0.5)));
            colorMap[j][2] = static_cast<uint8_t>(255. * gcem::sin(vPi<double>() * nval / 2.));
            if (nval < 0.02) {
                colorMap[j][0] = 255;
                colorMap[j][1] = 128;
                colorMap[j][2] = 128;
            }
            if (nval < 0.01) {
                colorMap[j][0] = 255;
                colorMap[j][1] = 255;
                colorMap[j][2] = 255;
            }
            break;
        case 3:
            nval = 255 * (j * oneOver255);
            colorMap[j][0] = static_cast<uint8_t>(255. *
                (0.998 * gcem::exp(-gcem::pow(7.7469e-03 * (nval - 160.), 6))
                    + 0.22 * gcem::exp(-gcem::pow(0.016818 * (nval - 305.), 4))));
            colorMap[j][1] = static_cast<uint8_t>(255. *
                (0.022 * gcem::exp(-gcem::pow(0.042045 * (nval - 25.), 4))
                    + 0.11 * gcem::exp(-gcem::pow(0.015289 * (nval - 120.), 4))
                    + 1 * gcem::exp(-gcem::pow(4.6889e-03 * (nval - 400.), 6))));
            colorMap[j][2] = static_cast<uint8_t>(255. *
                (gcem::exp(-gcem::pow(3.1101e-03 * (nval - 415), 10))));
            break;
        case 4:
            nval = j * oneOver255;
            colorMap[j][0] = static_cast<uint8_t>(255. * (1.0 * gcem::exp(-gcem::pow(4.5 * (nval - 0.05), 2))
                + 1.00 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][1] = static_cast<uint8_t>(255. * (0.95 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
            colorMap[j][2] = static_cast<uint8_t>(255. * (0.9 * gcem::exp(-gcem::pow(4.5 * (nval - 0.05), 2))
                + 0.2 * gcem::exp(-gcem::pow(3.5 * (nval - 1.05), 2))));
        }
    }
    return colorMap;
}

static constexpr std::array<std::array<std::array<uint8_t, 3>, 256>, 5> LweColorMaps{
        createColormap(0), 
        createColormap(1), 
        createColormap(2), 
        createColormap(3), 
        createColormap(4)};

class LweColor {
public:
    double r;
    double g;
    double b;
    double a;
    LweColor() : r(0.0), g(0.0), b(0.0), a(0.0) {
    }
    LweColor(double rIn, double gIn, double bIn, double aIn) {
        r = rIn;
        g = gIn;
        b = bIn;
        a = aIn;
    }
    ~LweColor() {};
    [[nodiscard]] constexpr int rHex() const noexcept { return (int)(15 * r); }
    [[nodiscard]] constexpr int gHex() const noexcept { return (int)(15 * g); }
    [[nodiscard]] constexpr int bHex() const noexcept { return (int)(15 * b); }
    [[nodiscard]] constexpr int aHex() const noexcept { return (int)(15 * a); }
    void setCairo(cairo_t* cr) { cairo_set_source_rgb(cr, r, g, b); }
    void setCairoA(cairo_t* cr) { cairo_set_source_rgba(cr, r, g, b, a); }
};

class LwePlot {
    bool makeSVG = false;
public:
    bool markers = true;
    double width = 0;
    double height = 0;
    double* data = nullptr;
    int ExtraLines = 0;
    double* data2 = nullptr;
    double* data3 = nullptr;
    double* data4 = nullptr;
    double* dataX = nullptr;
    const char* xLabel = nullptr;
    const char* yLabel = nullptr;
    bool hasDataX = false;
    bool logScale = false;
    double logMin = 0;
    int dataType = 0;
    double dx = 1.0;
    double x0 = 0.0;
    int64_t Npts = 0;
    double unitY = 1.0;
    bool forceYmin = false;
    double forcedYmin = 0.0;
    bool forceYmax = false;
    double forcedYmax = 0.0;
    bool forceXmin = false;
    double forcedXmin = 0.0;
    bool forceXmax = false;
    double forcedXmax = 0.0;
    LweColor axisColor = LweColor(0.5, 0.5, 0.5, 0.5);
    LweColor textColor = LweColor(0.8, 0.8, 0.8, 0.8);
    LweColor backgroundColor = LweColor(0.0, 0.0, 0.0, 0.0);
    LweColor color = LweColor(1, 1, 1, 1);
    LweColor color2 = LweColor(1, 1, 1, 1);
    LweColor color3 = LweColor(1, 1, 1, 1);
    LweColor color4 = LweColor(1, 1, 1, 1);
    std::string SVGString;
    std::string SVGPath;

    int plot(cairo_t* cr) {
        if (Npts == 0) return 1;
        std::unique_lock GTKlock(GTKmutex);
        if (SVGPath.length() > 5) makeSVG = true;
        int64_t iMin = 0;
        int64_t iMax = Npts;
        cairo_font_extents_t fe{};
        cairo_text_extents_t te{};
        double fontSize = 14.0;
        cairo_set_font_size(cr, fontSize);
        cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        cairo_font_extents(cr, &fe);
        double x1, y1, x2, y2;

        double layoutTop = 0.0;
        double layoutBottom = 0.0;
        double layoutLeft = 0.0;
        double layoutRight = 0.0;

        //Fixed parameters affecting aesthetics
        double radius = 2;
        double lineWidth = 1.5;
        double axisSpaceX = 75.0;
        double axisSpaceY = 35.0;
        double axisLabelSpaceX = 21.0;

        //get limits and make the plotting arrays
        double maxY = -1.0e300;
        double minY = 0.0;
        double maxX = 0.0;
        double minX = 0.0;
        double currentY;
        double currentX;
        std::vector<double> xValues(Npts + 2, 0.0);
        for (int i = 0; i < Npts; ++i) {
            if (hasDataX) currentX = (double)dataX[i];
            else { currentX = (double)(i * dx + x0); }
            if (i == 0) {
                minX = currentX;
                maxX = currentX;
            }
            xValues[i] = currentX;
            if (forceXmin && (currentX < forcedXmin)) {
                iMin = i + 1;
            }
            if (forceXmax && (currentX > forcedXmax)) {
                iMax = i;
                break;
            }
            maxX = maxN(currentX, maxX);
            minX = minN(currentX, minX);
        }
        if (iMin >= iMax || iMin >= Npts) return -1;
        for (int64_t i = iMin; i < iMax; ++i) {
            if (logScale) { currentY = (double)log10(data[i]); }
            else { currentY = (double)data[i]; }
            maxY = maxN(currentY, maxY);
            minY = minN(currentY, minY);
            if (ExtraLines > 0) {
                if (logScale) { currentY = (double)log10(data2[i]); }
                else { currentY = (double)data2[i]; }
                maxY = maxN(currentY, maxY);
                minY = minN(currentY, minY);
            }
        }

        if (minY == maxY) {
            minY = -1;
            maxY = 1;
        }
        if (forceYmin) {
            minY = (double)forcedYmin * unitY;
            if (logScale) minY = log10(minY);
        }
        if (forceYmax) {
            maxY = (double)forcedYmax * unitY;
            if (logScale) maxY = log10(maxY);
        }
        if (forceXmin) {
            minX = (double)forcedXmin;
        }
        if (forceXmax) {
            maxX = (double)forcedXmax;
        }

        //Tickmark labels
        int NyTicks = 3;
        std::string messageBuffer;
        double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
        double xTicks1[3] = { 
            minX + 0.25 * (maxX - minX), 
            minX + 0.5 * (maxX - minX), 
            minX + 0.75 * (maxX - minX) };

        //Start SVG file if building one
        auto SVGh = [&](double x) {
            return (int)(15 * x);
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
                SVGh(0.0f), SVGh(0.0f), SVGh(0.0f), width, height));
        }

        cairo_rectangle(cr, 0, 0, width, height);
        backgroundColor.setCairo(cr);
        cairo_fill(cr);
        width -= axisSpaceX;
        height -= axisSpaceY;
        double scaleX = width / ((double)(maxX - minX));
        double scaleY = height / ((double)(maxY - minY));
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

        auto SVGaddXYtoPolyLine = [&](double& a, double& b) {
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

        auto SVGlefttext = [&]() {
            if (makeSVG)SVGString.append(Sformat("<text font-family=\"Arial\" "
                "font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" "
                "y=\"{}\">\n{}\n</text>\n", 
                fontSize - 1, 
                currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 
                layoutLeft, layoutTop + fontSize, 
                messageBuffer));
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
            cairo_move_to(cr, layoutRight - te.width - 3, 
                0.5 * (layoutBottom + layoutTop - te.height));
            cairo_show_text(cr, messageBuffer.c_str());
        };
        auto cairoCenterText = [&]() {
            currentColor.setCairo(cr);
            cairo_text_extents(cr, messageBuffer.c_str(), &te);
            cairo_move_to(cr, 
                0.5 * (layoutLeft + layoutRight - te.width), 
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
                messageBuffer = Sformat("{:.1e}", ytVal);
            }
            else {
                messageBuffer = Sformat("{:4.4f}", ytVal);
            }

            layoutLeft = axisLabelSpaceX;
            layoutTop = (double)(i * (0.5 * (height)));
            if (i == 2) layoutTop -= 8.0f;
            if (i == 1) layoutTop -= 6.0f;
            layoutBottom = layoutTop + axisSpaceY;
            layoutRight = axisSpaceX;


            cairoRightText();
            SVGlefttext();
        }
        //y-axis name
        if (yLabel != NULL) {
            messageBuffer = std::string(yLabel);
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
        if (xLabel != NULL) {
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
            layoutLeft = (double)(axisSpaceX + 0.25 * width * ((int64_t)(i)+1) - axisSpaceX / 2);
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
            y1 = (double)(height - scaleY * (yTicks1[i] - minY));
            y2 = y1;
            x1 = axisSpaceX;
            x2 = x1 + 10;
            cairoLine();
            SVGstdline();
        }
        for (int i = 0; i < 3; ++i) {
            y1 = height;
            y2 = y1 - 10;
            x1 = (double)(axisSpaceX + scaleX * (xTicks1[i] - minX));
            x2 = x1;
            cairoLine();
            SVGstdline();
        }
        SVGendgroup();

        //Lambdas for plotting a line
        //currently dots are always there but has been refactored to make 
        //it easier to turn them off if I want.
        //I would think it would be faster to only call cairo_fill() 
        //at the end, but this requires calling cairo_cloase_path()
        //in the loop, which seems to be even slower....

        std::vector<double> scaledX(Npts);
        std::vector<double> scaledY(Npts);
        auto getNewScaledXY = [&](const std::vector<double>& xValues, const double* y) {
            if (logScale) {
#pragma omp parallel for num_threads(interfaceThreads)
                for (int i = 0; i < Npts; i++) {
                    scaledX[i] = scaleX * (xValues[i] - minX) + axisSpaceX;
                    scaledY[i] = height - scaleY * (static_cast<double>(log10(y[i])) - static_cast<double>(minY));
                    if (isnan(scaledY[i]) || isinf(scaledY[i])) scaledY[i] = maxX * 2;
                }
            }
            else {
#pragma omp parallel for num_threads(interfaceThreads)
                for (int i = 0; i < Npts; i++) {
                    scaledX[i] = scaleX * (xValues[i] - minX) + axisSpaceX;
                    scaledY[i] = height - scaleY * (static_cast<double>(y[i]) - static_cast<double>(minY));
                    if (isnan(scaledY[i]) || isinf(scaledY[i])) scaledY[i] = maxX * 2;
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
                if (y1 <= height && !isnan(y1) && !isinf(y1)) {
                    SVGstdcircle();
                }
            }
            SVGendgroup();
        };

        //Optional overlay curves
        if (ExtraLines > 0 && data2 != nullptr) {
            getNewScaledXY(xValues, data2);
            currentColor = color2;
            plotCairoPolyline();
            if (markers)plotCairoDots();
            if (makeSVG) {
                plotSVGPolyline();
                if (markers)plotSVGDots();
            }
        }
        if (ExtraLines > 1 && data2 != nullptr) {
            getNewScaledXY(xValues, data3);
            currentColor = color3;
            plotCairoPolyline();
            if (markers)plotCairoDots();
            if (makeSVG) {
                plotSVGPolyline();
                if (markers)plotSVGDots();
            }
        }
        if (ExtraLines > 2 && data2 != nullptr) {
            getNewScaledXY(xValues, data4);
            currentColor = color4;
            plotCairoPolyline();
            if (markers)plotCairoDots();
            if (makeSVG) {
                plotSVGPolyline();
                if (markers)plotSVGDots();
            }
        }

        //Plot the main line
        getNewScaledXY(xValues, data);
        currentColor = color;
        plotCairoPolyline();
        if (markers)plotCairoDots();
        if (makeSVG) {
            plotSVGPolyline();
            if (markers)plotSVGDots();
        }

        if (makeSVG) {
            SVGString.append("</svg>");
            std::ofstream fs(SVGPath);
            fs << SVGString;
        }
        return 0;
    }
};

class LweImage {
public:
    int width = 0;
    int height = 0;
    double* data = nullptr;
    int64_t dataXdim = 0;
    int64_t dataYdim = 0;
    std::complex<double>* complexData = nullptr;
    int colorMap = 4;
    bool logScale = false;
    double logMin = 0;
    int dataType = 0;

    void imagePlot(cairo_t* cr) {

        int dx = width;
        int dy = height;

        int64_t plotSize = (int64_t)dx * (int64_t)dy;
        float* plotarr2 = new float[plotSize];
        switch (dataType) {
        case 0:
            if (data == nullptr) break;
            linearRemapDoubleToFloat(
                data, 
                (int)dataYdim, 
                (int)dataXdim, 
                plotarr2, 
                (int)dy, 
                (int)dx);
            drawArrayAsBitmap(cr, width, height, plotarr2, colorMap);
            break;
        case 1:
            if (complexData == nullptr) break;
            linearRemapZToLogFloatShift(
                complexData, 
                (int)dataYdim, 
                (int)dataXdim, 
                plotarr2, 
                (int)dy, 
                (int)dx, 
                logMin);
            drawArrayAsBitmap(cr, width, height, plotarr2, colorMap);
            break;
        }

        delete[] plotarr2;
    }

    [[nodiscard]] constexpr double cModulusSquared(const std::complex<double>& x) {
        return x.real() * x.real() + x.imag() * x.imag();
    }

    void linearRemapZToLogFloatShift(
        const std::complex<double>* A, 
        const int nax, 
        const int nay, 
        float* B, 
        const int nbx, 
        const int nby, 
        const double logMin) {
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
        const std::complex<double>* A, 
        const int nax, 
        const int nay, 
        float* B, 
        const int nbx, 
        const int nby, 
        const double logMin) {
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
        const double* A, 
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

    void drawArrayAsBitmap(
        cairo_t* cr, 
        const int Nx, 
        const int Ny, 
        const float* data, 
        const int cm) {
        if (Nx * Ny == 0) return;
        std::unique_lock GTKlock(GTKmutex);
        // creating input
        const int64_t Ntot = Nx * Ny;
        uint8_t* pixels = new uint8_t[4 * Ntot]();
        if (pixels == nullptr) return;

        const std::array<std::array<uint8_t, 3>, 256> colorMap = LweColorMaps[cm];
        const int stride = 4;
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
                    (uint8_t)(255 * (data[p] - imin) / (imax - imin));
                pixels[stride * p + 0] = colorMap[currentValue][0];
                pixels[stride * p + 1] = colorMap[currentValue][1];
                pixels[stride * p + 2] = colorMap[currentValue][2];
            }
        }
        const int caiStride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, Nx);
        cairo_surface_t* cSurface = cairo_image_surface_create_for_data(
            pixels, 
            CAIRO_FORMAT_RGB24, 
            Nx, Ny, 
            caiStride);
        cairo_set_source_surface(cr, cSurface, 0, 0);
        cairo_paint(cr);
        cairo_surface_finish(cSurface);
        cairo_surface_destroy(cSurface);
        delete[] pixels;
    }
};

class LweGuiElement {
public:
    GtkWidget* label{};
    GtkWidget* elementHandle{};
    int _x{};
    int _y{};
    int _width{};
    int _height{};
    bool isAttached=false;
    GtkWidget* _grid{};

    void setPosition(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        if (_grid)gtk_grid_remove(GTK_GRID(_grid), elementHandle);
        _grid = grid;
        _x = x;
        _y = y;
        _width = width;
        _height = height;
        gtk_grid_attach(GTK_GRID(_grid), elementHandle, _x, _y, _width, _height);
    }
    void setLabel(int x, int y, const char* labelText) {
        std::unique_lock GTKlock(GTKmutex);
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), 45);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, _x + x, _y + y, 6, 1);
    }
    void setLabel(int x, int y, const char* labelText, int characters, int grids) {
        std::unique_lock GTKlock(GTKmutex);
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), characters);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, _x + x, _y + y, grids, 1);
    }
    void squeeze() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_halign(elementHandle, GTK_ALIGN_END);
    }
    void setTooltip(const char* tooltipText) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_tooltip_text(elementHandle, tooltipText);
    }
};

class LweTextBox : public LweGuiElement {
    int getNumberOfDecimalsToDisplay(double in, bool isExponential) {
        if (in == 0) return 0;
        in = abs(in);
        int digits = -1;
        int logValue = (int)floor(log10(in));
        in /= pow((double)10.0, logValue);
        while (digits < 15 && in > 1e-3 && in < 9.999) {
            in -= (int)in;
            in *= 10;
            digits++;
        }
        if (isExponential) {
            return maxN(0, digits);
        }
        return maxN(0, digits - logValue);
    }
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_entry_new();
        gtk_widget_set_halign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), 7);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }

    double valueDouble() {
        std::unique_lock GTKlock(GTKmutex);
        double sdata;
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> sdata;
            return sdata;
        }
        else {
            return 0.;
        }
    }

    int valueInt() {
        std::unique_lock GTKlock(GTKmutex);
        int sdata;
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> sdata;
            return sdata;
        }
        else {
            return 0;
        }
    }

    void valueToPointer(int* sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(int64_t* sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(double* sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToTwoPointers(double* sdata, double* sdata2) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        char c;
        *sdata2 = 0.0;
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata >> c >> *sdata2;
        }
    }

    void valueToTwoPointers(double multiplier, double* sdata, double* sdata2) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        char c;
        *sdata2 = 0.0;
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata >> c >> *sdata2;
            *sdata *= multiplier;
            *sdata2 *= multiplier;
        }
    }

    void valueToPointer(double multiplier, double* sdata) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
            *sdata *= multiplier;
        }
    }

    void setMaxCharacters(int charLimit) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), charLimit);
    }

    void setToDouble(double in) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = Sformat(std::string_view("{:g}"), in);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }
    template<typename... Args> void overwritePrint(std::string_view format, Args&&... args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = Svformat(format, Smake_format_args(args...));
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }

    void copyBuffer(char* destination, int64_t maxLength) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        std::string s(gtk_entry_buffer_get_text(buf));
        if (s.length() > 0) {
            s.append("\0");
            s.copy(destination, maxLength);
        }
        else {
            s.assign("None\0");
            s.copy(destination, maxLength);
        }
    }

    void copyBuffer(std::string& destination) {
        std::unique_lock GTKlock(GTKmutex);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));

        if (gtk_entry_buffer_get_length(buf) > 0) {
            std::string s(gtk_entry_buffer_get_text(buf));
            destination = s;
        }
        else {
            std::string s("None.");
            destination = s;
        }
    }

};

gboolean scrollTextViewToEndHandler(gpointer data) {
    std::unique_lock GTKlock(GTKmutex);
    GtkAdjustment* adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(data));
    gtk_adjustment_set_value(adjustment, gtk_adjustment_get_upper(adjustment));
    return false;
}

void formatSequenceEscapeAngleBrackets(std::string& s) {
    for (size_t i = 0; i < s.length(); ++i) {
        //find angle brackets signifying comments and escape
        //with &lt; or &gt; so they're not interpreted as
        //pango markup
        if (s[i] == '<') {
            s.erase(i, 1);
            s.insert(i, "&lt;");
            i += 3;
        }
        if (s[i] == '>') {
            s.erase(i, 1);
            s.insert(i, "&gt;");
            i += 3;
        }
    }
}

gboolean formatSequenceBuffer(gpointer data) {
    GtkTextBuffer* buf = GTK_TEXT_BUFFER(data);
    GtkTextIter start;
    GtkTextIter stop;
    GtkTextIter current;
    GtkTextIter currentTarget;
    gtk_text_buffer_get_start_iter(buf, &start);
    gtk_text_buffer_get_end_iter(buf, &stop);
    gtk_text_buffer_remove_all_tags(buf, &start, &stop);
    std::string s(gtk_text_buffer_get_text(buf, &start, &stop, false));
    std::vector<std::string> functionList{ "for",
        "plasma",
        "nonlinear",
        "linear",
        "sphericalMirror",
        "parabolicMirror",
        "init",
        "default",
        "rotate",
        "set",
        "plasmaReinject",
        "save",
        "savePlasma",
        "fresnelLoss",
        "aperture",
        "farFieldAperture",
        "energy",
        "filter",
        "lorentzian",
        "addPulse",
        "fdtd2d",
        "fdtd"
    };

    auto applyTag = [&](const char* tag, size_t a, size_t b) {
        current = start;
        currentTarget = start;
        gtk_text_iter_forward_chars(&current, (int)a);
        gtk_text_iter_forward_chars(&currentTarget, (int)b);
        gtk_text_buffer_apply_tag_by_name(buf, tag, &current, &currentTarget);
    };

    for (size_t i = 0; i < s.length(); ++i) {

        //anything enclosed between angle brackets is commented
        if (s[i] == '<') {
            size_t b = s.find_first_of('>', i);
            if (b == std::string::npos)break;
            applyTag("comment", i, b + 1);
            i = b + 1;
            if (i >= s.length())break;
        }

        //color function names and arguments
        if (s[i] == '(') {
            size_t nameStart = 0;
            size_t close = 0;
            //go backwards from (, amd upon encountering a newline, space, end of
            //another function, comment, or beginning of the buffer, check if the string
            //spanning that is in the functions list.
            //color it if it is.
            for (auto j = i; j > 0; --j) {
                if (j - 1 == 0 
                    || s[j - 1] == ' ' 
                    || s[j - 1] == '\n' 
                    || s[j - 1] == ')' 
                    || s[j - 1] == '>') {
                    nameStart = j - ((j - 1) == 0);
                    if (std::find(
                        std::begin(functionList),
                        std::end(functionList),
                        s.substr(nameStart, i - nameStart))
                        != std::end(functionList))
                        applyTag("function", j - 1, i);
                    break;
                }
            }
            //argument

            //if the ( isn't closed properly, paint it as an error
            close = s.find_first_of(')', i);
            nameStart = s.find_first_of('(', i + 1);
            if (close == std::string::npos || close > nameStart) {
                applyTag("error", i, i + 1);
            }
            //if it's closed, paint ( and ) 
            //and paint special variables in argument
            else {
                applyTag("parenthesis", i, i + 1);
                applyTag("parenthesis", close, close + 1);
                for (; i < close; ++i) {
                    if (s[i] == ' ') ++i;
                    if (s[i] == 'd') {
                        applyTag("delegate", i, i + 1);
                    }
                    if (s[i] == 'v') {
                        applyTag("variable", i, i + 3);
                    }
                    if (s[i] == 'i') {
                        applyTag("interface", i, i + 3);
                    }
                    if (s[i] == ',') {
                        applyTag("parenthesis", i, i + 1);
                    }
                }
            }

        }
    }
    return false;
}

class LweConsole : public LweGuiElement {
    GtkWidget* consoleText{};
    bool hasNewText{};
    int previousBufferSize{};
    GtkTextBuffer* buf{};
public:
    std::string textBuffer;
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        consoleText = gtk_text_view_new();
        gtk_text_view_set_accepts_tab(GTK_TEXT_VIEW(consoleText), false);
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        gtk_widget_set_vexpand(elementHandle, true);
        gtk_widget_set_vexpand(consoleText, true);
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_hexpand(consoleText, true);
        setPosition(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#9900CCFF", NULL);
    }
    void init(
        GtkWidget* grid, 
        int x, 
        int y, 
        int width, 
        int height, 
        int minWidth, 
        int minHeight) {
        consoleText = gtk_text_view_new();
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(elementHandle), minHeight);
        gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(elementHandle), minWidth);
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        setPosition(grid, x, y, width, height);
        buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_create_tag(buf, "function", "foreground", "#00FFFFFF", NULL);
        gtk_text_buffer_create_tag(buf, "comment", "foreground", "#006600FF", NULL);
        gtk_text_buffer_create_tag(buf, "variable", "foreground", "#FF00FFFF", NULL);
        gtk_text_buffer_create_tag(buf, "delegate", "foreground", "#FF8800FF", NULL);
        gtk_text_buffer_create_tag(buf, "interface", "foreground", "#FF0088FF", NULL);
        gtk_text_buffer_create_tag(buf, "error", "foreground", "#FF0000FF", NULL);
        gtk_text_buffer_create_tag(buf, "parenthesis", "foreground", "#9900CCFF", NULL);
    }

    template<typename... Args> void cPrint(std::string_view format, Args&&... args) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        std::string s = Svformat(format, Smake_format_args(args...));
        gtk_text_buffer_insert_markup(buf, &stop, s.c_str(), -1);
        GTKlock.unlock();
        scrollToEnd();
    }

    template<typename... Args> void tPrint(std::string_view format, Args&&... args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = Svformat(format, Smake_format_args(args...));
        textBuffer.append(s);
        hasNewText = true;
    }
    void scrollToEnd() {
        std::unique_lock GTKlock(GTKmutex);
        g_idle_add_full(
            G_PRIORITY_DEFAULT_IDLE, 
            scrollTextViewToEndHandler, 
            elementHandle, 
            NULL);
    }
    void updateFromBuffer() {
        if (hasNewText) {
            std::unique_lock GTKlock(GTKmutex);
            hasNewText = false;
            GtkTextIter end;
            gtk_text_buffer_get_end_iter(buf, &end);
            gtk_text_buffer_insert_markup(buf, &end, textBuffer.c_str(), -1);
            textBuffer.clear();
            GTKlock.unlock();
            scrollToEnd();
        }
    }

    void directOverwritePrint(const char* sIn) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_text_buffer_set_text(buf, sIn, -1);
        textBuffer.clear();
        GTKlock.unlock();
        scrollToEnd();
    }

    void directOverwritePrintSequence(const char* sIn) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_text_buffer_set_text(buf, sIn, -1);
        g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, formatSequenceBuffer, buf, NULL);
        textBuffer.clear();
        GTKlock.unlock();
        scrollToEnd();
    }

    void paintSequenceText() {
        if (previousBufferSize != gtk_text_buffer_get_char_count(buf)) {
            std::unique_lock GTKlock(GTKmutex);
            previousBufferSize = gtk_text_buffer_get_char_count(buf);
            g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, formatSequenceBuffer, buf, NULL);
        }

    }

    template<typename... Args> void overwritePrint(std::string_view format, Args&&... args) {
        std::unique_lock GTKlock(GTKmutex);
        std::string s = Svformat(format, Smake_format_args(args...));
        gtk_text_buffer_set_text(buf, s.c_str(), (int)s.length());
        textBuffer.clear();
        GTKlock.unlock();
        scrollToEnd();
    }

    void copyBuffer(char* destination, int64_t maxLength) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char* realBuf = gtk_text_buffer_get_text(buf, &start, &stop, false);
        std::string s(realBuf);
        s.copy(destination, maxLength);
    }

    void copyBuffer(std::string& s) {
        std::unique_lock GTKlock(GTKmutex);
        GtkTextIter start;
        GtkTextIter stop;
        if (gtk_text_buffer_get_char_count(buf) > 0) {
            gtk_text_buffer_get_start_iter(buf, &start);
            gtk_text_buffer_get_end_iter(buf, &stop);
            char* realBuf = gtk_text_buffer_get_text(buf, &start, &stop, false);
            std::string c(realBuf);
            s = c;
        }
        else {
            std::string c("None.");
            s = c;
        }
    }

    void clear() {
        std::unique_lock GTKlock(GTKmutex);
        char emptyBuffer[] = "";
        textBuffer.clear();
        GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_set_text(buf, emptyBuffer, 0);
        GTKlock.unlock();
        scrollToEnd();
    }
};

class LweButton : public LweGuiElement {
public:
    void init(
        const char* buttonName, 
        GtkWidget* grid, 
        int x, 
        int y, 
        int width, 
        int height, 
        auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction);
    }
    void init(
        const char* buttonName, 
        GtkWidget* grid, 
        int x, 
        int y, 
        int width, 
        int height, 
        auto buttonFunction, 
        gpointer functionData) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction, functionData);
    }
    void setFunction(auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), NULL);
    }
    void setFunction(auto buttonFunction, gpointer param) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), param);
    }
};

class LweCheckBox : public LweGuiElement {
public:
    void init(
        const char* buttonName, 
        GtkWidget* grid, 
        int x, 
        int y, 
        int width, 
        int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_check_button_new_with_label(buttonName);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    bool isChecked() {
        std::unique_lock GTKlock(GTKmutex);
        return (bool)gtk_check_button_get_active(GTK_CHECK_BUTTON(elementHandle));
    }
    void setFunction(auto buttonFunction) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect_after(elementHandle, "toggled", G_CALLBACK(buttonFunction), NULL);
    }
};

class LweProgressBar : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_progress_bar_new();
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_CENTER);
        gtk_widget_set_hexpand(elementHandle, true);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    void setValue(double fraction) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(elementHandle), fraction);
    }
};

class LwePulldown : public LweGuiElement {
    char** strArray;
    char* strings;
    int Nelements;
    int strLength;
public:
    LwePulldown() : Nelements(0), strLength(128) {
        strArray = new char* [100]();
        strings = new char[strLength * 100]();
    }
    ~LwePulldown() {
        delete[] strArray;
        delete[] strings;
    }
    void addElement(const char* newelement) {
        if (Nelements == 99) return;
        std::string s(newelement);
        stripLineBreaks(s);
        s.append("\0");
        strArray[Nelements] = &strings[Nelements * strLength];
        s.copy(strArray[Nelements], 127);
        ++Nelements;
    }
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_drop_down_new_from_strings(strArray);
        gtk_widget_set_hexpand(elementHandle, false);
        gtk_widget_set_vexpand(elementHandle, false);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    int getValue() {
        std::unique_lock GTKlock(GTKmutex);
        return (int)gtk_drop_down_get_selected(GTK_DROP_DOWN(elementHandle));
    }
    void setValue(int target) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_drop_down_set_selected(GTK_DROP_DOWN(elementHandle), target);
    }

    inline void removeCharacterFromString(std::string& s, char removedChar) {
        std::erase(s, removedChar);
    }

    void stripLineBreaks(std::string& s) {
        removeCharacterFromString(s, '\r');
        removeCharacterFromString(s, '\n');
    }
};

class LweWindow {
    GtkWidget* grid;
    GtkWidget* bigGrid;
    GtkWidget* consoleGrid;
    GtkWidget* plotGrid;
    GtkWidget* plotControlsGrid;
    GtkWidget* consoleControlsGrid;
    GtkWidget* consoleControlsSubgrid1;
    GtkWidget* consoleControlsSubgrid2;
    GtkWidget* plotControlsSubgrid1;
    GtkWidget* plotControlsSubgrid2;
public:
    GtkWidget* window;
    LweWindow() :grid(0), bigGrid(0), consoleGrid(0), plotGrid(0),
        plotControlsGrid(0),
        consoleControlsGrid(0),
        consoleControlsSubgrid1(0),
        consoleControlsSubgrid2(0),
        plotControlsSubgrid1(0),
        plotControlsSubgrid2(0),
        window(0) {}
    ~LweWindow() {};
    void init(
        GtkApplication* appHandle, 
        const char* windowName, 
        int width, 
        int height) {
        std::unique_lock GTKlock(GTKmutex);
        window = gtk_application_window_new(appHandle);
        gtk_window_set_title(GTK_WINDOW(window), windowName);
        gtk_window_set_default_size(GTK_WINDOW(window), width, height);
        bigGrid = gtk_grid_new();
        consoleGrid = gtk_grid_new();
        consoleControlsGrid = gtk_grid_new();
        consoleControlsSubgrid1 = gtk_grid_new();
        consoleControlsSubgrid2 = gtk_grid_new();
        plotGrid = gtk_grid_new();
        plotControlsGrid = gtk_grid_new();
        plotControlsSubgrid1 = gtk_grid_new();
        plotControlsSubgrid2 = gtk_grid_new();
        grid = gtk_grid_new();
        gtk_grid_set_row_spacing(GTK_GRID(grid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(grid), 1);

        gtk_grid_set_row_spacing(GTK_GRID(bigGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(bigGrid), 1);
        gtk_grid_set_row_spacing(GTK_GRID(plotGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(plotGrid), 1);
        gtk_grid_set_row_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(consoleControlsGrid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(grid), true);
        gtk_grid_set_column_homogeneous(GTK_GRID(plotGrid), true);
        gtk_grid_set_row_homogeneous(GTK_GRID(plotGrid), true);
        gtk_grid_set_row_spacing(GTK_GRID(consoleGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(consoleControlsGrid), 8);
        gtk_grid_set_column_spacing(GTK_GRID(consoleGrid), 1);
        gtk_window_set_child(GTK_WINDOW(window), bigGrid);
        gtk_widget_set_hexpand(grid, false);
        gtk_widget_set_halign(grid, GTK_ALIGN_START);
        gtk_widget_set_vexpand(grid, false);
        gtk_widget_set_valign(grid, GTK_ALIGN_START);
        gtk_widget_set_hexpand(consoleGrid, false);
        gtk_widget_set_hexpand(consoleControlsGrid, false);
        gtk_widget_set_vexpand(consoleGrid, true);
        gtk_widget_set_halign(consoleGrid, GTK_ALIGN_FILL);
        gtk_widget_set_halign(consoleControlsGrid, GTK_ALIGN_FILL);
        gtk_widget_set_valign(consoleGrid, GTK_ALIGN_FILL);
        gtk_widget_set_valign(consoleControlsSubgrid1, GTK_ALIGN_CENTER);
        gtk_widget_set_halign(consoleControlsSubgrid2, GTK_ALIGN_END);
        gtk_grid_attach(GTK_GRID(bigGrid), consoleGrid, 0, 1, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), consoleControlsGrid, 0, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(consoleControlsGrid), consoleControlsSubgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(consoleControlsGrid), consoleControlsSubgrid2, 1, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), plotGrid, 1, 0, 1, 2);
        gtk_grid_attach(GTK_GRID(bigGrid), grid, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(bigGrid), plotControlsGrid, 1, 2, 1, 1);
        gtk_grid_attach(GTK_GRID(plotControlsGrid), plotControlsSubgrid1, 0, 0, 1, 1);
        gtk_grid_attach(GTK_GRID(plotControlsGrid), plotControlsSubgrid2, 1, 0, 1, 1);
    }
    void present() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_window_present(GTK_WINDOW(window));
    }

    GtkWidget* parentHandle() {
        return grid;
    }
    GtkWidget* parentHandle(int index) {
        switch (index) {
        case 0:
            return grid;
        case 1:
            return consoleGrid;
        case 2:
            return plotGrid;
        case 3:
            return plotControlsSubgrid1;
        case 4:
            return plotControlsSubgrid2;
        case 5:
            return consoleControlsSubgrid1;
        case 6:
            return consoleControlsSubgrid2;
        }
        return grid;
    }

    GtkWindow* windowHandle() {
        return GTK_WINDOW(window);
    }
};

class LweDrawBox : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_drawing_area_new();
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_vexpand(elementHandle, true);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
        GTKlock.lock();
        gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(elementHandle), 12);
        gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(elementHandle), 12);
    }
    void setDrawingFunction(GtkDrawingAreaDrawFunc theFunction) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(elementHandle),
            theFunction,
            NULL, NULL);
    }
    void queueDraw() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_queue_draw(elementHandle);
    }
    void noVerticalExpantion() {
        std::unique_lock GTKlock(GTKmutex);
        gtk_widget_set_vexpand(elementHandle, false);
    }
};
class LweSpacer : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height, int spacing) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, spacing);
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_vexpand(elementHandle, true);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
};

class LweSlider : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        std::unique_lock GTKlock(GTKmutex);
        elementHandle = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, NULL);
        gtk_scale_set_draw_value(GTK_SCALE(elementHandle), true);
        gtk_scale_set_value_pos(GTK_SCALE(elementHandle), GTK_POS_LEFT);
        gtk_widget_set_hexpand(elementHandle, true);
        gtk_widget_set_margin_top(elementHandle, 0);
        gtk_widget_set_margin_bottom(elementHandle, 0);
        GTKlock.unlock();
        setPosition(grid, x, y, width, height);
    }
    int getIntValue() {
        std::unique_lock GTKlock(GTKmutex);
        return (int)gtk_range_get_value(GTK_RANGE(elementHandle));
    }
    double getDoubleValue() {
        std::unique_lock GTKlock(GTKmutex);
        return gtk_range_get_value(GTK_RANGE(elementHandle));
    }
    void setRange(double minVal, double maxVal) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_range_set_range(GTK_RANGE(elementHandle), minVal, maxVal);
    }
    void setDigits(int digits) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_scale_set_digits(GTK_SCALE(elementHandle), digits);
    }
    void setFunction(auto sliderFunction) {
        std::unique_lock GTKlock(GTKmutex);
        g_signal_connect_after(elementHandle, "change-value", G_CALLBACK(sliderFunction), NULL);
    }
    void setValue(int value) {
        std::unique_lock GTKlock(GTKmutex);
        gtk_range_set_value(GTK_RANGE(elementHandle), (double)value);
    }
};
