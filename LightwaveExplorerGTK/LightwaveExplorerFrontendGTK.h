#pragma once
#include <sstream>
#include <complex>
#include "../LightwaveExplorerUtilities.h"
#undef __noinline__
#include <gtk/gtk.h>

#define _T(x) x

#if defined __linux__ || defined __APPLE__
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
    int rHex() { return (int)(15 * r); }
    int gHex() { return (int)(15 * g); }
    int bHex() { return (int)(15 * b); }
    int aHex() { return (int)(15 * a); }
    void setCairo(cairo_t* cr) { cairo_set_source_rgb(cr, r, g, b); }
    void setCairoA(cairo_t* cr) { cairo_set_source_rgba(cr, r, g, b, a); }
};

typedef struct plotStruct {
    GtkDrawingArea* area = NULL;
    cairo_t* cr = NULL;
    int width = 0;
    int height = 0;
    double* data = NULL;
    int ExtraLines = 0;
    double* data2 = NULL;
    double* data3 = NULL;
    double* data4 = NULL;
    double* dataX = NULL;
    const char* xLabel = NULL;
    const char* yLabel = NULL;
    bool hasDataX = FALSE;
    std::complex<double>* complexData = NULL;
    bool logScale = FALSE;
    double logMin = 0;
    int dataType = 0;
    double dx = 1.0;
    double x0 = 0.0;
    size_t Npts = 0;
    double unitY = 1.0;
    bool forceYmin = FALSE;
    double forcedYmin = 0.0;
    bool forceYmax = FALSE;
    double forcedYmax = 0.0;
    bool forceXmin = FALSE;
    double forcedXmin = 0.0;
    bool forceXmax = FALSE;
    double forcedXmax = 0.0;
    LweColor axisColor = LweColor(0.5, 0.5, 0.5, 0.5);
    LweColor textColor = LweColor(0.8, 0.8, 0.8, 0.8);
    LweColor color = LweColor(1, 1, 1, 1);
    LweColor color2 = LweColor(1, 1, 1, 1);
    LweColor color3 = LweColor(1, 1, 1, 1);
    LweColor color4 = LweColor(1, 1, 1, 1);
    std::string SVG;
    bool makeSVG = FALSE;

} plotStruct;

typedef struct imagePlotStruct {
    GtkDrawingArea* area = NULL;
    cairo_t* cr = NULL;
    int width = 0;
    int height = 0;
    double* data = NULL;
    std::complex<double>* complexData = NULL;
    int colorMap = 4;
    bool logScale = FALSE;
    double logMin = 0;
    int dataType = 0;
} imagePlotStruct;

class LweGuiElement {
public:
    GtkWidget* label;
    GtkWidget* elementHandle;
    int _x;
    int _y;
    int _width;
    int _height;
    bool isAttached;
    GtkWidget* _grid;
    LweGuiElement() : label(0), elementHandle(0), _x(0), _y(0), _width(0), _height(0), isAttached(FALSE), _grid(0) {
    }
    ~LweGuiElement() {}
    void setPosition(GtkWidget* grid, int x, int y, int width, int height) {
        if (_grid)gtk_grid_remove(GTK_GRID(_grid), elementHandle);
        _grid = grid;
        _x = x;
        _y = y;
        _width = width;
        _height = height;
        gtk_grid_attach(GTK_GRID(_grid), elementHandle, _x, _y, _width, _height);
    }
    void setLabel(int x, int y, const char* labelText) {
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), 45);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, _x + x, _y + y, 6, 1);
    }
    void setLabel(int x, int y, const char* labelText, int characters, int grids) {
        label = gtk_label_new(labelText);
        gtk_label_set_xalign(GTK_LABEL(label), 0.0);
        gtk_label_set_max_width_chars(GTK_LABEL(label), characters);
        gtk_label_set_ellipsize(GTK_LABEL(label), PANGO_ELLIPSIZE_END);
        gtk_grid_attach(GTK_GRID(_grid), label, _x + x, _y + y, grids, 1);
    }
    void squeeze() {
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_halign(elementHandle, GTK_ALIGN_END);
    }
    void setTooltip(const char* tooltipText) {
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
    LweTextBox() {}
    ~LweTextBox() {}
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        elementHandle = gtk_entry_new();
        gtk_widget_set_halign(elementHandle, GTK_ALIGN_START);
        gtk_widget_set_hexpand(elementHandle, FALSE);
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), 7);
        setPosition(grid, x, y, width, height);
    }

    double valueDouble() {
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
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(size_t* sdata) {
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToPointer(double* sdata) {
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        int len = gtk_entry_buffer_get_length(buf);
        if (len > 0) {
            char* cbuf = (char*)gtk_entry_buffer_get_text(buf);
            std::stringstream s(cbuf);
            s >> *sdata;
        }
    }

    void valueToTwoPointers(double* sdata, double* sdata2) {
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
        gtk_editable_set_max_width_chars(GTK_EDITABLE(elementHandle), charLimit);
    }

    void setToDouble(double in) {
        std::string s = Sformat(std::string_view("{:g}"), in);
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }
    template<typename... Args> void overwritePrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        GtkEntryBuffer* buf = gtk_entry_get_buffer(GTK_ENTRY(elementHandle));
        gtk_entry_buffer_set_text(buf, s.c_str(), (int)s.length());
    }

    void copyBuffer(char* destination, size_t maxLength) {
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

};

gboolean scrollTextViewToEndHandler(gpointer data) {
    GtkAdjustment* adjustment = gtk_scrolled_window_get_vadjustment(GTK_SCROLLED_WINDOW(data));
    gtk_adjustment_set_value(adjustment, gtk_adjustment_get_upper(adjustment));
    return FALSE;
}

class LweConsole : public LweGuiElement {
    GtkWidget* consoleText;
    bool hasNewText;
public:
    std::string textBuffer;
    LweConsole() :
        consoleText(0), hasNewText(0) {
        _grid = NULL;
    }
    ~LweConsole() {
    }
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        consoleText = gtk_text_view_new();
        gtk_text_view_set_accepts_tab(GTK_TEXT_VIEW(consoleText), FALSE);
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        gtk_widget_set_vexpand(elementHandle, TRUE);
        gtk_widget_set_vexpand(consoleText, TRUE);
        gtk_widget_set_hexpand(elementHandle, TRUE);
        gtk_widget_set_hexpand(consoleText, TRUE);
        setPosition(grid, x, y, width, height);
    }
    void init(GtkWidget* grid, int x, int y, int width, int height, int minWidth, int minHeight) {
        consoleText = gtk_text_view_new();
        elementHandle = gtk_scrolled_window_new();
        gtk_scrolled_window_set_min_content_height(GTK_SCROLLED_WINDOW(elementHandle), minHeight);
        gtk_scrolled_window_set_min_content_width(GTK_SCROLLED_WINDOW(elementHandle), minWidth);
        gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(elementHandle), consoleText);
        setPosition(grid, x, y, width, height);
    }

    template<typename... Args> void cPrint(std::string_view format, Args&&... args) {
        GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        textBuffer.assign(gtk_text_buffer_get_text(buf, &start, &stop, FALSE));
        
        std::string s = Svformat(format, Smake_format_args(args...));
        textBuffer.append(s);

        gtk_text_buffer_set_text(buf, textBuffer.c_str(), (int)textBuffer.length());
        
        scrollToEnd();
    }

    template<typename... Args> void tPrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        textBuffer.append(s);
        hasNewText = TRUE;
    }
    void scrollToEnd() {
        g_idle_add_full(G_PRIORITY_DEFAULT_IDLE, scrollTextViewToEndHandler, elementHandle, NULL);
    }
    void updateFromBuffer() {
        if (hasNewText) {
            hasNewText = FALSE;
            GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
            gtk_text_buffer_set_text(buf, textBuffer.c_str(), (int)textBuffer.length());
            scrollToEnd();
        }
        

    }

    void directOverwritePrint(const char* sIn) {
        std::string s(sIn);
        textBuffer.clear();
        textBuffer.assign(s);
        GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_set_text(buf, textBuffer.c_str(), (int)textBuffer.length());
        scrollToEnd();
    }

    template<typename... Args> void overwritePrint(std::string_view format, Args&&... args) {
        std::string s = Svformat(format, Smake_format_args(args...));
        textBuffer.assign(s);
        GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        gtk_text_buffer_set_text(buf, textBuffer.c_str(), (int)textBuffer.length());
        scrollToEnd();
    }

    void copyBuffer(char* destination, size_t maxLength) {
        GtkTextBuffer* buf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(consoleText));
        GtkTextIter start;
        GtkTextIter stop;
        gtk_text_buffer_get_start_iter(buf, &start);
        gtk_text_buffer_get_end_iter(buf, &stop);
        char* realBuf = gtk_text_buffer_get_text(buf, &start, &stop, FALSE);
        std::string s(realBuf);
        s.copy(destination, maxLength);
    }
};

class LweButton : public LweGuiElement {
public:
    void init(const char* buttonName, GtkWidget* grid, int x, int y, int width, int height, auto buttonFunction) {
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, FALSE);
        gtk_widget_set_vexpand(elementHandle, FALSE);
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_START);
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction);
    }
    void init(const char* buttonName, GtkWidget* grid, int x, int y, int width, int height, auto buttonFunction, gpointer functionData) {
        elementHandle = gtk_button_new_with_label(buttonName);
        gtk_widget_set_hexpand(elementHandle, FALSE);
        gtk_widget_set_vexpand(elementHandle, FALSE);
        setPosition(grid, x, y, width, height);
        setFunction(buttonFunction, functionData);
    }
    void setFunction(auto buttonFunction) {
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), NULL);
    }
    void setFunction(auto buttonFunction, gpointer param) {
        g_signal_connect(elementHandle, "clicked", G_CALLBACK(buttonFunction), param);
    }
};

class LweCheckBox : public LweGuiElement {
public:
    void init(const char* buttonName, GtkWidget* grid, int x, int y, int width, int height) {
        elementHandle = gtk_check_button_new_with_label(buttonName);
        setPosition(grid, x, y, width, height);
    }
    bool isChecked() {
        return (bool)gtk_check_button_get_active(GTK_CHECK_BUTTON(elementHandle));
    }
    void setFunction(auto buttonFunction) {
        g_signal_connect_after(elementHandle, "toggled", G_CALLBACK(buttonFunction), NULL);
    }
};

class LweProgressBar : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        elementHandle = gtk_progress_bar_new();
        gtk_widget_set_valign(elementHandle, GTK_ALIGN_CENTER);
        gtk_widget_set_hexpand(elementHandle, TRUE);
        setPosition(grid, x, y, width, height);
    }
    void setValue(double fraction) {
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
        s.append("\0");
        strArray[Nelements] = &strings[Nelements * strLength];
        s.copy(strArray[Nelements], 127);
        stripLineBreaks(strArray[Nelements]);
        ++Nelements;
    }
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        elementHandle = gtk_drop_down_new_from_strings(strArray);
        gtk_widget_set_hexpand(elementHandle, FALSE);
        gtk_widget_set_vexpand(elementHandle, FALSE);
        setPosition(grid, x, y, width, height);
    }
    int getValue() {
        return (int)gtk_drop_down_get_selected(GTK_DROP_DOWN(elementHandle));
    }
    void setValue(int target) {
        gtk_drop_down_set_selected(GTK_DROP_DOWN(elementHandle), target);
    }
};

class LweWindow {
    GtkWidget* window;
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
    LweWindow() :window(0), grid(0), bigGrid(0), consoleGrid(0), plotGrid(0),
        plotControlsGrid(0),
        consoleControlsGrid(0), 
        consoleControlsSubgrid1(0), 
        consoleControlsSubgrid2(0), 
        plotControlsSubgrid1(0), 
        plotControlsSubgrid2(0) {}
    ~LweWindow() {};
    void init(GtkApplication* appHandle, const char* windowName, int width, int height) {
        window = gtk_application_window_new(appHandle);
        gtk_window_set_title(GTK_WINDOW(window), windowName);
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
        gtk_grid_set_row_homogeneous(GTK_GRID(grid), TRUE);
        gtk_grid_set_column_homogeneous(GTK_GRID(consoleControlsGrid), TRUE);
        gtk_grid_set_column_homogeneous(GTK_GRID(grid), TRUE);
        gtk_grid_set_column_homogeneous(GTK_GRID(plotGrid), TRUE);
        gtk_grid_set_row_homogeneous(GTK_GRID(plotGrid), TRUE);
        gtk_grid_set_row_spacing(GTK_GRID(consoleGrid), 1);
        gtk_grid_set_column_spacing(GTK_GRID(consoleControlsGrid), 8);
        gtk_grid_set_column_spacing(GTK_GRID(consoleGrid), 1);
        gtk_window_set_child(GTK_WINDOW(window), bigGrid);
        gtk_widget_set_hexpand(grid, FALSE);
        gtk_widget_set_halign(grid, GTK_ALIGN_START);
        gtk_widget_set_vexpand(grid, FALSE);
        gtk_widget_set_valign(grid, GTK_ALIGN_START);
        gtk_widget_set_hexpand(consoleGrid, FALSE);
        gtk_widget_set_hexpand(consoleControlsGrid, FALSE);
        gtk_widget_set_vexpand(consoleGrid, TRUE);
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
        elementHandle = gtk_drawing_area_new();
        gtk_widget_set_hexpand(elementHandle, TRUE);
        gtk_widget_set_vexpand(elementHandle, TRUE);
        setPosition(grid, x, y, width, height);
        gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(elementHandle), 12);
        gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(elementHandle), 12);
    }
    void setDrawingFunction(GtkDrawingAreaDrawFunc theFunction) {
        gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(elementHandle),
            theFunction,
            NULL, NULL);
    }
    void queueDraw() {
        gtk_widget_queue_draw(elementHandle);
    }
    void noVerticalExpantion() {
        gtk_widget_set_vexpand(elementHandle, FALSE);
    }
};
class LweSpacer : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height, int spacing) {
        elementHandle = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, spacing);
        gtk_widget_set_hexpand(elementHandle, TRUE);
        gtk_widget_set_vexpand(elementHandle, TRUE);
        setPosition(grid, x, y, width, height);
    }
};

class LweSlider : public LweGuiElement {
public:
    void init(GtkWidget* grid, int x, int y, int width, int height) {
        elementHandle = gtk_scale_new(GTK_ORIENTATION_HORIZONTAL, NULL);
        gtk_scale_set_draw_value(GTK_SCALE(elementHandle), TRUE);
        gtk_scale_set_value_pos(GTK_SCALE(elementHandle), GTK_POS_LEFT);
        gtk_widget_set_hexpand(elementHandle, TRUE);
        gtk_widget_set_margin_top(elementHandle, 0);
        gtk_widget_set_margin_bottom(elementHandle, 0);
        setPosition(grid, x, y, width, height);
    }
    int getIntValue() {
        return (int)gtk_range_get_value(GTK_RANGE(elementHandle));
    }
    double getDoubleValue() {
        return gtk_range_get_value(GTK_RANGE(elementHandle));
    }
    void setRange(double minVal, double maxVal) {
        gtk_range_set_range(GTK_RANGE(elementHandle), minVal, maxVal);
    }
    void setDigits(int digits) {
        gtk_scale_set_digits(GTK_SCALE(elementHandle), digits);
    }
    void setFunction(auto sliderFunction) {
        g_signal_connect_after(elementHandle, "change-value", G_CALLBACK(sliderFunction), NULL);
    }
    void setValue(int value) {
        gtk_range_set_value(GTK_RANGE(elementHandle), (double)value);
    }
};

int LwePlot2d(plotStruct* inputStruct) {
    plotStruct* s = (plotStruct*)inputStruct;
    if ((*s).Npts == 0) return 1;
    size_t iMin = 0;
    size_t iMax = (*s).Npts;
    std::string svgString;
    cairo_t* cr = (*s).cr;
    cairo_font_extents_t fe;
    memset(&fe, 0, sizeof(cairo_font_extents_t));
    cairo_text_extents_t te;
    memset(&te, 0, sizeof(cairo_text_extents_t));
    double fontSize = 14.0;
    cairo_set_font_size(cr, fontSize);
    cairo_select_font_face(cr, "Arial", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_font_extents(cr, &fe);
    double x1, y1, x2, y2;
    double width;
    double height;

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
    double* xValues = new double[(*s).Npts + 2]();
    for (int i = 0; i < (*s).Npts; ++i) {
        if ((*s).hasDataX) currentX = (double)(*s).dataX[i];
        else { currentX = (double)(i * (*s).dx + (*s).x0); }
        if (i == 0) {
            minX = currentX;
            maxX = currentX;
        }
        xValues[i] = currentX;
        if ((*s).forceXmin && (currentX < (*s).forcedXmin)) {
            iMin = i + 1;
        }
        if ((*s).forceXmax && (currentX > (*s).forcedXmax)) {
            iMax = i;
            break;
        }
        maxX = maxN(currentX, maxX);
        minX = minN(currentX, minX);
    }

    for (size_t i = iMin; i < iMax; ++i) {
        if ((*s).logScale) { currentY = (double)log10((*s).data[i]); }
        else { currentY = (double)(*s).data[i]; }
        maxY = maxN(currentY, maxY);
        minY = minN(currentY, minY);
        if ((*s).ExtraLines > 0) {
            if ((*s).logScale) { currentY = (double)log10((*s).data2[i]); }
            else { currentY = (double)(*s).data2[i]; }
            maxY = maxN(currentY, maxY);
            minY = minN(currentY, minY);
        }
    }

    if (minY == maxY) {
        minY = -1;
        maxY = 1;
    }
    if ((*s).forceYmin) {
        minY = (double)(*s).forcedYmin * (*s).unitY;
        if ((*s).logScale) minY = log10(minY);
    }
    if ((*s).forceYmax) {
        maxY = (double)(*s).forcedYmax * (*s).unitY;
        if ((*s).logScale) maxY = log10(maxY);
    }
    if ((*s).forceXmin) {
        minX = (double)(*s).forcedXmin;
    }
    if ((*s).forceXmax) {
        maxX = (double)(*s).forcedXmax;
    }

    //Tickmark labels
    int NyTicks = 3;
    std::string messageBuffer;
    double yTicks1[3] = { maxY, 0.5 * (maxY + minY), minY };
    double xTicks1[3] = { minX + 0.25 * (maxX - minX), minX + 0.5 * (maxX - minX), minX + 0.75 * (maxX - minX) };
    height = (double)(*s).height;
    width = (double)(*s).width;

    //Start SVG file if building one
    auto SVGh = [&](double x) {
        return (int)(15 * x);
    };
    if ((*s).makeSVG) {
        svgString.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        svgString.append(Sformat("<svg width=\"{}\" height=\"{}\" viewBox=\"0 0 {} {}\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
            width, height, width, height));
        svgString.append(Sformat("<rect fill=\"#{:x}{:x}{:x}\" stroke=\"#000\" x=\"0\" y=\"0\" width=\"{}\" height=\"{}\"/>\n",
            SVGh(0.0f), SVGh(0.0f), SVGh(0.0f), width, height));
    }
    LweColor black(0, 0, 0, 0);
    cairo_rectangle(cr, 0, 0, width, height);
    black.setCairo(cr);
    cairo_fill(cr);
    width -= axisSpaceX;
    height -= axisSpaceY;
    double scaleX = width / ((double)(maxX - minX));
    double scaleY = height / ((double)(maxY - minY));
    LweColor currentColor = (*s).textColor;


    //make the paintbrush
    //hr = pRenderTarget->CreateSolidColorBrush((*s).textColor, &pBrush);
    currentColor = (*s).textColor;
    //lambdas for writing components of SVG file
    auto SVGstdline = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#{:x}{:x}{:x}\" stroke-width=\"{}\"/>\n", x1, y1, x2, y2, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth));
    };

    auto SVGstartPolyLine = [&]() {
        svgString.append(Sformat("<polyline points=\""));
    };

    auto SVGendPolyLine = [&]() {
        svgString.append(Sformat("\" stroke=\"#{:x}{:x}{:x}\" stroke-width=\"{}\"/>\n", currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), lineWidth));
    };

    auto SVGaddXYtoPolyLine = [&](double& a, double& b) {
        svgString.append(Sformat("{},{} ", a, b));
    };

    auto SVGstdcircle = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<circle cx=\"{}\" cy=\"{}\" r=\"{}\" stroke=\"none\" fill=\"#{:x}{:x}{:x}\" />\n", x1, y1, radius, currentColor.rHex(), currentColor.gHex(), currentColor.bHex()));
    };

    auto SVGstartgroup = [&]() {
        if ((*s).makeSVG)svgString.append("<g>\n");
    };

    auto SVGendgroup = [&]() {
        if ((*s).makeSVG)svgString.append("</g>\n");
    };

    auto SVGcentertext = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" text-anchor=\"middle\">\n{}\n</text>\n", fontSize - 1, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), 0.5 * (layoutBottom + layoutTop - te.height), messageBuffer));
    };

    auto SVGlefttext = [&]() {
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\">\n{}\n</text>\n", fontSize - 1, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), layoutLeft, layoutTop + fontSize, std::string(messageBuffer)));
    };


    cairo_set_line_width(cr, lineWidth);
    auto cairoLine = [&]() {
        currentColor.setCairo(cr);
        cairo_move_to(cr, x1, y1);
        cairo_line_to(cr, x2, y2);
        cairo_stroke(cr);
    };
    auto cairoCircle = [&]() {
        currentColor.setCairo(cr);
        cairo_arc(cr, x1, y1, radius, 0, 6.2831853071795862);
        cairo_fill(cr);
    };
    auto cairoRightText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, layoutRight - te.width - 3, 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer.c_str());
    };
    auto cairoCenterText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), 0.5 * (layoutBottom + layoutTop - te.height));
        cairo_show_text(cr, messageBuffer.c_str());
    };

    auto cairoVerticalText = [&]() {
        currentColor.setCairo(cr);
        cairo_text_extents(cr, messageBuffer.c_str(), &te);
        cairo_move_to(cr, 0.0, height);
        cairo_rotate(cr, -3.1415926535897931 / 2);
        cairo_rel_move_to(cr, 0.5 * (layoutLeft + layoutRight - te.width), fontSize);
        cairo_show_text(cr, messageBuffer.c_str());
        cairo_rotate(cr, 3.1415926535897931 / 2);
    };

    currentColor = (*s).textColor;
    //y-tick text labels
    for (int i = 0; i < NyTicks; ++i) {
        double ytVal = yTicks1[i] / (*s).unitY;
        if ((*s).logScale) {
            switch (i) {
            case 0:
                ytVal = pow(10.0, maxY) / (*s).unitY;
                break;
            case 1:
                ytVal = pow(10.0, minY + 0.5 * (maxY - minY)) / (*s).unitY;
                break;
            case 2:
                ytVal = pow(10.0, minY) / (*s).unitY;
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
    if ((*s).yLabel != NULL) {
        messageBuffer = std::string((*s).yLabel);
        layoutLeft = 0;
        layoutTop = height;
        layoutBottom = height + axisSpaceY;
        layoutRight = height;

        cairoVerticalText();
        if ((*s).makeSVG)svgString.append(Sformat("<text font-family=\"Arial\" font-size=\"{}\" fill=\"#{:x}{:x}{:x}\" x=\"{}\" y=\"{}\" text-anchor=\"middle\" transform=\"translate({}, {}) rotate(-90)\">\n{}\n</text>\n", fontSize, currentColor.rHex(), currentColor.gHex(), currentColor.bHex(), 0.5 * (layoutLeft + layoutRight), layoutTop + fontSize, -(layoutLeft + layoutRight), height, messageBuffer));
    }

    //x-axis name
    if ((*s).xLabel != NULL) {
        layoutLeft = axisSpaceX;
        layoutTop = height + 2.8 * fontSize;
        layoutBottom = height + axisSpaceY;
        layoutRight = axisSpaceX + width;
        messageBuffer.assign((*s).xLabel);
        cairoCenterText();
        SVGcentertext();
    }

    //x-axis tick labels
    for (int i = 0; i < 3; ++i) {
        messageBuffer.assign(Sformat("{}", (int)round(xTicks1[i])));
        layoutLeft = (double)(axisSpaceX + 0.25 * width * ((size_t)(i)+1) - axisSpaceX / 2);
        layoutTop = height + 3;
        layoutBottom = height + axisSpaceY;
        layoutRight = layoutLeft + axisSpaceX;

        cairoCenterText();
        SVGcentertext();
    }

    //Draw axes and tickmarks
    SVGstartgroup();
    currentColor = (*s).axisColor;
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
    //currently dots are always there but has been refactored to make it easier to turn the off if I want.
    auto plotCairoDots = [&](double* y) {
        currentColor.setCairo(cr);
        for (size_t i = iMin; i < iMax - 1; ++i) {
            x1 = scaleX * (xValues[i] - minX) + axisSpaceX;
            if ((*s).logScale) {
                y1 = height - scaleY * ((double)log10(y[i]) - (double)minY);
            }
            else {
                y1 = height - scaleY * ((double)y[i] - (double)minY);
            }
            if (y1 <= height) {
                cairo_arc(cr, x1, y1, radius, 0, 6.2831853071795862);
                cairo_fill(cr);
            }
        }
    };

    auto plotCairoPolyline = [&](double* y) {
        currentColor.setCairo(cr);
        x2 = scaleX * (xValues[iMin] - minX) + axisSpaceX;
        if ((*s).logScale) {
            y2 = height - scaleY * (log10(y[iMin]) - minY);
        }
        else {
            y2 = height - scaleY * (y[iMin] - minY);
        }

        for (size_t i = iMin+1; i < iMax; ++i) {
            x1 = x2;
            x2 = scaleX * (xValues[i] - minX);
            y1 = y2;
            if ((*s).logScale) {
                y2 = height - scaleY * (log10(y[i]) - minY);
            }
            else {
                y2 = height - scaleY * (y[i] - minY);
            }
            x2 += axisSpaceX;

            if (y1 <= height) {
                if (y2 <= height) {
                    cairo_move_to(cr, x1, y1);
                    cairo_line_to(cr, x2, y2);
                }
                else {
                    x2 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                    y2 = height;
                    cairo_move_to(cr, x1, y1);
                    cairo_line_to(cr, x2, y2);
                }
            }
            else if (y2 <= height) {
                x1 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                y1 = height;
                cairo_move_to(cr, x1, y1);
                cairo_line_to(cr, x2, y2);
            }
        }
        cairo_stroke(cr);
    };

    auto plotSVGPolyline = [&](double* y) {
        bool lineOn = FALSE;
        x2 = scaleX * (xValues[iMin] - minX) + axisSpaceX;
        if ((*s).logScale) {
            y2 = height - scaleY * (log10(y[iMin]) - minY);
        }
        else {
            y2 = height - scaleY * (y[iMin] - minY);
        }
        if (y2 <= height) {
            SVGstartPolyLine();
            SVGaddXYtoPolyLine(x2, y2);
            lineOn = TRUE;
        }
        for (size_t i = iMin + 1; i < iMax; ++i) {
            x1 = x2;
            x2 = scaleX * (xValues[i] - minX);
            y1 = y2;
            if ((*s).logScale) {
                y2 = height - scaleY * (log10(y[i]) - minY);
            }
            else {
                y2 = height - scaleY * (y[i] - minY);
            }
            x2 += axisSpaceX;

            if (y1 <= height) {
                if (y2 <= height) {
                    if (!lineOn) {
                        SVGstartPolyLine();
                    }
                    SVGaddXYtoPolyLine(x2, y2);
                }
                else {
                    x2 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                    y2 = height;
                    SVGaddXYtoPolyLine(x2, y2);
                    SVGendPolyLine();
                }
            }
            else if (y2 <= height) {
                x1 = x1 + (height - y1) / ((y2 - y1) / (x2 - x1));
                y1 = height;
                if (!lineOn) {
                    SVGstartPolyLine();
                }
                SVGaddXYtoPolyLine(x2, y2);
            }
        }
        if (lineOn) {
            SVGendPolyLine();
        }
    };

    auto plotSVGDots = [&](double* y) {
        SVGstartgroup();
        for (size_t i = iMin; i < iMax - 1; ++i) {
            x1 = scaleX * (xValues[i] - minX) + axisSpaceX;
            if ((*s).logScale) {
                y1 = height - scaleY * ((double)log10(y[i]) - (double)minY);
            }
            else {
                y1 = height - scaleY * ((double)y[i] - (double)minY);
            }
            if (y1 <= height) {
                SVGstdcircle();
            }
        }
        SVGendgroup();
    };

    //Plot the main line
    currentColor = (*s).color;
    plotCairoPolyline((*s).data);
    plotCairoDots((*s).data);
    if ((*s).makeSVG) {
        plotSVGPolyline((*s).data);
        plotSVGDots((*s).data);
    }
    //Optional overlay curves
    if ((*s).ExtraLines > 0) {
        currentColor = (*s).color2;
        plotCairoPolyline((*s).data2);
        plotCairoDots((*s).data2);
        if ((*s).makeSVG) {
            plotSVGPolyline((*s).data2);
            plotSVGDots((*s).data2);
        }
    }
    if ((*s).ExtraLines > 1) {
        currentColor = (*s).color3;
        plotCairoPolyline((*s).data3);
        plotCairoDots((*s).data3);
        if ((*s).makeSVG) {
            plotSVGPolyline((*s).data3);
            plotSVGDots((*s).data3);
        }
    }
    if ((*s).ExtraLines > 2) {
        currentColor = (*s).color4;
        plotCairoPolyline((*s).data4);
        plotCairoDots((*s).data4);
        if ((*s).makeSVG) {
            plotSVGPolyline((*s).data4);
            plotSVGDots((*s).data4);
        }
    }

    delete[] xValues;
    if ((*s).makeSVG) {
        svgString.append("</svg>");
        (*s).SVG.assign(svgString);
    }
    return 0;
}

void openFileDialogCallback(GtkWidget* widget, gpointer pathTarget);
void saveFileDialogCallback(GtkWidget* widget, gpointer pathTarget);
void drawFourierImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawFourierImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawTimeImage1(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawTimeImage2(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawField1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawField2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrum1Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawSpectrum2Plot(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
void drawProgress(GtkDrawingArea* area, cairo_t* cr, int width, int height, gpointer data);
int drawArrayAsBitmap(cairo_t* cr, int Nx, int Ny, float* data, int cm);
int linearRemapZToLogFloatShift(std::complex<double>* A, int nax, int nay, float* B, int nbx, int nby, double logMin);
int LwePlot2d(plotStruct* inputStruct);
void checkLibraryAvailability();
void setInterfaceValuesToActiveValues();
int formatSequence(char* cString, size_t N);
int insertLineBreaksAfterSemicolons(char* cString, size_t N);
int freeSemipermanentGrids();
void mainSimThread(int pulldownSelection, int pulldownSelection2);
void launchRunThread();
void independentPlotQueue();
void loadCallback(GtkWidget* widget, gpointer pathTarget);
void launchFitThread();
void fittingThread(int pulldownSelection);
void stopButtonCallback();
void svgCallback();
void createRunFile();
static void buttonAddSameCrystal();
static void buttonAddDefault();
static void buttonAddRotation();
static void buttonAddPulse();
static void buttonAddMirror();
static void buttonAddFilter();
static void buttonAddLinear();
static void buttonAddAperture();
static void buttonAddFarFieldAperture();
static void buttonAddForLoop();

