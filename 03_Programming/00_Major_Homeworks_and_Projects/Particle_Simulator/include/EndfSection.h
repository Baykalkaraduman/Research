#pragma once

#include <cstddef>
#include <string>
#include <vector>

// ==============================================================================
// Class: EndfSection
// Sequential record reader (CONT, LIST, TAB1) for one (MF, MT) section of an
// ENDF-6 formatted data file.
// ==============================================================================
class EndfSection {
private:
    // --------------------------------------------------------------------------
    // Member Data & Internal Helpers
    // --------------------------------------------------------------------------
    std::vector<std::string> lines;
    std::size_t line = 0;
    std::size_t field = 6;

    double nextValue();
    void startNewLine() {field = 6;}

public:
    // --------------------------------------------------------------------------
    // Record Types
    // --------------------------------------------------------------------------
    struct Cont {
        double c1 = 0.0, c2 = 0.0;
        long l1 = 0, l2 = 0, n1 = 0, n2 = 0;
    };

    struct List {
        Cont head;
        std::vector<double> values;
    };

    struct Tab1 {
        Cont head;
        std::vector<long> nbt, interpolation;
        std::vector<double> x, y;
    };

    // --------------------------------------------------------------------------
    // Construction & Utilities
    // --------------------------------------------------------------------------
    EndfSection(const std::string& path, int mf, int mt);

    [[nodiscard]] static std::vector<int> findSections(const std::string& path, int mf);
    [[nodiscard]] static double parseNumber(const std::string& text);

    [[nodiscard]] bool atEnd() const {return line >= lines.size();}

    // --------------------------------------------------------------------------
    // Sequential Record Readers
    // --------------------------------------------------------------------------
    [[nodiscard]] Cont readCont();
    [[nodiscard]] List readList();
    [[nodiscard]] Tab1 readTab1();
};