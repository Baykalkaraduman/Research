#include "EndfSection.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string_view>

// ==============================================================================
// Internal Helpers
// ==============================================================================
namespace {

// Parses a fixed-width integer field without allocating a substring.
int parseIntField(std::string_view text, std::size_t pos, std::size_t len) noexcept {
    int value = 0;
    int sign = 1;
    for (std::size_t i = 0; i < len; ++i) {
        const char c = text[pos + i];
        if (c == '-') {
            sign = -1;
        } else if (c >= '0' && c <= '9') {
            value = value * 10 + (c - '0');
        }
    }
    return sign * value;
}

bool isSection(std::string_view text, int mf, int mt) noexcept {
    if (text.size() < 75) {
        return false;
    }
    return parseIntField(text, 70, 2) == mf && parseIntField(text, 72, 3) == mt;
}

// Parses an ENDF numeric field ("1.234567+3", "-2.5-4", "3.8D-4") using a stack buffer.
double parseNumberView(std::string_view text) {
    char buf[64];
    std::size_t len = 0;

    for (char c : text) {
        if (c == ' ') {
            continue;
        }
        if (c == 'D' || c == 'd') {
            c = 'e';
        }
        if ((c == '+' || c == '-') && len > 0 && buf[len - 1] != 'e' && buf[len - 1] != 'E') {
            if (len + 2 >= sizeof(buf)) {
                throw std::runtime_error("EndfSection: number too long '" + std::string(text) + "'");
            }
            buf[len++] = 'e';
        }
        if (len + 1 >= sizeof(buf)) {
            throw std::runtime_error("EndfSection: number too long '" + std::string(text) + "'");
        }
        buf[len++] = c;
    }

    if (len == 0) {
        return 0.0;
    }
    buf[len] = '\0';

    char* end = nullptr;
    const double value = std::strtod(buf, &end);
    if (end == buf) {
        throw std::runtime_error("EndfSection: cannot read number '" + std::string(text) + "'");
    }
    return value;
}

} // namespace

// ==============================================================================
// Construction & Section Lookup
// ==============================================================================
EndfSection::EndfSection(const std::string& path, int mf, int mt) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("EndfSection: cannot open " + path);
    }

    std::string text;
    while (std::getline(in, text)) {
        if (isSection(text, mf, mt)) {
            lines.emplace_back(text, 0, 66);
        }
    }

    if (lines.empty()) {
        throw std::runtime_error("EndfSection: no MF=" + std::to_string(mf) + " MT=" + std::to_string(mt) + " in " + path);
    }
}

std::vector<int> EndfSection::findSections(const std::string& path, int mf) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("EndfSection: cannot open " + path);
    }

    std::vector<int> mts;
    std::string text;
    while (std::getline(in, text)) {
        if (text.size() < 75 || parseIntField(text, 70, 2) != mf) {
            continue;
        }
        const int mt = parseIntField(text, 72, 3);
        if (mt > 0 && (mts.empty() || mts.back() != mt)) {
            mts.push_back(mt);
        }
    }
    return mts;
}

double EndfSection::parseNumber(const std::string& text) {
    return parseNumberView(text);
}

double EndfSection::nextValue() {
    if (field >= 6) {
        if (line >= lines.size()) {
            throw std::runtime_error("EndfSection: read past the end of the section");
        }
        ++line;
        field = 0;
    }

    const std::string& text = lines[line - 1];
    const std::size_t start = 11 * field++;
    if (start >= text.size()) {
        return 0.0;
    }
    return parseNumberView(std::string_view(text).substr(start, 11));
}

// ==============================================================================
// Sequential Record Readers
// ==============================================================================
EndfSection::Cont EndfSection::readCont() {
    startNewLine();
    Cont c;
    c.c1 = nextValue();
    c.c2 = nextValue();
    c.l1 = std::lround(nextValue());
    c.l2 = std::lround(nextValue());
    c.n1 = std::lround(nextValue());
    c.n2 = std::lround(nextValue());
    return c;
}

EndfSection::List EndfSection::readList() {
    List list;
    list.head = readCont();
    startNewLine();

    const auto count = static_cast<std::size_t>(std::max(0L, list.head.n1));
    list.values.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        list.values.push_back(nextValue());
    }
    return list;
}

EndfSection::Tab1 EndfSection::readTab1() {
    Tab1 tab;
    tab.head = readCont();

    startNewLine();
    const auto nRegions = static_cast<std::size_t>(std::max(0L, tab.head.n1));
    tab.nbt.reserve(nRegions);
    tab.interpolation.reserve(nRegions);
    for (std::size_t i = 0; i < nRegions; ++i) {
        tab.nbt.push_back(std::lround(nextValue()));
        tab.interpolation.push_back(std::lround(nextValue()));
    }

    startNewLine();
    const auto nPoints = static_cast<std::size_t>(std::max(0L, tab.head.n2));
    tab.x.reserve(nPoints);
    tab.y.reserve(nPoints);
    for (std::size_t i = 0; i < nPoints; ++i) {
        tab.x.push_back(nextValue());
        tab.y.push_back(nextValue());
    }
    return tab;
}