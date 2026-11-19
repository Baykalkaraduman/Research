#pragma once

namespace Units {

    // Base Units
    constexpr double cm = 1;               // centimeter
    constexpr double GeV = 1;              // Giga Electron Volt
    constexpr double ns = 1;               // nanosecond

    // Other Length Units
    constexpr double fm = cm*1e-13;        // femtometer
    constexpr double pm = cm*1e-10;        // picometer
    constexpr double nm = cm*1e-7;         // nanometer
    constexpr double um = cm*1e-4;         // micrometer
    constexpr double mm = cm*0.1;          // millimeter
    constexpr double m = cm*100;           // meter
    constexpr double km = cm*1e+5;         // kilometer

    // Other Energy Units
    constexpr double MeV = GeV*1e-3;       // Mega Electron Volt
    constexpr double keV = GeV*1e-6;       // kilo Electron Volt
    constexpr double eV = GeV*1e-9;        // Electron Volt
    constexpr double TeV = GeV*1e+3;       // Tera Electron Volt

    // Other Time Units
    constexpr double fs = ns*1e-6;         // femtoseconds
    constexpr double ps = ns*1e-3;         // picoseconds
    constexpr double ms = ns*1e+6;         // milliseconds
    constexpr double us = ns*1e+3;         // microseconds
    constexpr double s = ns*1e+9;          // seconds

    // Area Units for Cross Section
    constexpr double barn = 1e-24*cm*cm;   // barn
    constexpr double mbarn = barn*1e-3;    // millibarn

}