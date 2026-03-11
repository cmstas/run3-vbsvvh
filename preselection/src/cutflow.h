#ifndef CUTFLOW_H
#define CUTFLOW_H

#pragma once

#include <string>
#include <vector>
#include <iostream>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

using RNode = ROOT::RDF::RNode;

class Cutflow
{
public:
    static void Enable(bool enabled = true);
    static bool IsEnabled();
    static void SetWeightCol(const std::string &col);
    static void Add(RNode df, const std::string &label);
    static void Print(std::ostream &out = std::cout);
    static void Clear();
    static std::size_t Size();

private:
    struct Entry {
        std::string label;
        mutable ROOT::RDF::RResultPtr<double> wResult;
        mutable ROOT::RDF::RResultPtr<ULong64_t> cResult;
        bool isCount{false};

        double value() const {
            return isCount ? static_cast<double>(*cResult) : *wResult;
        }
    };

    static std::string &weightCol();
    static std::vector<Entry> &entries();
    static bool &enabled();
};

#endif // CUTFLOW_H