#ifndef UTILS_H
#define UTILS_H

#pragma once

#include <limits>
#include <filesystem>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "TLorentzVector.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "TString.h"

using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;
using ROOT::RDF::RSampleInfo;


/*
############################################
SNAPSHOT
############################################
*/
std::string setOutputDirectory(const std::string &outdir);
void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool dumpInput);

#endif
