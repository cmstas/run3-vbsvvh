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
DEFINE METADATA
############################################
*/

RNode defineMetadata(RNode df);

/*
############################################
LUMIMASK
############################################
*/

class lumiMask {
public:
    using Run = unsigned int;
    using LumiBlock = unsigned int;

    class LumiBlockRange {
    public:
        LumiBlockRange(Run run, LumiBlock firstLumi, LumiBlock lastLumi) : m_run(run), m_firstLumi(firstLumi), m_lastLumi(lastLumi ? lastLumi : std::numeric_limits<LumiBlock>::max()) {}
        Run run() const { return m_run; }
        LumiBlock firstLumi() const { return m_firstLumi; }
        LumiBlock lastLumi () const { return m_lastLumi ; }
    private:
        Run m_run;
        LumiBlock m_firstLumi;
        LumiBlock m_lastLumi;
    };

    explicit lumiMask(const std::vector<LumiBlockRange>& accept) : m_accept(accept) {
        std::sort(m_accept.begin(), m_accept.end());
    }

    bool accept(Run run, LumiBlock lumi) const { 
        return std::binary_search(m_accept.begin(), m_accept.end(), LumiBlockRange(run, lumi, lumi)); 
    }

    static lumiMask fromJSON(const std::string& fileName, lumiMask::Run firstRun=0, lumiMask::Run lastRun=0);

private:
    std::vector<LumiBlockRange> m_accept;
};

bool operator< ( const lumiMask::LumiBlockRange& lh, const lumiMask::LumiBlockRange& rh );

/*
############################################
REMOVE DUPLICATES
############################################
*/

class FilterOnePerKind {
    std::unordered_set<size_t> _seenCategories;
public:
    bool operator()(unsigned int run, unsigned int luminosityBlock, unsigned long long event) {
        std::hash<std::string> categoryHasher;
        std::string eventStr = std::to_string(run) + "," + std::to_string(luminosityBlock) + "," + std::to_string(event);
        size_t category = categoryHasher(eventStr);
        {
        // build char category from run, luminosityBlock, event
        R__READ_LOCKGUARD(ROOT::gCoreMutex); // many threads can take a read lock concurrently
        if (_seenCategories.count(category) == 1)
            return false;
        }
        // if we are here, `category` was not already in _seenCategories
        R__WRITE_LOCKGUARD(ROOT::gCoreMutex); // only one thread at a time can take the write lock
        _seenCategories.insert(category);
        return true;
    }
};

RNode removeDuplicates(RNode df);

/*
############################################
CUTFLOW
############################################
*/

class Cutflow {
public:
    Cutflow(RNode df, const std::vector<std::string>& cuts);
    void Print(std::string output_file);
private:
    RNode _df;
    std::vector<std::string> _cuts;
    std::vector<std::pair<ROOT::RDF::RResultPtr<double>, ROOT::RDF::RResultPtr<double>>> _cutflow;
};

/*
############################################
SELECTION UTILS
############################################
*/

RVec<float> VfDeltaR (RVec<float> vec_eta, RVec<float> vec_phi, float obj_eta, float obj_phi);
RVec<float> VfInvariantMass(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass);
RVec<float> VfInvariantPt(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass);
RVec<float> VfInvariantPhi(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass);
RVec<float> MT(RVec<float> vec_pt, RVec<float> vec_phi, float obj_pt, float obj_phi);
float fInvariantMass(float obj1_pt, float obj1_eta, float obj1_phi, float obj1_mass, float obj2_pt, float obj2_eta, float obj2_phi, float obj2_mass);
RVec<int> VBS_MaxE( RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass );
RVec<int> VBS_MaxEtaJJ(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass);
RVec<float> VfdRfromClosestJet(const ROOT::RVecF &ak4_eta, const ROOT::RVecF &ak4_phi, const ROOT::RVecF &ak8_eta, const ROOT::RVecF &ak8_phi);

/*
############################################
SNAPSHOT
############################################
*/
std::string setOutputDirectory(const std::string &dir);
void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool isData = false);

#endif