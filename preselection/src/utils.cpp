#include "utils.h"

#include <algorithm>
#include "TFile.h"
#include "TTree.h"

/*
############################################
RDF UTILS
############################################
*/

RNode defineMetadata(RNode df, bool isData = false) {
    if (isData) {df = df.Define("genWeight", []() { return 1.f; }, {});}
    return df.DefinePerSample("xsec", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("xsec");})
        .DefinePerSample("lumi", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("lumi");})
        .DefinePerSample("sumw", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("sumw");})
        .DefinePerSample("kind", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("kind");})
        .DefinePerSample("year", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("year");})
        .DefinePerSample("shortname", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("shortname");})
        .DefinePerSample("name", [](unsigned int slot, const RSampleInfo &id) { return id.GetSampleName();})
        .DefinePerSample("do_ewk_corr", [](unsigned int slot, const RSampleInfo &id) { return id.GetI("do_ewk_corr");})
        .Define("isData", "kind == \"data\"")
        .Define("is2016", "year == \"2016preVFP\" || year == \"2016postVFP\"")
        .Define("is2017", "year == \"2017\"")
        .Define("is2018", "year == \"2018\"")
        .Define("is2022", "year == \"2022Re-recoBCD\" || year == \"2022Re-recoE+PromptFG\"")
        .Define("is2023", "year == \"2023PromptC\" || year == \"2023PromptD\"")
        .Define("is2024", "year == \"2024Prompt\"")
        .Define("isRun2", "is2016 || is2017 || is2018")
        .Define("isRun3", "is2022 || is2023 || is2024")
        .Define("xsecweight", "isData ? 1 : 1000 * xsec * lumi / sumw")
        .Define("baseweight", "xsecweight * genWeight")
        .Define("weight", "baseweight");

}

// Extract sample kind from the JSON config file
// Returns the kind of the first sample found (assumes all samples in a config have the same kind)
std::string getCategoryFromConfig(const std::string& config_path) {
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(config_path, pt);

    // Navigate to samples and get the first sample's kind
    for (const auto& sample : pt.get_child("samples")) {
        return sample.second.get<std::string>("metadata.kind");
    }
    return "";
}

RNode removeDuplicates(RNode df){
    return df.Filter(FilterOnePerKind(), {"run", "luminosityBlock", "event"}, "REMOVED DUPLICATES");
}

RNode applyObjectMask(RNode df, const std::string& maskName, const std::string& objectName) {
    auto columnNames = df.GetColumnNames();
    for (const auto& colName : columnNames) {
        if (colName.starts_with(objectName + "_")) {
            df = df.Redefine(colName, colName + "[" + maskName + "]");
        }
    }
    df = df.Redefine("n" + objectName, "Sum(" + maskName + ")");
    return df;
}

RNode applyObjectMaskNewAffix(RNode df, const std::string &maskName, const std::string &objectName, const std::string &newAffix)
{
    auto columnNames = df.GetColumnNames();
    for (const auto &colName : columnNames)
    {
        if (colName.starts_with(objectName + "_"))
        {
            std::string suffix = colName.substr(objectName.size() + 1);
            std::string newCol = newAffix + "_" + suffix;
            df = df.Define(newCol, colName + "[" + maskName + "]");
        }
    }
    df = df.Define("n" + newAffix, "Sum(" + maskName + ")");
    return df;
}

/*
############################################
LUMIMASK - GOLDEN JSON
############################################
*/

bool operator< ( const lumiMask::LumiBlockRange& lh, const lumiMask::LumiBlockRange& rh )
{
    return ( lh.run() == rh.run() ) ? ( lh.lastLumi() < rh.firstLumi() ) : lh.run() < rh.run();
}

lumiMask lumiMask::fromJSON(const std::vector<std::string>& files, lumiMask::Run firstRun, lumiMask::Run lastRun)
{
  const bool noRunFilter = ( firstRun == 0 ) && ( lastRun == 0 );

  std::vector<lumiMask::LumiBlockRange> accept;

  for ( const auto& file : files ) {
    boost::property_tree::ptree ptree;
    boost::property_tree::read_json(file, ptree);
    for ( const auto& runEntry : ptree ) {
        const lumiMask::Run run = std::stoul(runEntry.first);
        if ( noRunFilter || ( ( firstRun <= run ) && ( run <= lastRun ) ) ) {
        for ( const auto& lrEntry : runEntry.second ) {
            const auto lrNd = lrEntry.second;
            const lumiMask::LumiBlock firstLumi = std::stoul(lrNd.begin()->second.data());
            const lumiMask::LumiBlock lastLumi  = std::stoul((++lrNd.begin())->second.data());
            accept.emplace_back(run, firstLumi, lastLumi);
        }
        }
    }
  }
  
  return lumiMask(accept);
}

/*
############################################
SELECTION UTILS
############################################
*/

float fdR(float eta1, float phi1, float eta2, float phi2) {
    return ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
}

RVec<float> VdR(const RVec<float>& vec_eta, const RVec<float>& vec_phi, float obj_eta, float obj_phi) {
    RVec<float> out(vec_eta.size());
    if (obj_eta == -999 || obj_phi == -999) {
        std::fill(out.begin(), out.end(), 1.0f);
        return out;
    }
    for (size_t i = 0; i < vec_eta.size(); i++) {
        out[i] = ROOT::VecOps::DeltaR(vec_eta[i], obj_eta, vec_phi[i], obj_phi);
    }
    return out;
}

RVec<float> VVdR(const RVec<float>& vec_eta1, const RVec<float>& vec_phi1, const RVec<float>& vec_eta2, const RVec<float>& vec_phi2) {
    if (vec_eta1.empty()) {
        return RVec<float>();
    }
    if (vec_eta2.empty()) {
        return RVec<float>(vec_eta1.size(), 999.0f);
    }
    RVec<float> out(vec_eta1.size());
    for (size_t i = 0; i < vec_eta1.size(); i++) {
        float mindR = 999.;
        for (size_t j = 0; j < vec_eta2.size(); j++) {
            float dR = ROOT::VecOps::DeltaR(vec_eta1[i], vec_eta2[j], vec_phi1[i], vec_phi2[j]);
            if (dR < mindR) {
                mindR = dR;
            }
        }
        out[i] = mindR;
    }
    return out;
}

float fInvariantMass(float obj1_pt, float obj1_eta, float obj1_phi, float obj1_mass, 
                    float obj2_pt, float obj2_eta, float obj2_phi, float obj2_mass) {
    TLorentzVector obj1, obj2;
    obj1.SetPtEtaPhiM(obj1_pt, obj1_eta, obj1_phi, obj1_mass);
    obj2.SetPtEtaPhiM(obj2_pt, obj2_eta, obj2_phi, obj2_mass);
    return (obj1 + obj2).M();
}

RVec<float> VInvariantMass(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                          const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass[i] = (obj1 + obj2).M();
    }
    return invMass;
}

RVec<float> VInvariantPt(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                        const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPt(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPt[i] = (obj1 + obj2).Pt();
    }
    return invPt;
}

RVec<float> VInvariantPhi(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                         const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPhi(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPhi[i] = (obj1 + obj2).Phi();
    }
    return invPhi;
}

RVec<float> VTransverseMass(const RVec<float>& vec_pt, const RVec<float>& vec_phi, float obj_pt, float obj_phi) {
    RVec<float> mt(vec_pt.size());
    for (size_t i = 0; i < vec_pt.size(); i++) {
        mt[i] = std::sqrt(2 * vec_pt[i] * obj_pt * (1 - std::cos(ROOT::VecOps::DeltaPhi(vec_phi[i], obj_phi))));
    }
    return mt;
}

RVec<float> VVInvariantPt(const RVec<float>& pt1, const RVec<float>& eta1, const RVec<float>& phi1, const RVec<float>& m1,
                         const RVec<float>& pt2, const RVec<float>& eta2, const RVec<float>& phi2, const RVec<float>& m2) {
    RVec<float> ptjj;
    for (size_t i = 0; i < pt1.size(); ++i) {
        if (pt1[i] < 0 || pt2[i] < 0) {
            ptjj.push_back(-999.0f);
            continue;
        }
        auto vec1 = ROOT::Math::PtEtaPhiMVector(pt1[i], eta1[i], phi1[i], m1[i]);
        auto vec2 = ROOT::Math::PtEtaPhiMVector(pt2[i], eta2[i], phi2[i], m2[i]);
        ptjj.push_back((vec1 + vec2).Pt());
    }
    return ptjj;
}

RVec<float> VVInvariantMass(const RVec<float>& pt1, const RVec<float>& eta1, const RVec<float>& phi1, const RVec<float>& m1,
                                   const RVec<float>& pt2, const RVec<float>& eta2, const RVec<float>& phi2, const RVec<float>& m2) {
    RVec<float> invariant_mass;
    for (size_t i = 0; i < pt1.size(); ++i) {
        if (pt1[i] < 0 || pt2[i] < 0) {
            invariant_mass.push_back(-999.0f);
            continue;
        }
        auto vec1 = ROOT::Math::PtEtaPhiMVector(pt1[i], eta1[i], phi1[i], m1[i]);
        auto vec2 = ROOT::Math::PtEtaPhiMVector(pt2[i], eta2[i], phi2[i], m2[i]);
        invariant_mass.push_back((vec1 + vec2).M());
    }
    return invariant_mass;
}

RVec<float> VVDeltaR(const RVec<float>& eta1, const RVec<float>& phi1, const RVec<float>& eta2, const RVec<float>& phi2) {
    RVec<float> dR;
    for (size_t i = 0; i < eta1.size(); ++i) {
        if (eta1[i] < -900 || eta2[i] < -900) {
            dR.push_back(999.0f);
            continue;
        }
        float dphi = std::abs(phi1[i] - phi2[i]);
        if (dphi > M_PI) dphi = 2 * M_PI - dphi;
        float deta = eta1[i] - eta2[i];
        dR.push_back(std::sqrt(deta * deta + dphi * dphi));
    }
    return dR;
}

// Return for each ak4 jet, the dR from the closest ak8 jet
RVec<float> dRfromClosestJet(const RVec<float>& ak4_eta, const RVec<float>& ak4_phi, const RVec<float>& ak8_eta, const RVec<float>& ak8_phi) {
    RVec<float> vec_minDR = {};
    for (size_t i = 0; i < ak4_eta.size(); i++)
    {
        float mindR = 999.;
        for (size_t j = 0; j < ak8_eta.size(); j++)
        {
            float dR = ROOT::VecOps::DeltaR(ak4_eta.at(i), ak8_eta.at(j), ak4_phi.at(i), ak8_phi.at(j));
            if (dR < mindR) {
                mindR = dR;
            }
        }
        vec_minDR.push_back(mindR);
    }
    return vec_minDR;
}

RVec<RVec<int>> getJetPairs(const RVec<float>& goodJets) {
    if (goodJets.size() >= 2) {
        return ROOT::VecOps::Combinations(goodJets, 2);
    } else {
        RVec<RVec<int>> result;
        result.emplace_back(RVec<int>{999});
        result.emplace_back(RVec<int>{999});
        return result;
    }
}

RVec<int> findJetPairWithMaxDeltaEta(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass) {
    // find pair of jets with max delta eta
    RVec<int> good_jet_idx = {};
    RVec<float> Jet_Pt = {};
    for (size_t i = 0; i < Jet_pt.size(); i++) {
        TLorentzVector jet;
        jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        Jet_Pt.push_back(jet.Pt());
    }
    int Nvbfjet1 = -1;
    int Nvbfjet2 = -1;
    float maxvbfjetdeta = 0;
    for (size_t i = 0; i < Jet_eta.size(); i++) {
        for (size_t j = i+1; j < Jet_eta.size(); j++) {
            float deta = std::abs(Jet_eta[i] - Jet_eta[j]);
            if (deta > maxvbfjetdeta) {
                maxvbfjetdeta = deta;
                Nvbfjet1 = i;
                Nvbfjet2 = j;
            }
        }
    }
    if (Jet_Pt[Nvbfjet1] > Jet_Pt[Nvbfjet2]) {
        good_jet_idx.push_back(Nvbfjet1);
        good_jet_idx.push_back(Nvbfjet2);
    }
    else {
        good_jet_idx.push_back(Nvbfjet2);
        good_jet_idx.push_back(Nvbfjet1);
    }
    return good_jet_idx;
}

RVec<float> VBSBDTInfer(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass, bool isRun2) {
    if (Jet_pt.size() < 2) {
        return RVec<float>{-1, -1, -1};
    }
    auto combination_idxs = ROOT::VecOps::Combinations(Jet_pt, 2);

    auto jet1_pt = ROOT::VecOps::Take(Jet_pt, combination_idxs[0]);
    auto jet1_eta = ROOT::VecOps::Take(Jet_eta, combination_idxs[0]);
    auto jet1_phi = ROOT::VecOps::Take(Jet_phi, combination_idxs[0]);
    auto jet1_mass = ROOT::VecOps::Take(Jet_mass, combination_idxs[0]);
    auto jet2_pt = ROOT::VecOps::Take(Jet_pt, combination_idxs[1]);
    auto jet2_eta = ROOT::VecOps::Take(Jet_eta, combination_idxs[1]);
    auto jet2_phi = ROOT::VecOps::Take(Jet_phi, combination_idxs[1]);
    auto jet2_mass = ROOT::VecOps::Take(Jet_mass, combination_idxs[1]);
    auto detajj = ROOT::VecOps::abs(jet1_eta - jet2_eta);

    auto pt_m_jj = [](const RVec<float>& jet1_pt, const RVec<float>& jet1_eta, const RVec<float>& jet1_phi, const RVec<float>& jet1_mass, 
                        const RVec<float>& jet2_pt, const RVec<float>& jet2_eta, const RVec<float>& jet2_phi, const RVec<float>& jet2_mass) {
        RVec<float> pt_jj;
        RVec<float> m_jj;
        for (size_t i = 0; i < jet1_pt.size(); ++i) {
            auto v_jj = ROOT::Math::PtEtaPhiMVector(jet1_pt[i], jet1_eta[i], jet1_phi[i], jet1_mass[i]) + ROOT::Math::PtEtaPhiMVector(jet2_pt[i], jet2_eta[i], jet2_phi[i], jet2_mass[i]);
            pt_jj.push_back(v_jj.Pt());
            m_jj.push_back(v_jj.M());
        }
        return std::make_pair(pt_jj, m_jj);
    };

    auto [ptjj, mjj] = pt_m_jj(jet1_pt, jet1_eta, jet1_phi, jet1_mass, jet2_pt, jet2_eta, jet2_phi, jet2_mass);
    auto dphijj = ROOT::VecOps::DeltaPhi(jet1_phi, jet2_phi);

    RVec<float> scores;
    float score;
    for (size_t i = 0; i < mjj.size(); i++) {
        score = bdt.Compute({
                    jet1_pt[i], jet2_pt[i],
                    jet1_eta[i], jet2_eta[i],
                    jet1_phi[i], jet2_phi[i],
                    jet1_mass[i], jet2_mass[i],
                    ptjj[i], detajj[i], 
                    dphijj[i], mjj[i]
                })[0];
        scores.push_back(score);
    }
    auto max_score_idx = std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));
    if (scores.size() > 0) {
        return RVec<float>{static_cast<float>(combination_idxs[0][max_score_idx]), 
                         static_cast<float>(combination_idxs[1][max_score_idx]),
                         scores[max_score_idx]};
    }
    return RVec<float>{-1, -1, -1};
}

/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &outdir, bool spanet_training) {
    std::string output_dir = "";
    if (spanet_training) {
        output_dir = outdir + "/spanet_training/";
    }
    else {
        output_dir = outdir;
    }

    std::filesystem::path directory_path(output_dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        std::cerr << "Output directory: " << directory_path << std::endl;
    }
    // Try to create the directory and any missing parent directories
    else if (std::filesystem::create_directories(directory_path)) {
        std::cout << "Created output directory : " << directory_path << std::endl;
    }
    else {
        std::cerr << "Failed to create output directory: " << directory_path << std::endl;
        std::exit(EXIT_FAILURE); 
    }

    return directory_path;
}

void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool isSig, bool dumpInput, bool storeHLT)
{
    auto ColNames = df.GetDefinedColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("run");
    final_variables.push_back("luminosityBlock");
    final_variables.push_back("event");

    // do not store branches that start with "_" nor raw NanoAOD collections
    for (auto &&ColName : ColNames) {
        if (ColName.starts_with("_") || ColName.starts_with("Jet") || ColName.starts_with("FatJet_") || ColName.starts_with("Electron_") || ColName.starts_with("Muon_"))
        {
            continue;
        }
        final_variables.push_back(ColName);
    }

    // Optionally store HLT branches from input NanoAOD, providing default values
    // for branches that may not exist in all files of a multi-file chain
    if (storeHLT) {
        auto allColNames = df.GetColumnNames();
        for (auto &&colName : allColNames) {
            if (colName.starts_with("HLT_") &&
                std::find(final_variables.begin(), final_variables.end(), colName) == final_variables.end()) {
                df = df.DefaultValueFor(colName, (bool)false);
                final_variables.push_back(colName);
            }
        }
    }

    if (isSig) {
        final_variables.push_back("LHEReweightingWeight");
        final_variables.push_back("nLHEReweightingWeight");
    }

    // store all columns from input nanoAOD tree
    if (dumpInput) {
        auto nanoColNames = df.GetColumnNames();
        for (auto &&colName : nanoColNames) {
            if ((std::find(final_variables.begin(), final_variables.end(), colName) == final_variables.end()) &&
                (colName.find("HLT") == std::string::npos) && (colName.find("L1") == std::string::npos)) {
                final_variables.push_back(colName);
            }
        }
    }

    std::string outputFile = outputDir + "/" + outputFileName + ".root";
    df.Snapshot("Events", outputFile, final_variables);

    // RDataFrame::Snapshot() in multi-threaded mode does not write a TTree when
    // 0 events pass the filters, producing a ROOT file with no keys.  Ensure the
    // output always contains an "Events" TTree so downstream code can open it.
    {
        TFile f(outputFile.c_str(), "UPDATE");
        if (!f.Get("Events")) {
            TTree t("Events", "Events");
            t.Write();
        }
        f.Close();
    }

    std::cout << " -> Stored output file: " << outputFile << std::endl;
}


void saveSpanetSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName)
{
    auto ColNames = df.GetColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");

    for (auto &&ColName : ColNames) {
        if (ColName.starts_with("jet_") ||
            ColName.starts_with("fatjet_") ||
            ColName.starts_with("PuppiMET_") ||
            ColName.starts_with("GenPart_") ||
            ColName.starts_with("lepton_") ||
            ColName.starts_with("gen_") ||
            ColName.starts_with("truth_")) {
                final_variables.push_back(ColName);
            }
    }

    std::string outputFile = outputDir + "/" + outputFileName + "_spanet_training_data.root";
    df.Snapshot("Events", outputFile, final_variables);

    {
        TFile f(outputFile.c_str(), "UPDATE");
        if (!f.Get("Events")) {
            TTree t("Events", "Events");
            t.Write();
        }
        f.Close();
    }

    std::cout << " -> Stored output file: " << outputFile << std::endl;
}
