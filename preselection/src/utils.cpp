#include "utils.h"

/*
############################################
RDF UTILS
############################################
*/

RNode defineMetadata(RNode df) {
    return df.DefinePerSample("xsec", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("xsec");})
        .DefinePerSample("lumi", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("lumi");})
        .DefinePerSample("nevents", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("nevents");})
        .DefinePerSample("sample_category", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_category");})
        .DefinePerSample("sample_type", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_type");})
        .DefinePerSample("sample_year", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_year");})
        .Define("xsec_weight", "1000 * xsec * lumi / nevents")
        .Define("isData", "sample_category == \"data\"");
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
CUTFLOW
############################################
*/

Cutflow::Cutflow(RNode df) : _df(df){
    for (auto colName : _df.GetDefinedColumnNames()) {
        if (colName.starts_with("_cut")) {
            _cuts.push_back(colName);
        }
    }

    if (_cuts.size() != 0) { 
        _df = _df.Define("weight2", "weight * weight");
        _cutflow.push_back(std::make_pair(_df.Sum<double>("weight"), _df.Sum<double>("weight2")));
        std::string cut_string = "";
        for (size_t i = 0; i < _cuts.size(); i++) {
            if (i == 0) {
                cut_string = _cuts[i];
            } else {
                cut_string += " && " + _cuts[i];
            }
            _cutflow.push_back(std::make_pair(_df.Filter(cut_string).Sum<double>("weight"), _df.Filter(cut_string).Sum<double>("weight2")));
        }        
    }
}

void Cutflow::Print(std::string output_file) {
    if (_cuts.size() == 0) {
        std::cout << "No cuts _cut_* defined" << std::endl;
        return;
    }
    tabulate::Table table;
    table.add_row({"Cut", "Count", "Percent Change"});
    {
        std::stringstream ss_count, ss_error;
        ss_count << std::fixed << std::setprecision(3) << _cutflow[0].first.GetValue();
        ss_error << std::fixed << std::setprecision(3) << std::sqrt(_cutflow[0].second.GetValue());
        table.add_row({
            "initialCount",
            ss_count.str() + " +/- " + ss_error.str(),
            std::string("")
        });
    }
    
    for (size_t i = 0; i < _cuts.size(); i++) {
        std::stringstream ss_count, ss_error, ss_percent;
        ss_count << std::fixed << std::setprecision(3) << _cutflow[i + 1].first.GetValue();
        ss_error << std::fixed << std::setprecision(3) << std::sqrt(_cutflow[i + 1].second.GetValue());
        ss_percent << std::fixed << std::setprecision(3) << 100 * (_cutflow[i].first.GetValue() - _cutflow[i + 1].first.GetValue()) / _cutflow[i].first.GetValue();
        
        table.add_row({
            _cuts[i],
            ss_count.str() + " +/- " + ss_error.str(),
            ss_percent.str()
        });
    }

    table.column(1).format().font_align(tabulate::FontAlign::center);
    table.column(2).format().font_align(tabulate::FontAlign::center);

    if (!output_file.empty()) {
        std::ofstream out(output_file);
        out << table;
        out.close();
    }
    else {
        std::cout << table << std::endl;
    }
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

RVec<RVec<int>> getVBSPairs(const RVec<int>& goodJets, const RVec<float>& jet_var) {
    if (Sum(goodJets) >= 2) {
        return ROOT::VecOps::Combinations(jet_var, 2);
    } else {
    // Create properly matched return type: vector of vector
        RVec<RVec<int>> result;
        // Add two empty vectors
        result.emplace_back(RVec<int>{-999});
        result.emplace_back(RVec<int>{-999});
        return result;
    }
}

RVec<int> VBS_MaxEtaJJ(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass) {
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

int num_hadronic_gauge_bosons(RVec<int> pdgId, RVec<short> motherIdx) {
    // particles at idx 2 and 3 are hadronic gauge bosons
    auto daughterIdx = [&] (int idx) {
        RVec<int> daughters;
        for (size_t i = 0; i < pdgId.size(); ++i) {
            if (motherIdx[i] == idx) {
                daughters.push_back(i);
            }
        }
        return daughters;
    };
                       
    std::function<bool(int)> is_hadronic = [&] (int idx) -> bool {
        if (idx < 0 || idx >= pdgId.size()) return false;
        if (pdgId[idx] == 23 || pdgId[idx] == 24 || pdgId[idx] == -24) {
            auto daughters = daughterIdx(idx);
            for (int daughter : daughters) {
                if (pdgId[daughter] == 23 || pdgId[daughter] == 24 || pdgId[daughter] == -24) {
                    if (is_hadronic(daughter)) {
                        return true;
                    }
                } else if (abs(pdgId[daughter]) < 6) {
                    return true;
                } else if (abs(pdgId[daughter]) == 11 || abs(pdgId[daughter]) == 13 || abs(pdgId[daughter]) == 15) {
                    return false;
                }
            }
        }
        return false;
    };

    // return 0 if both are leptonic, 1 if only genpart 2 is hadronic, 2 if only genpart 3 is hadronic, 3 if both are hadronic
    int count = 0;
    if (is_hadronic(2)) count += 1;
    if (is_hadronic(3)) count += 2;
    return count;
}

int get_higgs_boson_idx(RVec<int>& pdgId, RVec<short>& motherIdx) {
    const int bId = 5;
    const int hId = 25;
    RVec<int> bIndices, bBarIndices;

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] == bId) bIndices.push_back(i);
        else if (pdgId[i] == -bId) bBarIndices.push_back(i);
    }

    auto traceToTopHiggs = [&](int idx) {
        while (idx >= 0 && pdgId[idx] == hId && motherIdx[idx] >= 0 && pdgId[motherIdx[idx]] == hId) {
            idx = motherIdx[idx];
        }
        return (pdgId[idx] == hId) ? idx : -1;
    };

    for (int bIdx : bIndices) {
        for (int bBarIdx : bBarIndices) {
            int motherB = motherIdx[bIdx];
            int motherBBar = motherIdx[bBarIdx];

            if (motherB == -1 || motherBBar == -1) continue;
            if (motherB != motherBBar) continue;
            if (pdgId[motherB] != hId) continue;

            int topHiggsIdx = traceToTopHiggs(motherB);
            if (topHiggsIdx == 4) {
                return motherB;
            }
        }
    }
    return -1;
}

int find_matching_jet(RVec<float> dR_values, RVec<int> excluded_indices) {
    int max_jets = 10;
    const float dR_cut = 0.4f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break;
        if (idx >= max_jets) continue;
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet(RVec<float> dR_values, RVec<int> excluded_indices) {
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break;
        if (idx >= max_fatjets) continue;
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}

int find_matching_jet_conditional(int check_idx, RVec<float> dR_values, RVec<int> excluded_indices) {
    int max_jets = 10;
    if (check_idx == -1 || dR_values.empty()) return -1;

    const float dR_cut = 0.4f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break;
        if (idx >= max_jets) continue;
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet_conditional(int check_idx, RVec<float> dR_values, RVec<int> excluded_indices) {
    int max_fatjets = 3;
    if (check_idx == -1 || dR_values.empty()) return -1;
    
    const float dR_cut = 0.8f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break;
        if (idx >= max_fatjets) continue;
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}

RVec<float> get_dR(float eta1, float phi1, RVec<float> eta2, RVec<float> phi2) {
    RVec<float> dR;
    for (size_t i = 0; i < eta2.size(); ++i) {
        dR.push_back(ROOT::VecOps::DeltaR(eta1, eta2[i], phi1, phi2[i]));
    }
    return dR;
}

RVec<float> get_dR_conditional(int idx, float eta1, float phi1, RVec<float> eta2, RVec<float> phi2) {
    if (idx == -1) {
        return RVec<float>();
    }
    RVec<float> dR;
    for (size_t i = 0; i < eta2.size(); ++i) {
        dR.push_back(ROOT::VecOps::DeltaR(eta1, eta2[i], phi1, phi2[i]));
    }
    return dR;
}

std::vector<int> assign_all_objects(
    std::vector<std::vector<float>> vbs_assignment,
    std::vector<std::vector<float>> h_assignment, 
    std::vector<std::vector<float>> bh_assignment,
    std::vector<std::vector<float>> v1_assignment,
    std::vector<std::vector<float>> v2_assignment,
    std::vector<std::vector<float>> bv1_assignment,
    std::vector<std::vector<float>> bv2_assignment,
    float vbs_detection,
    float h_detection,
    float bh_detection,
    float v1_detection,
    float v2_detection,
    float bv1_detection,
    float bv2_detection,
    RVec<float> Jet_eta,
    RVec<float> Jet_phi,
    RVec<float> FatJet_eta,
    RVec<float> FatJet_phi
) {
    std::vector<int> result(11, -1);
    std::vector<int> assigned_jets;
    std::vector<int> assigned_fatjets;
    
    std::vector<std::pair<float, int>> detection_order = {
        {vbs_detection, 0},
        {h_detection, 1},
        {bh_detection, 2},
        {v1_detection, 3},
        {v2_detection, 4},
        {bv1_detection, 5},
        {bv2_detection, 6} 
    };
    
    std::sort(detection_order.begin(), detection_order.end(), 
        [](const auto& a, const auto& b) {
            return a.first > b.first;
        });
    
    auto checkJetOverlap = [&](int jet_idx, const std::vector<int>& assigned) -> bool {
        if (jet_idx < 0 || jet_idx >= Jet_eta.size()) return true;
        return std::find(assigned.begin(), assigned.end(), jet_idx) != assigned.end();
    };
    
    auto checkFatJetOverlap = [&](int fatjet_idx, const std::vector<int>& assigned) -> bool {
        if (fatjet_idx < 0 || fatjet_idx >= FatJet_eta.size()) return true;
        return std::find(assigned.begin(), assigned.end(), fatjet_idx) != assigned.end();
    };
    
    auto checkDeltaR = [&](int idx1, int idx2, bool is_fatjet1, bool is_fatjet2) -> bool {
        if (idx1 < 0 || idx2 < 0) return true;
        
        float eta1 = is_fatjet1 ? FatJet_eta[idx1] : Jet_eta[idx1];
        float phi1 = is_fatjet1 ? FatJet_phi[idx1] : Jet_phi[idx1];
        float eta2 = is_fatjet2 ? FatJet_eta[idx2] : Jet_eta[idx2];
        float phi2 = is_fatjet2 ? FatJet_phi[idx2] : Jet_phi[idx2];
        
        float dR = ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
        if (is_fatjet1 || is_fatjet2) {
            return dR >= 0.8;
        }
        return dR >= 0.4;
    };
    
    auto checkAllOverlaps = [&](int candidate_idx, bool is_fatjet) -> bool {
        for (int assigned_jet : assigned_jets) {
            if (!checkDeltaR(candidate_idx, assigned_jet, is_fatjet, false)) {
                return false;
            }
        }
        for (int assigned_fatjet : assigned_fatjets) {
            if (!checkDeltaR(candidate_idx, assigned_fatjet, is_fatjet, true)) {
                return false;
            }
        }
        return true;
    };
    
    for (const auto& [prob, obj_type] : detection_order) {
        switch (obj_type) {
            case 0: {
                for (size_t i = 0; i < vbs_assignment.size(); i++) {
                    int j1_candidate = static_cast<int>(vbs_assignment[i][1]);
                    int j2_candidate = static_cast<int>(vbs_assignment[i][2]);
                    
                    if (!checkJetOverlap(j1_candidate, assigned_jets) && 
                        !checkJetOverlap(j2_candidate, assigned_jets) &&
                        j1_candidate != j2_candidate &&
                        checkDeltaR(j1_candidate, j2_candidate, false, false) &&
                        checkAllOverlaps(j1_candidate, false) &&
                        checkAllOverlaps(j2_candidate, false)) {
                        
                        result[0] = j1_candidate;
                        result[1] = j2_candidate;
                        assigned_jets.push_back(j1_candidate);
                        assigned_jets.push_back(j2_candidate);
                        break;
                    }
                }
                break;
            }
            case 1: {
                for (size_t i = 0; i < h_assignment.size(); i++) {
                    int j1_candidate = static_cast<int>(h_assignment[i][1]);
                    int j2_candidate = static_cast<int>(h_assignment[i][2]);
                    
                    if (!checkJetOverlap(j1_candidate, assigned_jets) && 
                        !checkJetOverlap(j2_candidate, assigned_jets) &&
                        j1_candidate != j2_candidate &&
                        checkDeltaR(j1_candidate, j2_candidate, false, false) &&
                        checkAllOverlaps(j1_candidate, false) &&
                        checkAllOverlaps(j2_candidate, false)) {
                        
                        result[2] = j1_candidate;  // h1_idx
                        result[3] = j2_candidate;  // h2_idx
                        assigned_jets.push_back(j1_candidate);
                        assigned_jets.push_back(j2_candidate);
                        break;
                    }
                }
                break;
            }
            case 2: {
                for (size_t i = 0; i < bh_assignment.size(); i++) {
                    int fatjet_candidate = static_cast<int>(bh_assignment[i][1]);
                    
                    if (!checkFatJetOverlap(fatjet_candidate, assigned_fatjets) &&
                        checkAllOverlaps(fatjet_candidate, true)) {
                        
                        result[4] = fatjet_candidate;  // bh_idx
                        assigned_fatjets.push_back(fatjet_candidate);
                        break;
                    }
                }
                break;
            }
            case 3: {
                for (size_t i = 0; i < v1_assignment.size(); i++) {
                    int j1_candidate = static_cast<int>(v1_assignment[i][1]);
                    int j2_candidate = static_cast<int>(v1_assignment[i][2]);
                    
                    if (!checkJetOverlap(j1_candidate, assigned_jets) && 
                        !checkJetOverlap(j2_candidate, assigned_jets) &&
                        j1_candidate != j2_candidate &&
                        checkDeltaR(j1_candidate, j2_candidate, false, false) &&
                        checkAllOverlaps(j1_candidate, false) &&
                        checkAllOverlaps(j2_candidate, false)) {
                        
                        result[5] = j1_candidate;  // v1_j1_idx
                        result[6] = j2_candidate;  // v1_j2_idx
                        assigned_jets.push_back(j1_candidate);
                        assigned_jets.push_back(j2_candidate);
                        break;
                    }
                }
                break;
            }
            case 4: {
                for (size_t i = 0; i < v2_assignment.size(); i++) {
                    int j1_candidate = static_cast<int>(v2_assignment[i][1]);
                    int j2_candidate = static_cast<int>(v2_assignment[i][2]);
                    
                    if (!checkJetOverlap(j1_candidate, assigned_jets) && 
                        !checkJetOverlap(j2_candidate, assigned_jets) &&
                        j1_candidate != j2_candidate &&
                        checkDeltaR(j1_candidate, j2_candidate, false, false) &&
                        checkAllOverlaps(j1_candidate, false) &&
                        checkAllOverlaps(j2_candidate, false)) {
                        
                        result[7] = j1_candidate;  // v2_j1_idx
                        result[8] = j2_candidate;  // v2_j2_idx
                        assigned_jets.push_back(j1_candidate);
                        assigned_jets.push_back(j2_candidate);
                        break;
                    }
                }
                break;
            }
            case 5: {
                for (size_t i = 0; i < bv1_assignment.size(); i++) {
                    int fatjet_candidate = static_cast<int>(bv1_assignment[i][1]);
                    
                    if (!checkFatJetOverlap(fatjet_candidate, assigned_fatjets) &&
                        checkAllOverlaps(fatjet_candidate, true)) {
                        
                        result[9] = fatjet_candidate;  // bv1_idx
                        assigned_fatjets.push_back(fatjet_candidate);
                        break;
                    }
                }
                break;
            }
            case 6: {
                for (size_t i = 0; i < bv2_assignment.size(); i++) {
                    int fatjet_candidate = static_cast<int>(bv2_assignment[i][1]);
                    
                    if (!checkFatJetOverlap(fatjet_candidate, assigned_fatjets) &&
                        checkAllOverlaps(fatjet_candidate, true)) {
                        
                        result[10] = fatjet_candidate;  // bv2_idx
                        assigned_fatjets.push_back(fatjet_candidate);
                        break;
                    }
                }
                break;
            }
        }
    }
    
    return result;
}

/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &ana, const std::string &output_subdir) {
    // If on UAF and the USER environment variable is defined, store output on ceph
    const char* userEnv = getenv("USER");
    std::string storage_dir = "/data/userdata/";
    std::string output_dir = "./";
    if (userEnv != nullptr && std::filesystem::exists(storage_dir) && std::filesystem::is_directory(storage_dir)) {
        output_dir = storage_dir + std::string(userEnv) + "/vbsvvhAnalysis/preselection/" + ana + "/";
    }

    if (!output_subdir.empty()) {
        output_dir += "/" + output_subdir + "/";
    }

    std::filesystem::path directory_path(output_dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        std::cerr << "Output directory already exists: " << directory_path << std::endl;
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

void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool isData, bool dumpInput)
{
    auto ColNames = df.GetDefinedColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");

    for (auto &&ColName : ColNames) {
        if (ColName.starts_with("_")) {
            if (ColName.starts_with("_cut")) {
            } else
                continue;
        }
        final_variables.push_back(ColName);
    }

    // add LHE info
    if (!isData) {
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

    std::string outputFile = outputDir + "/" + outputFileName + ".root" df.Snapshot("Events", output_file, final_variables);
    std::cout << " -> Stored output file: " << outputFile << std::endl;
}