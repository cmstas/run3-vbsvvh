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

RVec<int> get_higgs_boson_idx(RVec<int>& pdgId, RVec<int>& status, RVec<short>& motherIdx) {
    // output [h1, b1, b2]
    RVec<int> result = {-1, -1, -1};
    
    int higgs_idx = -1;
    RVec<int> hdecay_idx;
    
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && part_pdgId == 25) {
            higgs_idx = igen;
            break;
        }
    }
    
    if (higgs_idx == -1) return result;
    
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int mother_idx = motherIdx[igen];
        if (mother_idx > 0 && mother_idx < pdgId.size()) {
            int mother_pdgId = pdgId[mother_idx];
            int part_pdgId = pdgId[igen];
            
            if (mother_pdgId == 25 && part_pdgId != 25) {
                hdecay_idx.push_back(igen);
                    }
        }
    }
    
    result[0] = higgs_idx;
    if (hdecay_idx.size() >= 1) result[1] = hdecay_idx[0];
    if (hdecay_idx.size() >= 2) result[2] = hdecay_idx[1];
    
    return result;
                }

int findLastIndex(int current_idx, int current_pdgId, RVec<int>& pdgId, RVec<short>& motherIdx) {
    int outIdx = current_idx;
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == current_idx) {
            if (part_pdgId == current_pdgId) {
                outIdx = findLastIndex(igen, part_pdgId, pdgId, motherIdx);
            }
        }
    }
    return outIdx;
}

RVec<int> get_v_boson_idx(RVec<int>& pdgId, RVec<int>& status, RVec<short>& motherIdx) {
    // output [v1, v1d1, v1d2, v2, v2d1, v2d2]
    RVec<int> result = {-1, -1, -1, -1, -1, -1};
    
    RVec<int> firstVs_idx;

    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && (part_pdgId == 23 || abs(part_pdgId) == 24)) {
            firstVs_idx.push_back(igen);
        }
    }
    
    if (firstVs_idx.size() < 2) return result;
    
    RVec<int> hadronic_V_indices;
    RVec<int> leptonic_V_indices;
    
    for (size_t iV = 0; iV < firstVs_idx.size() && iV < 2; ++iV) {
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        int hadronic_daughters = 0;
        int leptonic_daughters = 0;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                    hadronic_daughters++;
                } else {
                    leptonic_daughters++;
                }
            }
        }
        
        if (hadronic_daughters > 0) {
            hadronic_V_indices.push_back(iV);
        } else if (leptonic_daughters > 0) {
            leptonic_V_indices.push_back(iV);
        }
    }
    
    if (hadronic_V_indices.size() > 0) {
        int iV = hadronic_V_indices[0];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) { // quarks only
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[0] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[1] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[2] = vdecays_idx[1];
    }
    
    if (hadronic_V_indices.size() >= 2) {
        int iV = hadronic_V_indices[1];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[3] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[4] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[5] = vdecays_idx[1];
    }
    
    return result;
}

RVec<int> get_vbs_quarks_idxs(RVec<int>& pdgId, RVec<int>& status, RVec<short>& motherIdx) {
    RVec<int> result = {-1, -1};
    RVec<int> vbsquarks_idx;
    
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];

        if (mother_idx == 0 && part_status == 23 && abs(part_pdgId) <= 6 && abs(part_pdgId) >= 1) {
            vbsquarks_idx.push_back(igen);
        }
    }
    
    if (vbsquarks_idx.size() >= 1) result[0] = vbsquarks_idx[0];
    if (vbsquarks_idx.size() >= 2) result[1] = vbsquarks_idx[1];
    
    return result;
}

int find_matching_jet(int target_idx, float target_eta, float target_phi, RVec<int> already_matched_jet_indices, RVec<int> already_matched_fatjet_indices, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> fatjet_eta, RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.4f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.4f;
                       
    if (target_idx < 0) {
        return -1;
    }
    
    // Calculate dR values for all jets
    RVec<float> dR_values;
    for (size_t i = 0; i < jet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, jet_eta[i], target_phi, jet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_jets) continue; // Skip indices beyond padding limit
        
        // Check if this jet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_jet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;
        
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(jet_eta[idx], jet_eta[matched_j_idx], jet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }

        // Check overlap with matched fatjets only
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(jet_eta[idx], fatjet_eta[matched_fj_idx], jet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }
        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet(int target_idx, float target_eta, float target_phi, RVec<int> already_matched_jet_indices, RVec<int> already_matched_fatjet_indices, RVec<float> jet_eta, RVec<float> jet_phi, RVec<float> fatjet_eta, RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.8f;
    
    if (target_idx < 0) {
        return -1;
    }
                       
    // Calculate dR values for all fatjets
    RVec<float> dR_values;
    for (size_t i = 0; i < fatjet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, fatjet_eta[i], target_phi, fatjet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_fatjets) continue; // Skip indices beyond padding limit
        
        // Check if this fatjet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_fatjet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;

        // Check overlap with matched jets only
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(fatjet_eta[idx], jet_eta[matched_j_idx], fatjet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }
        
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(fatjet_eta[idx], fatjet_eta[matched_fj_idx], fatjet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }

        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
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

    std::string outputFile = outputDir + "/" + outputFileName + ".root";
    df.Snapshot("Events", outputFile, final_variables);
    std::cout << " -> Stored output file: " << outputFile << std::endl;
}