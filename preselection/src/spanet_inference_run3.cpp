#include "spanet_inference_run3.h"
namespace SPANet {

SPANetInferenceRun3::SPANetInferenceRun3(const std::string &model_path, size_t batch_size)
    : SPANetInferenceBase(model_path, batch_size, MAX_AK4_JETS, AK4_FEATURES, MAX_AK8_JETS, AK8_FEATURES, MET_FEATURES, MAX_LEPTONS, LEPTON_FEATURES) {}

std::vector<std::shared_ptr<EventDataBase>> SPANetInferenceRun3::extractEventsFromDataFrame(RNode df) {
    auto ak4_pt_vec = df.Take<RVec<float>>("Jet_pt").GetValue();
    auto ak4_eta_vec = df.Take<RVec<float>>("Jet_eta").GetValue();
    auto ak4_phi_vec = df.Take<RVec<float>>("Jet_phi").GetValue();
    auto ak4_mass_vec = df.Take<RVec<float>>("Jet_mass").GetValue();
    auto ak4_isTightBTag_vec = df.Take<RVec<int>>("Jet_isTightBTag").GetValue();
    auto ak4_isMediumBTag_vec = df.Take<RVec<int>>("Jet_isMediumBTag").GetValue();
    auto ak4_isLooseBTag_vec = df.Take<RVec<int>>("Jet_isLooseBTag").GetValue();

    auto ak8_pt_vec = df.Take<RVec<float>>("FatJet_pt").GetValue();
    auto ak8_eta_vec = df.Take<RVec<float>>("FatJet_eta").GetValue();
    auto ak8_phi_vec = df.Take<RVec<float>>("FatJet_phi").GetValue();
    auto ak8_mass_vec = df.Take<RVec<float>>("FatJet_mass").GetValue();

    auto ak8_nConstituents_vec = df.Take<RVec<unsigned char>>("FatJet_nConstituents").GetValue();
    auto ak8_XbbScore_vec = df.Take<RVec<float>>("FatJet_globalParT3_Xbb").GetValue();
    auto ak8_XqqScore_vec = df.Take<RVec<float>>("FatJet_globalParT3_Xqq").GetValue();
    auto ak8_XccScore_vec = df.Take<RVec<float>>("FatJet_globalParT3_Xcc").GetValue();
    auto ak8_XcsScore_vec = df.Take<RVec<float>>("FatJet_globalParT3_Xcs").GetValue();
    auto ak8_XqcdScore_vec = df.Take<RVec<float>>("FatJet_globalParT3_QCD").GetValue();

    auto met_pt_vec = df.Take<float>("PuppiMET_pt").GetValue();
    auto met_phi_vec = df.Take<float>("PuppiMET_phi").GetValue();

    auto lepton_pt_vec = df.Take<RVec<float>>("Lepton_pt").GetValue();
    auto lepton_eta_vec = df.Take<RVec<float>>("Lepton_eta").GetValue();
    auto lepton_phi_vec = df.Take<RVec<float>>("Lepton_phi").GetValue();
    auto lepton_mass_vec = df.Take<RVec<float>>("Lepton_mass").GetValue();
    auto lepton_charge_vec = df.Take<RVec<int>>("Lepton_charge").GetValue();

    auto rdf_entry_vec = df.Take<ULong64_t>("rdfentry_").GetValue();
    auto rdf_slot_vec = df.Take<unsigned int>("rdfslot_").GetValue();

    std::vector<std::shared_ptr<EventDataBase>> events;
    size_t n_events = ak4_pt_vec.size();
    events.reserve(n_events);

    for (size_t i = 0; i < n_events; ++i) {
        auto event = std::make_shared<EventData>();

        event->ak4_pt.assign(ak4_pt_vec[i].begin(), ak4_pt_vec[i].end());
        event->ak4_eta.assign(ak4_eta_vec[i].begin(), ak4_eta_vec[i].end());
        event->ak4_phi.assign(ak4_phi_vec[i].begin(), ak4_phi_vec[i].end());
        event->ak4_mass.assign(ak4_mass_vec[i].begin(), ak4_mass_vec[i].end());
        event->ak4_isTightBTag.assign(ak4_isTightBTag_vec[i].begin(), ak4_isTightBTag_vec[i].end());
        event->ak4_isMediumBTag.assign(ak4_isMediumBTag_vec[i].begin(), ak4_isMediumBTag_vec[i].end());
        event->ak4_isLooseBTag.assign(ak4_isLooseBTag_vec[i].begin(), ak4_isLooseBTag_vec[i].end());

        event->ak8_pt.assign(ak8_pt_vec[i].begin(), ak8_pt_vec[i].end());
        event->ak8_eta.assign(ak8_eta_vec[i].begin(), ak8_eta_vec[i].end());
        event->ak8_phi.assign(ak8_phi_vec[i].begin(), ak8_phi_vec[i].end());
        event->ak8_mass.assign(ak8_mass_vec[i].begin(), ak8_mass_vec[i].end());
        event->ak8_nConstituents.assign(ak8_nConstituents_vec[i].begin(), ak8_nConstituents_vec[i].end());
        event->ak8_XbbScore.assign(ak8_XbbScore_vec[i].begin(), ak8_XbbScore_vec[i].end());
        event->ak8_XqqScore.assign(ak8_XqqScore_vec[i].begin(), ak8_XqqScore_vec[i].end());
        event->ak8_XccScore.assign(ak8_XccScore_vec[i].begin(), ak8_XccScore_vec[i].end());
        event->ak8_XcsScore.assign(ak8_XcsScore_vec[i].begin(), ak8_XcsScore_vec[i].end());
        event->ak8_XqcdScore.assign(ak8_XqcdScore_vec[i].begin(), ak8_XqcdScore_vec[i].end());

        event->met_pt = met_pt_vec[i];
        event->met_phi = met_phi_vec[i];

        event->lepton1_pt = lepton_pt_vec[i].empty() ? 0.0f : lepton_pt_vec[i][0];
        event->lepton1_eta = lepton_eta_vec[i].empty() ? 0.0f : lepton_eta_vec[i][0];
        event->lepton1_phi = lepton_phi_vec[i].empty() ? 0.0f : lepton_phi_vec[i][0];
        event->lepton1_mass = lepton_mass_vec[i].empty() ? 0.0f : lepton_mass_vec[i][0];
        event->lepton1_charge = lepton_charge_vec[i].empty() ? 0 : lepton_charge_vec[i][0];

        event->lepton2_pt = lepton_pt_vec[i].size() < 2 ? 0.0f : lepton_pt_vec[i][1];
        event->lepton2_eta = lepton_eta_vec[i].size() < 2 ? 0.0f : lepton_eta_vec[i][1];
        event->lepton2_phi = lepton_phi_vec[i].size() < 2 ? 0.0f : lepton_phi_vec[i][1];
        event->lepton2_mass = lepton_mass_vec[i].size() < 2 ? 0.0f : lepton_mass_vec[i][1];
        event->lepton2_charge = lepton_charge_vec[i].size() < 2 ? 0 : lepton_charge_vec[i][1];

        event->rdf_entry = rdf_entry_vec[i];
        event->rdf_slot = rdf_slot_vec[i];

        events.push_back(event);
    }

    return events;
}

void SPANetInferenceRun3::fillBatchTensors(const std::vector<std::shared_ptr<EventDataBase>>& events_base, size_t actual_batch_size) {
    std::fill(ak4_flat_jets_.begin(), ak4_flat_jets_.begin() + actual_batch_size * MAX_AK4_JETS * AK4_FEATURES, 0.0f);
    std::fill(ak8_flat_jets_.begin(), ak8_flat_jets_.begin() + actual_batch_size * MAX_AK8_JETS * AK8_FEATURES, 0.0f);
    std::fill(ak4_mask_char_.begin(), ak4_mask_char_.begin() + actual_batch_size * MAX_AK4_JETS, 0);
    std::fill(ak8_mask_char_.begin(), ak8_mask_char_.begin() + actual_batch_size * MAX_AK8_JETS, 0);

    for (size_t batch_idx = 0; batch_idx < actual_batch_size; ++batch_idx) {
        auto event_base = events_base[batch_idx];
        auto event = std::dynamic_pointer_cast<EventData>(event_base);

        const size_t max_ak4 = std::min<size_t>(MAX_AK4_JETS, event->ak4_pt.size());
        for (size_t i = 0; i < max_ak4; ++i) {
            ak4_mask_char_[batch_idx * MAX_AK4_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK4_JETS * AK4_FEATURES + i * AK4_FEATURES;
            ak4_flat_jets_[base_idx]     = std::log(event->ak4_mass[i] + 1.0f);
            ak4_flat_jets_[base_idx + 1] = std::log(event->ak4_pt[i] + 1.0f);
            ak4_flat_jets_[base_idx + 2] = event->ak4_eta[i];
            ak4_flat_jets_[base_idx + 3] = std::sin(event->ak4_phi[i]);
            ak4_flat_jets_[base_idx + 4] = std::cos(event->ak4_phi[i]);
            ak4_flat_jets_[base_idx + 5] = static_cast<float>(event->ak4_isTightBTag[i]);
            ak4_flat_jets_[base_idx + 6] = static_cast<float>(event->ak4_isMediumBTag[i]);
            ak4_flat_jets_[base_idx + 7] = static_cast<float>(event->ak4_isLooseBTag[i]);
        }

        const size_t max_ak8 = std::min<size_t>(MAX_AK8_JETS, event->ak8_pt.size());
        for (size_t i = 0; i < max_ak8; ++i) {
            ak8_mask_char_[batch_idx * MAX_AK8_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK8_JETS * AK8_FEATURES + i * AK8_FEATURES;
            ak8_flat_jets_[base_idx]     = std::log(event->ak8_mass[i] + 1.0f);
            ak8_flat_jets_[base_idx + 1] = std::log(event->ak8_pt[i] + 1.0f);
            ak8_flat_jets_[base_idx + 2] = event->ak8_eta[i];
            ak8_flat_jets_[base_idx + 3] = std::sin(event->ak8_phi[i]);
            ak8_flat_jets_[base_idx + 4] = std::cos(event->ak8_phi[i]);
            ak8_flat_jets_[base_idx + 5] = static_cast<float>(event->ak8_nConstituents[i]);
            ak8_flat_jets_[base_idx + 6] = event->ak8_XbbScore[i];
            ak8_flat_jets_[base_idx + 7] = event->ak8_XqqScore[i];
            ak8_flat_jets_[base_idx + 8] = event->ak8_XccScore[i];
            ak8_flat_jets_[base_idx + 9] = event->ak8_XcsScore[i];
            ak8_flat_jets_[base_idx + 10] = event->ak8_XqcdScore[i];
        }

        const size_t event_base_idx = batch_idx * MET_FEATURES;
        met_inputs_[event_base_idx]     = std::log(event->met_pt + 1.0f);
        met_inputs_[event_base_idx + 1] = std::sin(event->met_phi);
        met_inputs_[event_base_idx + 2] = std::cos(event->met_phi);

        met_mask_char_[batch_idx] = 1;

        const size_t lepton1_base_idx = batch_idx * LEPTON_FEATURES;
        lepton1_inputs_[lepton1_base_idx]     = std::log(event->lepton1_pt + 1.0f);
        lepton1_inputs_[lepton1_base_idx + 1] = event->lepton1_eta;
        lepton1_inputs_[lepton1_base_idx + 2] = std::sin(event->lepton1_phi);
        lepton1_inputs_[lepton1_base_idx + 3] = std::cos(event->lepton1_phi);
        lepton1_inputs_[lepton1_base_idx + 4] = std::log(event->lepton1_mass + 1.0f);
        lepton1_inputs_[lepton1_base_idx + 5] = static_cast<float>(event->lepton1_charge);

        lepton1_mask_char_[batch_idx] = 1;

        const size_t lepton2_base_idx = batch_idx * LEPTON_FEATURES;
        lepton2_inputs_[lepton2_base_idx]     = std::log(event->lepton2_pt + 1.0f);
        lepton2_inputs_[lepton2_base_idx + 1] = event->lepton2_eta;
        lepton2_inputs_[lepton2_base_idx + 2] = std::sin(event->lepton2_phi);
        lepton2_inputs_[lepton2_base_idx + 3] = std::cos(event->lepton2_phi);
        lepton2_inputs_[lepton2_base_idx + 4] = std::log(event->lepton2_mass + 1.0f);
        lepton2_inputs_[lepton2_base_idx + 5] = static_cast<float>(event->lepton2_charge);

        lepton2_mask_char_[batch_idx] = 1;
    }
}

RNode SPANetInferenceRun3::ParseSpanetInference(RNode df_) {
    auto df = df_.Define("_all_assignments", SPANetInferenceBase::assign_all_objects, {
        "spanet_vbs_assignment", "spanet_h_assignment", "spanet_bh_assignment",
        "spanet_v1_assignment", "spanet_v2_assignment", "spanet_bv1_assignment", "spanet_bv2_assignment",
        "spanet_vbs_detection", "spanet_h_detection", "spanet_bh_detection",
        "spanet_v1_detection", "spanet_v2_detection", "spanet_bv1_detection", "spanet_bv2_detection",
        "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"
    })
    .Define("vbs1_idx", "_all_assignments[0]")
    .Define("vbs2_idx", "_all_assignments[1]")
    .Define("h1_idx", "_all_assignments[2]")
    .Define("h2_idx", "_all_assignments[3]")
    .Define("bh_idx", "_all_assignments[4]")
    .Define("v1_j1_idx", "_all_assignments[5]")
    .Define("v1_j2_idx", "_all_assignments[6]")
    .Define("v2_j1_idx", "_all_assignments[7]")
    .Define("v2_j2_idx", "_all_assignments[8]")
    .Define("bv1_idx", "_all_assignments[9]")
    .Define("bv2_idx", "_all_assignments[10]");

    // VBS jet variables
    df = df.Define("spanet_vbs1_pt", "vbs1_idx >= 0 ? Jet_pt[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_eta", "vbs1_idx >= 0 ? Jet_eta[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_phi", "vbs1_idx >= 0 ? Jet_phi[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_mass", "vbs1_idx >= 0 ? Jet_mass[vbs1_idx] : -999.0f")
           .Define("spanet_vbs2_pt", "vbs2_idx >= 0 ? Jet_pt[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_eta", "vbs2_idx >= 0 ? Jet_eta[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_phi", "vbs2_idx >= 0 ? Jet_phi[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_mass", "vbs2_idx >= 0 ? Jet_mass[vbs2_idx] : -999.0f")
           .Define("spanet_vbs_detajj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? abs(spanet_vbs1_eta - spanet_vbs2_eta) : -999.0f")
           .Define("spanet_vbs_mjj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_vbs1_pt, spanet_vbs1_eta, spanet_vbs1_phi, spanet_vbs1_mass) + "
                                     "ROOT::Math::PtEtaPhiMVector(spanet_vbs2_pt, spanet_vbs2_eta, spanet_vbs2_phi, spanet_vbs2_mass)).M() : -999.0f");

    // Resolved Higgs variables
    df = df.Define("spanet_h1_pt", "h1_idx >= 0 ? Jet_pt[h1_idx] : -999.0f")
           .Define("spanet_h1_eta", "h1_idx >= 0 ? Jet_eta[h1_idx] : -999.0f")
           .Define("spanet_h1_phi", "h1_idx >= 0 ? Jet_phi[h1_idx] : -999.0f")
           .Define("spanet_h1_mass", "h1_idx >= 0 ? Jet_mass[h1_idx] : -999.0f")
           .Define("spanet_h2_pt", "h2_idx >= 0 ? Jet_pt[h2_idx] : -999.0f")
           .Define("spanet_h2_eta", "h2_idx >= 0 ? Jet_eta[h2_idx] : -999.0f")
           .Define("spanet_h2_phi", "h2_idx >= 0 ? Jet_phi[h2_idx] : -999.0f")
           .Define("spanet_h2_mass", "h2_idx >= 0 ? Jet_mass[h2_idx] : -999.0f")
           .Define("spanet_h_mjj", "h1_idx >= 0 && h2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_h1_pt, spanet_h1_eta, spanet_h1_phi, spanet_h1_mass) + "
                                   "ROOT::Math::PtEtaPhiMVector(spanet_h2_pt, spanet_h2_eta, spanet_h2_phi, spanet_h2_mass)).M() : -999.0f");

    // Boosted Higgs variables
    df = df.Define("spanet_bh_eta", "bh_idx >= 0 ? FatJet_eta[bh_idx] : -999.0f")
           .Define("spanet_bh_phi", "bh_idx >= 0 ? FatJet_phi[bh_idx] : -999.0f")
           .Define("spanet_bh_msoftdrop", "bh_idx >= 0 ? FatJet_msoftdrop[bh_idx] : -999.0f")
           .Define("spanet_bh_pt", "bh_idx >= 0 ? FatJet_pt[bh_idx] : -999.0f")
           .Define("spanet_bh_score", "bh_idx >= 0 ? FatJet_globalParT3_Xbb[bh_idx] / (FatJet_globalParT3_Xbb[bh_idx] + FatJet_globalParT3_QCD[bh_idx]) : -999.0f");

    // Resolved V1 variables
    df = df.Define("spanet_v1_j1_pt", "v1_j1_idx >= 0 ? Jet_pt[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_eta", "v1_j1_idx >= 0 ? Jet_eta[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_phi", "v1_j1_idx >= 0 ? Jet_phi[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_mass", "v1_j1_idx >= 0 ? Jet_mass[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j2_pt", "v1_j2_idx >= 0 ? Jet_pt[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_eta", "v1_j2_idx >= 0 ? Jet_eta[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_phi", "v1_j2_idx >= 0 ? Jet_phi[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_mass", "v1_j2_idx >= 0 ? Jet_mass[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_mjj", "v1_j1_idx >= 0 && v1_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v1_j1_pt, spanet_v1_j1_eta, spanet_v1_j1_phi, spanet_v1_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v1_j2_pt, spanet_v1_j2_eta, spanet_v1_j2_phi, spanet_v1_j2_mass)).M() : -999.0f");

    // Resolved V2 variables
    df = df.Define("spanet_v2_j1_pt", "v2_j1_idx >= 0 ? Jet_pt[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_eta", "v2_j1_idx >= 0 ? Jet_eta[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_phi", "v2_j1_idx >= 0 ? Jet_phi[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_mass", "v2_j1_idx >= 0 ? Jet_mass[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j2_pt", "v2_j2_idx >= 0 ? Jet_pt[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_eta", "v2_j2_idx >= 0 ? Jet_eta[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_phi", "v2_j2_idx >= 0 ? Jet_phi[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_mass", "v2_j2_idx >= 0 ? Jet_mass[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_mjj", "v2_j1_idx >= 0 && v2_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v2_j1_pt, spanet_v2_j1_eta, spanet_v2_j1_phi, spanet_v2_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v2_j2_pt, spanet_v2_j2_eta, spanet_v2_j2_phi, spanet_v2_j2_mass)).M() : -999.0f");

    // Boosted V1 variables
    df = df.Define("spanet_bv1_eta", "bv1_idx >= 0 ? FatJet_eta[bv1_idx] : -999.0f")
           .Define("spanet_bv1_phi", "bv1_idx >= 0 ? FatJet_phi[bv1_idx] : -999.0f")
           .Define("spanet_bv1_msoftdrop", "bv1_idx >= 0 ? FatJet_msoftdrop[bv1_idx] : -999.0f")
           .Define("spanet_bv1_pt", "bv1_idx >= 0 ? FatJet_pt[bv1_idx] : -999.0f")
           .Define("spanet_bv1_score", "bv1_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv1_idx] : -999.0f")
           .Define("spanet_bv1_w_score", "bv1_idx >= 0 ? FatJet_globalParT3_Xqq[bv1_idx] / (FatJet_globalParT3_Xqq[bv1_idx] + FatJet_globalParT3_Xcs[bv1_idx] + FatJet_globalParT3_QCD[bv1_idx]) : -999.0f");

    // Boosted V2 variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? FatJet_eta[bv2_idx] : -999.0f")
           .Define("spanet_bv2_phi", "bv2_idx >= 0 ? FatJet_phi[bv2_idx] : -999.0f")
           .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? FatJet_msoftdrop[bv2_idx] : -999.0f")
           .Define("spanet_bv2_pt", "bv2_idx >= 0 ? FatJet_pt[bv2_idx] : -999.0f")
           .Define("spanet_bv2_score", "bv2_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv2_idx] : -999.0f")
           .Define("spanet_bv2_w_score", "bv2_idx >= 0 ? FatJet_globalParT3_Xqq[bv2_idx] / (FatJet_globalParT3_Xqq[bv2_idx] + FatJet_globalParT3_Xcs[bv2_idx] + FatJet_globalParT3_QCD[bv2_idx]) : -999.0f");

    return df;
}

}