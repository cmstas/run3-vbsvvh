#include "spanet_run2.h"

std::vector<SPANetRun2::SPANetInference::EventData> SPANetRun2::SPANetInference::extractEventsFromDataFrame(RNode df) {
    std::cout << " -> SPANetInferenceRun2::extractEventsFromDataFrame()" << std::endl;

    auto ak4_pt_vec = df.Take<RVec<float>>("jet_pt").GetValue();
    auto ak4_eta_vec = df.Take<RVec<float>>("jet_eta").GetValue();
    auto ak4_phi_vec = df.Take<RVec<float>>("jet_phi").GetValue();
    auto ak4_mass_vec = df.Take<RVec<float>>("jet_mass").GetValue();
    auto ak4_isTightBTag_vec = df.Take<RVec<bool>>("jet_isTightBTag").GetValue();
    auto ak4_isMediumBTag_vec = df.Take<RVec<bool>>("jet_isMediumBTag").GetValue();
    auto ak4_isLooseBTag_vec = df.Take<RVec<bool>>("jet_isLooseBTag").GetValue();

    auto ak8_pt_vec = df.Take<RVec<float>>("fatjet_pt").GetValue();
    auto ak8_eta_vec = df.Take<RVec<float>>("fatjet_eta").GetValue();
    auto ak8_phi_vec = df.Take<RVec<float>>("fatjet_phi").GetValue();
    auto ak8_mass_vec = df.Take<RVec<float>>("fatjet_msoftdrop").GetValue();

    auto ak8_nConstituents_vec = df.Take<RVec<unsigned char>>("fatjet_nConstituents").GetValue();
    auto ak8_HbbScore_vec = df.Take<RVec<float>>("fatjet_HbbScore").GetValue();
    auto ak8_WqqScore_vec = df.Take<RVec<float>>("fatjet_WqqScore").GetValue();

    auto met_pt_vec = df.Take<float>("MET_pt").GetValue();
    auto ht_ak4_vec = df.Take<float>("ht_jets").GetValue();
    auto ht_ak8_vec = df.Take<float>("ht_fatjets").GetValue();

    auto rdf_entry_vec = df.Take<ULong64_t>("rdfentry_").GetValue();
    auto rdf_slot_vec = df.Take<unsigned int>("rdfslot_").GetValue();

    std::vector<EventData> events;
    size_t n_events = ak4_pt_vec.size();
    events.reserve(n_events);

    for (size_t i = 0; i < n_events; ++i) {
        EventData event;

        event.ak4_pt.assign(ak4_pt_vec[i].begin(), ak4_pt_vec[i].end());
        event.ak4_eta.assign(ak4_eta_vec[i].begin(), ak4_eta_vec[i].end());
        event.ak4_phi.assign(ak4_phi_vec[i].begin(), ak4_phi_vec[i].end());
        event.ak4_mass.assign(ak4_mass_vec[i].begin(), ak4_mass_vec[i].end());
        event.ak4_isTightBTag.assign(ak4_isTightBTag_vec[i].begin(), ak4_isTightBTag_vec[i].end());
        event.ak4_isMediumBTag.assign(ak4_isMediumBTag_vec[i].begin(), ak4_isMediumBTag_vec[i].end());
        event.ak4_isLooseBTag.assign(ak4_isLooseBTag_vec[i].begin(), ak4_isLooseBTag_vec[i].end());

        event.ak8_pt.assign(ak8_pt_vec[i].begin(), ak8_pt_vec[i].end());
        event.ak8_eta.assign(ak8_eta_vec[i].begin(), ak8_eta_vec[i].end());
        event.ak8_phi.assign(ak8_phi_vec[i].begin(), ak8_phi_vec[i].end());
        event.ak8_mass.assign(ak8_mass_vec[i].begin(), ak8_mass_vec[i].end());
        event.ak8_nConstituents.assign(ak8_nConstituents_vec[i].begin(), ak8_nConstituents_vec[i].end());
        event.ak8_HbbScore.assign(ak8_HbbScore_vec[i].begin(), ak8_HbbScore_vec[i].end());
        event.ak8_WqqScore.assign(ak8_WqqScore_vec[i].begin(), ak8_WqqScore_vec[i].end());

        event.met_pt = met_pt_vec[i];
        event.ht_ak4 = ht_ak4_vec[i];
        event.ht_ak8 = ht_ak8_vec[i];

        event.rdf_entry = rdf_entry_vec[i];
        event.rdf_slot = rdf_slot_vec[i];

        events.push_back(event);
    }

    return events;
}

std::vector<std::vector<std::vector<std::vector<float>>>> SPANetRun2::SPANetInference::runBatchInference(const std::vector<EventData>& events) {
    std::vector<std::vector<std::vector<std::vector<float>>>> all_results;
    all_results.reserve(events.size());
    
    auto extract_idxs_from_inference = [this](const Ort::Value& output_tensor, size_t batch_idx, size_t tensor_idx, 
                                         std::vector<std::vector<std::vector<std::vector<float>>>>& batch_results) {
        const Ort::TensorTypeAndShapeInfo& shape_info = output_tensor.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> output_shape = shape_info.GetShape();
        const float* output_data = output_tensor.GetTensorData<float>();
        

        if (tensor_idx < 7) { 
            static thread_local std::vector<std::pair<float, size_t>> value_idx_pairs;
            
            if (output_shape.size() == 3) { // [batch, jets, jets] - always jet pair assignment
                size_t jets_dim1 = output_shape[1];
                size_t jets_dim2 = output_shape[2];
                size_t offset = batch_idx * jets_dim1 * jets_dim2;
                
                size_t tri_size = (jets_dim1 * (jets_dim1 - 1)) / 2;
                value_idx_pairs.clear();
                value_idx_pairs.reserve(tri_size);
                
                for (size_t i = 0; i < jets_dim1; ++i) {
                    size_t row_offset = offset + i * jets_dim2;
                    for (size_t j = i + 1; j < jets_dim2; ++j) { 
                        value_idx_pairs.emplace_back(output_data[row_offset + j], i * jets_dim2 + j);
                    }
                }
                
                const size_t _top_k = std::min(static_cast<size_t>(this->top_k), value_idx_pairs.size());
                std::partial_sort(value_idx_pairs.begin(), value_idx_pairs.begin() + _top_k, value_idx_pairs.end(),
                                [](const auto& a, const auto& b) { return a.first > b.first; });
                
                auto& result_vec = batch_results[batch_idx][tensor_idx];
                result_vec.clear();
                result_vec.reserve(_top_k);
                
                for (size_t k = 0; k < _top_k; ++k) {
                    float value = value_idx_pairs[k].first;
                    size_t linear_idx = value_idx_pairs[k].second;
                    int jet1_idx = static_cast<int>(linear_idx / jets_dim2);
                    int jet2_idx = static_cast<int>(linear_idx % jets_dim2);
                    result_vec.push_back({value, static_cast<float>(jet1_idx), static_cast<float>(jet2_idx)});
                }
            } // [batch, jets, jets]

            else if (output_shape.size() == 2) {
                size_t jets_dim = output_shape[1];
                size_t offset = batch_idx * jets_dim;
                
                value_idx_pairs.clear();
                value_idx_pairs.reserve(jets_dim);
                
                for (size_t j = MAX_AK4_JETS; j < jets_dim; ++j) { // start after AK4 jets
                    value_idx_pairs.emplace_back(output_data[offset + j], j - MAX_AK4_JETS); // store value and index relative to AK8 jets
                }
                
                const size_t _top_k = std::min(static_cast<size_t>(this->top_k), value_idx_pairs.size());
                std::partial_sort(value_idx_pairs.begin(), value_idx_pairs.begin() + _top_k, value_idx_pairs.end(),
                                [](const auto& a, const auto& b) { return a.first > b.first; });
                
                auto& result_vec = batch_results[batch_idx][tensor_idx];
                result_vec.clear();
                result_vec.reserve(_top_k);
                
                for (size_t k = 0; k < _top_k; ++k) {
                    result_vec.push_back({value_idx_pairs[k].first, static_cast<float>(value_idx_pairs[k].second)});
                }
            } // [batch, jets]
        } 
        
        else if (tensor_idx < 14) { 
            batch_results[batch_idx][tensor_idx] = {{output_data[batch_idx]}};
        } 
            
        else {
            size_t offset = batch_idx * 2;
            batch_results[batch_idx][tensor_idx] = {{output_data[offset + 1]}}; 
        }
    };

    for (size_t start_idx = 0; start_idx < events.size(); start_idx += batch_size) {
        size_t end_idx = std::min(start_idx + batch_size, events.size());
        size_t current_batch_size = end_idx - start_idx;
        
        std::vector<EventData> batch_events(events.begin() + start_idx, events.begin() + end_idx);
        
        std::vector<int64_t> current_ak4_shape = {static_cast<int64_t>(current_batch_size), MAX_AK4_JETS, AK4_FEATURES};
        std::vector<int64_t> current_ak8_shape = {static_cast<int64_t>(current_batch_size), MAX_AK8_JETS, AK8_FEATURES};
        std::vector<int64_t> current_met_shape = {static_cast<int64_t>(current_batch_size), 1, MET_FEATURES};
        std::vector<int64_t> current_ak4_mask_shape = {static_cast<int64_t>(current_batch_size), MAX_AK4_JETS};
        std::vector<int64_t> current_ak8_mask_shape = {static_cast<int64_t>(current_batch_size), MAX_AK8_JETS};
        std::vector<int64_t> current_met_mask_shape = {static_cast<int64_t>(current_batch_size), 1};
        std::vector<int64_t> current_lepton1_shape = {static_cast<int64_t>(current_batch_size), 1, LEPTON_FEATURES};
        std::vector<int64_t> current_lepton2_shape = {static_cast<int64_t>(current_batch_size), 1, LEPTON_FEATURES};
        std::vector<int64_t> current_lepton1_mask_shape = {static_cast<int64_t>(current_batch_size), 1};
        std::vector<int64_t> current_lepton2_mask_shape = {static_cast<int64_t>(current_batch_size), 1};

        fillBatchTensors(batch_events, current_batch_size);
        
        std::vector<Ort::Value> input_tensors;
        input_tensors.reserve(10);
        
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak4_flat_jets.data(), current_batch_size * MAX_AK4_JETS * AK4_FEATURES,
            current_ak4_shape.data(), current_ak4_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak4_mask_char.data()), current_batch_size * MAX_AK4_JETS,
            current_ak4_mask_shape.data(), current_ak4_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak8_flat_jets.data(), current_batch_size * MAX_AK8_JETS * AK8_FEATURES,
            current_ak8_shape.data(), current_ak8_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak8_mask_char.data()), current_batch_size * MAX_AK8_JETS,
            current_ak8_mask_shape.data(), current_ak8_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, met_inputs.data(), current_batch_size * MET_FEATURES,
            current_met_shape.data(), current_met_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(met_mask_char.data()), current_batch_size * 1,
            current_met_mask_shape.data(), current_met_mask_shape.size()));


        if (lepton1_inputs.size()>0) {
            input_tensors.push_back(Ort::Value::CreateTensor<float>(
                memory_info, lepton1_inputs.data(), current_batch_size * lepton1_input_shape[2], current_lepton1_shape.data(), current_lepton1_shape.size()));
            input_tensors.push_back(Ort::Value::CreateTensor<bool>(
                memory_info, reinterpret_cast<bool*>(lepton1_mask_char.data()), current_batch_size * 1, current_lepton1_mask_shape.data(), current_lepton1_mask_shape.size()));
        }
        if (lepton2_inputs.size()>0) {
            input_tensors.push_back(Ort::Value::CreateTensor<float>(
                memory_info, lepton2_inputs.data(), current_batch_size * lepton2_input_shape[2], current_lepton2_shape.data(), current_lepton2_shape.size()));
            input_tensors.push_back(Ort::Value::CreateTensor<bool>(
                memory_info, reinterpret_cast<bool*>(lepton2_mask_char.data()), current_batch_size * 1, current_lepton2_mask_shape.data(), current_lepton2_mask_shape.size()));
        }

        auto output_tensors = session->Run(
            Ort::RunOptions{nullptr},
            input_names.data(), input_tensors.data(), input_tensors.size(),
            output_names.data(), output_names.size());
        
        std::vector<std::vector<std::vector<std::vector<float>>>> batch_results(current_batch_size);

        for (size_t tensor_idx = 0; tensor_idx < output_tensors.size(); ++tensor_idx) {
            Ort::Value& output_tensor = output_tensors[tensor_idx];
            for (size_t batch_idx = 0; batch_idx < current_batch_size; ++batch_idx) {
                if (batch_results[batch_idx].empty()) {
                    batch_results[batch_idx].resize(output_tensors.size());
                }
                extract_idxs_from_inference(output_tensor, batch_idx, tensor_idx, batch_results);
            }
        }
        
        for (size_t batch_idx = 0; batch_idx < current_batch_size; ++batch_idx) {
            all_results.push_back(std::move(batch_results[batch_idx]));
        }
    }
    
    return all_results;
}


void SPANetRun2::SPANetInference::fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size){
    //std::cout << "  -> SPANetInference::fillBatchTensors() " << std::endl;

    std::fill(ak4_flat_jets.begin(), ak4_flat_jets.begin() + actual_batch_size * MAX_AK4_JETS * AK4_FEATURES, 0.0f);
    std::fill(ak8_flat_jets.begin(), ak8_flat_jets.begin() + actual_batch_size * MAX_AK8_JETS * AK8_FEATURES, 0.0f);
    std::fill(ak4_mask_char.begin(), ak4_mask_char.begin() + actual_batch_size * MAX_AK4_JETS, 0);
    std::fill(ak8_mask_char.begin(), ak8_mask_char.begin() + actual_batch_size * MAX_AK8_JETS, 0);

    for (size_t batch_idx = 0; batch_idx < actual_batch_size; ++batch_idx) {
        const EventData& event = events[batch_idx];

        const size_t max_ak4 = std::min<size_t>(MAX_AK4_JETS, event.ak4_pt.size());
        for (size_t i = 0; i < max_ak4; ++i) {
            ak4_mask_char[batch_idx * MAX_AK4_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK4_JETS * AK4_FEATURES + i * AK4_FEATURES;
            ak4_flat_jets[base_idx]     = std::log(event.ak4_mass[i] + 1.0f);
            ak4_flat_jets[base_idx + 1] = std::log(event.ak4_pt[i] + 1.0f);
            ak4_flat_jets[base_idx + 2] = event.ak4_eta[i];
            ak4_flat_jets[base_idx + 3] = std::sin(event.ak4_phi[i]);
            ak4_flat_jets[base_idx + 4] = std::cos(event.ak4_phi[i]);
            ak4_flat_jets[base_idx + 5] = static_cast<float>(event.ak4_isTightBTag[i]);
            ak4_flat_jets[base_idx + 6] = static_cast<float>(event.ak4_isMediumBTag[i]);
            ak4_flat_jets[base_idx + 7] = static_cast<float>(event.ak4_isLooseBTag[i]);
        }

        const size_t max_ak8 = std::min<size_t>(MAX_AK8_JETS, event.ak8_pt.size());
        for (size_t i = 0; i < max_ak8; ++i) {
            ak8_mask_char[batch_idx * MAX_AK8_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK8_JETS * AK8_FEATURES + i * AK8_FEATURES;
            ak8_flat_jets[base_idx]     = std::log(event.ak8_mass[i] + 1.0f);
            ak8_flat_jets[base_idx + 1] = std::log(event.ak8_pt[i] + 1.0f);
            ak8_flat_jets[base_idx + 2] = event.ak8_eta[i];
            ak8_flat_jets[base_idx + 3] = std::sin(event.ak8_phi[i]);
            ak8_flat_jets[base_idx + 4] = std::cos(event.ak8_phi[i]);
            ak8_flat_jets[base_idx + 5] = static_cast<float>(event.ak8_nConstituents[i]);
            ak8_flat_jets[base_idx + 6] = event.ak8_HbbScore[i];
            ak8_flat_jets[base_idx + 7] = event.ak8_WqqScore[i];
        }

        const size_t event_base_idx = batch_idx * MET_FEATURES;
        met_inputs[event_base_idx]     = std::log(event.met_pt + 1.0f);
        met_inputs[event_base_idx + 1] = std::log(event.ht_ak4 + 1.0f);
        met_inputs[event_base_idx + 2] = std::log(event.ht_ak8 + 1.0f);

        met_mask_char[batch_idx] = 1;
    }
    
}

RNode SPANetRun2::SPANetInference::addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs, const std::vector<EventData>& events) {
    // Create maps for fast lookup using entry+slot as key
    std::map<std::pair<ULong64_t, unsigned int>, size_t> entry_slot_to_output_idx;
    
    for (size_t i = 0; i < events.size(); ++i) {
        auto key = std::make_pair(events[i].rdf_entry, events[i].rdf_slot);
        entry_slot_to_output_idx[key] = i;
    }
    
    // Create column-wise accessors that use entry+slot matching
    auto getHAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][0];
        }
        return {};
    };
    
    auto getBHAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][1];
        }
        return {};
    };
    
    auto getV1Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][2];
        }
        return {};
    };
    
    auto getV2Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][3];
        }
        return {};
    };
    
    auto getBV1Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][4];
        }
        return {};
    };
    
    auto getBV2Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][5];
        }
        return {};
    };
    
    auto getVBSAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>> {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][6];
        }
        return {};
    };
    
    auto getHDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][7][0][0];
        }
        return -999.0f;
    };
    
    auto getBHDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][8][0][0];
        }
        return -999.0f;
    };
    
    auto getV1Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][9][0][0];
        }
        return -999.0f;
    };
    
    auto getV2Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][10][0][0];
        }
        return -999.0f;
    };
    
    auto getBV1Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][11][0][0];
        }
        return -999.0f;
    };
    
    auto getBV2Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][12][0][0];
        }
        return -999.0f;
    };
    
    auto getVBSDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][13][0][0];
        }
        return -999.0f;
    };
    
    auto getEventSignal = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end()) {
            return all_outputs[it->second][14][0][0];
        }
        return -999.0f;
    };
    
    return df.DefineSlotEntry("spanet_h_assignment", getHAssignment, {})
             .DefineSlotEntry("spanet_bh_assignment", getBHAssignment, {})
             .DefineSlotEntry("spanet_v1_assignment", getV1Assignment, {})
             .DefineSlotEntry("spanet_v2_assignment", getV2Assignment, {})
             .DefineSlotEntry("spanet_bv1_assignment", getBV1Assignment, {})
             .DefineSlotEntry("spanet_bv2_assignment", getBV2Assignment, {})
             .DefineSlotEntry("spanet_vbs_assignment", getVBSAssignment, {})
             .DefineSlotEntry("spanet_h_detection", getHDetection, {})
             .DefineSlotEntry("spanet_bh_detection", getBHDetection, {})
             .DefineSlotEntry("spanet_v1_detection", getV1Detection, {})
             .DefineSlotEntry("spanet_v2_detection", getV2Detection, {})
             .DefineSlotEntry("spanet_bv1_detection", getBV1Detection, {})
             .DefineSlotEntry("spanet_bv2_detection", getBV2Detection, {})
             .DefineSlotEntry("spanet_vbs_detection", getVBSDetection, {})
             .DefineSlotEntry("spanet_event_signal", getEventSignal, {});
}

RNode SPANetRun2::SPANetInference::RunSPANetInference(RNode df) {
    auto start_time = std::chrono::high_resolution_clock::now();
    auto events = extractEventsFromDataFrame(df);
    auto extract_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "Running batch inference on " << events.size() << " events..." << std::endl;
    auto all_outputs = runBatchInference(events);
    auto inference_time = std::chrono::high_resolution_clock::now();
    
    auto result = addSPANetOutputsToDataFrame(df, all_outputs, events);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    auto extract_ms = std::chrono::duration_cast<std::chrono::milliseconds>(extract_time - start_time).count();
    auto inference_ms = std::chrono::duration_cast<std::chrono::milliseconds>(inference_time - extract_time).count();
    auto output_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - inference_time).count();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    
    std::cout << "SPANet processing times: " 
              << "Data extraction: " << extract_ms << "ms, "
              << "Inference: " << inference_ms << "ms, "
              << "Output creation: " << output_ms << "ms, "
              << "Total: " << total_ms << "ms" << std::endl;
    
    return result;
}

std::vector<int> SPANetRun2::SPANetInference::assign_all_objects(
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


std::vector<int> SPANetRun2::SPANetInference::assign_all_objects_maxprob(
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
    bool DEBUG = false;

    if (DEBUG) {
        std::cout << "\n\n EVENT STARTS HERE" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "\n N fat jets = " << FatJet_phi.size() << ", N jets = " << Jet_phi.size() << std::endl;

        std::cout << "\nStarting Boosted Decays Selection" << std::endl;
        std::cout << "====================" << std::endl;
        std::cout << " Boosted assignments:" << std::endl;
        std::cout << " bh = {";
        for (const auto& v : bh_assignment) {
            if (!v.empty() && v.size() > 1) std::cout << "[" << v[0] << ", " << static_cast<int>(v[1]) << "], ";
        }
        std::cout << "}" << std::endl;
        std::cout << " bv1 = {";
        for (const auto& v : bv1_assignment) {
            if (!v.empty() && v.size() > 1) std::cout << "[" << v[0] << ", " << static_cast<int>(v[1]) << "], ";
        }
        std::cout << "}" << std::endl;
        std::cout << " bv2 = {";
        for (const auto& v : bv2_assignment) {
            if (!v.empty() && v.size() > 1) std::cout << "[" << v[0] << ", " << static_cast<int>(v[1]) << "], ";
        }
        std::cout << "}" << std::endl;
    }

    std::vector<int> result(11, -1);
    std::set<int> used_jets;
    std::set<int> used_fatjets;

    // Use inputs directly as candidates (already flattened and sorted descending by prob)
    auto cand_vbs = vbs_assignment;
    auto cand_h = h_assignment;
    auto cand_v1 = v1_assignment;
    auto cand_v2 = v2_assignment;
    auto cand_bh = bh_assignment;
    auto cand_bv1 = bv1_assignment;
    auto cand_bv2 = bv2_assignment;

    // Filter for boosted
    auto filter_candidates = [](const std::vector<std::vector<float>>& cand, float th) {
        std::vector<std::vector<float>> filtered;
        for (const auto& c : cand) {
            if (!c.empty() && c.size() > 1 && c[0] > th) {
                filtered.push_back(c);
            }
        }
        return filtered;
    };

    float threshold = 0.5f;
    std::vector<std::vector<std::vector<float>>> boosted_cands = {
        filter_candidates(cand_bh, threshold),
        filter_candidates(cand_bv1, threshold),
        filter_candidates(cand_bv2, threshold)
    };

    // Helper to check jet idx
    auto passJetIdxCheck = [&](int jet_idx) -> bool {
        if (jet_idx < 0 || jet_idx >= Jet_eta.size()) return false;
        else if (used_jets.count(jet_idx) > 0) return false;
        return true;
    };

    // Helper to check fatjet idx
    auto passFatJetIdxCheck = [&](int fatjet_idx) -> bool {
        if (fatjet_idx < 0 || fatjet_idx >= FatJet_eta.size()) return false;
        else if (used_fatjets.count(fatjet_idx) > 0) return false;
        return true;
    };

    // Helper to check DeltaR
    auto passDeltaRCheck = [&](int idx1, int idx2, bool is_fatjet1, bool is_fatjet2) -> bool {
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

    // Helper to check all overlaps
    auto checkAllOverlaps = [&](int candidate_idx, bool is_fatjet) -> bool {
        for (int assigned_jet : used_jets) {
            if (!passDeltaRCheck(candidate_idx, assigned_jet, is_fatjet, false)) {
                if (DEBUG) { std::cout << "Failed deltaR(cand, assigned_jet) check with assigned jet idx " << assigned_jet << std::endl; }
                return false;
            }
        }
        for (int assigned_fatjet : used_fatjets) {
            if (!passDeltaRCheck(candidate_idx, assigned_fatjet, is_fatjet, true)) {
                if (DEBUG) {std::cout << "Failed deltaR(cand, assigned_fatjet) check with assigned fatjet idx " << assigned_fatjet << std::endl;}
                return false;
            }
        }
        return true;
    };

    // Helper to find best fatjet
    auto get_best_fatjet = [&](const std::vector<std::vector<float>>& candidates, const std::set<int>& used) -> std::pair<int, float> {
        for (const auto& cand : candidates) {
            int fj = static_cast<int>(cand[1]);
            if (fj >= 0 && fj < FatJet_eta.size() && passFatJetIdxCheck(fj) && checkAllOverlaps(fj, true)) {
                return {fj, cand[0]};
            }
        }
        return {-1, -std::numeric_limits<float>::infinity()};
    };

    // Helper to find best pair
    auto get_best_pair = [&](const std::vector<std::vector<float>>& candidates, const std::set<int>& used) -> std::tuple<int, int, float> {
        for (const auto& cand : candidates) {
            int j1 = static_cast<int>(cand[1]);
            int j2 = static_cast<int>(cand[2]);
            if (j1 != j2 && j1 >= 0 && j1 < Jet_eta.size() && j2 >= 0 && j2 < Jet_eta.size() && 
                passJetIdxCheck(j1) && passJetIdxCheck(j2) && 
                passDeltaRCheck(j1, j2, false, false) && 
                checkAllOverlaps(j1, false) && checkAllOverlaps(j2, false)) {
                return {j1, j2, cand[0]};
            }
        }
        return {-1, -1, -std::numeric_limits<float>::infinity()};
    };

    // Boosted assignment logic
    std::vector<int> boosted_result_pos = {4, 9, 10}; // bh, bv1, bv2

    std::vector<float> top_probs(3, -std::numeric_limits<float>::infinity());
    for (int i = 0; i < 3; ++i) {
        if (!boosted_cands[i].empty()) {
            top_probs[i] = boosted_cands[i][0][0];
        }
    }

    int first_idx = -1;
    float max_prob = -std::numeric_limits<float>::infinity();
    for (int i = 0; i < 3; ++i) {
        if (top_probs[i] > max_prob) {
            max_prob = top_probs[i];
            first_idx = i;
        }
    }

    if (DEBUG && boosted_cands[0].empty() && boosted_cands[1].empty() && boosted_cands[2].empty()) {
        std::cout << "No boosted boson found with any p > 0.5" << std::endl;
    }

    if (first_idx != -1) {
        auto [fj, prob] = get_best_fatjet(boosted_cands[first_idx], used_fatjets);
        if (fj != -1) {
            if (DEBUG) {
                std::cout << "First matched boosted boson: jet idx = " << fj << ", original idx = " << first_idx << std::endl;
            }
            result[boosted_result_pos[first_idx]] = fj;
            used_fatjets.insert(fj);
        }

        std::vector<int> remaining;
        for (int i = 0; i < 3; ++i) {
            if (i != first_idx && !boosted_cands[i].empty()) {
                remaining.push_back(i);
            }
        }

        if (!remaining.empty()) {
            if (remaining.size() == 1) {
                int idx = remaining[0];
                auto [fj, prob] = get_best_fatjet(boosted_cands[idx], used_fatjets);
                if (fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found second matched boosted boson: idx = " << fj << ", original idx = " << idx << std::endl;
                    }
                    result[boosted_result_pos[idx]] = fj;
                    used_fatjets.insert(fj);
                } else if (DEBUG) {
                    std::cout << "No boson left after overlap removal." << std::endl;
                }
            } else {
                std::vector<std::tuple<int, float, int>> best_probs;
                for (int idx : remaining) {
                    auto [fj, prob] = get_best_fatjet(boosted_cands[idx], used_fatjets);
                    best_probs.push_back({idx, prob, fj});
                }

                if (DEBUG) {
                    std::cout << "Potential second boosted bosons (boson idx, prob, fj): {";
                    for (const auto& [idx, prob, fj] : best_probs) {
                        std::cout << "(" << idx << ", " << prob << ", " << fj << "), ";
                    }
                    std::cout << "}" << std::endl;
                }

                // Sort descending by prob
                std::sort(best_probs.begin(), best_probs.end(),
                          [](const std::tuple<int, float, int>& a, const std::tuple<int, float, int>& b) {
                              return std::get<1>(a) > std::get<1>(b);
                          });

                // Assign the best one
                int best_idx = std::get<0>(best_probs[0]);
                int best_fj = std::get<2>(best_probs[0]);
                float best_p = std::get<1>(best_probs[0]);
                if (best_fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found second matched boosted boson: idx = " << best_fj << ", original idx = " << best_idx << std::endl;
                    }
                    result[boosted_result_pos[best_idx]] = best_fj;
                    used_fatjets.insert(best_fj);
                } else if (DEBUG) {
                    std::cout << "No jets left after overlap removal for either boson." << std::endl;
                }

                // Last one: pick the remaining index
                int last_idx = remaining[best_idx == remaining[0] ? 1 : 0];
                auto [last_fj, last_prob] = get_best_fatjet(boosted_cands[last_idx], used_fatjets);
                if (last_fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found third matched boosted boson: idx = " << last_fj << ", original idx = " << last_idx << std::endl;
                    }
                    result[boosted_result_pos[last_idx]] = last_fj;
                    used_fatjets.insert(last_fj);
                }
            }
        } else if (DEBUG) {
            std::cout << "No boosted boson left to match." << std::endl;
        }
    }

    // Resolved
    std::vector<bool> is_boosted(3, false);
    is_boosted[0] = (result[4] != -1);
    is_boosted[1] = (result[9] != -1);
    is_boosted[2] = (result[10] != -1);

    std::vector<std::vector<std::vector<float>>> resolved_cands;
    std::vector<std::pair<int, int>> resolved_pos;

    if (!is_boosted[0]) {
        resolved_cands.push_back(cand_h);
        resolved_pos.emplace_back(2, 3);
    }
    if (!is_boosted[1]) {
        resolved_cands.push_back(cand_v1);
        resolved_pos.emplace_back(5, 6);
    }
    if (!is_boosted[2]) {
        resolved_cands.push_back(cand_v2);
        resolved_pos.emplace_back(7, 8);
    }

    if (DEBUG) {
        std::cout << "\nStarting Resolved Decays Selection" << std::endl;
        std::cout << "====================" << std::endl;
        std::cout << " Resolved assignments:" << std::endl;
        std::vector<std::string> res_labels = {"h", "v1", "v2"};
        for (size_t ri = 0; ri < 3; ++ri) {
            const auto& cand = (ri == 0 ? cand_h : (ri == 1 ? cand_v1 : cand_v2));
            std::cout << " " << res_labels[ri] << " = {";
            for (const auto& c : cand) {
                if (c[0] > 0 && c[0] <= 1) { // Skip garbage probs
                    std::cout << "[" << c[0] << ", " << static_cast<int>(c[1]) << ", " << static_cast<int>(c[2]) << "], ";
                }
            }
            std::cout << "}" << std::endl;
        }
        std::cout << "Resolved bosons to match: {";
        for (size_t i = 0; i < resolved_pos.size(); ++i) {
            std::string label = (resolved_pos[i].first == 2 ? "h" : (resolved_pos[i].first == 5 ? "v1" : "v2"));
            std::cout << label << ", ";
        }
        std::cout << "}" << std::endl;
    }

    if (!resolved_cands.empty()) {
        std::vector<float> res_top_probs(resolved_cands.size(), -std::numeric_limits<float>::infinity());
        for (size_t i = 0; i < resolved_cands.size(); ++i) {
            if (!resolved_cands[i].empty()) {
                res_top_probs[i] = resolved_cands[i][0][0];
            }
        }

        int res_first_idx = -1;
        float res_max_prob = -std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < resolved_cands.size(); ++i) {
            if (res_top_probs[i] > res_max_prob) {
                res_max_prob = res_top_probs[i];
                res_first_idx = static_cast<int>(i);
            }
        }

        if (DEBUG && res_first_idx != -1) {
            std::cout << "Boson with max reconstruction prob: index = " << res_first_idx << std::endl;
        }

        if (res_first_idx != -1) {
            auto [j1, j2, prob] = get_best_pair(resolved_cands[res_first_idx], used_jets);
            if (j1 != -1) {
                if (DEBUG) {
                    std::cout << "Selected pair: (" << j1 << ", " << j2 << ")" << std::endl;
                }
                result[resolved_pos[res_first_idx].first] = j1;
                result[resolved_pos[res_first_idx].second] = j2;
                used_jets.insert(j1);
                used_jets.insert(j2);
            }

            std::vector<int> res_remaining;
            for (size_t i = 0; i < resolved_cands.size(); ++i) {
                if (static_cast<int>(i) != res_first_idx && !resolved_cands[i].empty()) {
                    res_remaining.push_back(static_cast<int>(i));
                }
            }

            if (!res_remaining.empty()) {
                if (res_remaining.size() == 1) {
                    if (DEBUG) {
                        std::cout << "Only one remaining boson." << std::endl;
                    }
                    int idx = res_remaining[0];
                    auto [j1, j2, prob] = get_best_pair(resolved_cands[idx], used_jets);
                    if (j1 != -1) {
                        result[resolved_pos[idx].first] = j1;
                        result[resolved_pos[idx].second] = j2;
                        used_jets.insert(j1);
                        used_jets.insert(j2);
                    }
                } else {
                    if (DEBUG) {
                        std::cout << "\nAfter first selection:" << std::endl;
                    }
                    std::vector<std::tuple<int, float, int, int>> res_best_probs;
                    for (int idx : res_remaining) {
                        auto [j1, j2, prob] = get_best_pair(resolved_cands[idx], used_jets);
                        res_best_probs.push_back({idx, prob, j1, j2});
                    }

                    if (DEBUG) {
                        std::cout << "Best probs for remaining (idx, prob, j1, j2): {";
                        for (const auto& [idx, prob, j1, j2] : res_best_probs) {
                            std::cout << "(" << idx << ", " << prob << ", " << j1 << ", " << j2 << "), ";
                        }
                        std::cout << "}" << std::endl;
                    }

                    // Sort descending by prob
                    std::sort(res_best_probs.begin(), res_best_probs.end(),
                              [](const std::tuple<int, float, int, int>& a, const std::tuple<int, float, int, int>& b) {
                                  return std::get<1>(a) > std::get<1>(b);
                              });

                    int best_idx = std::get<0>(res_best_probs[0]);
                    int best_j1 = std::get<2>(res_best_probs[0]);
                    int best_j2 = std::get<3>(res_best_probs[0]);
                    float best_p = std::get<1>(res_best_probs[0]);
                    if (best_j1 != -1) {
                        result[resolved_pos[best_idx].first] = best_j1;
                        result[resolved_pos[best_idx].second] = best_j2;
                        used_jets.insert(best_j1);
                        used_jets.insert(best_j2);
                    } else if (DEBUG) {
                        std::cout << "No valid pairs for remaining bosons." << std::endl;
                    }

                    int last_idx = res_remaining[best_idx == res_remaining[0] ? 1 : 0];
                    auto [last_j1, last_j2, last_prob] = get_best_pair(resolved_cands[last_idx], used_jets);
                    if (last_j1 != -1) {
                        if (DEBUG) {
                            std::cout << "\nAfter second selection:" << std::endl;
                            std::cout << "Last pair: (" << last_j1 << ", " << last_j2 << ")" << std::endl;
                        }
                        result[resolved_pos[last_idx].first] = last_j1;
                        result[resolved_pos[last_idx].second] = last_j2;
                        used_jets.insert(last_j1);
                        used_jets.insert(last_j2);
                    }
                }
            } else if (DEBUG) {
                std::cout << "No remaining bosons, returning." << std::endl;
            }
        }
    }

    if (DEBUG) {
        std::cout << "\nStarting VBS Jets Selection" << std::endl;
        std::cout << "====================" << std::endl;
        std::cout << "vbs = {";
        for (const auto& c : cand_vbs) {
            std::cout << "[" << c[0] << ", " << static_cast<int>(c[1]) << ", " << static_cast<int>(c[2]) << "], ";
        }
        std::cout << "}" << std::endl;
    }

    auto [vbs_j1, vbs_j2, vbs_prob] = get_best_pair(cand_vbs, used_jets);
    if (vbs_j1 != -1 && vbs_prob >= 0.0f) {
        if (DEBUG) {
            std::cout << "vbs_best_pair = (" << vbs_j1 << ", " << vbs_j2 << "), vbs_best_p = " << vbs_prob << std::endl;
        }
        result[0] = vbs_j1;
        result[1] = vbs_j2;
    }

    if (DEBUG) {
        std::cout << "\n Selected boosted bosons = {";
        std::vector<std::string> bb_labels = {"bh", "bv1", "bv2"};
        for (int i = 0; i < 3; ++i) {
            int pos = boosted_result_pos[i];
            if (result[pos] != -1) {
                std::cout << bb_labels[i] << ": " << result[pos] << ", ";
            }
        }
        std::cout << "}" << std::endl;
        std::cout << "Selected resolved jet pairs: {";
        for (size_t i = 0; i < resolved_pos.size(); ++i) {
            int j1 = result[resolved_pos[i].first];
            int j2 = result[resolved_pos[i].second];
            if (j1 != -1 && j2 != -1) {
                std::string label = (resolved_pos[i].first == 2 ? "h" : (resolved_pos[i].first == 5 ? "v1" : "v2"));
                std::cout << label << ": (" << j1 << ", " << j2 << "), ";
            }
        }
        std::cout << "}" << std::endl;
        if (result[0] != -1 && result[1] != -1) {
            std::cout << "Selcted vbs pair = (" << result[0] << ", " << result[1] << ")" << std::endl;
        }
    }

    return result;
}


RNode SPANetRun2::SPANetInference::ParseSpanetInference(RNode df_){
    std::cout << " -> SPANetInferenceRun2::ParseSpanetInference()" << std::endl;

    auto df = df_.Define("_all_assignments", assign_all_objects_maxprob, {
        "spanet_vbs_assignment", "spanet_h_assignment", "spanet_bh_assignment",
        "spanet_v1_assignment", "spanet_v2_assignment", "spanet_bv1_assignment", "spanet_bv2_assignment",
        "spanet_vbs_detection", "spanet_h_detection", "spanet_bh_detection",
        "spanet_v1_detection", "spanet_v2_detection", "spanet_bv1_detection", "spanet_bv2_detection",
        "jet_eta", "jet_phi", "fatjet_eta", "fatjet_phi"
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
    df = df.Define("spanet_vbs1_pt", "vbs1_idx >= 0 ? jet_pt[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_eta", "vbs1_idx >= 0 ? jet_eta[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_phi", "vbs1_idx >= 0 ? jet_phi[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_mass", "vbs1_idx >= 0 ? jet_mass[vbs1_idx] : -999.0f")
           .Define("spanet_vbs2_pt", "vbs2_idx >= 0 ? jet_pt[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_eta", "vbs2_idx >= 0 ? jet_eta[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_phi", "vbs2_idx >= 0 ? jet_phi[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_mass", "vbs2_idx >= 0 ? jet_mass[vbs2_idx] : -999.0f")
           .Define("spanet_vbs_detajj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? abs(spanet_vbs1_eta - spanet_vbs2_eta) : -999.0f")
           .Define("spanet_vbs_mjj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_vbs1_pt, spanet_vbs1_eta, spanet_vbs1_phi, spanet_vbs1_mass) + "
                                     "ROOT::Math::PtEtaPhiMVector(spanet_vbs2_pt, spanet_vbs2_eta, spanet_vbs2_phi, spanet_vbs2_mass)).M() : -999.0f");

    // Resolved Higgs variables
    df = df.Define("spanet_h1_pt", "h1_idx >= 0 ? jet_pt[h1_idx] : -999.0f")
           .Define("spanet_h1_eta", "h1_idx >= 0 ? jet_eta[h1_idx] : -999.0f")
           .Define("spanet_h1_phi", "h1_idx >= 0 ? jet_phi[h1_idx] : -999.0f")
           .Define("spanet_h1_mass", "h1_idx >= 0 ? jet_mass[h1_idx] : -999.0f")
           .Define("spanet_h2_pt", "h2_idx >= 0 ? jet_pt[h2_idx] : -999.0f")
           .Define("spanet_h2_eta", "h2_idx >= 0 ? jet_eta[h2_idx] : -999.0f")
           .Define("spanet_h2_phi", "h2_idx >= 0 ? jet_phi[h2_idx] : -999.0f")
           .Define("spanet_h2_mass", "h2_idx >= 0 ? jet_mass[h2_idx] : -999.0f")
           .Define("spanet_h_mjj", "h1_idx >= 0 && h2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_h1_pt, spanet_h1_eta, spanet_h1_phi, spanet_h1_mass) + "
                                   "ROOT::Math::PtEtaPhiMVector(spanet_h2_pt, spanet_h2_eta, spanet_h2_phi, spanet_h2_mass)).M() : -999.0f");

    // Boosted Higgs variables
    df = df.Define("spanet_bh_eta", "bh_idx >= 0 ? fatjet_eta[bh_idx] : -999.0f")
             .Define("spanet_bh_phi", "bh_idx >= 0 ? fatjet_phi[bh_idx] : -999.0f")
             .Define("spanet_bh_msoftdrop", "bh_idx >= 0 ? fatjet_msoftdrop[bh_idx] : -999.0f")
             .Define("spanet_bh_pt", "bh_idx >= 0 ? fatjet_pt[bh_idx] : -999.0f")
             .Define("spanet_bh_HbbScore", "bh_idx >= 0 ? fatjet_HbbScore[bh_idx] : -999.0f")
             .Define("spanet_bh_WqqScore", "bh_idx >= 0 ? fatjet_WqqScore[bh_idx] : -999.0f");

    // Resolved V1 variables
    df = df.Define("spanet_v1_j1_pt", "v1_j1_idx >= 0 ? jet_pt[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_eta", "v1_j1_idx >= 0 ? jet_eta[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_phi", "v1_j1_idx >= 0 ? jet_phi[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_mass", "v1_j1_idx >= 0 ? jet_mass[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j2_pt", "v1_j2_idx >= 0 ? jet_pt[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_eta", "v1_j2_idx >= 0 ? jet_eta[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_phi", "v1_j2_idx >= 0 ? jet_phi[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_mass", "v1_j2_idx >= 0 ? jet_mass[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_mjj", "v1_j1_idx >= 0 && v1_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v1_j1_pt, spanet_v1_j1_eta, spanet_v1_j1_phi, spanet_v1_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v1_j2_pt, spanet_v1_j2_eta, spanet_v1_j2_phi, spanet_v1_j2_mass)).M() : -999.0f");

    // Resolved V2 variables
    df = df.Define("spanet_v2_j1_pt", "v2_j1_idx >= 0 ? jet_pt[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_eta", "v2_j1_idx >= 0 ? jet_eta[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_phi", "v2_j1_idx >= 0 ? jet_phi[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_mass", "v2_j1_idx >= 0 ? jet_mass[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j2_pt", "v2_j2_idx >= 0 ? jet_pt[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_eta", "v2_j2_idx >= 0 ? jet_eta[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_phi", "v2_j2_idx >= 0 ? jet_phi[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_mass", "v2_j2_idx >= 0 ? jet_mass[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_mjj", "v2_j1_idx >= 0 && v2_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v2_j1_pt, spanet_v2_j1_eta, spanet_v2_j1_phi, spanet_v2_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v2_j2_pt, spanet_v2_j2_eta, spanet_v2_j2_phi, spanet_v2_j2_mass)).M() : -999.0f");

    // Boosted V1 variables
    df = df.Define("spanet_bv1_eta", "bv1_idx >= 0 ? fatjet_eta[bv1_idx] : -999.0f")
             .Define("spanet_bv1_phi", "bv1_idx >= 0 ? fatjet_phi[bv1_idx] : -999.0f")
             .Define("spanet_bv1_msoftdrop", "bv1_idx >= 0 ? fatjet_msoftdrop[bv1_idx] : -999.0f")
             .Define("spanet_bv1_pt", "bv1_idx >= 0 ? fatjet_pt[bv1_idx] : -999.0f")
             .Define("spanet_bv1_HbbScore", "bv1_idx >= 0 ? fatjet_HbbScore[bv1_idx] : -999.0f")
             .Define("spanet_bv1_WqqScore", "bv1_idx >= 0 ? fatjet_WqqScore[bv1_idx] : -999.0f");

    // Boosted V2 variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? fatjet_eta[bv2_idx] : -999.0f")
             .Define("spanet_bv2_phi", "bv2_idx >= 0 ? fatjet_phi[bv2_idx] : -999.0f")
             .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? fatjet_msoftdrop[bv2_idx] : -999.0f")
             .Define("spanet_bv2_pt", "bv2_idx >= 0 ? fatjet_pt[bv2_idx] : -999.0f")
             .Define("spanet_bv2_HbbScore", "bv1_idx >= 0 ? fatjet_HbbScore[bv1_idx] : -999.0f")
             .Define("spanet_bv2_WqqScore", "bv1_idx >= 0 ? fatjet_WqqScore[bv1_idx] : -999.0f");

    return df;
}
