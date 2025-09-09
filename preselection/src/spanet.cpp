#include "spanet.h"

std::vector<SPANet::SPANetInference::EventData> SPANet::SPANetInference::extractEventsFromDataFrame(RNode df) {
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
    
    // Extract RDF entry and slot information for proper matching
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
        event.ak8_XbbScore.assign(ak8_XbbScore_vec[i].begin(), ak8_XbbScore_vec[i].end());
        event.ak8_XqqScore.assign(ak8_XqqScore_vec[i].begin(), ak8_XqqScore_vec[i].end());
        event.ak8_XccScore.assign(ak8_XccScore_vec[i].begin(), ak8_XccScore_vec[i].end());
        event.ak8_XcsScore.assign(ak8_XcsScore_vec[i].begin(), ak8_XcsScore_vec[i].end());
        event.ak8_XqcdScore.assign(ak8_XqcdScore_vec[i].begin(), ak8_XqcdScore_vec[i].end());
        
        event.met_pt = met_pt_vec[i];
        event.met_phi = met_phi_vec[i];
        event.rdf_entry = rdf_entry_vec[i];
        event.rdf_slot = rdf_slot_vec[i];
        
        events.push_back(std::move(event));
    }
    
    return events;
}

std::vector<std::vector<std::vector<std::vector<float>>>> SPANet::SPANetInference::runBatchInference(const std::vector<EventData>& events) {
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
                
                for (size_t j = 10; j < jets_dim; ++j) { // start at 10 since the first 10 are AK4 jets
                    value_idx_pairs.emplace_back(output_data[offset + j], j - 10); // store value and index relative to AK8 jets
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
        
        std::vector<int64_t> current_ak4_shape = {static_cast<int64_t>(current_batch_size), 10, 8};
        std::vector<int64_t> current_ak8_shape = {static_cast<int64_t>(current_batch_size), 3, 11};
        std::vector<int64_t> current_met_shape = {static_cast<int64_t>(current_batch_size), 1, 3};
        std::vector<int64_t> current_ak4_mask_shape = {static_cast<int64_t>(current_batch_size), 10};
        std::vector<int64_t> current_ak8_mask_shape = {static_cast<int64_t>(current_batch_size), 3};
        std::vector<int64_t> current_met_mask_shape = {static_cast<int64_t>(current_batch_size), 1};

        fillBatchTensors(batch_events, current_batch_size);
        
        std::vector<Ort::Value> input_tensors;
        input_tensors.reserve(8);
        
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak4_flat_jets.data(), current_batch_size * 10 * 8,
            current_ak4_shape.data(), current_ak4_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak4_mask_char.data()), current_batch_size * 10,
            current_ak4_mask_shape.data(), current_ak4_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak8_flat_jets.data(), current_batch_size * 3 * 11,
            current_ak8_shape.data(), current_ak8_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak8_mask_char.data()), current_batch_size * 3,
            current_ak8_mask_shape.data(), current_ak8_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, met_inputs.data(), current_batch_size * 3,
            current_met_shape.data(), current_met_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(met_mask_char.data()), current_batch_size * 1,
            current_met_mask_shape.data(), current_met_mask_shape.size()));

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

void SPANet::SPANetInference::fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size) {
    std::fill(ak4_flat_jets.begin(), ak4_flat_jets.begin() + actual_batch_size * 10 * 8, 0.0f);
    std::fill(ak8_flat_jets.begin(), ak8_flat_jets.begin() + actual_batch_size * 3 * 11, 0.0f);
    std::fill(ak4_mask_char.begin(), ak4_mask_char.begin() + actual_batch_size * 10, 0);
    std::fill(ak8_mask_char.begin(), ak8_mask_char.begin() + actual_batch_size * 3, 0);
    
    for (size_t batch_idx = 0; batch_idx < actual_batch_size; ++batch_idx) {
        const EventData& event = events[batch_idx];
        
        const size_t max_ak4 = std::min<size_t>(10, event.ak4_pt.size());
        for (size_t i = 0; i < max_ak4; ++i) {
            ak4_mask_char[batch_idx * 10 + i] = 1;
            
            const size_t base_idx = batch_idx * 10 * 8 + i * 8;
            ak4_flat_jets[base_idx]     = std::log(event.ak4_mass[i] + 1.0f);
            ak4_flat_jets[base_idx + 1] = std::log(event.ak4_pt[i] + 1.0f);
            ak4_flat_jets[base_idx + 2] = event.ak4_eta[i];
            ak4_flat_jets[base_idx + 3] = std::sin(event.ak4_phi[i]);
            ak4_flat_jets[base_idx + 4] = std::cos(event.ak4_phi[i]);
            ak4_flat_jets[base_idx + 5] = static_cast<float>(event.ak4_isTightBTag[i]);
            ak4_flat_jets[base_idx + 6] = static_cast<float>(event.ak4_isMediumBTag[i]);
            ak4_flat_jets[base_idx + 7] = static_cast<float>(event.ak4_isLooseBTag[i]);
        }
    
        const size_t max_ak8 = std::min<size_t>(3, event.ak8_pt.size());
        for (size_t i = 0; i < max_ak8; ++i) {
            ak8_mask_char[batch_idx * 3 + i] = 1;
            
            const size_t base_idx = batch_idx * 3 * 11 + i * 11;
            ak8_flat_jets[base_idx]     = std::log(event.ak8_mass[i] + 1.0f);
            ak8_flat_jets[base_idx + 1] = std::log(event.ak8_pt[i] + 1.0f);
            ak8_flat_jets[base_idx + 2] = event.ak8_eta[i];
            ak8_flat_jets[base_idx + 3] = std::sin(event.ak8_phi[i]);
            ak8_flat_jets[base_idx + 4] = std::cos(event.ak8_phi[i]);
            ak8_flat_jets[base_idx + 5] = static_cast<float>(event.ak8_nConstituents[i]);
            ak8_flat_jets[base_idx + 6] = event.ak8_XbbScore[i];
            ak8_flat_jets[base_idx + 7] = event.ak8_XqqScore[i];
            ak8_flat_jets[base_idx + 8] = event.ak8_XccScore[i];
            ak8_flat_jets[base_idx + 9] = event.ak8_XcsScore[i];
            ak8_flat_jets[base_idx + 10] = event.ak8_XqcdScore[i];
        }
        
        const size_t event_base_idx = batch_idx * 3;
        met_inputs[event_base_idx]     = std::log(event.met_pt + 1.0f);
        met_inputs[event_base_idx + 1] = std::sin(event.met_phi);
        met_inputs[event_base_idx + 2] = std::cos(event.met_phi);
        
        met_mask_char[batch_idx] = 1;
    }
}

RNode SPANet::SPANetInference::addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs, const std::vector<EventData>& events) {
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

RNode SPANet::SPANetInference::RunSPANetInference(RNode df) {
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

std::vector<int> SPANet::SPANetInference::assign_all_objects(
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

RNode SPANet::SPANetInference::ParseSpanetInference(RNode df_){
    auto df = df_.Define("_all_assignments", assign_all_objects, {
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

    // Boosted V2 final_variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? FatJet_eta[bv2_idx] : -999.0f")
            .Define("spanet_bv2_phi", "bv2_idx >= 0 ? FatJet_phi[bv2_idx] : -999.0f")
            .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? FatJet_msoftdrop[bv2_idx] : -999.0f")
            .Define("spanet_bv2_pt", "bv2_idx >= 0 ? FatJet_pt[bv2_idx] : -999.0f")
            .Define("spanet_bv2_score", "bv2_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv2_idx] : -999.0f")
            .Define("spanet_bv2_w_score", "bv2_idx >= 0 ? FatJet_globalParT3_Xqq[bv2_idx] / (FatJet_globalParT3_Xqq[bv2_idx] + FatJet_globalParT3_Xcs[bv2_idx] + FatJet_globalParT3_QCD[bv2_idx]) : -999.0f");

    return df;
}
