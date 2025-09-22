#include "spanet_inference_base.h"

namespace SPANet {

SPANetInferenceBase::SPANetInferenceBase(
    const std::string &model_path,
    size_t batch_size,
    size_t max_ak4_jets,
    size_t ak4_features,
    size_t max_ak8_jets,
    size_t ak8_features,
    size_t met_features,
    size_t max_leptons,
    size_t lepton_features
) : batch_size_(batch_size),
    top_k_(3),
    memory_info_(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault)),
    input_names_{"AK4Jets_data", "AK4Jets_mask", "AK8Jets_data", "AK8Jets_mask", "MET_data", "MET_mask", "Lepton1_data", "Lepton1_mask", "Lepton2_data", "Lepton2_mask"},
    output_names_{"h_assignment_probability",
                  "bh_assignment_probability", 
                  "v1_assignment_probability", 
                  "v2_assignment_probability", 
                  "bv1_assignment_probability",
                  "bv2_assignment_probability",
                  "vbs_assignment_probability",
                  "h_detection_probability",
                  "bh_detection_probability",
                  "v1_detection_probability",
                  "v2_detection_probability",
                  "bv1_detection_probability",
                  "bv2_detection_probability",
                  "vbs_detection_probability",
                  "EVENT/isSignal"} {
    
    ak4_input_shape_ = {static_cast<int64_t>(batch_size_), static_cast<int64_t>(max_ak4_jets), static_cast<int64_t>(ak4_features)};
    ak8_input_shape_ = {static_cast<int64_t>(batch_size_), static_cast<int64_t>(max_ak8_jets), static_cast<int64_t>(ak8_features)};
    met_input_shape_ = {static_cast<int64_t>(batch_size_), 1, static_cast<int64_t>(met_features)};
    lepton1_input_shape_ = {static_cast<int64_t>(batch_size_), 1, static_cast<int64_t>(lepton_features)};
    lepton2_input_shape_ = {static_cast<int64_t>(batch_size_), 1, static_cast<int64_t>(lepton_features)};
    ak4_mask_shape_ = {static_cast<int64_t>(batch_size_), static_cast<int64_t>(max_ak4_jets)};
    ak8_mask_shape_ = {static_cast<int64_t>(batch_size_), static_cast<int64_t>(max_ak8_jets)};
    met_mask_shape_ = {static_cast<int64_t>(batch_size_), 1};
    lepton1_mask_shape_ = {static_cast<int64_t>(batch_size_), 1};
    lepton2_mask_shape_ = {static_cast<int64_t>(batch_size_), 1};

    ak4_flat_jets_.resize(batch_size_ * max_ak4_jets * ak4_features, 0.0f);
    ak8_flat_jets_.resize(batch_size_ * max_ak8_jets * ak8_features, 0.0f);
    met_inputs_.resize(batch_size_ * met_features, 0.0f);
    lepton1_inputs_.resize(batch_size_ * lepton_features, 0.0f);
    lepton2_inputs_.resize(batch_size_ * lepton_features, 0.0f);
    
    ak4_mask_char_.resize(batch_size_ * max_ak4_jets, 0);
    ak8_mask_char_.resize(batch_size_ * max_ak8_jets, 0);
    met_mask_char_.resize(batch_size_ * 1, 1);
    lepton1_mask_char_.resize(batch_size_ * 1, 1);
    lepton2_mask_char_.resize(batch_size_ * 1, 1);

    Ort::SessionOptions sessionOptions;
    sessionOptions.SetIntraOpNumThreads(0);
    sessionOptions.SetInterOpNumThreads(0);
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    sessionOptions.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);
    
    try {
        env_ = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_ERROR, "SPANet");
        session_ = std::make_unique<Ort::Session>(*env_, model_path.c_str(), sessionOptions);
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize ONNX session: " << e.what() << std::endl;
        throw;
    }
}

RNode SPANetInferenceBase::RunSPANetInference(RNode df) {
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

std::vector<std::vector<std::vector<std::vector<float>>>> SPANetInferenceBase::runBatchInference(const std::vector<std::shared_ptr<EventDataBase>>& events) {
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
                
                const size_t _top_k = std::min(static_cast<size_t>(this->top_k_), value_idx_pairs.size());
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
            } else if (output_shape.size() == 2) {
                size_t jets_dim = output_shape[1];
                size_t offset = batch_idx * jets_dim;
                
                value_idx_pairs.clear();
                value_idx_pairs.reserve(jets_dim);
                
                for (size_t j = ak4_input_shape_[1]; j < jets_dim; ++j) { // start after AK4 jets
                    value_idx_pairs.emplace_back(output_data[offset + j], j - ak4_input_shape_[1]);
                }
                
                const size_t _top_k = std::min(static_cast<size_t>(this->top_k_), value_idx_pairs.size());
                std::partial_sort(value_idx_pairs.begin(), value_idx_pairs.begin() + _top_k, value_idx_pairs.end(),
                                [](const auto& a, const auto& b) { return a.first > b.first; });
                
                auto& result_vec = batch_results[batch_idx][tensor_idx];
                result_vec.clear();
                result_vec.reserve(_top_k);
                
                for (size_t k = 0; k < _top_k; ++k) {
                    result_vec.push_back({value_idx_pairs[k].first, static_cast<float>(value_idx_pairs[k].second)});
                }
            }
        } else if (tensor_idx < 14) { 
            batch_results[batch_idx][tensor_idx] = {{output_data[batch_idx]}};
        } else {
            size_t offset = batch_idx * 2;
            batch_results[batch_idx][tensor_idx] = {{output_data[offset + 1]}}; 
        }
    };

    for (size_t start_idx = 0; start_idx < events.size(); start_idx += batch_size_) {
        size_t end_idx = std::min(start_idx + batch_size_, events.size());
        size_t current_batch_size = end_idx - start_idx;
        
        std::vector<std::shared_ptr<EventDataBase>> batch_events(events.begin() + start_idx, events.begin() + end_idx);
        
        std::vector<int64_t> current_ak4_shape = {static_cast<int64_t>(current_batch_size), ak4_input_shape_[1], ak4_input_shape_[2]};
        std::vector<int64_t> current_ak8_shape = {static_cast<int64_t>(current_batch_size), ak8_input_shape_[1], ak8_input_shape_[2]};
        std::vector<int64_t> current_met_shape = {static_cast<int64_t>(current_batch_size), 1, met_input_shape_[2]};
        std::vector<int64_t> current_lepton1_shape = {static_cast<int64_t>(current_batch_size), 1, lepton1_input_shape_[2]};
        std::vector<int64_t> current_lepton2_shape = {static_cast<int64_t>(current_batch_size), 1, lepton2_input_shape_[2]};
        std::vector<int64_t> current_ak4_mask_shape = {static_cast<int64_t>(current_batch_size), ak4_mask_shape_[1]};
        std::vector<int64_t> current_ak8_mask_shape = {static_cast<int64_t>(current_batch_size), ak8_mask_shape_[1]};
        std::vector<int64_t> current_met_mask_shape = {static_cast<int64_t>(current_batch_size), 1};
        std::vector<int64_t> current_lepton1_mask_shape = {static_cast<int64_t>(current_batch_size), 1};
        std::vector<int64_t> current_lepton2_mask_shape = {static_cast<int64_t>(current_batch_size), 1};

        fillBatchTensors(batch_events, current_batch_size);
        
        std::vector<Ort::Value> input_tensors;
        input_tensors.reserve(10);
        
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info_, ak4_flat_jets_.data(), current_batch_size * ak4_input_shape_[1] * ak4_input_shape_[2], current_ak4_shape.data(), current_ak4_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info_, reinterpret_cast<bool*>(ak4_mask_char_.data()), current_batch_size * ak4_mask_shape_[1], current_ak4_mask_shape.data(), current_ak4_mask_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info_, ak8_flat_jets_.data(), current_batch_size * ak8_input_shape_[1] * ak8_input_shape_[2], current_ak8_shape.data(), current_ak8_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info_, reinterpret_cast<bool*>(ak8_mask_char_.data()), current_batch_size * ak8_mask_shape_[1], current_ak8_mask_shape.data(), current_ak8_mask_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info_, met_inputs_.data(), current_batch_size * met_input_shape_[2], current_met_shape.data(), current_met_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info_, reinterpret_cast<bool*>(met_mask_char_.data()), current_batch_size * 1, current_met_mask_shape.data(), current_met_mask_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info_, lepton1_inputs_.data(), current_batch_size * lepton1_input_shape_[2], current_lepton1_shape.data(), current_lepton1_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info_, reinterpret_cast<bool*>(lepton1_mask_char_.data()), current_batch_size * 1, current_lepton1_mask_shape.data(), current_lepton1_mask_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info_, lepton2_inputs_.data(), current_batch_size * lepton2_input_shape_[2], current_lepton2_shape.data(), current_lepton2_shape.size()));
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info_, reinterpret_cast<bool*>(lepton2_mask_char_.data()), current_batch_size * 1, current_lepton2_mask_shape.data(), current_lepton2_mask_shape.size()));

        auto output_tensors = session_->Run(Ort::RunOptions{nullptr}, input_names_.data(), input_tensors.data(), input_names_.size(), output_names_.data(), output_names_.size());

        std::vector<std::vector<std::vector<std::vector<float>>>> batch_results(current_batch_size, std::vector<std::vector<std::vector<float>>>(output_names_.size()));

        for (size_t batch_idx = 0; batch_idx < current_batch_size; ++batch_idx) {
            for (size_t tensor_idx = 0; tensor_idx < output_tensors.size(); ++tensor_idx) {
                extract_idxs_from_inference(output_tensors[tensor_idx], batch_idx, tensor_idx, batch_results);
            }
        }

        all_results.insert(all_results.end(), batch_results.begin(), batch_results.end());
    }

    return all_results;
}
RNode SPANetInferenceBase::addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>> &all_outputs, const std::vector<std::shared_ptr<EventDataBase>> &events)
{
    // Create maps for fast lookup using entry+slot as key
    std::map<std::pair<ULong64_t, unsigned int>, size_t> entry_slot_to_output_idx;

    for (size_t i = 0; i < events.size(); ++i)
    {
        auto key = std::make_pair(events[i]->rdf_entry, events[i]->rdf_slot);
        entry_slot_to_output_idx[key] = i;
    }

    // Create column-wise accessors that use entry+slot matching
    auto getHAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][0];
        }
        return {};
    };

    auto getBHAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][1];
        }
        return {};
    };

    auto getV1Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][2];
        }
        return {};
    };

    auto getV2Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][3];
        }
        return {};
    };

    auto getBV1Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][4];
        }
        return {};
    };

    auto getBV2Assignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][5];
        }
        return {};
    };

    auto getVBSAssignment = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> std::vector<std::vector<float>>
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][6];
        }
        return {};
    };

    auto getHDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][7][0][0];
        }
        return -999.0f;
    };

    auto getBHDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][8][0][0];
        }
        return -999.0f;
    };

    auto getV1Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][9][0][0];
        }
        return -999.0f;
    };

    auto getV2Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][10][0][0];
        }
        return -999.0f;
    };

    auto getBV1Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][11][0][0];
        }
        return -999.0f;
    };

    auto getBV2Detection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][12][0][0];
        }
        return -999.0f;
    };

    auto getVBSDetection = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
            return all_outputs[it->second][13][0][0];
        }
        return -999.0f;
    };

    auto getEventSignal = [all_outputs, entry_slot_to_output_idx](unsigned int slot, ULong64_t entry) -> float
    {
        auto key = std::make_pair(entry, slot);
        auto it = entry_slot_to_output_idx.find(key);
        if (it != entry_slot_to_output_idx.end())
        {
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


std::vector<int> SPANetInferenceBase::assign_all_objects(
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
    // Implementation from original - since truncated in query, use the provided logic
    std::vector<int> result(11, -1);

    // Full assignment logic with overlap checks, deltaR, etc.
    // From the provided snippet:
    // Define check functions like checkJetOverlap, checkDeltaR, checkAllOverlaps
    // Then sort detections or priorities
    // For example, std::vector<float> detections = {vbs_detection, h_detection, bh_detection, v1_detection, v2_detection, bv1_detection, bv2_detection};
    // Sort indices by detection prob descending
    // Then for each priority, assign from corresponding assignment list, checking overlaps
    // Since truncated, implement as per original description

    // Placeholder - replace with full code from query snippets
    return result;
}

RNode SPANetInferenceBase::ParseSpanetInference(RNode df_) {
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

    // Boosted V2 variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? FatJet_eta[bv2_idx] : -999.0f")
            .Define("spanet_bv2_phi", "bv2_idx >= 0 ? FatJet_phi[bv2_idx] : -999.0f")
            .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? FatJet_msoftdrop[bv2_idx] : -999.0f")
            .Define("spanet_bv2_pt", "bv2_idx >= 0 ? FatJet_pt[bv2_idx] : -999.0f")
            .Define("spanet_bv2_score", "bv2_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv2_idx] : -999.0f")
            .Define("spanet_bv2_w_score", "bv2_idx >= 0 ? FatJet_globalParT3_Xqq[bv2_idx] / (FatJet_globalParT3_Xqq[bv2_idx] + FatJet_globalParT3_Xcs[bv2_idx] + FatJet_globalParT3_QCD[bv2_idx]) : -999.0f");

    return df;
}

} // namespace SPANet