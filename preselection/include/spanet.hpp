#pragma once

#ifndef SPANET_H
#define SPANET_H

#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;

namespace SPANet {
    class SPANetInference {
    public:
        SPANetInference(const std::string &model_path, size_t batch_size = 64) 
            : batch_size(batch_size),
              memory_info(Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault)),
              input_names{"AK4Jets_data", "AK4Jets_mask", "AK8Jets_data", "AK8Jets_mask", "MET_data", "MET_mask"},
              output_names{"h_assignment_probability",
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
                          "EVENT/isSignal"},
              ak4_input_shape{static_cast<int64_t>(batch_size), 10, 8},
              ak8_input_shape{static_cast<int64_t>(batch_size), 3, 11},
              met_input_shape{static_cast<int64_t>(batch_size), 1, 3},
              ak4_mask_shape{static_cast<int64_t>(batch_size), 10},
              ak8_mask_shape{static_cast<int64_t>(batch_size), 3},
              met_mask_shape{static_cast<int64_t>(batch_size), 1},
              top_k(3) {
            
            Ort::SessionOptions sessionOptions;
            sessionOptions.SetIntraOpNumThreads(0);
            sessionOptions.SetInterOpNumThreads(0);
            sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
            sessionOptions.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);
            
            try {
                env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_ERROR, "SPANet");
                session = std::make_unique<Ort::Session>(*env, model_path.c_str(), sessionOptions);
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize ONNX session: " << e.what() << std::endl;
                throw;
            }
            
            ak4_flat_jets.resize(batch_size * 10 * 8, 0.0f);
            ak8_flat_jets.resize(batch_size * 3 * 11, 0.0f);
            met_inputs.resize(batch_size * 3, 0.0f);
            
            ak4_mask_char.resize(batch_size * 10, 0);
            ak8_mask_char.resize(batch_size * 3, 0);
            met_mask_char.resize(batch_size * 1, 1);
        }
        
        RNode RunSPANetInference(RNode df);

    private:
        std::unique_ptr<Ort::Env> env;
        std::unique_ptr<Ort::Session> session;
        size_t batch_size;
        int64_t top_k;
        Ort::MemoryInfo memory_info;
        
        std::vector<const char*> input_names;
        std::vector<const char*> output_names;
        std::vector<int64_t> ak4_input_shape;
        std::vector<int64_t> ak8_input_shape;
        std::vector<int64_t> met_input_shape;
        std::vector<int64_t> ak4_mask_shape;
        std::vector<int64_t> ak8_mask_shape;
        std::vector<int64_t> met_mask_shape;
        
        std::vector<float> ak4_flat_jets;
        std::vector<float> ak8_flat_jets;
        std::vector<float> met_inputs;
        std::vector<char> ak4_mask_char;
        std::vector<char> ak8_mask_char;
        std::vector<char> met_mask_char;
        
        struct EventData {
            std::vector<float> ak4_pt, ak4_eta, ak4_phi, ak4_mass;
            std::vector<int> ak4_isTightBTag, ak4_isMediumBTag, ak4_isLooseBTag;
            std::vector<float> ak8_pt, ak8_eta, ak8_phi, ak8_mass;
            std::vector<float> ak8_XbbScore, ak8_XqqScore, ak8_XccScore, ak8_XcsScore, ak8_XqcdScore;
            std::vector<unsigned char> ak8_nConstituents;
            float met_pt, met_phi;
        };
        
        std::vector<std::vector<std::vector<std::vector<float>>>> runBatchInference(const std::vector<EventData>& events);
        void fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size);
        
        std::vector<EventData> extractEventsFromDataFrame(RNode df);
        
        RNode addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs);
    };
}

inline std::vector<SPANet::SPANetInference::EventData> SPANet::SPANetInference::extractEventsFromDataFrame(RNode df) {
    using ROOT::VecOps::RVec;
    using RVecF = RVec<float>;
    using RVecI = RVec<int>;
    using RVecUC = RVec<unsigned char>;
    
    auto ak4_pt_vec = df.Take<RVecF>("Jet_pt").GetValue();
    auto ak4_eta_vec = df.Take<RVecF>("Jet_eta").GetValue();
    auto ak4_phi_vec = df.Take<RVecF>("Jet_phi").GetValue();
    auto ak4_mass_vec = df.Take<RVecF>("Jet_mass").GetValue();
    auto ak4_isTightBTag_vec = df.Take<RVecI>("Jet_isTightBTag").GetValue();
    auto ak4_isMediumBTag_vec = df.Take<RVecI>("Jet_isMediumBTag").GetValue();
    auto ak4_isLooseBTag_vec = df.Take<RVecI>("Jet_isLooseBTag").GetValue();
    
    auto ak8_pt_vec = df.Take<RVecF>("FatJet_pt").GetValue();
    auto ak8_eta_vec = df.Take<RVecF>("FatJet_eta").GetValue();
    auto ak8_phi_vec = df.Take<RVecF>("FatJet_phi").GetValue();
    auto ak8_mass_vec = df.Take<RVecF>("FatJet_mass").GetValue();
    
    auto ak8_nConstituents_vec = df.Take<RVecUC>("FatJet_nConstituents").GetValue();
    auto ak8_XbbScore_vec = df.Take<RVecF>("FatJet_globalParT3_Xbb").GetValue();
    auto ak8_XqqScore_vec = df.Take<RVecF>("FatJet_globalParT3_Xqq").GetValue();
    auto ak8_XccScore_vec = df.Take<RVecF>("FatJet_globalParT3_Xcc").GetValue();
    auto ak8_XcsScore_vec = df.Take<RVecF>("FatJet_globalParT3_Xcs").GetValue();
    auto ak8_XqcdScore_vec = df.Take<RVecF>("FatJet_globalParT3_QCD").GetValue();

    auto met_pt_vec = df.Take<float>("PuppiMET_pt").GetValue();
    auto met_phi_vec = df.Take<float>("PuppiMET_phi").GetValue();
    
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
        
        events.push_back(std::move(event));
    }
    
    return events;
}

inline std::vector<std::vector<std::vector<std::vector<float>>>> SPANet::SPANetInference::runBatchInference(const std::vector<EventData>& events) {
    std::vector<std::vector<std::vector<std::vector<float>>>> all_results;
    all_results.reserve(events.size());
    
    auto extract_idxs_from_inference = [this](const Ort::Value& output_tensor, size_t batch_idx, size_t tensor_idx, 
                                         std::vector<std::vector<std::vector<std::vector<float>>>>& batch_results) {
        const Ort::TensorTypeAndShapeInfo& shape_info = output_tensor.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> output_shape = shape_info.GetShape();
        const float* output_data = output_tensor.GetTensorData<float>();
        

        // ASSIGNMENT PROBABILITIES
        if (tensor_idx < 7) { 
            static thread_local std::vector<std::pair<float, size_t>> value_idx_pairs;
            
            if (output_shape.size() == 3) { // [batch, jets, jets] - always jet pair assignment
                size_t jets_dim1 = output_shape[1];
                size_t jets_dim2 = output_shape[2];
                size_t offset = batch_idx * jets_dim1 * jets_dim2;
                
                // Triangular size = n*(n-1)/2
                size_t tri_size = (jets_dim1 * (jets_dim1 - 1)) / 2;
                value_idx_pairs.clear();
                value_idx_pairs.reserve(tri_size);
                
                // Find values and their indices from upper triangular matrix only
                for (size_t i = 0; i < jets_dim1; ++i) {
                    size_t row_offset = offset + i * jets_dim2;
                    for (size_t j = i + 1; j < jets_dim2; ++j) { // Only upper triangular: i < j
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
        
            // DETECTION PROBABILITIES
        else if (tensor_idx < 14) { 
            // Detection probabilities have shape [batch_size], one value per event
            batch_results[batch_idx][tensor_idx] = {{output_data[batch_idx]}};
        } 
            
            // EVENT-LEVEL OUTPUT
        else {
            // Event-level output has shape [batch_size, 2] for [background, signal]
            size_t offset = batch_idx * 2;
            batch_results[batch_idx][tensor_idx] = {{output_data[offset + 1]}}; // Take signal probability
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

inline void SPANet::SPANetInference::fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size) {
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
        
        // Fill event inputs
        const size_t event_base_idx = batch_idx * 3;
        met_inputs[event_base_idx]     = std::log(event.met_pt + 1.0f);
        met_inputs[event_base_idx + 1] = std::sin(event.met_phi);
        met_inputs[event_base_idx + 2] = std::cos(event.met_phi);
        
        met_mask_char[batch_idx] = 1;
    }
}

inline RNode SPANet::SPANetInference::addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs) {
    // Create separate vectors for each output type to match expected format
    std::vector<std::vector<std::vector<float>>> h_assignment_prob, bh_assignment_prob, v1_assignment_prob, v2_assignment_prob, bv1_assignment_prob, bv2_assignment_prob, vbs_assignment_prob;
    std::vector<float> h_detection_prob, bh_detection_prob, v1_detection_prob, v2_detection_prob, bv1_detection_prob, bv2_detection_prob, vbs_detection_prob, event_signal_prob;

    // Reserve space for efficiency
    size_t n_events = all_outputs.size();
    h_assignment_prob.reserve(n_events);
    bh_assignment_prob.reserve(n_events);
    v1_assignment_prob.reserve(n_events);
    v2_assignment_prob.reserve(n_events);
    bv1_assignment_prob.reserve(n_events);
    bv2_assignment_prob.reserve(n_events);
    vbs_assignment_prob.reserve(n_events);
    h_detection_prob.reserve(n_events);
    bh_detection_prob.reserve(n_events);
    v1_detection_prob.reserve(n_events);
    v2_detection_prob.reserve(n_events);
    bv1_detection_prob.reserve(n_events);
    bv2_detection_prob.reserve(n_events);
    vbs_detection_prob.reserve(n_events);
    event_signal_prob.reserve(n_events);
    
    // Extract outputs for each event
    for (const auto& event_outputs : all_outputs) {
        h_assignment_prob.push_back(event_outputs[0]);
        bh_assignment_prob.push_back(event_outputs[1]);
        v1_assignment_prob.push_back(event_outputs[2]);
        v2_assignment_prob.push_back(event_outputs[3]);
        bv1_assignment_prob.push_back(event_outputs[4]);
        bv2_assignment_prob.push_back(event_outputs[5]);
        vbs_assignment_prob.push_back(event_outputs[6]);

        h_detection_prob.push_back(event_outputs[7][0][0]);
        bh_detection_prob.push_back(event_outputs[8][0][0]);
        v1_detection_prob.push_back(event_outputs[9][0][0]);
        v2_detection_prob.push_back(event_outputs[10][0][0]);
        bv1_detection_prob.push_back(event_outputs[11][0][0]);
        bv2_detection_prob.push_back(event_outputs[12][0][0]);
        vbs_detection_prob.push_back(event_outputs[13][0][0]);

        event_signal_prob.push_back(event_outputs[14][0][0]);
    }
    
    auto getHAssignment = [h_assignment_prob](unsigned int slot, ULong64_t entry) {
        return h_assignment_prob[entry];
    };
    auto getBHAssignment = [bh_assignment_prob](unsigned int slot, ULong64_t entry) {
        return bh_assignment_prob[entry];
    };
    auto getV1Assignment = [v1_assignment_prob](unsigned int slot, ULong64_t entry) {
        return v1_assignment_prob[entry];
    };
    auto getV2Assignment = [v2_assignment_prob](unsigned int slot, ULong64_t entry) {
        return v2_assignment_prob[entry];
    };
    auto getBV1Assignment = [bv1_assignment_prob](unsigned int slot, ULong64_t entry) {
        return bv1_assignment_prob[entry];
    };
    auto getBV2Assignment = [bv2_assignment_prob](unsigned int slot, ULong64_t entry) {
        return bv2_assignment_prob[entry];
    };
    auto getVBSAssignment = [vbs_assignment_prob](unsigned int slot, ULong64_t entry) {
        return vbs_assignment_prob[entry];
    };
    auto getHDetection = [h_detection_prob](unsigned int slot, ULong64_t entry) {
        return h_detection_prob[entry];
    };
    auto getBHDetection = [bh_detection_prob](unsigned int slot, ULong64_t entry) {
        return bh_detection_prob[entry];
    };
    auto getV1Detection = [v1_detection_prob](unsigned int slot, ULong64_t entry) {
        return v1_detection_prob[entry];
    };
    auto getV2Detection = [v2_detection_prob](unsigned int slot, ULong64_t entry) {
        return v2_detection_prob[entry];
    };
    auto getBV1Detection = [bv1_detection_prob](unsigned int slot, ULong64_t entry) {
        return bv1_detection_prob[entry];
    };
    auto getBV2Detection = [bv2_detection_prob](unsigned int slot, ULong64_t entry) {
        return bv2_detection_prob[entry];
    };
    auto getVBSDetection = [vbs_detection_prob](unsigned int slot, ULong64_t entry) {
        return vbs_detection_prob[entry];
    };
    auto getEventSignal = [event_signal_prob](unsigned int slot, ULong64_t entry) {
        return event_signal_prob[entry];
    };
    
    // Add all SPANet outputs using DefineSlotEntry for proper entry mapping
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

inline RNode SPANet::SPANetInference::RunSPANetInference(RNode df) {
    // Efficient data extraction with timing
    auto start_time = std::chrono::high_resolution_clock::now();
    auto events = extractEventsFromDataFrame(df);
    auto extract_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "Running batch inference on " << events.size() << " events..." << std::endl;
    auto all_outputs = runBatchInference(events);
    auto inference_time = std::chrono::high_resolution_clock::now();
    
    // Free event data memory as soon as possible
    events.clear();
    events.shrink_to_fit();
    
    auto result = addSPANetOutputsToDataFrame(df, all_outputs);
    auto end_time = std::chrono::high_resolution_clock::now();
    
    // Report timing for profiling
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
#endif
