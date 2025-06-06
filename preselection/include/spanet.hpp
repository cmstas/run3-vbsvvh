#pragma once

#ifndef SPANET_H
#define SPANET_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <random>

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
                          "v_assignment_probability", 
                          "bh_assignment_probability", 
                          "bv_assignment_probability",
                          "vbs_assignment_probability",
                          "h_detection_probability",
                          "v_detection_probability",
                          "bh_detection_probability",
                          "bv_detection_probability",
                          "vbs_detection_probability",
                          "EVENT/isSignal"},
              ak4_input_shape{static_cast<int64_t>(batch_size), 10, 8},
              ak8_input_shape{static_cast<int64_t>(batch_size), 3, 8},
              event_input_shape{static_cast<int64_t>(batch_size), 1, 3},
              ak4_mask_shape{static_cast<int64_t>(batch_size), 10},
              ak8_mask_shape{static_cast<int64_t>(batch_size), 3},
              event_mask_shape{static_cast<int64_t>(batch_size), 1} {
            
            Ort::SessionOptions sessionOptions;
            // Configure for maximum CPU core usage
            sessionOptions.SetIntraOpNumThreads(0);  // 0 means use all available cores for an operation
            sessionOptions.SetInterOpNumThreads(0);  // 0 means use all available cores between operations
            sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
            sessionOptions.SetExecutionMode(ExecutionMode::ORT_SEQUENTIAL);
            
            try {
                env = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_ERROR, "SPANet");
                session = std::make_unique<Ort::Session>(*env, model_path.c_str(), sessionOptions);
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize ONNX session: " << e.what() << std::endl;
                throw;
            }
            
            // Pre-allocate buffers for batch processing
            ak4_flat_jets.resize(batch_size * 10 * 8, 0.0f);
            ak8_flat_jets.resize(batch_size * 3 * 8, 0.0f);
            event_inputs.resize(batch_size * 3, 0.0f);
            
            ak4_mask_char.resize(batch_size * 10, 0);
            ak8_mask_char.resize(batch_size * 3, 0);
            event_mask_char.resize(batch_size * 1, 1);
        }
        
        RNode RunSPANetInference(RNode df);

    private:
        std::unique_ptr<Ort::Env> env;
        std::unique_ptr<Ort::Session> session;
        size_t batch_size;
        Ort::MemoryInfo memory_info;
        
        // Pre-defined names and shapes
        std::vector<const char*> input_names;
        std::vector<const char*> output_names;
        std::vector<int64_t> ak4_input_shape;
        std::vector<int64_t> ak8_input_shape;
        std::vector<int64_t> event_input_shape;
        std::vector<int64_t> ak4_mask_shape;
        std::vector<int64_t> ak8_mask_shape;
        std::vector<int64_t> event_mask_shape;
        
        // Pre-allocated buffers
        std::vector<float> ak4_flat_jets;
        std::vector<float> ak8_flat_jets;
        std::vector<float> event_inputs;
        std::vector<char> ak4_mask_char;
        std::vector<char> ak8_mask_char;
        std::vector<char> event_mask_char;
        
        // Event data structure
        struct EventData {
            std::vector<float> ak4_pt, ak4_eta, ak4_phi, ak4_mass;
            std::vector<int> ak4_isTightBTag, ak4_isMediumBTag, ak4_isLooseBTag;
            std::vector<float> ak8_pt, ak8_eta, ak8_phi, ak8_mass;
            std::vector<float> ak8_HbbScore, ak8_WqqScore;
            std::vector<unsigned char> ak8_nConstituents;
            float met_pt, met_phi;
        };
        
        // Batch processing functions
        std::vector<std::vector<std::vector<float>>> runBatchInference(const std::vector<EventData>& events);
        void fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size);
        
        // Extract data from RDataFrame
        std::vector<EventData> extractEventsFromDataFrame(RNode df);
        
        // Create new dataframe with SPANet outputs
        RNode addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<float>>>& all_outputs);
    };
}

inline std::vector<SPANet::SPANetInference::EventData> SPANet::SPANetInference::extractEventsFromDataFrame(RNode df) {
    using ROOT::VecOps::RVec;
    using RVecF = RVec<float>;
    using RVecI = RVec<int>;
    using RVecUC = RVec<unsigned char>;
    
    // Take all the required columns from the dataframe
    auto ak4_pt_vec = df.Take<RVecF>("ak4jet_pt").GetValue();
    auto ak4_eta_vec = df.Take<RVecF>("ak4jet_eta").GetValue();
    auto ak4_phi_vec = df.Take<RVecF>("ak4jet_phi").GetValue();
    auto ak4_mass_vec = df.Take<RVecF>("ak4jet_mass").GetValue();
    auto ak4_isTightBTag_vec = df.Take<RVecI>("ak4jet_isTightBTag").GetValue();
    auto ak4_isMediumBTag_vec = df.Take<RVecI>("ak4jet_isMediumBTag").GetValue();
    auto ak4_isLooseBTag_vec = df.Take<RVecI>("ak4jet_isLooseBTag").GetValue();
    
    auto ak8_pt_vec = df.Take<RVecF>("ak8jet_pt").GetValue();
    auto ak8_eta_vec = df.Take<RVecF>("ak8jet_eta").GetValue();
    auto ak8_phi_vec = df.Take<RVecF>("ak8jet_phi").GetValue();
    auto ak8_mass_vec = df.Take<RVecF>("ak8jet_mass").GetValue();
    auto ak8_HbbScore_vec = df.Take<RVecF>("ak8jet_xbbvsqcd").GetValue();
    auto ak8_WqqScore_vec = df.Take<RVecF>("ak8jet_xqqvsqcd").GetValue();
    auto ak8_nConstituents_vec = df.Take<RVecUC>("ak8jet_nConstituents").GetValue();
    
    auto met_pt_vec = df.Take<float>("MET_pt").GetValue();
    auto met_phi_vec = df.Take<float>("MET_phi").GetValue();
    
    // Convert to EventData vector
    std::vector<EventData> events;
    size_t n_events = ak4_pt_vec.size();
    events.reserve(n_events);
    
    for (size_t i = 0; i < n_events; ++i) {
        EventData event;
        
        // Convert RVec to std::vector
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
        event.ak8_HbbScore.assign(ak8_HbbScore_vec[i].begin(), ak8_HbbScore_vec[i].end());
        event.ak8_WqqScore.assign(ak8_WqqScore_vec[i].begin(), ak8_WqqScore_vec[i].end());
        event.ak8_nConstituents.assign(ak8_nConstituents_vec[i].begin(), ak8_nConstituents_vec[i].end());
        
        event.met_pt = met_pt_vec[i];
        event.met_phi = met_phi_vec[i];
        
        events.push_back(std::move(event));
    }
    
    return events;
}

inline std::vector<std::vector<std::vector<float>>> SPANet::SPANetInference::runBatchInference(const std::vector<EventData>& events) {
    std::vector<std::vector<std::vector<float>>> all_results;
    all_results.reserve(events.size());
    
    // Process events in batches
    for (size_t start_idx = 0; start_idx < events.size(); start_idx += batch_size) {
        size_t end_idx = std::min(start_idx + batch_size, events.size());
        size_t current_batch_size = end_idx - start_idx;
        
        // Extract batch
        std::vector<EventData> batch_events(events.begin() + start_idx, events.begin() + end_idx);
        
        // Update shapes for current batch size
        std::vector<int64_t> current_ak4_shape = {static_cast<int64_t>(current_batch_size), 10, 8};
        std::vector<int64_t> current_ak8_shape = {static_cast<int64_t>(current_batch_size), 3, 8};
        std::vector<int64_t> current_event_shape = {static_cast<int64_t>(current_batch_size), 1, 3};
        std::vector<int64_t> current_ak4_mask_shape = {static_cast<int64_t>(current_batch_size), 10};
        std::vector<int64_t> current_ak8_mask_shape = {static_cast<int64_t>(current_batch_size), 3};
        std::vector<int64_t> current_event_mask_shape = {static_cast<int64_t>(current_batch_size), 1};
        
        fillBatchTensors(batch_events, current_batch_size);
        
        // Create input tensors
        std::vector<Ort::Value> input_tensors;
        input_tensors.reserve(6);
        
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak4_flat_jets.data(), current_batch_size * 10 * 8,
            current_ak4_shape.data(), current_ak4_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak4_mask_char.data()), current_batch_size * 10,
            current_ak4_mask_shape.data(), current_ak4_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, ak8_flat_jets.data(), current_batch_size * 3 * 8,
            current_ak8_shape.data(), current_ak8_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(ak8_mask_char.data()), current_batch_size * 3,
            current_ak8_mask_shape.data(), current_ak8_mask_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<float>(
            memory_info, event_inputs.data(), current_batch_size * 3,
            current_event_shape.data(), current_event_shape.size()));
            
        input_tensors.push_back(Ort::Value::CreateTensor<bool>(
            memory_info, reinterpret_cast<bool*>(event_mask_char.data()), current_batch_size * 1,
            current_event_mask_shape.data(), current_event_mask_shape.size()));
        
        // Run batch inference
        auto output_tensors = session->Run(
            Ort::RunOptions{nullptr},
            input_names.data(), input_tensors.data(), input_tensors.size(),
            output_names.data(), output_names.size());
        
        // Process batch outputs
        std::vector<std::vector<std::vector<float>>> batch_results(current_batch_size);

        for (size_t tensor_idx = 0; tensor_idx < output_tensors.size(); ++tensor_idx) {
            Ort::Value& output_tensor = output_tensors[tensor_idx];
            const Ort::TensorTypeAndShapeInfo& shape_info = output_tensor.GetTensorTypeAndShapeInfo();
            std::vector<int64_t> output_shape = shape_info.GetShape();
            const float* output_data = output_tensor.GetTensorData<float>();
            
            for (size_t batch_idx = 0; batch_idx < current_batch_size; ++batch_idx) {
                if (batch_results[batch_idx].size() <= tensor_idx) {
                    batch_results[batch_idx].resize(output_tensors.size());
                }
                
                if (tensor_idx < 5) { // Assignment probabilities
                    if (output_shape.size() == 3) { // [batch, jets, jets] - always jet pair assignment
                        size_t jets_dim1 = output_shape[1];
                        size_t jets_dim2 = output_shape[2];
                        size_t offset = batch_idx * jets_dim1 * jets_dim2;
                        
                        float max_value = output_data[offset];
                        size_t max_idx = 0;
                        for (size_t j = 1; j < jets_dim1 * jets_dim2; ++j) {
                            if (output_data[offset + j] > max_value) {
                                max_value = output_data[offset + j];
                                max_idx = j;
                            }
                        }
                        
                        int jet1_idx = max_idx / jets_dim2;
                        int jet2_idx = max_idx % jets_dim2;
                        batch_results[batch_idx][tensor_idx] = {max_value, static_cast<float>(jet1_idx), static_cast<float>(jet2_idx)};
                    } else if (output_shape.size() == 2) { // [batch, jets] - single jet assignment
                        size_t jets_dim = output_shape[1];
                        size_t offset = batch_idx * jets_dim;
                        
                        float max_value = output_data[offset];
                        size_t max_idx = 0;
                        for (size_t j = 1; j < jets_dim; ++j) {
                            if (output_data[offset + j] > max_value) {
                                max_value = output_data[offset + j];
                                max_idx = j;
                            }
                        }
                        batch_results[batch_idx][tensor_idx] = {max_value, static_cast<float>(max_idx)};
                    } else {
                        // Fallback for unexpected dimensions
                        batch_results[batch_idx][tensor_idx] = {0.0f, 0.0f};
                    }
                } else if (tensor_idx < 10) { // Detection probabilities
                    if (output_shape.size() == 1) { // [batch]
                        batch_results[batch_idx][tensor_idx] = {output_data[batch_idx]};
                    } else if (output_shape.size() == 2) { // [batch, 1] 
                        size_t offset = batch_idx * output_shape[1];
                        batch_results[batch_idx][tensor_idx] = {output_data[offset]};
                    } else {
                        batch_results[batch_idx][tensor_idx] = {0.0f};
                    }
                } else { // Event-level output
                    if (output_shape.size() == 1) { // [batch]
                        batch_results[batch_idx][tensor_idx] = {output_data[batch_idx]};
                    } else if (output_shape.size() == 2) { // [batch, 2] - binary classification
                        size_t offset = batch_idx * output_shape[1];
                        if (output_shape[1] >= 2) {
                            batch_results[batch_idx][tensor_idx] = {output_data[offset + 1]}; // Take signal probability
                        } else {
                            batch_results[batch_idx][tensor_idx] = {output_data[offset]};
                        }
                    } else {
                        batch_results[batch_idx][tensor_idx] = {0.0f};
                    }
                }
            }
        }
        
        // Add batch results to overall results
        all_results.insert(all_results.end(), batch_results.begin(), batch_results.end());
    }
    
    return all_results;
}

inline void SPANet::SPANetInference::fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size) {
    // Clear buffers
    std::fill(ak4_flat_jets.begin(), ak4_flat_jets.begin() + actual_batch_size * 10 * 8, 0.0f);
    std::fill(ak8_flat_jets.begin(), ak8_flat_jets.begin() + actual_batch_size * 3 * 8, 0.0f);
    std::fill(ak4_mask_char.begin(), ak4_mask_char.begin() + actual_batch_size * 10, 0);
    std::fill(ak8_mask_char.begin(), ak8_mask_char.begin() + actual_batch_size * 3, 0);
    
    for (size_t batch_idx = 0; batch_idx < actual_batch_size; ++batch_idx) {
        const EventData& event = events[batch_idx];
        
        // Fill AK4 jets - only if vectors are not empty
        if (!event.ak4_pt.empty()) {
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
        }
        
        // Fill AK8 jets - only if vectors are not empty
        if (!event.ak8_pt.empty()) {
            const size_t max_ak8 = std::min<size_t>(3, event.ak8_pt.size());
            for (size_t i = 0; i < max_ak8; ++i) {
                ak8_mask_char[batch_idx * 3 + i] = 1;
                
                const size_t base_idx = batch_idx * 3 * 8 + i * 8;
                ak8_flat_jets[base_idx]     = std::log(event.ak8_mass[i] + 1.0f);
                ak8_flat_jets[base_idx + 1] = std::log(event.ak8_pt[i] + 1.0f);
                ak8_flat_jets[base_idx + 2] = event.ak8_eta[i];
                ak8_flat_jets[base_idx + 3] = std::sin(event.ak8_phi[i]);
                ak8_flat_jets[base_idx + 4] = std::cos(event.ak8_phi[i]);
                ak8_flat_jets[base_idx + 5] = static_cast<float>(event.ak8_nConstituents[i]);
                ak8_flat_jets[base_idx + 6] = event.ak8_HbbScore[i];
                ak8_flat_jets[base_idx + 7] = event.ak8_WqqScore[i];
            }
        }
        
        // Fill event inputs
        const size_t event_base_idx = batch_idx * 3;
        event_inputs[event_base_idx]     = std::log(event.met_pt + 1.0f);
        event_inputs[event_base_idx + 1] = std::sin(event.met_phi);
        event_inputs[event_base_idx + 2] = std::cos(event.met_phi);
        
        event_mask_char[batch_idx] = 1;
    }
}

inline RNode SPANet::SPANetInference::addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<float>>>& all_outputs) {
    // Create separate vectors for each output type to match expected format
    std::vector<std::vector<float>> h_assignment_prob, v_assignment_prob, bh_assignment_prob, bv_assignment_prob, vbs_assignment_prob;
    std::vector<float> h_detection_prob, v_detection_prob, bh_detection_prob, bv_detection_prob, vbs_detection_prob, event_signal_prob;
    
    // Reserve space for efficiency
    size_t n_events = all_outputs.size();
    h_assignment_prob.reserve(n_events);
    v_assignment_prob.reserve(n_events);
    bh_assignment_prob.reserve(n_events);
    bv_assignment_prob.reserve(n_events);
    vbs_assignment_prob.reserve(n_events);
    h_detection_prob.reserve(n_events);
    v_detection_prob.reserve(n_events);
    bh_detection_prob.reserve(n_events);
    bv_detection_prob.reserve(n_events);
    vbs_detection_prob.reserve(n_events);
    event_signal_prob.reserve(n_events);
    
    // Extract outputs for each event
    for (const auto& event_outputs : all_outputs) {
        if (event_outputs.size() >= 11) {
            // Assignment probabilities (0-4): [max_prob, jet1_idx, jet2_idx] or [max_prob, jet_idx]
            h_assignment_prob.push_back(event_outputs[0]);
            v_assignment_prob.push_back(event_outputs[1]);
            bh_assignment_prob.push_back(event_outputs[2]);
            bv_assignment_prob.push_back(event_outputs[3]);
            vbs_assignment_prob.push_back(event_outputs[4]);
            
            // Detection probabilities (5-9): single values
            h_detection_prob.push_back(event_outputs[5].empty() ? 0.0f : event_outputs[5][0]);
            v_detection_prob.push_back(event_outputs[6].empty() ? 0.0f : event_outputs[6][0]);
            bh_detection_prob.push_back(event_outputs[7].empty() ? 0.0f : event_outputs[7][0]);
            bv_detection_prob.push_back(event_outputs[8].empty() ? 0.0f : event_outputs[8][0]);
            vbs_detection_prob.push_back(event_outputs[9].empty() ? 0.0f : event_outputs[9][0]);
            
            // Event signal probability (10)
            event_signal_prob.push_back(event_outputs[10].empty() ? 0.0f : event_outputs[10][0]);
        } else {
            // Default values for incomplete outputs
            h_assignment_prob.push_back({0.0f, 0.0f, 0.0f});
            v_assignment_prob.push_back({0.0f, 0.0f, 0.0f});
            bh_assignment_prob.push_back({0.0f, 0.0f});
            bv_assignment_prob.push_back({0.0f, 0.0f});
            vbs_assignment_prob.push_back({0.0f, 0.0f, 0.0f});
            h_detection_prob.push_back(0.0f);
            v_detection_prob.push_back(0.0f);
            bh_detection_prob.push_back(0.0f);
            bv_detection_prob.push_back(0.0f);
            vbs_detection_prob.push_back(0.0f);
            event_signal_prob.push_back(0.0f);
        }
    }
    
    // Create lambdas that use slot and entry for proper mapping
    auto getHAssignment = [h_assignment_prob](unsigned int slot, ULong64_t entry) {
        return entry < h_assignment_prob.size() ? h_assignment_prob[entry] : std::vector<float>{0.0f, 0.0f, 0.0f};
    };
    
    auto getVAssignment = [v_assignment_prob](unsigned int slot, ULong64_t entry) {
        return entry < v_assignment_prob.size() ? v_assignment_prob[entry] : std::vector<float>{0.0f, 0.0f, 0.0f};
    };
    
    auto getBHAssignment = [bh_assignment_prob](unsigned int slot, ULong64_t entry) {
        return entry < bh_assignment_prob.size() ? bh_assignment_prob[entry] : std::vector<float>{0.0f, 0.0f};
    };
    
    auto getBVAssignment = [bv_assignment_prob](unsigned int slot, ULong64_t entry) {
        return entry < bv_assignment_prob.size() ? bv_assignment_prob[entry] : std::vector<float>{0.0f, 0.0f};
    };
    
    auto getVBSAssignment = [vbs_assignment_prob](unsigned int slot, ULong64_t entry) {
        return entry < vbs_assignment_prob.size() ? vbs_assignment_prob[entry] : std::vector<float>{0.0f, 0.0f, 0.0f};
    };
    
    auto getHDetection = [h_detection_prob](unsigned int slot, ULong64_t entry) {
        return entry < h_detection_prob.size() ? h_detection_prob[entry] : 0.0f;
    };
    
    auto getVDetection = [v_detection_prob](unsigned int slot, ULong64_t entry) {
        return entry < v_detection_prob.size() ? v_detection_prob[entry] : 0.0f;
    };
    
    auto getBHDetection = [bh_detection_prob](unsigned int slot, ULong64_t entry) {
        return entry < bh_detection_prob.size() ? bh_detection_prob[entry] : 0.0f;
    };
    
    auto getBVDetection = [bv_detection_prob](unsigned int slot, ULong64_t entry) {
        return entry < bv_detection_prob.size() ? bv_detection_prob[entry] : 0.0f;
    };
    
    auto getVBSDetection = [vbs_detection_prob](unsigned int slot, ULong64_t entry) {
        return entry < vbs_detection_prob.size() ? vbs_detection_prob[entry] : 0.0f;
    };
    
    auto getEventSignal = [event_signal_prob](unsigned int slot, ULong64_t entry) {
        return entry < event_signal_prob.size() ? event_signal_prob[entry] : 0.0f;
    };
    
    // Add all SPANet outputs using DefineSlotEntry for proper entry mapping
    return df.DefineSlotEntry("spanet_h_assignment", getHAssignment, {})
             .DefineSlotEntry("spanet_v_assignment", getVAssignment, {})
             .DefineSlotEntry("spanet_bh_assignment", getBHAssignment, {})
             .DefineSlotEntry("spanet_bv_assignment", getBVAssignment, {})
             .DefineSlotEntry("spanet_vbs_assignment", getVBSAssignment, {})
             .DefineSlotEntry("spanet_h_detection", getHDetection, {})
             .DefineSlotEntry("spanet_v_detection", getVDetection, {})
             .DefineSlotEntry("spanet_bh_detection", getBHDetection, {})
             .DefineSlotEntry("spanet_bv_detection", getBVDetection, {})
             .DefineSlotEntry("spanet_vbs_detection", getVBSDetection, {})
             .DefineSlotEntry("spanet_event_signal", getEventSignal, {});
}

inline RNode SPANet::SPANetInference::RunSPANetInference(RNode df) {
    auto events = extractEventsFromDataFrame(df);
    
    std::cout << "Running batch inference on " << events.size() << " events..." << std::endl;
    auto all_outputs = runBatchInference(events);
    
    return addSPANetOutputsToDataFrame(df, all_outputs);
}
#endif
