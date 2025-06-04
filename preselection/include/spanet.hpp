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
        SPANetInference(const std::string &model_path) 
            : session(nullptr) { // Initialize with nullptr first
            Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "SPANet");
            session = Ort::Session(env, model_path.c_str(), Ort::SessionOptions{});
        }
        RNode RunSPANetInference(RNode df);

    private:
        Ort::Session session;
        std::vector<std::vector<float>> runSPANetInference(Ort::Session& session,
            std::vector<float> ak4_pt, std::vector<float> ak4_eta, std::vector<float> ak4_phi, std::vector<float> ak4_mass,
            std::vector<int> ak4_isTightBTag, std::vector<int> ak4_isMediumBTag, std::vector<int> ak4_isLooseBTag,
            std::vector<float> ak8_pt, std::vector<float> ak8_eta, std::vector<float> ak8_phi, std::vector<float> ak8_msoftdrop,
            std::vector<float> ak8_HbbScore, std::vector<float> ak8_WqqScore, std::vector<unsigned char> ak8_nConstituents,
            float met_pt, float met_phi);
    };
}

inline std::vector<std::vector<float>> SPANet::SPANetInference::runSPANetInference(
    Ort::Session& session,
    std::vector<float> ak4_pt, std::vector<float> ak4_eta, std::vector<float> ak4_phi, std::vector<float> ak4_mass,
    std::vector<int> ak4_isTightBTag, std::vector<int> ak4_isMediumBTag, std::vector<int> ak4_isLooseBTag,
    std::vector<float> ak8_pt, std::vector<float> ak8_eta, std::vector<float> ak8_phi, std::vector<float> ak8_msoftdrop,
    std::vector<float> ak8_HbbScore, std::vector<float> ak8_WqqScore, std::vector<unsigned char> ak8_nConstituents,
    float met_pt, float met_phi)
{
    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    
    // Prepare AK4 jet inputs (up to 10 jets, 8 features each)
    std::vector<std::vector<float>> ak4_jets(10, std::vector<float>(8, 0.0f));
    std::vector<bool> ak4_mask(10, false);
    
    for (size_t i = 0; i < std::min<size_t>(10, ak4_pt.size()); ++i)
    {
        ak4_mask[i] = true;
        ak4_jets[i] = {
            std::log(ak4_mass[i] + 1.0f),            // Feature 5: log(mass+1)
            std::log(ak4_pt[i] + 1.0f),              // Feature 1: log(pT+1)
            ak4_eta[i],                              // Feature 2: eta
            std::sin(ak4_phi[i]),                    // Feature 4: sin(phi)
            std::cos(ak4_phi[i]),                    // Feature 3: cos(phi)
            static_cast<float>(ak4_isTightBTag[i]),  // Feature 6: tight b-tag
            static_cast<float>(ak4_isMediumBTag[i]), // Feature 7: medium b-tag
            static_cast<float>(ak4_isLooseBTag[i])   // Feature 8: loose b-tag
        };
    }

    // Prepare AK8 jet inputs (exactly 3 jets, 8 features each)
    std::vector<std::vector<float>> ak8_jets(3, std::vector<float>(8, 0.0f));
    std::vector<bool> ak8_mask(3, false);
    
    for (size_t i = 0; i < std::min<size_t>(3, ak8_pt.size()); ++i)
    {
        ak8_mask[i] = true;
        ak8_jets[i] = {
            std::log(ak8_msoftdrop[i] +  1.0f),       // Feature 5: log(msoftdrop+1)
            std::log(ak8_pt[i] +  1.0f),              // Feature 1: log(pT+1)
            ak8_eta[i],                               // Feature 2: eta
            std::sin(ak8_phi[i]),                     // Feature 4: sin(phi)
            std::cos(ak8_phi[i]),                     // Feature 3: cos(phi)
            static_cast<float>(ak8_nConstituents[i]), // Feature 8: nConstituents
            ak8_HbbScore[i],                          // Feature 6: Hbb score
            ak8_WqqScore[i]                           // Feature 7: Wqq score
        };
    }

    // Event-level inputs
    float event_input_1 = std::log(met_pt + 1.0f); // log(MET)
    float event_input_2 = std::sin(met_phi);
    float event_input_3 = std::cos(met_phi);
    bool event_mask = true;                        // Event mask is always True

    // Combine AK4 jets into a flattened vector
    std::vector<float> ak4_flat_jets;
    for (const auto &jet : ak4_jets)
    {
        ak4_flat_jets.insert(ak4_flat_jets.end(), jet.begin(), jet.end());
    }

    // Combine AK8 jets into a flattened vector
    std::vector<float> ak8_flat_jets;
    for (const auto &jet : ak8_jets)
    {
        ak8_flat_jets.insert(ak8_flat_jets.end(), jet.begin(), jet.end());
    }

    std::vector<float> event_inputs = {event_input_1, event_input_2, event_input_3};
    std::vector<bool> event_mask_bool = {event_mask};

    // Prepare input tensor shapes
    std::vector<int64_t> ak4_input_shape = {1, 10, 8};  // [batch, num_jets, features]
    std::vector<int64_t> ak8_input_shape = {1, 3, 8};   // [batch, num_jets, features]
    std::vector<int64_t> event_input_shape = {1, 1, 3}; // [batch, 1, features]

    // Create AK4 input tensor
    Ort::Value ak4_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        ak4_flat_jets.data(),
        ak4_flat_jets.size(),
        ak4_input_shape.data(),
        ak4_input_shape.size());

    // Create AK8 input tensor
    Ort::Value ak8_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        ak8_flat_jets.data(),
        ak8_flat_jets.size(),
        ak8_input_shape.data(),
        ak8_input_shape.size());

    // Create event input tensor
    Ort::Value event_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        event_inputs.data(),
        event_inputs.size(),
        event_input_shape.data(),
        event_input_shape.size());

    // Prepare mask tensors
    std::vector<int64_t> ak4_mask_shape = {1, 10};
    std::vector<int64_t> ak8_mask_shape = {1, 3};
    std::vector<int64_t> event_mask_shape = {1, 1};

    // Convert bool masks to char for tensor creation (workaround)
    std::vector<char> ak4_mask_char(ak4_mask.begin(), ak4_mask.end());
    std::vector<char> ak8_mask_char(ak8_mask.begin(), ak8_mask.end());
    std::vector<char> event_mask_char(event_mask_bool.begin(), event_mask_bool.end());

    Ort::Value ak4_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool *>(ak4_mask_char.data()),
        ak4_mask_char.size(),
        ak4_mask_shape.data(),
        ak4_mask_shape.size());

    Ort::Value ak8_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool *>(ak8_mask_char.data()),
        ak8_mask_char.size(),
        ak8_mask_shape.data(),
        ak8_mask_shape.size());

    Ort::Value event_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool *>(event_mask_char.data()),
        event_mask_char.size(),
        event_mask_shape.data(),
        event_mask_shape.size());

    // Prepare input tensors vector
    std::vector<Ort::Value> input_tensors;
    input_tensors.push_back(std::move(ak4_input_tensor));
    input_tensors.push_back(std::move(ak4_mask_tensor));
    input_tensors.push_back(std::move(ak8_input_tensor));
    input_tensors.push_back(std::move(ak8_mask_tensor));
    input_tensors.push_back(std::move(event_input_tensor));
    input_tensors.push_back(std::move(event_mask_tensor));

    // Get model input/output names
    std::vector<const char *> input_names = {
        "AK4Jets_data", "AK4Jets_mask", "AK8Jets_data", "AK8Jets_mask", "MET_data", "MET_mask"
    };

    std::vector<const char *> output_names = {
        "h_assignment_probability", 
        "v_assignment_probability", 
        "bh_assignment_probability", 
        "bv_assignment_probability",
        "vbs_assignment_probability",
        "h_detection_probability",
        "v_detection_probability",
        "bh_detection_probability",
        "bv_detection_probability",
        "vbs_detection_probability",
        "EVENT/isSignal"
    };

    // Run inference
    auto output_tensors = session.Run(
        Ort::RunOptions{nullptr},
        input_names.data(),
        input_tensors.data(),
        input_tensors.size(),
        output_names.data(),
        output_names.size());

    // Process outputs
    std::vector<std::vector<float>> outputs;
    for (size_t i = 0; i < output_tensors.size(); ++i)
    {
        Ort::Value &output_tensor = output_tensors[i];
        const Ort::TensorTypeAndShapeInfo &shape_info = output_tensor.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> output_shape = shape_info.GetShape();
        size_t num_elements = shape_info.GetElementCount();

        if (output_tensor.IsTensor())
        {
            const float *output_data = output_tensor.GetTensorData<float>();
            std::vector<float> output_vector(output_data, output_data + num_elements);
            // get argmax for the vector and resolve the index in the original two dimensional tensor
            // assignment probability outputs are (1, num_jets, num_jets) if object is composed of two jets, (1, num_jets) if one jet
            if (i < 5) { // Assignment probabilities
               if (output_shape.size() == 3) {
                    auto max_iter = std::max_element(output_vector.begin(), output_vector.end());
                    float max_value = *max_iter;
                    int max_idx = std::distance(output_vector.begin(), max_iter);
                    int num_jets = output_shape[1];
                    int jet1_idx = max_idx / num_jets;
                    int jet2_idx = max_idx % num_jets;
                    outputs.push_back({max_value, static_cast<float>(jet1_idx), static_cast<float>(jet2_idx)});
                } else if (output_shape.size() == 2) {
                    auto max_iter = std::max_element(output_vector.begin(), output_vector.end());
                    float max_value = *max_iter;
                    int max_idx = std::distance(output_vector.begin(), max_iter);
                    outputs.push_back({max_value, static_cast<float>(max_idx)});
                }
            }
            else if (i < 10) { // Detection probabilities
                outputs.push_back({output_vector[0]});
            } else { // Event-level output
                outputs.push_back({output_vector[1]});
            }

        }
    }
    
    return outputs;
}

inline RNode SPANet::SPANetInference::RunSPANetInference(RNode df)
{
    using ROOT::VecOps::RVec;
    using RVecF = RVec<float>;
    using RVecI = RVec<int>;
    using RVecUC = RVec<unsigned char>;
    
    auto runInference = [this](const RVecF &ak4_pt, const RVecF &ak4_eta, const RVecF &ak4_phi, const RVecF &ak4_mass,
                               const RVecI &ak4_isTightBTag, const RVecI &ak4_isMediumBTag, const RVecI &ak4_isLooseBTag,
                               const RVecF &ak8_pt, const RVecF &ak8_eta, const RVecF &ak8_phi, const RVecF &ak8_msoftdrop,
                               const RVecF &ak8_HbbScore, const RVecF &ak8_WqqScore, const RVecUC &ak8_nConstituents,
                               float met_pt, float met_phi)
    {
        // Convert RVec to std::vector
        std::vector<float> std_ak4_pt(ak4_pt.begin(), ak4_pt.end());
        std::vector<float> std_ak4_eta(ak4_eta.begin(), ak4_eta.end());
        std::vector<float> std_ak4_phi(ak4_phi.begin(), ak4_phi.end());
        std::vector<float> std_ak4_mass(ak4_mass.begin(), ak4_mass.end());
        std::vector<int> std_ak4_isTightBTag(ak4_isTightBTag.begin(), ak4_isTightBTag.end());
        std::vector<int> std_ak4_isMediumBTag(ak4_isMediumBTag.begin(), ak4_isMediumBTag.end());
        std::vector<int> std_ak4_isLooseBTag(ak4_isLooseBTag.begin(), ak4_isLooseBTag.end());
        std::vector<float> std_ak8_pt(ak8_pt.begin(), ak8_pt.end());
        std::vector<float> std_ak8_eta(ak8_eta.begin(), ak8_eta.end());
        std::vector<float> std_ak8_phi(ak8_phi.begin(), ak8_phi.end());
        std::vector<float> std_ak8_msoftdrop(ak8_msoftdrop.begin(), ak8_msoftdrop.end());
        std::vector<float> std_ak8_HbbScore(ak8_HbbScore.begin(), ak8_HbbScore.end());
        std::vector<float> std_ak8_WqqScore(ak8_WqqScore.begin(), ak8_WqqScore.end());
        std::vector<unsigned char> std_ak8_nConstituents(ak8_nConstituents.begin(), ak8_nConstituents.end());
        // std::vector<unsigned char> std_ak8_nConstituents;
        // std_ak8_nConstituents.reserve(ak8_nConstituents.size());
        // for (auto val : ak8_nConstituents) {
        //     std_ak8_nConstituents.push_back(static_cast<unsigned char>(val));
        // }
        

        // Run SPANet inference
        auto outputs = runSPANetInference(
            this->session,
            std_ak4_pt, std_ak4_eta, std_ak4_phi, std_ak4_mass,
            std_ak4_isTightBTag, std_ak4_isMediumBTag, std_ak4_isLooseBTag,
            std_ak8_pt, std_ak8_eta, std_ak8_phi, std_ak8_msoftdrop,
            std_ak8_HbbScore, std_ak8_WqqScore, std_ak8_nConstituents,
            met_pt, met_phi);

        return outputs;

    };
 
    const std::vector<std::string> input_columns = {
        "ak4jet_pt", "ak4jet_eta", "ak4jet_phi", "ak4jet_mass",
        "ak4jet_isTightBTag", "ak4jet_isMediumBTag", "ak4jet_isLooseBTag",
        "ak8jet_pt", "ak8jet_eta", "ak8jet_phi", "ak8jet_msoftdrop",
        "ak8jet_xbbvsqcd", "ak8jet_xqqvsqcd", "ak8jet_nConstituents",
        "MET_pt", "MET_phi"};

    return df.Define("spanet_outputs", runInference, input_columns);
}


#endif
