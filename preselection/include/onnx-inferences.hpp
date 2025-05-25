#pragma once
#include <onnxruntime/onnxruntime_cxx_api.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <random>

std::vector<std::vector<float>> runSPANetInference(Ort::Session& session,
    std::vector<float> ak4_jet_1, std::vector<float> ak4_jet_2, std::vector<float> ak4_jet_3,
    std::vector<float> ak4_jet_4, std::vector<float> ak4_jet_5, std::vector<float> ak4_jet_6,
    std::vector<float> ak4_jet_7, std::vector<float> ak4_jet_8, std::vector<float> ak4_jet_9,
    std::vector<float> ak4_jet_10, std::vector<bool> ak4_mask,
    std::vector<float> ak8_jet_1, std::vector<float> ak8_jet_2, std::vector<float> ak8_jet_3,
    std::vector<bool> ak8_mask,
    float event_input_1, float event_input_2, float event_input_3, bool event_mask)
{

    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    
    // Combine AK4 jets into a flattened vector
    std::vector<std::vector<float>> ak4_jets = {
        ak4_jet_1, ak4_jet_2, ak4_jet_3, ak4_jet_4, ak4_jet_5,
        ak4_jet_6, ak4_jet_7, ak4_jet_8, ak4_jet_9, ak4_jet_10};
    std::vector<float> ak4_flat_jets;
    for (const auto &jet : ak4_jets)
    {
        ak4_flat_jets.insert(ak4_flat_jets.end(), jet.begin(), jet.end());
    }

    // Combine AK8 jets into a flattened vector
    std::vector<std::vector<float>> ak8_jets = {ak8_jet_1, ak8_jet_2, ak8_jet_3};
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

    // Get model input/output names # FIXME: repeated in commonSelections.cpp. Take as input?
    std::vector<const char *> input_names = {
        "AK4Jets_data", "AK4Jets_mask",
        "AK8Jets_data", "AK8Jets_mask",
        "HadronicActivity_data", "HadronicActivity_mask"};

    std::vector<const char *> output_names = {
        "h_assignment_probability",
        "v1_assignment_probability",
        "v2_assignment_probability",
        "bh_assignment_probability",
        "bv1_assignment_probability",
        "bv2_assignment_probability",
        "vbs_assignment_probability",
        "h_detection_probability",
        "v1_detection_probability",
        "v2_detection_probability",
        "bh_detection_probability",
        "bv1_detection_probability",
        "bv2_detection_probability",
        "vbs_detection_probability",
        "EVENT/isSignal"};

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
            outputs.push_back(output_vector);
        }
    }
    
    return outputs;
}
