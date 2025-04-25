// filepath: /home/users/aaarora/phys/run3/onnx.cpp
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <random>

std::vector<std::vector<float>> runSPANetInference(
    std::vector<float> ak4_input_1,
    std::vector<float> ak4_input_2,
    std::vector<float> ak4_input_3,
    std::vector<float> ak4_input_4,
    std::vector<float> ak4_input_5,
    std::vector<float> ak4_input_6,
    std::vector<float> ak4_input_7,
    std::vector<float> ak4_input_8,
    std::vector<bool> ak4_mask,
    std::vector<float> ak8_input_1,
    std::vector<float> ak8_input_2,
    std::vector<float> ak8_input_3,
    std::vector<float> ak8_input_4,
    std::vector<float> ak8_input_5,
    std::vector<float> ak8_input_6,
    std::vector<float> ak8_input_7,
    std::vector<float> ak8_input_8,
    std::vector<bool> ak8_mask,
    float event_input_2,
    float event_input_3,
    float event_input_1,
    bool event_mask
) {
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNXInference");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);

    const char* model_path = "/home/users/mmazza/public/forAashay/spanet_v31_model/spanet.onnx";
    Ort::Session session(env, model_path, session_options);
    Ort::AllocatorWithDefaultOptions allocator;

    std::vector<std::vector<float>> ak4_inputs = {
        ak4_input_1, ak4_input_2, ak4_input_3, ak4_input_4,
        ak4_input_5, ak4_input_6, ak4_input_7, ak4_input_8
    };

    std::vector<std::vector<float>> ak8_inputs = {
        ak8_input_1, ak8_input_2, ak8_input_3, ak8_input_4,
        ak8_input_5, ak8_input_6, ak8_input_7, ak8_input_8
    };

    std::vector<float> event_inputs = {event_input_1, event_input_2, event_input_3};
    std::vector<bool> event_mask_bool = {event_mask};
    
    // Prepare input tensor
    std::vector<int64_t> ak4_input_shape = {1, 10, 8};
    std::vector<int64_t> ak8_input_shape = {1, 3, 8};
    std::vector<int64_t> event_input_shape = {1, 1, 3};


    Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
    Ort::Value ak4_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        ak4_inputs.data(),
        ak4_inputs.size(),
        ak4_input_shape.data(),
        ak4_input_shape.size()
    );
    Ort::Value ak8_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        ak8_inputs.data(),
        ak8_inputs.size(),
        ak8_input_shape.data(),
        ak8_input_shape.size()
    );

    Ort::Value event_input_tensor = Ort::Value::CreateTensor<float>(
        memory_info,
        event_inputs.data(),
        event_inputs.size(),
        event_input_shape.data(),
        event_input_shape.size()
    );


    std::vector<int64_t> ak4_mask_shape = {1, 10};  // Match the first dimension of the data tensor
    std::vector<int64_t> ak8_mask_shape = {1, 3};  // Match the first dimension of the data tensor
    std::vector<int64_t> event_mask_shape = {1, 1};

    std::vector<char> ak4_mask_char(ak4_mask.begin(), ak4_mask.end());
    std::vector<char> ak8_mask_char(ak8_mask.begin(), ak8_mask.end());
    std::vector<char> event_mask_char(event_mask_bool.begin(), event_mask_bool.end());

    Ort::Value ak4_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool*>(ak4_mask_char.data()),
        ak4_mask_char.size(),
        ak4_mask_shape.data(),
        ak4_mask_shape.size()
    );

    Ort::Value ak8_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool*>(ak8_mask_char.data()),
        ak8_mask_char.size(),
        ak8_mask_shape.data(),
        ak8_mask_shape.size()
    );

    Ort::Value event_mask_tensor = Ort::Value::CreateTensor<bool>(
        memory_info,
        reinterpret_cast<bool*>(event_mask_char.data()),
        event_mask_char.size(),
        event_mask_shape.data(),
        event_mask_shape.size()
    );

    // Prepare input tensors vector
    std::vector<Ort::Value> input_tensors;
    input_tensors.push_back(std::move(ak4_input_tensor));
    input_tensors.push_back(std::move(ak4_mask_tensor));
    input_tensors.push_back(std::move(ak8_input_tensor));
    input_tensors.push_back(std::move(ak8_mask_tensor));
    input_tensors.push_back(std::move(event_input_tensor));
    input_tensors.push_back(std::move(event_mask_tensor));

    // Get model input/output names
    size_t num_input_nodes = session.GetInputCount();
    size_t num_output_nodes = session.GetOutputCount();
    
    std::vector<const char*> input_names = {
        "AK4Jets_data", "AK4Jets_mask", 
        "AK8Jets_data", "AK8Jets_mask",
        "HadronicActivity_data", "HadronicActivity_mask"
    };
    
    std::vector<const char*> output_names = {
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
        "EVENT/isSignal"
    };

    auto output_tensors = session.Run(
        Ort::RunOptions{nullptr},
        input_names.data(),
        input_tensors.data(),
        input_tensors.size(),
        output_names.data(),
        output_names.size()
    );


    std::vector<std::vector<float>> outputs;
    for (size_t i = 0; i < output_tensors.size(); ++i) {
        Ort::Value& output_tensor = output_tensors[i];
        const Ort::TensorTypeAndShapeInfo& shape_info = output_tensor.GetTensorTypeAndShapeInfo();
        std::vector<int64_t> output_shape = shape_info.GetShape();
        size_t num_elements = shape_info.GetElementCount();

        if (output_tensor.IsTensor()) {
            const float* output_data = output_tensor.GetTensorData<float>();
            std::vector<float> output_vector(output_data, output_data + num_elements);
            outputs.push_back(output_vector);
        }
    }
    
    return outputs;
}