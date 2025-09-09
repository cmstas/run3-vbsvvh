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
using ROOT::VecOps::RVec;

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
        RNode ParseSpanetInference(RNode df_);

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
            ULong64_t rdf_entry;
            unsigned int rdf_slot;
        };
        
        std::vector<std::vector<std::vector<std::vector<float>>>> runBatchInference(const std::vector<EventData>& events);
        void fillBatchTensors(const std::vector<EventData>& events, size_t actual_batch_size);
        
        std::vector<EventData> extractEventsFromDataFrame(RNode df);
        
        RNode addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs, const std::vector<EventData>& events);

        static std::vector<int> assign_all_objects(
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
        );
    };
}

#endif