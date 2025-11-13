#pragma once

#ifndef SPANETINFERENCEBASE_H
#define SPANETINFERENCEBASE_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <memory>

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

namespace SPANet {

    class EventDataBase
    {
    public:
        virtual ~EventDataBase() = default;
        ULong64_t rdf_entry;
        unsigned int rdf_slot;
    };

    class SPANetInferenceBase {
    public:
        SPANetInferenceBase(
            const std::string &model_path,
            size_t batch_size,
             std::vector<const char*> input_names,
            size_t max_ak4_jets,
            size_t ak4_features,
            size_t max_ak8_jets,
            size_t ak8_features,
            size_t met_features,
            size_t max_leptons,
            size_t lepton_features
        );

        virtual ~SPANetInferenceBase() = default;

        // Functions that depend on input features are defined as pure virtual and have to be defined in derived class
        RNode RunSPANetInference(RNode df);
        virtual RNode ParseSpanetInference(RNode df_) = 0;
        

    protected:
        std::unique_ptr<Ort::Env> env_;
        std::unique_ptr<Ort::Session> session_;
        size_t batch_size_;
        int64_t top_k_;
        Ort::MemoryInfo memory_info_;

        std::vector<const char*> input_names_;
        std::vector<const char*> output_names_;

        std::vector<int64_t> ak4_input_shape_;
        std::vector<int64_t> ak8_input_shape_;
        std::vector<int64_t> met_input_shape_;
        std::vector<int64_t> lepton1_input_shape_;
        std::vector<int64_t> lepton2_input_shape_;
        std::vector<int64_t> ak4_mask_shape_;
        std::vector<int64_t> ak8_mask_shape_;
        std::vector<int64_t> met_mask_shape_;
        std::vector<int64_t> lepton1_mask_shape_;
        std::vector<int64_t> lepton2_mask_shape_;

        std::vector<float> ak4_flat_jets_;
        std::vector<float> ak8_flat_jets_;
        std::vector<float> met_inputs_;
        std::vector<float> lepton1_inputs_;
        std::vector<float> lepton2_inputs_;
        std::vector<char> ak4_mask_char_;
        std::vector<char> ak8_mask_char_;
        std::vector<char> met_mask_char_;
        std::vector<char> lepton1_mask_char_;
        std::vector<char> lepton2_mask_char_;

        // Functions that depend on input features are defined as pure virtual and have to be defined in derived class
        virtual std::vector<std::shared_ptr<EventDataBase>> extractEventsFromDataFrame(RNode df) = 0;
        virtual void fillBatchTensors(const std::vector<std::shared_ptr<EventDataBase>>& events, size_t actual_batch_size) = 0;
        
        
        // Functions that should be identical for all derived classes are defined in the base class
        std::vector<std::vector<std::vector<std::vector<float>>>> runBatchInference(const std::vector<std::shared_ptr<EventDataBase>>& events);
        RNode addSPANetOutputsToDataFrame(RNode df, const std::vector<std::vector<std::vector<std::vector<float>>>>& all_outputs, const std::vector<std::shared_ptr<EventDataBase>>& events);

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