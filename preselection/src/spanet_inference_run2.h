#pragma once

#ifndef SPANETINFERENCERUN2_H
#define SPANETINFERENCERUN2_H

#include "spanet_inference_base.h"

namespace SPANet {

    class SPANetInferenceRun2 : public SPANetInferenceBase {
    public:
        static const std::vector<const char*> INPUT_NAMES;

        static constexpr size_t MAX_AK4_JETS = 10;
        static constexpr size_t MAX_AK8_JETS = 3;
        static constexpr size_t AK4_FEATURES = 8;
        static constexpr size_t AK8_FEATURES = 8;
        static constexpr size_t MET_FEATURES = 3;
        static constexpr size_t MAX_LEPTONS = 0;
        static constexpr size_t LEPTON_FEATURES = 0;

        struct EventData : public EventDataBase
        {
          std::vector<float> ak4_pt, ak4_eta, ak4_phi, ak4_mass;
          std::vector<int> ak4_isTightBTag, ak4_isMediumBTag, ak4_isLooseBTag;
          std::vector<float> ak8_pt, ak8_eta, ak8_phi, ak8_mass;
          std::vector<float> ak8_HbbScore, ak8_WqqScore;
          std::vector<unsigned char> ak8_nConstituents;
          float met_pt, ht_ak4, ht_ak8;
        };

        SPANetInferenceRun2(const std::string &model_path, size_t batch_size = 64);

        RNode ParseSpanetInference(RNode df_);

    private:
        std::vector<std::shared_ptr<EventDataBase>> extractEventsFromDataFrame(RNode df) override;
        void fillBatchTensors(const std::vector<std::shared_ptr<EventDataBase>>& events_base, size_t actual_batch_size) override;

        static std::vector<int> assign_all_objects_maxprob(
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