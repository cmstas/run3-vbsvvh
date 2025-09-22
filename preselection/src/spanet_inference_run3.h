#pragma once

#ifndef SPANETINFERENCERUN3_H
#define SPANETINFERENCERUN3_H

#include "spanet_inference_base.h"

namespace SPANet {

    class SPANetInferenceRun3 : public SPANetInferenceBase {
    public:
        static constexpr size_t MAX_AK4_JETS = 10;
        static constexpr size_t MAX_AK8_JETS = 3;
        static constexpr size_t AK4_FEATURES = 8;
        static constexpr size_t AK8_FEATURES = 11;
        static constexpr size_t MET_FEATURES = 3;
        static constexpr size_t MAX_LEPTONS = 2;
        static constexpr size_t LEPTON_FEATURES = 6;

        struct EventData : public EventDataBase
        {
          std::vector<float> ak4_pt, ak4_eta, ak4_phi, ak4_mass;
          std::vector<int> ak4_isTightBTag, ak4_isMediumBTag, ak4_isLooseBTag;
          std::vector<float> ak8_pt, ak8_eta, ak8_phi, ak8_mass;
          std::vector<float> ak8_XbbScore, ak8_XqqScore, ak8_XccScore, ak8_XcsScore, ak8_XqcdScore;
          std::vector<unsigned char> ak8_nConstituents;
          float met_pt, met_phi;
          float lepton1_pt, lepton1_eta, lepton1_phi, lepton1_mass, lepton1_charge;
          float lepton2_pt, lepton2_eta, lepton2_phi, lepton2_mass, lepton2_charge;
        };

        SPANetInferenceRun3(const std::string &model_path, size_t batch_size = 64);

        RNode ParseSpanetInference(RNode df_);

    private:
        std::vector<std::shared_ptr<EventDataBase>> extractEventsFromDataFrame(RNode df) override;
        void fillBatchTensors(const std::vector<std::shared_ptr<EventDataBase>>& events_base, size_t actual_batch_size) override;
    };

}

#endif