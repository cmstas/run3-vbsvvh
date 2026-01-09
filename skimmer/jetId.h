#include <ROOT/RVec.hxx>
#include "correction.h"
#include <string>
using namespace ROOT::VecOps;

#define JET_ID_JSON_2024 "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2024_Summer24/jetid.json.gz"

RVec<float> evalJetID(const std::string& jet_id_json, const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                        const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                        const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                        const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    auto cset_jetId = correction::CorrectionSet::from_file(jet_id_json);
    RVec<float> jetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        jetId[i] += 2 * cset_jetId->at("AK4PUPPI_Tight")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
        jetId[i] += 4 * cset_jetId->at("AK4PUPPI_TightLeptonVeto")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
    }
    return jetId;
}

RVec<float> evalFatJetID(const std::string& jet_id_json, const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                            const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                            const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                            const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    auto cset_fatJetId = correction::CorrectionSet::from_file(jet_id_json);
    RVec<float> fatJetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        fatJetId[i] += 2 * cset_fatJetId->at("AK8PUPPI_Tight")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
        fatJetId[i] += 4 * cset_fatJetId->at("AK8PUPPI_TightLeptonVeto")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
    }
    return fatJetId;
}

RVec<float> evalJetID2024(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                        const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                        const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                        const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    return evalJetID(JET_ID_JSON_2024, eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity);
}

RVec<float> evalFatJetID2024(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                            const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                            const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                            const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    return evalFatJetID(JET_ID_JSON_2024, eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity);
}