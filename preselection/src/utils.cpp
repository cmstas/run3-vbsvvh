#include "utils.h"


/*
############################################
DEFINE METADATA
############################################
*/

RNode defineMetadata(RNode df){
    return df.DefinePerSample("xsec", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("xsec");})
        .DefinePerSample("lumi", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("lumi");})
        .DefinePerSample("nevents", [](unsigned int slot, const RSampleInfo &id) { return id.GetD("nevents");})
        .DefinePerSample("sample_category", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_category");})
        .DefinePerSample("sample_type", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_type");})
        .DefinePerSample("sample_year", [](unsigned int slot, const RSampleInfo &id) { return id.GetS("sample_year");})
        .Define("xsec_weight", "1000 * xsec * lumi / nevents")
        .Define("isData", "sample_category == \"data\"");
}

/*
############################################
LUMIMASK
############################################
*/

bool operator< ( const lumiMask::LumiBlockRange& lh, const lumiMask::LumiBlockRange& rh )
{
    return ( lh.run() == rh.run() ) ? ( lh.lastLumi() < rh.firstLumi() ) : lh.run() < rh.run();
}

lumiMask lumiMask::fromJSON(const std::string& file, lumiMask::Run firstRun, lumiMask::Run lastRun)
{
  const bool noRunFilter = ( firstRun == 0 ) && ( lastRun == 0 );
  boost::property_tree::ptree ptree;
  boost::property_tree::read_json(file, ptree);

  std::vector<lumiMask::LumiBlockRange> accept;
  for ( const auto& runEntry : ptree ) {
    const lumiMask::Run run = std::stoul(runEntry.first);
    if ( noRunFilter || ( ( firstRun <= run ) && ( run <= lastRun ) ) ) {
      for ( const auto& lrEntry : runEntry.second ) {
        const auto lrNd = lrEntry.second;
        const lumiMask::LumiBlock firstLumi = std::stoul(lrNd.begin()->second.data());
        const lumiMask::LumiBlock lastLumi  = std::stoul((++lrNd.begin())->second.data());
        accept.emplace_back(run, firstLumi, lastLumi);
      }
    }
  }
  return lumiMask(accept);
}

/*
############################################
REMOVE DUPLICATES
############################################
*/

RNode removeDuplicates(RNode df){
    return df.Filter(FilterOnePerKind(), {"run", "luminosityBlock", "event"}, "REMOVED DUPLICATES");
}

/*
############################################
CUTFLOW
############################################
*/

Cutflow::Cutflow(RNode df, const std::vector<std::string>& cuts) : _df(df), _cuts(cuts){
    _df = _df.Define("weight2", "weight * weight");
    _cutflow.push_back(std::make_pair(_df.Sum<double>("weight"), _df.Sum<double>("weight2")));
    for (size_t i = 0; i < _cuts.size(); i++) {
        _cutflow.push_back(std::make_pair(_df.Filter(_cuts[i]).Sum<double>("weight"), _df.Filter(_cuts[i]).Sum<double>("weight2")));
    }
}

void Cutflow::Print(std::string output_file) {
    if (!output_file.empty()) {
        std::ofstream out(output_file);
        out << "Cutflow\n";
        out << "initialCount: " << '[' << _cutflow[0].first.GetValue() << ',' << std::sqrt(_cutflow[0].second.GetValue()) << ']' << "\n";
        for (size_t i = 0; i < _cuts.size(); i++) {
            out << _cuts[i] << ": " << '[' << _cutflow[i+1].first.GetValue() << ',' << std::sqrt(_cutflow[i+1].second.GetValue()) << ']' << "\n";
        }
        out.close();
    }
    else {
        printf("%-15s: %10.2f \\pm %10.2f\n", "initialCount", _cutflow[0].first.GetValue(), std::sqrt(_cutflow[0].second.GetValue()));
        for (size_t i = 0; i < _cuts.size(); i++) {
            printf("%-15s: %10.2f \\pm %10.2f\n", _cuts[i].c_str(), _cutflow[i+1].first.GetValue(), std::sqrt(_cutflow[i+1].second.GetValue()));
        }
    }
}

/*
############################################
SELECTION UTILS
############################################
*/

RVec<float> VfDeltaR (RVec<float> vec_eta, RVec<float> vec_phi, float obj_eta, float obj_phi) { 
    RVec<float> deltaR = {};
    if (obj_eta == -999 || obj_phi == -999) {
        for (size_t i = 0; i < vec_eta.size(); i++) {
            deltaR.push_back(1);
        }
        return deltaR;
    }
    else {
        for (size_t i = 0; i < vec_eta.size(); i++) {
            deltaR.push_back(ROOT::VecOps::DeltaR(vec_eta[i], obj_eta, vec_phi[i], obj_phi));
        }
    }
    return deltaR;
}

RVec<float> VfInvariantMass(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass = {};
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);
    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass.push_back((obj1 + obj2).M());
    }
    return invMass;
}

RVec<float> VfInvariantPt(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass = {};
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);
    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass.push_back((obj1 + obj2).Pt());
    }
    return invMass;
}

RVec<float> VfInvariantPhi(RVec<float> vec_pt, RVec<float> vec_eta, RVec<float> vec_phi, RVec<float> vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass = {};
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);
    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass.push_back((obj1 + obj2).Phi());
    }
    return invMass;
}

RVec<float> MT(RVec<float> vec_pt, RVec<float> vec_phi, float obj_pt, float obj_phi) {
    RVec<float> mt = {};
    for (size_t i = 0; i < vec_pt.size(); i++) {
        mt.push_back(std::sqrt(2 * vec_pt[i] * obj_pt * (1 - std::cos(ROOT::VecOps::DeltaPhi(vec_phi[i], obj_phi)))));
    }
    return mt;
}

float fInvariantMass(float obj1_pt, float obj1_eta, float obj1_phi, float obj1_mass, float obj2_pt, float obj2_eta, float obj2_phi, float obj2_mass) {
    RVec<float> invMass = {};
    TLorentzVector obj1, obj2;
    obj1.SetPtEtaPhiM(obj1_pt, obj1_eta, obj1_phi, obj1_mass);
    obj2.SetPtEtaPhiM(obj2_pt, obj2_eta, obj2_phi, obj2_mass);
    return (obj1 + obj2).M();
}

RVec<int> VBS_MaxE( RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass) {
    RVec<int> good_jet_idx = {};
    int Nvbfjet1 = -1;
    int Nvbfjet2 = -1;
    float maxvbfjetp = 0;
    float secvbfjetp = 0;
    RVec<float> Jet_P = {};
    for (size_t i = 0; i < Jet_pt.size(); i++) {
        TLorentzVector jet;
        jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        Jet_P.push_back(jet.P());
    }
    for (size_t i = 0; i < Jet_eta.size(); i++) {
        if (Jet_eta[i] > 0 && maxvbfjetp < Jet_P[i]){
            maxvbfjetp = Jet_P[i];
            Nvbfjet1 = i;
        } else if(Jet_eta[i] < 0 && secvbfjetp < Jet_P[i]){
            secvbfjetp = Jet_P[i];
            Nvbfjet2 = i;
        }
    }
    if (Nvbfjet1 != -1 && Nvbfjet2 != -1) {
        if (maxvbfjetp > secvbfjetp) {
            good_jet_idx.push_back(Nvbfjet1);
            good_jet_idx.push_back(Nvbfjet2);
        }
        else {
            good_jet_idx.push_back(Nvbfjet2);
            good_jet_idx.push_back(Nvbfjet1);
        }
    }
    else {
        auto sorted_jet_p = Argsort(-Jet_P);
        good_jet_idx.push_back(sorted_jet_p[0]);
        good_jet_idx.push_back(sorted_jet_p[1]);        
    }
    return good_jet_idx;
}

RVec<int> VBS_MaxEtaJJ(RVec<float> Jet_pt, RVec<float> Jet_eta, RVec<float> Jet_phi, RVec<float> Jet_mass) {
    // find pair of jets with max delta eta
    RVec<int> good_jet_idx = {};
    RVec<float> Jet_P = {};
    for (size_t i = 0; i < Jet_pt.size(); i++) {
        TLorentzVector jet;
        jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
        Jet_P.push_back(jet.P());
    }
    int Nvbfjet1 = -1;
    int Nvbfjet2 = -1;
    float maxvbfjetdeta = 0;
    for (size_t i = 0; i < Jet_eta.size(); i++) {
        for (size_t j = i+1; j < Jet_eta.size(); j++) {
            float deta = std::abs(Jet_eta[i] - Jet_eta[j]);
            if (deta > maxvbfjetdeta) {
                maxvbfjetdeta = deta;
                Nvbfjet1 = i;
                Nvbfjet2 = j;
            }
        }
    }
    if (Jet_P[Nvbfjet1] > Jet_P[Nvbfjet2]) {
        good_jet_idx.push_back(Nvbfjet1);
        good_jet_idx.push_back(Nvbfjet2);
    }
    else {
        good_jet_idx.push_back(Nvbfjet2);
        good_jet_idx.push_back(Nvbfjet1);
    }
    return good_jet_idx;
}

// Return for each ak4 jet, the dR from the closest ak8 jet
RVec<float> VfdRfromClosestJet(const ROOT::RVecF &ak4_eta, const ROOT::RVecF &ak4_phi, const ROOT::RVecF &ak8_eta, const ROOT::RVecF &ak8_phi)
{
    RVec<float> vec_minDR = {};
    for (size_t i = 0; i < ak4_eta.size(); i++)
    {
        float mindR = 999.;
        for (size_t j = 0; j < ak8_eta.size(); j++)
        {
            float dR = ROOT::VecOps::DeltaR(ak4_eta.at(i), ak8_eta.at(j), ak4_phi.at(i), ak8_phi.at(j));
            if (dR < mindR) {
                mindR = dR;
            }
        }
        vec_minDR.push_back(mindR);
    }
    return vec_minDR;
}

/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &dir) {
    // If on UAF and the USER environment variable is defined, store output on ceph
    const char* userEnv = getenv("USER");
    std::string ceph_dir = "/ceph/cms/store/user/";
    std::string output_dir = "./";
    if (userEnv != nullptr && std::filesystem::exists(ceph_dir) && std::filesystem::is_directory(ceph_dir)) {
        output_dir = ceph_dir + std::string(userEnv) + "/vbsvvhAnalysis/preselection/";
    }

    std::filesystem::path directory_path(output_dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        std::cerr << "Output directory already exists: " << directory_path << std::endl;
    }
    // Try to create the directory and any missing parent directories
    else if (std::filesystem::create_directories(directory_path)) {
        std::cout << "Created output directory : " << directory_path << std::endl;
    }
    else {
        std::cerr << "Failed to create output directory: " << directory_path << std::endl;
        std::exit(EXIT_FAILURE); 
    }

    return directory_path;
}

void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool isData)
{
    auto ColNames = df.GetDefinedColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");
    
    for (auto &&ColName : ColNames) {
        TString colName = ColName;
        std::string name = colName.Data();
        final_variables.push_back(name);
    }

    if (!isData)
        final_variables.push_back("LHEReweightingWeight");

    df.Snapshot("Events", outputDir + "/" + outputFileName + ".root", final_variables);
}