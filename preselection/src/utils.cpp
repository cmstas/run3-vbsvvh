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

Cutflow::Cutflow(RNode df) : _df(df){
    for (auto colName : _df.GetDefinedColumnNames()) {
        if (colName.starts_with("_cut")) {
            _cuts.push_back(colName);
        }
    }

    if (_cuts.size() != 0) { 
        _df = _df.Define("weight2", "weight * weight");
        _cutflow.push_back(std::make_pair(_df.Sum<double>("weight"), _df.Sum<double>("weight2")));
        std::string cut_string = "";
        for (size_t i = 0; i < _cuts.size(); i++) {
            if (i == 0) {
                cut_string = _cuts[i];
            } else {
                cut_string += " && " + _cuts[i];
            }
            _cutflow.push_back(std::make_pair(_df.Filter(cut_string).Sum<double>("weight"), _df.Filter(cut_string).Sum<double>("weight2")));
        }        
    }
}

void Cutflow::Print(std::string output_file) {
    if (_cuts.size() == 0) {
        std::cout << "No cuts _cut_* defined" << std::endl;
        return;
    }
    tabulate::Table table;
    table.add_row({"Cut", "Count", "Percent Change"});
    {
        std::stringstream ss_count, ss_error;
        ss_count << std::fixed << std::setprecision(3) << _cutflow[0].first.GetValue();
        ss_error << std::fixed << std::setprecision(3) << std::sqrt(_cutflow[0].second.GetValue());
        table.add_row({
            "initialCount",
            ss_count.str() + " +/- " + ss_error.str(),
            std::string("")
        });
    }
    
    for (size_t i = 0; i < _cuts.size(); i++) {
        std::stringstream ss_count, ss_error, ss_percent;
        ss_count << std::fixed << std::setprecision(3) << _cutflow[i + 1].first.GetValue();
        ss_error << std::fixed << std::setprecision(3) << std::sqrt(_cutflow[i + 1].second.GetValue());
        ss_percent << std::fixed << std::setprecision(3) << 100 * (_cutflow[i].first.GetValue() - _cutflow[i + 1].first.GetValue()) / _cutflow[i].first.GetValue();
        
        table.add_row({
            _cuts[i],
            ss_count.str() + " +/- " + ss_error.str(),
            ss_percent.str()
        });
    }

    table.column(1).format().font_align(tabulate::FontAlign::center);
    table.column(2).format().font_align(tabulate::FontAlign::center);

    if (!output_file.empty()) {
        std::ofstream out(output_file);
        out << table;
        out.close();
    }
    else {
        std::cout << table << std::endl;
    }
}

/*
############################################
SELECTION UTILS
############################################
*/

float fdR(float eta1, float phi1, float eta2, float phi2) {
    return ROOT::VecOps::DeltaR(eta1, eta2, phi1, phi2);
}

RVec<float> VdR(const RVec<float>& vec_eta, const RVec<float>& vec_phi, float obj_eta, float obj_phi) {
    RVec<float> out(vec_eta.size());
    if (obj_eta == -999 || obj_phi == -999) {
        std::fill(out.begin(), out.end(), 1.0f);
        return out;
    }
    for (size_t i = 0; i < vec_eta.size(); i++) {
        out[i] = ROOT::VecOps::DeltaR(vec_eta[i], obj_eta, vec_phi[i], obj_phi);
    }
    return out;
}

float fInvariantMass(float obj1_pt, float obj1_eta, float obj1_phi, float obj1_mass, 
                    float obj2_pt, float obj2_eta, float obj2_phi, float obj2_mass) {
    TLorentzVector obj1, obj2;
    obj1.SetPtEtaPhiM(obj1_pt, obj1_eta, obj1_phi, obj1_mass);
    obj2.SetPtEtaPhiM(obj2_pt, obj2_eta, obj2_phi, obj2_mass);
    return (obj1 + obj2).M();
}

RVec<float> VInvariantMass(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                          const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invMass(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invMass[i] = (obj1 + obj2).M();
    }
    return invMass;
}

RVec<float> VInvariantPt(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                        const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPt(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPt[i] = (obj1 + obj2).Pt();
    }
    return invPt;
}

RVec<float> VInvariantPhi(const RVec<float>& vec_pt, const RVec<float>& vec_eta, const RVec<float>& vec_phi, 
                         const RVec<float>& vec_mass, float obj_pt, float obj_eta, float obj_phi, float obj_mass) {
    RVec<float> invPhi(vec_pt.size());
    TLorentzVector obj1;
    obj1.SetPtEtaPhiM(obj_pt, obj_eta, obj_phi, obj_mass);

    for (size_t i = 0; i < vec_pt.size(); i++) {
        TLorentzVector obj2;
        obj2.SetPtEtaPhiM(vec_pt[i], vec_eta[i], vec_phi[i], vec_mass[i]);
        invPhi[i] = (obj1 + obj2).Phi();
    }
    return invPhi;
}

RVec<float> VTransverseMass(const RVec<float>& vec_pt, const RVec<float>& vec_phi, float obj_pt, float obj_phi) {
    RVec<float> mt(vec_pt.size());
    for (size_t i = 0; i < vec_pt.size(); i++) {
        mt[i] = std::sqrt(2 * vec_pt[i] * obj_pt * (1 - std::cos(ROOT::VecOps::DeltaPhi(vec_phi[i], obj_phi))));
    }
    return mt;
}

// Return for each ak4 jet, the dR from the closest ak8 jet
RVec<float> dRfromClosestJet(const RVec<float>& ak4_eta, const RVec<float>& ak4_phi, const RVec<float>& ak8_eta, const RVec<float>& ak8_phi) {
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

RVec<RVec<int>> getVBSPairs(const RVec<int>& goodJets, const RVec<float>& jet_var) {
    if (Sum(goodJets) >= 2) {
        return ROOT::VecOps::Combinations(jet_var, 2);
    } else {
    // Create properly matched return type: vector of vector
        RVec<RVec<int>> result;
        // Add two empty vectors
        result.emplace_back(RVec<int>{-990});
        result.emplace_back(RVec<int>{-999});
        return result;
    }
}

/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &ana) {
    // If on UAF and the USER environment variable is defined, store output on ceph
    const char* userEnv = getenv("USER");
    std::string storage_dir = "/data/userdata/";
    std::string output_dir = "./";
    if (userEnv != nullptr && std::filesystem::exists(storage_dir) && std::filesystem::is_directory(storage_dir)) {
        output_dir = storage_dir + std::string(userEnv) + "/vbsvvhAnalysis/preselection/" + ana + "/";
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
        if (ColName.starts_with("_")) {
            continue;
        }
        final_variables.push_back(ColName);
    }

    df.Snapshot("Events", outputDir + "/" + outputFileName + ".root", final_variables);
}