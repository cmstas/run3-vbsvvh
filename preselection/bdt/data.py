import ROOT as r
import awkward as ak
import pyarrow as pa
import pyarrow.parquet as pq

r.EnableImplicitMT(96)

r.gInterpreter.Declare("""
auto pt_m_jj = [](const ROOT::RVec<float>& jet1_pt, const ROOT::RVec<float>& jet1_eta, const ROOT::RVec<float>& jet1_phi, const ROOT::RVec<float>& jet1_mass, const ROOT::RVec<float>& jet2_pt, const ROOT::RVec<float>& jet2_eta, const ROOT::RVec<float>& jet2_phi, const ROOT::RVec<float>& jet2_mass) {
        ROOT::VecOps::RVec<float> pt_jj;
        ROOT::VecOps::RVec<float> m_jj;
        for (size_t i = 0; i < jet1_pt.size(); ++i) {
            auto v_jj = ROOT::Math::PtEtaPhiMVector(jet1_pt[i], jet1_eta[i], jet1_phi[i], jet1_mass[i]) + ROOT::Math::PtEtaPhiMVector(jet2_pt[i], jet2_eta[i], jet2_phi[i], jet2_mass[i]);
            pt_jj.push_back(v_jj.Pt());
            m_jj.push_back(v_jj.M());
        }
    return std::make_pair(pt_jj, m_jj);
};
""")

def process_signal():
    df = r.RDataFrame("Events", ["/data/userdata/aaarora/spanet_training/run2.root", "/data/userdata/aaarora/spanet_training/run3.root"])
    r.RDF.Experimental.AddProgressBar(df)

    df = df.Filter("jet_pt.size() >= 2", "At least 2 jets")
    df = df.Filter("truth_vbs1_idx >= 0 && truth_vbs2_idx >= 0", "Valid truth VBS jet indices")

    df = df.Define("jet_pair_idx", "ROOT::VecOps::Combinations(jet_pt, 2)") \
        .Define("jet1_pt", "ROOT::VecOps::Take(jet_pt, jet_pair_idx[0])") \
        .Define("jet2_pt", "ROOT::VecOps::Take(jet_pt, jet_pair_idx[1])") \
        .Define("jet1_eta", "ROOT::VecOps::Take(jet_eta, jet_pair_idx[0])") \
        .Define("jet2_eta", "ROOT::VecOps::Take(jet_eta, jet_pair_idx[1])") \
        .Define("jet1_phi", "ROOT::VecOps::Take(jet_phi, jet_pair_idx[0])") \
        .Define("jet2_phi", "ROOT::VecOps::Take(jet_phi, jet_pair_idx[1])") \
        .Define("jet1_mass", "ROOT::VecOps::Take(jet_mass, jet_pair_idx[0])") \
        .Define("jet2_mass", "ROOT::VecOps::Take(jet_mass, jet_pair_idx[1])") \
        .Define("pt_m_jj", "pt_m_jj(jet1_pt, jet1_eta, jet1_phi, jet1_mass, jet2_pt, jet2_eta, jet2_phi, jet2_mass)") \
        .Define("pt_jj", "pt_m_jj.first") \
        .Define("m_jj", "pt_m_jj.second") \
        .Define("deta_jj", "abs(jet1_eta - jet2_eta)") \
        .Define("dphi_jj", "ROOT::VecOps::DeltaPhi(jet1_phi, jet2_phi)") \
        .Define("labels", "(jet_pair_idx[0] == truth_vbs1_idx && jet_pair_idx[1] == truth_vbs2_idx || jet_pair_idx[1] == truth_vbs1_idx && jet_pair_idx[0] == truth_vbs2_idx)")
    
    return df

if __name__ == "__main__":
    output_cols = ["jet1_pt", "jet2_pt", "jet1_eta", "jet2_eta", "jet1_phi", "jet2_phi", 
                "jet1_mass", "jet2_mass", "pt_jj", "deta_jj", "dphi_jj", "m_jj", "labels"]
    
    df_sig = ak.to_parquet(ak.from_rdataframe(process_signal(), output_cols), "sig.parquet")