#!/usr/bin/python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import time
import shutil
import subprocess
from argparse import ArgumentParser

import ROOT as r

r.gInterpreter.Declare("""
ROOT::RVec<int> findHiggsAndDaughters(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    ROOT::RVec<int> result = {-1, -1, -1};
    
    int higgs_idx = -1;
    ROOT::RVec<int> hdecay_idx;
    
    // Find intermediate Higgs with status 22 and mother_idx 0
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && part_pdgId == 25) {
            higgs_idx = igen;
            break;
        }
    }
    
    if (higgs_idx == -1) return result;
    
    // Find first daughters of Higgs (non-Higgs particles with Higgs mother)
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int mother_idx = motherIdx[igen];
        if (mother_idx > 0 && mother_idx < pdgId.size()) {
            int mother_pdgId = pdgId[mother_idx];
            int part_pdgId = pdgId[igen];
            
            if (mother_pdgId == 25 && part_pdgId != 25) {
                hdecay_idx.push_back(igen);
            }
        }
    }
    
    result[0] = higgs_idx;
    if (hdecay_idx.size() >= 1) result[1] = hdecay_idx[0];
    if (hdecay_idx.size() >= 2) result[2] = hdecay_idx[1];
    
    return result;
}

int findLastIndex(int current_idx, int current_pdgId, ROOT::RVec<int>& pdgId, ROOT::RVec<short>& motherIdx) {
    int outIdx = current_idx;
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == current_idx) {
            if (part_pdgId == current_pdgId) {
                outIdx = findLastIndex(igen, part_pdgId, pdgId, motherIdx);
            }
        }
    }
    return outIdx;
}

ROOT::RVec<int> findVBosonsAndDaughters(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    // Returns: [v1_idx, v1_d1_idx, v1_d2_idx, v2_idx, v2_d1_idx, v2_d2_idx] or [-1, -1, -1, -1, -1, -1]
    ROOT::RVec<int> result = {-1, -1, -1, -1, -1, -1};
    
    ROOT::RVec<int> firstVs_idx;
    
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && (part_pdgId == 23 || abs(part_pdgId) == 24)) {
            firstVs_idx.push_back(igen);
        }
    }
    
    if (firstVs_idx.size() < 2) return result;
    
    ROOT::RVec<int> hadronic_V_indices;
    ROOT::RVec<int> leptonic_V_indices;
    
    for (size_t iV = 0; iV < firstVs_idx.size() && iV < 2; ++iV) {
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        int hadronic_daughters = 0;
        int leptonic_daughters = 0;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                    hadronic_daughters++;
                } else {
                    leptonic_daughters++;
                }
            }
        }
        
        if (hadronic_daughters > 0) {
            hadronic_V_indices.push_back(iV);
        } else if (leptonic_daughters > 0) {
            leptonic_V_indices.push_back(iV);
        }
    }
    
    if (hadronic_V_indices.size() > 0) {
        int iV = hadronic_V_indices[0];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        ROOT::RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) { // quarks only
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[0] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[1] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[2] = vdecays_idx[1];
    }
    
    if (hadronic_V_indices.size() >= 2) {
        int iV = hadronic_V_indices[1];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        ROOT::RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[3] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[4] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[5] = vdecays_idx[1];
    }
    
    return result;
}

ROOT::RVec<int> findVBSQuarks(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    // Returns: [vbs1_idx, vbs2_idx] or [-1, -1]
    ROOT::RVec<int> result = {-1, -1};
    ROOT::RVec<int> vbsquarks_idx;
    
    // Find outgoing quarks with status 23 and mother_idx 0
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 23 && abs(part_pdgId) <= 6 && abs(part_pdgId) >= 1) {
            vbsquarks_idx.push_back(igen);
        }
    }
    
    if (vbsquarks_idx.size() >= 1) result[0] = vbsquarks_idx[0];
    if (vbsquarks_idx.size() >= 2) result[1] = vbsquarks_idx[1];
    
    return result;
}
""")

subprocess.run("python3 -m pip install --user --no-binary=correctionlib correctionlib", shell=True, check=True)
import importlib
correctionlib = importlib.import_module("correctionlib")
correctionlib.register_pyroot_binding()

# Constants
CONDOR_OUTPUT_DIR = "output"
XROOTD_REDIRECTOR = "root://xrootd-cms.infn.it/"
OUTPUT_XRD = "davs://redirector.t2.ucsd.edu:1095//store/user/aaarora/skims"
MAX_RETRIES = 10
SLEEP_DURATION = 60  # 1 minute in seconds

JET_ID_JSONS = {"2024": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2024_Summer24/jetid.json.gz"}

class Skimmer():
    def __init__(self, inFiles, outDir, keepDropFile):
        self.inFiles = inFiles
        self.outDir = outDir
        self.keepDropFile = keepDropFile

        self.df = r.RDataFrame("Events", self.inFiles)
        r.RDF.Experimental.AddProgressBar(self.df)
        columns = self.df.GetColumnNames()
        for col in columns:
            if col.startswith("Muon_") or col.startswith("Electron_") or col.startswith("Jet_") or col.startswith("FatJet_"):
                colType = self.df.GetColumnType(col)
                if colType == "ROOT::VecOps::RVec<Float_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Float_t')())
                elif colType == "ROOT::VecOps::RVec<Int_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Int_t')())
                elif colType == "ROOT::VecOps::RVec<UShort_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UShort_t')())
                elif colType == "ROOT::VecOps::RVec<Bool_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Bool_t')())
                elif colType == "ROOT::VecOps::RVec<UChar_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UChar_t')())
                elif colType == "ROOT::VecOps::RVec<Double_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Double_t')())
                elif colType == "ROOT::VecOps::RVec<Long64_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Long64_t')())
                elif colType == "ROOT::VecOps::RVec<Short_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Short_t')())
                else:
                    print(f"Unknown column type: {colType} for column {col}")

        self.df.Display(["Electron_pt"]).Print()

    @staticmethod
    def genSelection(df):
        # define bb, qq definition
        df = df.Define("higgs_info", "findHiggsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_h_idx", "higgs_info[0]") \
            .Define("gen_b1_idx", "higgs_info[1]") \
            .Define("gen_b2_idx", "higgs_info[2]")

        df = df.Define("vboson_info", "findVBosonsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_v1_idx", "vboson_info[0]") \
            .Define("gen_v1q1_idx", "vboson_info[1]") \
            .Define("gen_v1q2_idx", "vboson_info[2]") \
            .Define("gen_v2_idx", "vboson_info[3]") \
            .Define("gen_v2q1_idx", "vboson_info[4]") \
            .Define("gen_v2q2_idx", "vboson_info[5]")

        df = df.Define("vbs_info", "findVBSQuarks(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_vbs1_idx", "vbs_info[0]") \
            .Define("gen_vbs2_idx", "vbs_info[1]") 
        
        return df
        
    def analyze(self, is_signal):
        self.df = self.df.Define("tight_mu_mask", "Muon_pt > 35. && abs(Muon_eta) < 2.4 && Muon_tightId") \
            .Define("tight_ele_mask", "Electron_pt > 35. && abs(Electron_eta) < 2.5 && Electron_cutBased >= 4") \
            .Filter("Sum(tight_mu_mask) + Sum(tight_ele_mask) < 2") \
            .Define("fatjet_mask", "FatJet_pt > 200 && FatJet_msoftdrop > 10 && FatJet_mass > 10") \
            .Filter("Sum(fatjet_mask) >= 1")

        if self.sample_year in JET_ID_JSONS:
            jet_id_json = JET_ID_JSONS[self.sample_year]

        self.df = self.df.Define("Jet_multiplicity", "Jet_chMultiplicity + Jet_neMultiplicity") \
            .Define("FatJet_multiplicity", "FatJet_chMultiplicity + FatJet_neMultiplicity")

        r.gInterpreter.Declare("""
            #include <ROOT/RVec.hxx>
            using namespace ROOT::VecOps;
            
            RVec<float> evalJetID(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                                  const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                                  const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                                  const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
                auto cset_jetId = correction::CorrectionSet::from_file(\"""" + jet_id_json + """\");
                RVec<float> jetId(eta.size(), 0.0);
                for (size_t i = 0; i < eta.size(); ++i) {
                    jetId[i] += 2 * cset_jetId->at(\"AK4PUPPI_Tight\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                    jetId[i] += 4 * cset_jetId->at(\"AK4PUPPI_TightLeptonVeto\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                }
                return jetId;
            }

            RVec<float> evalFatJetID(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                                     const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                                     const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                                     const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
                auto cset_fatJetId = correction::CorrectionSet::from_file(\"""" + jet_id_json + """\");
                RVec<float> fatJetId(eta.size(), 0.0);
                for (size_t i = 0; i < eta.size(); ++i) {
                    fatJetId[i] += 2 * cset_fatJetId->at(\"AK8PUPPI_Tight\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                    fatJetId[i] += 4 * cset_fatJetId->at(\"AK8PUPPI_TightLeptonVeto\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                }
                return fatJetId;
            }

        """)

        self.df = self.df.Define("Jet_jetId", "evalJetID(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, Jet_multiplicity)") \
            .Define("FatJet_jetId", "evalFatJetID(FatJet_eta, FatJet_chHEF, FatJet_neHEF, FatJet_chEmEF, FatJet_neEmEF, FatJet_muEF, FatJet_chMultiplicity, FatJet_neMultiplicity, FatJet_multiplicity)")

        if is_signal:
            self.df = self.genSelection(self.df)

        return self.df.Count().GetValue()

    def Snapshot(self, tag):
        all_cols = [str(col) for col in self.df.GetColumnNames()]
        keep_cols = {col: 0 for col in all_cols}
        comment = re.compile(r"#.*")
        ops = []
        with open(self.keepDropFile, 'r') as f:
            for line in f:
                # convert to python regex
                if len(line) == 0 or line[0] == '#': 
                    continue
                line = re.sub(comment, "", line)
                (op, sel) = line.split()
                if op == "keep":
                    ops.append((sel, 1))
                elif op == "drop":
                    ops.append((sel, 0))
        
        for bre, stat in ops:
            try:
                re.compile(bre)
                for n in all_cols:
                    if re.match(bre, n):
                        keep_cols[n] = stat
            except re.error:
                keep_cols[bre] = stat

        keep_cols = [k for k, v in keep_cols.items() if v == 1]

        self.df.Snapshot("Events", self.outDir + "/" + tag + ".root", keep_cols)
        
        snap_opts = r.RDF.RSnapshotOptions()
        snap_opts.fMode = "UPDATE"

        runs_df = r.RDataFrame("Runs", self.inFiles)
        cols = [str(col) for col in runs_df.GetColumnNames()]
        runs_df.Snapshot("Runs", self.outDir + "/" + tag + ".root", cols, snap_opts)

    @property
    def sample_year(self):
        match = re.search(r'Run3Summer24|RunIII2024Summer24NanoAODv15', self.inFiles[0])
        if match:
            return "2024"
        else:
            return None

def run_skimmer(input_file, output_dir, is_signal):
    print(f"Running skimmer on {input_file}")
    os.makedirs(output_dir, exist_ok=True)
    
    inFiles = [XROOTD_REDIRECTOR + input_file if input_file.startswith('/store') else 'file://' + input_file]
    keepDropFile = "keep_and_drop_skim.txt"
    
    skimmer = Skimmer(inFiles, output_dir, keepDropFile)
    passed = skimmer.analyze(is_signal)
    if passed:
        skimmer.Snapshot("skim")
        return True
    else:
        print("No entries in output")
        return False


def merge_skims(output_dir):
    skim_files = glob.glob(f"{output_dir}/*")
    
    if len(skim_files) == 0:
        print("No output files to merge; exiting...")
        return True
    elif len(skim_files) == 1:
        shutil.move(skim_files[0], f"{output_dir}/output.root")
        return True
    else:
        merge_cmd = ["hadd", f"{output_dir}/output.root"] + skim_files
        print(" ".join(merge_cmd))
        result = subprocess.run(merge_cmd)
        return result.returncode == 0


def determine_output_paths(input_file, is_signal):
    if not is_signal:
        era = input_file.split('/')[3]
        sample_name = input_file.split('/')[4]
        campaign = input_file.split('/')[6]
    else:
        era = input_file.split('/')[6]
        sample_name = input_file.split('/')[7]
        campaign = "private"
        
    output_dir = f"{OUTPUT_XRD}/{era}/{campaign}/{sample_name}"
    return output_dir

def copy_output_file(source, destination):
    print(f"Copying {source} to {destination}")
    
    # Create destination directory
    subprocess.run(["gfal-mkdir", "-p", os.path.dirname(destination)])
    
    # Copy with retries
    for i in range(1, MAX_RETRIES + 1):
        print(f"Attempt {i}")
        result = subprocess.run(["gfal-copy", "-f", source, destination])
        if result.returncode == 0:
            return True
        
        print(f"Failed to copy {source} to {destination}; sleeping for 60s")
        time.sleep(SLEEP_DURATION)
    
    return False

if __name__ == "__main__":
    parser = ArgumentParser(description='Run the NanoAOD skimmer with file transfer.')
    parser.add_argument('proxy', help="Path to the X509 proxy")
    parser.add_argument('input_file', help="Input file path")
    parser.add_argument('job_id', help="Job ID")
    parser.add_argument('is_signal', help='Flag indicating if this is a signal sample', type=int)
    args = parser.parse_args()
    
    # Set up X509 proxy
    os.environ['X509_USER_PROXY'] = args.proxy

    # Run the skimmer
    success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    # Retry once if failed
    if not success:
        print("Skimmer failed; retrying one more time...")
        success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    # Merge results
    merge_skims(CONDOR_OUTPUT_DIR)
    
    # Determine output paths
    output_dir = determine_output_paths(args.input_file, args.is_signal)
    
    # Copy the output file
    copy_src = os.path.join(os.getcwd(), f"{CONDOR_OUTPUT_DIR}/output.root")
    copy_dest = f"{output_dir}/output_{args.job_id}.root"
    
    success = copy_output_file(copy_src, copy_dest)
    if not success:
        print(f"Failed to copy output file after {MAX_RETRIES} attempts")
        sys.exit(1)
    
    sys.exit(0)
