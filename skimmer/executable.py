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

r.gInterpreter.Declare('#include "truthSelections.h"')

subprocess.run("python3 -m pip install --user --no-binary=correctionlib correctionlib", shell=True, check=True)
import importlib
correctionlib = importlib.import_module("correctionlib")
correctionlib.register_pyroot_binding()

# Constants
CONDOR_OUTPUT_DIR = "output"
XROOTD_REDIRECTOR = "root://xrootd-cms.infn.it/"
OUTPUT_XRD = "davs://redirector.t2.ucsd.edu:1095//store/user/USER_UAF_DIR/skim/" # Change to user's skim directory on UAF
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
        # define truth event from gen particles
        df = df.Define("truth_indices", "getTruthEventInfo(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass, false)") \
            .Define("gen_h_idx", "truth_indices[0]") \
            .Define("gen_b1_idx", "truth_indices[1]") \
            .Define("gen_b2_idx", "truth_indices[2]") \
            .Define("gen_v1_idx", "truth_indices[3]") \
            .Define("gen_v1q1_idx", "truth_indices[4]") \
            .Define("gen_v1q2_idx", "truth_indices[5]") \
            .Define("gen_v2_idx", "truth_indices[6]") \
            .Define("gen_v2q1_idx", "truth_indices[7]") \
            .Define("gen_v2q2_idx", "truth_indices[8]") \
            .Define("gen_vbs1_idx", "truth_indices[9]") \
            .Define("gen_vbs2_idx", "truth_indices[10]") 
        
        return df
        
    def analyze(self, is_signal):
        self.df = self.df.Define("__tight_mu_mask", "Muon_pt > 35. && abs(Muon_eta) < 2.4 && Muon_tightId") \
            .Define("__tight_ele_mask", "Electron_pt > 35. && abs(Electron_eta) < 2.5 && Electron_cutBased >= 4") \
            .Define("__n_tight_leptons", "Sum(__tight_mu_mask) + Sum(__tight_ele_mask)") \
            .Define("__fatjet_mask", "FatJet_pt > 200 && FatJet_msoftdrop > 10") \
            .Define("__n_fatjets", "Sum(__fatjet_mask)") \
            .Filter("(__n_fatjets + __n_tight_leptons) >= 1")

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

        # Run3 event filters
        self.df = self.df.Filter("Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter")

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


def determine_output_paths(input_file, is_signal, output_tag):
    if not is_signal:
        era = input_file.split('/')[3]
        sample_name = input_file.split('/')[4]
        campaign = input_file.split('/')[6]
    else:
        era = input_file.split('/')[6]
        sample_name = input_file.split('/')[7]
        campaign = "private"

    output_dir = f"{OUTPUT_XRD}/skims_{output_tag}/{campaign}/{sample_name}"
    return output_dir

def check_output_liveness(file):
    with r.TFile.Open(file) as f:
        t = f.Get("Events")
        nevts = t.GetEntries()
        for i in range(0,t.GetEntries(),1):
            if t.GetEntry(i) < 0:
                return False
        return True

def copy_output_file(source, destination):
    print(f"Copying {source} to {destination}")
    
    subprocess.run(["gfal-mkdir", "-p", os.path.dirname(destination)])

    # check if output is good
    if not check_output_liveness(source):
        print(f"Output file {source} is corrupted; not copying")
        return False
    
    result = subprocess.run(["gfal-copy", "-f", source, destination])

    if result.returncode != 0:    
        print(f"Failed to copy {source} to {destination}; sleeping for 60s")
        return False

    # check copied file liveness
    if not check_output_liveness(destination):
        print(f"Copied file {destination} is corrupted")
        subprocess.run(["gfal-rm", destination])
        print(f"Removed corrupted file {destination}")
        return False

    return True

if __name__ == "__main__":
    parser = ArgumentParser(description='Run the NanoAOD skimmer with file transfer.')
    parser.add_argument('proxy', help="Path to the X509 proxy")
    parser.add_argument('input_file', help="Input file path")
    parser.add_argument('job_id', help="Job ID")
    parser.add_argument('is_signal', help='Flag indicating if this is a signal sample', type=int)
    parser.add_argument('output_tag', help='Output tag, including version of skims eg. v2', type=str)
    args = parser.parse_args()
    
    os.environ['X509_USER_PROXY'] = args.proxy
    
    success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    if not success:
        print("Skimmer failed; retrying one more time...")
        success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    merge_skims(CONDOR_OUTPUT_DIR)
    
    output_dir = determine_output_paths(args.input_file, args.is_signal, args.output_tag)
    
    copy_src = os.path.join(os.getcwd(), f"{CONDOR_OUTPUT_DIR}/output.root")
    copy_dest = f"{output_dir}/output_{args.job_id}.root"
    
    for attempt in range(MAX_RETRIES + 1):
        success = copy_output_file(copy_src, copy_dest)
        if success:
            break
        
        if attempt < MAX_RETRIES:
            print(f"Retrying copy attempt {attempt + 1} of {MAX_RETRIES}...")
            time.sleep(SLEEP_DURATION)

    if not success:
        print(f"Failed to copy output file after {MAX_RETRIES} attempts")
        sys.exit(1)
        
    sys.exit(0)
