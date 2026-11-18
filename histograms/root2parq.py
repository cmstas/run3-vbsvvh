"""
Skimmed ROOT files in /groups/cjessop/users/jkil/HVV_2L_RDF/run*_*FJ/{longname}/output_*.root takes long time to load.
This script takes multiple ROOT file per dataset into one parquet file, and save in /groups/cjessop/users/jkil/HVV_2L_RDF_parq/run*_*FJ/{shortname}.parq.

To run on condorfe:
ssh condorfe
cd ~/run3-vbsvvh/histograms && conda activate hvv
for ch in run2_1FJ run2_2FJ run3_1FJ run3_2FJ; do python root2parq.py -m condor -c $ch -f; [ $ch != run3_2FJ ] && sleep 5; done
"""
import os, sys, glob
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
from subprocess import run
import awkward as ak
import numpy as np
import uproot
from PlotterTools import SHORT_TO_FULL

RDF_PATH  = "/groups/cjessop/users/jkil/HVV_2L_RDF"
PARQ_PATH = "/groups/cjessop/users/jkil/HVV_2L_RDF_parq"
FULL_TO_SHORT = {full: short for short, full in SHORT_TO_FULL.items()}
branches_to_load = [
    'weight', 'passes_2lep_?FJ_nom', 'njet', 'nLep_Sel',
    'lepton_pt', 'lepton_eta', 'lepton_phi', 'lepton_mass',
    'electron_pt', 'electron_cutBased', 'muon_pt', # makeDataset.LooseLeptons
    'met_pt', 'met_phi',
    'fatjet_pt', 'fatjet_eta', 'fatjet_phi', 'fatjet_massGloParT3', 'fatjet_msoftdrop', 'fatjet_tau1', 'fatjet_tau2',
    'fatjet_Hbb', 'ht_fatjets',
    'vbs_score', 'vbs_mjj', 'vbs_detajj', 'vbs_jet1_idx', 'vbs_jet2_idx', 'jet_pt', 'jet_eta', 'jet_jetId', # makeDataset.VBFCut
    'jet_phi', 'jet_mass', # VBSJ1/VBSJ2 kinematics
    'jet_isMediumBTag', # CalcDerivedVars
]

def MergedChunks(files, max_bytes=5e8):
    # uproot.iterate never merges across files, so buffer the many small per-file chunks into ~max_bytes row groups
    buffer = [uproot.open(f"{files[0]}:Events").arrays(filter_name=branches_to_load, entry_stop=0)] # empty seed keeps the schema when no file has events
    for chunk in uproot.iterate([f"{f}:Events" for f in files], filter_name=branches_to_load, step_size="500 MB"):
        buffer.append(chunk)
        if sum(c.nbytes for c in buffer) > max_bytes:
            yield ak.concatenate(buffer)
            buffer = []
    if buffer: yield ak.concatenate(buffer)

def CalcDerivedVars(chunk, L1, L2):
    dEta, dPhi = L1.eta - L2.eta, abs((L1.phi - L2.phi + np.pi) % (2 * np.pi) - np.pi)
    return {
        'nbjet': ak.sum(chunk.jet_isMediumBTag, axis=1), # same as VBSVVH's: jet_* is already cleaned against fatjets
        'dilep_dR': np.hypot(dEta, dPhi),
        'dilep_dphi': dPhi,
        'dilep_mass': np.sqrt(2 * L1.pt * L2.pt * (np.cosh(dEta) - np.cos(dPhi))), # massless: lepton_mass is rounding noise for Run 2 electrons
    }

def Flatten(chunk):
    # Jagged collections -> 2 None-padded objects each (L1/L2 by pT, FJ1/FJ2, VBSJ1/VBSJ2), so all channels share one flat schema
    collection = lambda prefix: ak.zip({f.removeprefix(prefix): chunk[f] for f in chunk.fields if f.startswith(prefix)})
    ptOrder = ak.argsort(ak.concatenate([chunk.electron_pt, chunk.muon_pt], axis=1), ascending=False)
    isLoose = ak.concatenate([chunk.electron_cutBased >= 2, chunk.muon_pt > 0], axis=1)[ptOrder] # lepton_* = pT-sorted veto e + loose mu; keep the loose ones nLep_Sel counts
    
    leptons = ak.pad_none(collection('lepton_')[isLoose], 2)
    fatjets = ak.pad_none(collection('fatjet_'), 2)
    jets    = collection('jet_')
    jetIdx  = ak.local_index(chunk.jet_pt, axis=1) # vbs_jet*_idx = -1 (no pair) matches nothing -> None
    objects = {
        'L1': leptons[:, 0], 
        'L2': leptons[:, 1], 
        'FJ1': fatjets[:, 0], 
        'FJ2': fatjets[:, 1],
        'VBSJ1': ak.firsts(jets[jetIdx == chunk.vbs_jet1_idx]), 
        'VBSJ2': ak.firsts(jets[jetIdx == chunk.vbs_jet2_idx])
    }
    scalars = {'passes_2lep_nom' if f.startswith('passes_') else f: chunk[f] for f in chunk.fields if chunk[f].ndim == 1}
    flatBranches = scalars | CalcDerivedVars(chunk, objects['L1'], objects['L2']) | {f"{name}_{f}": obj[f] for name, obj in objects.items() for f in obj.fields} # combining dicts
    return ak.zip(flatBranches)

def Convert(sampleDir, outPath):
    files = sorted(glob.glob(f"{sampleDir}/*.root"))
    ak.to_parquet_row_groups(map(Flatten, MergedChunks(files)), f"{outPath}.tmp")
    os.replace(f"{outPath}.tmp", outPath) # an interrupted write never leaves a complete-looking parquet
    src, dst = Path(sampleDir), Path(outPath)
    print(f"\033[1;32m{src.parent.parent.name}/{src.parent.name}/{dst.stem} ({len(files)} files) -> {'/'.join(dst.parts[-3:])} done!\033[0m", flush=True)

def SamplesToConvert(channels, outdir, force):
    assert len(FULL_TO_SHORT) == len(SHORT_TO_FULL), "SHORT_TO_FULL is not bijective!"
    for ch in channels: os.makedirs(f"{outdir}/{ch}", exist_ok=True)
    todo = {f"{outdir}/{ch}/{FULL_TO_SHORT[full]}.parq": f"{RDF_PATH}/{ch}/{full}" for ch in channels for full in sorted(os.listdir(f"{RDF_PATH}/{ch}"))}
    return {out: src for out, src in todo.items() if force or not os.path.exists(out)}

def SubmitCondor(todo, memory, dryRun):
    # One job per sample on the CRC opportunistic pool; jobs run this python and read/write /groups directly, so nothing is transferred
    if not dryRun and not os.uname().nodename.startswith('condorfe'):
        sys.exit("Submit from condorfe.crc.nd.edu; condor_submit on glados goes to the ndcms pool.")
    logDir = f"{os.path.dirname(os.path.abspath(__file__))}/{datetime.now().strftime('%y%m%d')}_root2parq_condor"
    os.makedirs(logDir, exist_ok=True)
    with open(f"{logDir}/samples.txt", 'w') as f: f.writelines(f"{src} {out} {Path(out).parent.name}_{Path(out).stem}\n" for out, src in todo.items())
    with open(f"{logDir}/submit.jdl", 'w') as f: 
        f.write(
            f"""
            universe              = vanilla
            executable            = {sys.executable}
            arguments             = {os.path.abspath(__file__)} --sample $(src) $(out)
            transfer_executable   = False
            should_transfer_files = NO
            request_cpus          = 1
            request_memory        = {memory}
            max_retries           = 3
            log                   = {logDir}/$(name).log
            output                = {logDir}/$(name).out
            error                 = {logDir}/$(name).err
            +JobBatchName         = "root2parq"
            queue src, out, name from {logDir}/samples.txt
            """
        )
    print(f"submit file: {logDir}/submit.jdl")
    if not dryRun: run(["condor_submit", f"{logDir}/submit.jdl"], check=True)

if __name__ == "__main__":
    parser = ArgumentParser(description='Merge the RDF output ROOT files of each sample into one parquet.')
    parser.add_argument('-m', '--mode'    , type=str, default='local', choices=['local', 'condor'], help='condor: one job per sample on condorfe opportunistic resources.')
    parser.add_argument('-c', '--channels', type=str, nargs='+', default=['run2_1FJ', 'run2_2FJ', 'run3_1FJ', 'run3_2FJ'])
    parser.add_argument('-o', '--outdir'  , type=str, default=PARQ_PATH)
    parser.add_argument('-l', '--limit'   , type=int, default=None, help='Convert only the first N samples (for testing).')
    parser.add_argument('-f', '--force'   , action='store_true', help='If True, overwrites existing parquets. Default: skip them.')
    parser.add_argument('-j', '--jobs'    , type=int, default=1, help='local: number of samples converted in parallel.')
    parser.add_argument('--memory'        , type=str, default='4 GB', help='condor: request_memory per job.')
    parser.add_argument('-n', '--dry-run' , action='store_true', help='condor: writes the submit file without submitting.')
    parser.add_argument('--sample'        , type=str, nargs=2, metavar=('SRC', 'OUT'), help='Converts one sample; run by each condor job.')
    args = parser.parse_args()

    if args.sample:
        Convert(*args.sample)
        sys.exit()

    todo = dict(list(SamplesToConvert(args.channels, args.outdir, args.force).items())[:args.limit])
    print(f"{len(todo)} samples to convert")
    if args.mode == 'condor':
        SubmitCondor(todo, args.memory, args.dry_run)
    else:
        with ProcessPoolExecutor(args.jobs) as pool: list(pool.map(Convert, todo.values(), todo.keys()))

"""
Field names:
 1. weight           12. passes_2lep_nom  23. L2_phi           34. FJ2_eta          45. VBSJ1_phi
 2. met_pt           13. nbjet            24. L2_mass          35. FJ2_massGloParT3 46. VBSJ1_pt
 3. met_phi          14. dilep_dR         25. FJ1_Hbb          36. FJ2_msoftdrop    47. VBSJ2_eta
 4. nLep_Sel         15. dilep_dphi       26. FJ1_eta          37. FJ2_phi          48. VBSJ2_isMediumBTag
 5. ht_fatjets       16. dilep_mass       27. FJ1_massGloParT3 38. FJ2_pt           49. VBSJ2_jetId
 6. njet             17. L1_pt            28. FJ1_msoftdrop    39. FJ2_tau1         50. VBSJ2_mass
 7. vbs_jet1_idx     18. L1_eta           29. FJ1_phi          40. FJ2_tau2         51. VBSJ2_phi
 8. vbs_jet2_idx     19. L1_phi           30. FJ1_pt           41. VBSJ1_eta        52. VBSJ2_pt
 9. vbs_score        20. L1_mass          31. FJ1_tau1         42. VBSJ1_isMediumBTag
10. vbs_mjj          21. L2_pt            32. FJ1_tau2         43. VBSJ1_jetId
11. vbs_detajj       22. L2_eta           33. FJ2_Hbb          44. VBSJ1_mass
"""