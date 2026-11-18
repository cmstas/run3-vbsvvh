"""
Creates abcd/dataset/{signal}/{tag}_{train/valid}.pt, as well as .parqs in same dir.
Same as makeDataset.py, but reads the per-sample parquets from histograms/root2parq.py instead of the RDF ROOT files.
"""

from TrainingTools import *
import re, time, yaml, json, warnings
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "histograms"))
from functions import chooseOverlaps
from root2parq import PARQ_PATH, FULL_TO_SHORT
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
warnings.filterwarnings("ignore", message=".*reduce_op.*")

def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def extract_expression_variables(expression):
    tokens = re.findall(r"(?<!\.)\b[A-Za-z_][A-Za-z0-9_]*\b", expression)
    excluded = {"np", "abs", "and", "or", "not", "True", "False"}
    return [tok for tok in tokens if tok not in excluded]

def parse_training_features(cfg, flavor):
    # ——————— Fetching features/transforms ————————————————————————————
    transforms = {feature: transform for f in cfg['training_features'] for feature, transform in f.items()} # list of dicts
    features   = list(transforms)

    # ——————— Fetching constraint —————————————————————————————————————
    constraint = cfg['constraint_var']
    if flavor == "double" and constraint not in features:
        features.append(constraint)
        transforms[constraint] = cfg.get("constraint_as_feature_transform", "none")

    return features, transforms, constraint

def parse_derived_vars(cfg, features, transforms):
    # ——————— Fetching derived vars ———————————————————————————————————
    dv_dict = cfg['derived_vars']
    derivedVars = set(dv_dict.keys())
    input_for_derivedVars = []
    for variable in dv_dict.values():
        input_for_derivedVars.extend(var for var in extract_expression_variables(variable) if var not in derivedVars)
    input_for_derivedVars = list(dict.fromkeys(input_for_derivedVars))

    # ——————— Adding derived vars to features (optional) ——————————————
    include_dv = cfg.get("auto_include_derived_vars", False)
    dvt = cfg.get("auto_include_derived_vars_transform", "none")
    if include_dv:
        to_add = [name for name in derivedVars if name not in features]
        features.extend(to_add)
        transforms.update({name: dvt for name in to_add})

    return dv_dict, derivedVars, input_for_derivedVars, features, transforms

def auto_include_derived_vars(cfg, features, transforms, derived_vars_cfg):
    if not cfg.get("auto_include_derived_vars", False):
        return
    auto_tf = cfg.get("auto_include_derived_vars_transform", "none")
    added = [name for name in derived_vars_cfg if name not in features]
    features.extend(added)
    transforms.update({name: auto_tf for name in added})
    print(f"Auto-include derived vars enabled: added {len(added)} derived features with transform='{auto_tf}'")

def apply_derived_vars(events, DerivedVars):
    for name, expression in DerivedVars.items():
        events[name] = eval(expression, locals={f: events[f] for f in events.fields})
    return events

def evaluate_preselection(events, expression):
    cut = eval(expression, {"__builtins__": {}, "np": np, "ak": ak}, {f: events[f] for f in events.fields})
    return events[ak.fill_none(cut, False)] # cuts on None-padded objects (e.g. FJ2_*) give None, which would otherwise keep the event as None

VBF_BRANCHES = ["vbs_detajj", "vbs_mjj"] + [f"VBSJ{i}_{v}" for i in (1, 2) for v in ("pt", "eta", "jetId")]
def VBFCut(events, features):
    # Same as isVBF in VBSVVH/01_processing/utils/CategorizerTools.py:setVBFJets: deta > 2.5, mjj > 750, and both BDT-tagged VBS jets pass SelectGoodAK4Jets
    goodJet = lambda pt, eta, jetId: (pt > 30) & (jetId > 2) & (abs(eta) < 5.0) & ~((pt <= 50) & (abs(eta) > 2.5) & (abs(eta) < 3.0))
    goodVBS = goodJet(events.VBSJ1_pt, events.VBSJ1_eta, events.VBSJ1_jetId) & goodJet(events.VBSJ2_pt, events.VBSJ2_eta, events.VBSJ2_jetId)
    events  = events[(events.vbs_detajj > 2.5) & (events.vbs_mjj > 750) & ak.fill_none(goodVBS, False)] # VBSJ* are None when there is no VBS pair
    return events[[f for f in events.fields if f not in VBF_BRANCHES or f in features]] # drop the extra branches

# LEPTON_ID_BRANCHES = ["electron_pt", "muon_pt", "electron_cutBased"]
# def LooseLeptons(arr, lep):
#     # In preselection, 2lep selection counts Loose leptons (veto for muon). So it's possible that third veto electron enter events in RDF.
#     # When selecting leadlep and subleadlep in makeDataset, if that veto e has enough pT, it can enter field[:,0] or field[:,1].
#     # Thus, I am making electron-level cut to drop any veto electrons.
#     # lepton_* = pT-sorted veto electrons + loose muons; keep only the loose ones nLep_Sel counts (mirrors LeptonSelections in selections.cpp)
#     isLoose = ak.concatenate([lep.electron_cutBased >= 2, lep.muon_pt > 0], axis=1)[ak.argsort(ak.concatenate([lep.electron_pt, lep.muon_pt], axis=1), ascending=False)]
#     for field in [f for f in arr.fields if f.startswith("lepton_")]: arr[field] = arr[field][isLoose]
#     return arr

def local_paths_from_json(jsonpath, parq_dir, kind='sig', signal="c2v1p5_c3_1p0"):
    with open(jsonpath, 'r') as f: samples = json.load(f)["samples"]
    names = [name for name, s in samples.items() if s["metadata"]["kind"] == kind and (kind != "sig" or s["metadata"]["shortname"].endswith(signal))]
    paths = {name: f"{parq_dir}/{FULL_TO_SHORT[name]}.parq" for name in names}
    return {name: path for name, path in paths.items() if os.path.exists(path)} # samples absent from the RDF outputs have no parquet

def LoadEvents(paths, features, label):
    """Trying to migrate to awkward array"""
    rejected = chooseOverlaps()
    paths = {name: path for name, path in sorted(paths.items()) if name not in rejected}
    passFlag = "passes_2lep_nom"
    arrays, samples = [], {}
    for file_idx, (name, path) in enumerate(paths.items()):
        arr = ak.from_parquet(path, columns=features + [passFlag])
        assert len(arr.fields) == len(features) + 1, f"{path} lacks {set(features) - set(arr.fields)}" # from_parquet silently skips missing columns
        # arr = LooseLeptons(arr, ak.from_parquet(path, columns=LEPTON_ID_BRANCHES))
        arr["file_idx"] = file_idx
        arr["proc_idx"] = samples.setdefault(name, len(samples)) # increments per sample (e.g. DY=0, ttbar=1)
        arrays.append(arr)
    events = ak.concatenate(arrays)

    # Additional cutflow
    events = events[events.weight>0]
    events = events[events.vbs_score != -999] # TEMP: -999 means nominal njet < 2; mirrors the njet >= 2 cut now in selections.cpp
    events = events[events[passFlag]] # nominal channel cut; preselection keeps events passing under ANY JES/JER variation
    # events = VBFCut(events, features)

    # Giving sig/bkg label for training
    events["label"] = int(label)

    return events, list(samples) # samples indexed by proc_idx

def log_transform(quantity):
    # NOTE: Negatives are rounded-m² artifacts, physically ~0. Clipping at 0 instead of -0.999
    # keeps log1p(-0.999)=-6.9 from dominating the min-max range that follows.
    return np.log1p(np.clip(quantity, a_min=0, a_max=None))

def minmax_transform(quantity):
    scaler = MinMaxScaler(feature_range=(0,1))
    quantity = np.asarray(quantity).reshape(-1,1)
    quantity = scaler.fit_transform(quantity).ravel()
    return quantity, float(scaler.data_min_[0]), float(scaler.data_max_[0])

def preprocess(events, features, transforms, constraint):
    features = [f for f in features if f in events.fields]
    if constraint not in features: transforms[constraint] = 'minmax'
    features += [constraint]

    scaling = {}
    for feature in features:
        quantity = events[feature]
        if transforms[feature] == "log":
            quantity = log_transform(quantity)
        # NOTE: Every transform, incl. none, ends in minmax
        events[feature], lo, hi = minmax_transform(quantity)
        scaling[feature] = {"transform": transforms[feature], "min": lo, "max": hi}

    return events, scaling

def normalize_weights(events):
    events["weight"] = events.weight / sum(events.weight)
    return events

def makeDataLoaders(events, features, constraint, batch_size):
    labl_str = np.asarray(events.label).astype(str) # ['1' '1' '1' ... '0' '0' '0']
    proc_str = np.asarray(events.proc_idx).astype(str) # ['0' '0' '0' ... '1' '1' '1']
    stratify_key = np.char.add(np.char.add(labl_str, "_"), proc_str) # '{label}_{process}'
    unique_keys, unique_counts = np.unique(stratify_key, return_counts=True) # ['0_0'(bkg, DY) '0_1'(bkg, tt) '1_0'(sig, c2v1p5_c3_1p0)], [463 88812 410]
    top_key = {k[0]: k for k in unique_keys[np.argsort(unique_counts)]} # last wins: most populous key per label
    rare = np.isin(stratify_key, unique_keys[unique_counts < 2]) # single-event processes can't be split
    stratify_key[rare] = [top_key[l] for l in labl_str[rare]]
    all_indices = np.arange(len(events)) # [    0     1     2 ... 89682 89683 89684] len(events)=89685

    train_idx, valid_idx = train_test_split(
        all_indices, 
        test_size = 0.2,   # nValid / (nTrain + nValid)
        random_state = 42, # make it deterministic
        stratify = stratify_key
    )

    print(f"Feature column order: {dict(enumerate(features))}")
    feature_matrix = np.column_stack([ak.to_numpy(events[f]) for f in features]).astype(np.float32, copy=False) # row:events, col:features. Col order matches features order.
    constraint_values = np.asarray(events[constraint], dtype=np.float32).reshape(-1, 1)
    labels  = np.asarray(events.label , dtype=np.float32)
    weights = np.asarray(events.weight, dtype=np.float32)

    train_loader = get_dataloader(
        dnn_input_data  = torch.from_numpy(feature_matrix[train_idx]),
        constraint_data = torch.from_numpy(constraint_values[train_idx]),
        labels          = torch.from_numpy(labels[train_idx]),
        weights         = torch.from_numpy(weights[train_idx]),
        batch_size      = batch_size,
        use_sampler     = False
    )
    valid_loader = get_dataloader(
        dnn_input_data  = torch.from_numpy(feature_matrix[valid_idx]),
        constraint_data = torch.from_numpy(constraint_values[valid_idx]),
        labels          = torch.from_numpy(labels[valid_idx]),
        weights         = torch.from_numpy(weights[valid_idx]),
        batch_size      = batch_size,
        use_sampler     = False,
        is_validation   = True
    )
    return train_loader, valid_loader, train_idx, valid_idx

def saveParquet(events, procNames, train_idx, valid_idx, constraint, path):
    """ Rows ordered [train, valid] so they align with ConcatDataset(train.pt, valid.pt), letting the scan just assign ABCDscore. """
    order = np.concatenate([train_idx, valid_idx])
    labels = np.asarray(events.label)[order]
    procs  = np.asarray(events.proc_idx)[order]
    df = pd.DataFrame({
        "VBSscore"     : np.asarray(events.rawVBS   , dtype=np.float32)[order],
        "weight"       : np.asarray(events.rawWeight, dtype=np.float64)[order],
        "norm_VBSscore": np.asarray(events[constraint], dtype=np.float32)[order], # min-max scaled; matches disco in the .pt
        "norm_weight"  : np.asarray(events.weight   , dtype=np.float64)[order],   # class-normalized; sums to 1 per label
        "label"   : labels.astype(np.int8),
        "process" : [procNames[int(l)][p] for l, p in zip(labels, procs)], # proc_idx restarts at 0 per sig/bkg
        "split"   : np.where(np.arange(len(order)) < len(train_idx), "train", "valid"),
    })
    df.to_parquet(path, index=False)
    print(f"\033[1;32mSaved abcd/{path}!\033[0m")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(epilog="Ex) python makeDataset_parq.py --config single/run2_2L_1FJ.yaml --flavor single")
    parser.add_argument('-c', "--config"    , required=True      , help="Path to YAML config")
    parser.add_argument('-f', "--flavor"    , default="single"   , choices=["single", "double"], help="Training flavor: single (one output) or double (two outputs). Later prob should put this in config.")
    parser.add_argument('-s', "--signal"    , default="c2v1p5_c3_1p0", choices=["c2v1p0_c3_1p0", "c2v1p0_c3_10p0", "c2v1p5_c3_1p0"], help="Signal coupling point (c2v1p0_c3_1p0 is SM)")
    args = parser.parse_args()

    # —————————— Load config ————————————————————————————————————————————————————————————
    print(f"\n—————————— Processing {Path(args.config).stem} ——————————")
    cfg = load_config(args.config)
    TrainingFeatures, FeatureTransforms, constraint = parse_training_features(cfg, args.flavor)
    derived_vars_cfg, DerivedVars, input_for_derivedVars, TrainingFeatures, FeatureTransforms = parse_derived_vars(cfg, TrainingFeatures, FeatureTransforms)

    # —————————— Loading samples ———————————————————————————————————————————
    VarsToLoad = TrainingFeatures + [constraint] + input_for_derivedVars + cfg.get("extra_vars", [])
    VarsToLoad = list(dict.fromkeys(f for f in VarsToLoad if f not in DerivedVars))
    new_load_start = time.time()
    parq_dir  = f"{PARQ_PATH}/{Path(cfg['local_base_path']).name}" # e.g. run2_1FJ
    sig_paths = local_paths_from_json(cfg["sample_json"], parq_dir, "sig", args.signal)
    bkg_paths = local_paths_from_json(cfg["sample_json"], parq_dir, "bkg")
    sig_events, sig_procs = LoadEvents(sig_paths, VarsToLoad, 1)
    bkg_events, bkg_procs = LoadEvents(bkg_paths, VarsToLoad, 0)
    print(f"Loaded sig+bkg in {time.time()-new_load_start:.1f} s")

    # —————————— Preselection cut ———————————————————————————————————————————
    i_sigW = sum(sig_events.weight) ; i_bkgW = sum(bkg_events.weight)
    print(f"\033[1mBefore preselection\033[0m: sig:{len(sig_events.weight)} (weighted: {i_sigW:.3f})  bkg:{len(bkg_events.weight)} (weighted: {i_bkgW:.3f})")
    if cfg.get("preselection"):
        sig_events = evaluate_preselection(sig_events, cfg["preselection"])
        bkg_events = evaluate_preselection(bkg_events, cfg["preselection"])
    sig_percent = round(sum(sig_events.weight)/i_sigW*100, 1) ; bkg_percent = round(sum(bkg_events.weight)/i_bkgW*100, 1)
    print(f"\033[1mAfter preselection\033[0m : sig:{len(sig_events)} (weighted: {sum(sig_events.weight):.3f}, {sig_percent}% survived)  bkg:{len(bkg_events)} (weighted: {sum(bkg_events.weight):.3f}, {bkg_percent}% survived)")

    for events in (sig_events, bkg_events): # kept for the ABCD scan, which needs real yields and uncscaled cut values
        events["rawWeight"] = events.weight
        events["rawVBS"]    = events[constraint]
    sig_events = normalize_weights(sig_events)
    bkg_events = normalize_weights(bkg_events)

    # —————————— Preprocess data ———————————————————————————————————————————
    data = ak.concatenate([sig_events, bkg_events])
    data = apply_derived_vars(data, derived_vars_cfg)
    data, scaling = preprocess(data, TrainingFeatures, FeatureTransforms, constraint)
    inputs = {f: scaling[f] for f in TrainingFeatures} # DNN input columns in order; doTraining bundles it into the model checkpoint

    # —————————— Make datasets ———————————————————————————————————————————
    train_loader, valid_loader, train_idx, valid_idx = makeDataLoaders(
        events     = data,
        features   = TrainingFeatures,
        constraint = constraint,
        batch_size = cfg.get("batch_size", 4096),
    )
    dataset_dir = f"{time.strftime('%y%m%d')}_dataset/{args.signal}"
    os.makedirs(dataset_dir, exist_ok=True)
    for split, loader in (("train", train_loader), ("valid", valid_loader)):
        datasetPath = f"{dataset_dir}/{Path(args.config).stem}_{split}.pt"
        loader.dataset.inputs = inputs
        torch.save(loader.dataset, datasetPath)
        print(f"\033[1;32mSaved abcd/{datasetPath}!\033[0m")
    saveParquet(data, [bkg_procs, sig_procs], train_idx, valid_idx, constraint, f"{dataset_dir}/{Path(args.config).stem}.parq")