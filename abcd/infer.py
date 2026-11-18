"""
Assigns the ABCDNet score of each benchmark model to every event (sig, bkg and data) in PARQ_PATH/{tag}/*.parq,
and writes them back into the same parquet as ABCDscore_{signal}.
"""

import glob
from TrainingTools import *
from makeDataset_parq import log_transform, PARQ_PATH

SIGNALS = ["c2v1p0_c3_1p0", "c2v1p0_c3_10p0", "c2v1p5_c3_1p0"]

def loadModel(path):
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    model = makeModel(len(checkpoint["inputs"]))[0]
    model.load_state_dict(checkpoint["state_dict"])
    return model.eval(), checkpoint["inputs"] # eval: dropout off, batch-norm uses its stored running stats

@torch.no_grad()
def inferABCDscore(model, inputs, events):
    """ Scales each input like makeDataset_parq.preprocess; events with a None input (e.g. no FJ2) get NaN. """
    scale = lambda x, s: ((log_transform(x) if s["transform"] == "log" else x) - s["min"]) / (s["max"] - s["min"])
    features = np.column_stack([ak.to_numpy(ak.fill_none(scale(events[f], s), np.nan)) for f, s in inputs.items()]).astype(np.float32)
    device = next(model.parameters()).device
    return torch.cat([torch.sigmoid(model(batch.to(device))).cpu() for batch in torch.split(torch.from_numpy(features), 65536)]).reshape(-1).numpy()

def ScoredRowGroups(path, models):
    # one row group at a time keeps memory bounded
    for i in range(ak.metadata_from_parquet(path)["num_row_groups"]):
        events = ak.from_parquet(path, row_groups=[i])
        for signal, (model, inputs) in models.items():
            events[f"ABCDscore_{signal}"] = inferABCDscore(model, inputs, events)
        yield events

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tag", default='run2_1FJ', choices=["run2_1FJ", "run2_2FJ", "run3_1FJ", "run3_2FJ"])
    args = parser.parse_args()

    run, nFJ = args.tag.split("_")
    scriptPath = os.path.dirname(os.path.abspath(__file__))
    models = {signal: loadModel(f"{scriptPath}/models/{signal}/{run}_2L_{nFJ}_best_model.pt") for signal in SIGNALS}

    for path in sorted(glob.glob(f"{PARQ_PATH}/{args.tag}/*.parq")):
        ak.to_parquet_row_groups(ScoredRowGroups(path, models), f"{path}.tmp")
        os.replace(f"{path}.tmp", path) # an interrupted write never corrupts the original parquet
        print(f"\033[1;32mSaved ABCDscore_* into {path}\033[0m", flush=True)
