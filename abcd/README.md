# VBS ABCDNet

DNN (+ optional BDT constraint) training with DisCo decorrelation, followed by
ABCD cut optimization and closure tests.

## Layout

| File | Purpose |
| --- | --- |
| `main.py` | CLI entry point: train / infer / data-infer orchestration |
| `config.py` | YAML config parsing into a resolved `RunConfig` |
| `root_io.py` | Reading ROOT ntuples into flat numpy columns |
| `predictions.py` | Reading/writing the Parquet prediction files |
| `preprocessing.py` | Derived vars, min-max scaling (fit + frozen apply), class weights |
| `training.py` | Dataloaders, Lightning training run, run-artifact collection |
| `inference.py` | Checkpoint scoring, prediction files, permutation importance |
| `checkpoints.py` | Locating version dirs / checkpoints / saved scaler params |
| `plots.py` | All diagnostic plots |
| `model.py`, `dataloader.py`, `utilities.py` | Lightning module, dataset, DisCo loss |
| `bdt.py` | XGBoost BDT used as the decorrelation constraint |
| `channels.py` | Channel definitions (taggers + cut combination) for 1fj/2fj/3fj |
| `run_analysis.py` | ABCD cut optimization for one channel (replaces `{1,2,3}fj_analysis.py`) |
| `closure.py` | ABCD closure systematic from the data control regions |

## Run

```bash
# Train + run inference with the best and last checkpoints
python3 main.py --config single/config_1lep_boosted_run3.yaml --flavor single

# Inference only (latest checkpoint, or --checkpoint <path>)
python3 main.py --config ... --infer
python3 main.py --config ... --infer --data   # on real data

# Cut optimization on the prediction files
python3 run_analysis.py -c 2fj -i <predictions.parquet> -d <predictions_data.parquet> -n 5

# Closure systematic
python3 closure.py -i <predictions_data.parquet> --channel 2fj --regions <regions.yaml> --region all
```

The batch drivers `train.sh`, `infer.sh`, `analyze.sh`, and `closure.sh` run the
above over all configured channels. Monitor training with TensorBoard:
`tensorboard --logdir <output_dir>`.

## Config

Required keys:

- `sig_path`, `bkg_path` (and `data_path` for `--data`), plus optional
  `sig_base_path` / `bkg_base_path` / `data_base_path` (or a common
  `input_base_path`) that the path globs are resolved against
- `sig_infer_path` / `bkg_infer_path` (optional): files to score in `--infer`
  mode when they should differ from the training files (e.g. a larger skim).
  Each takes an optional `sig_infer_base_path` / `bkg_infer_base_path` (else the
  training base path). If unset, `--infer` reuses the training `sig_path` /
  `bkg_path`.
- `training_features`: list of branch names, or `name: transform` mappings
  (transforms: `none`, `log`; everything is min-max scaled after the transform)
- `constraint_var`: decorrelation variable (ignored when
  `use_bdt_as_constraint` selects the BDT score)
- `output`: output directory

Optional keys (defaults in parentheses):

- `extra_vars`: extra branches to carry into the prediction files (`weight` is
  always included)
- `preselection`: numpy expression over the loaded branches
- `derived_vars` (`{}`): mapping of `new_var: expression`, computed before
  preprocessing, e.g.
  `delta_eta_abs: "np.abs(boosted_h_candidate_eta - boosted_v_candidate_eta)"`
- `auto_include_derived_vars` (`false`) and
  `auto_include_derived_vars_transform` (`none`)
- `flavor` (`single`): `single` for classical DisCo against `constraint_var`,
  `double` for double DisCo (the constraint var becomes a regular input;
  `constraint_as_feature_transform` sets its transform, default `none`)

The data-driven QCD GloParT score correction now lives in the preselection
(`preselection/src/utils.cpp` `applyQCDScoreResampling`, applied in the 0lep_3FJ
channel), so the corrected scores are already baked into the ntuples the MVA
reads — there is no longer a `qcd_resampling` config block here.

Training hyperparameters:

- `learning_rate` (`1e-3`), `batch_size` (`4096`), `n_epochs` (`100`)
- `architecture` (`[64, 32, 16]`), `use_batchnorm` (`true`), `dropout` (`0.0`)
- `bce_weight` (`1.0`), `disco_lambda` (`0.0`)
- `weight_decay` (`1e-2`), `label_smoothing` (`0.0`)
- `use_lr_scheduler` (`true`), `lr_scheduler_patience` (`40`),
  `lr_scheduler_factor` (`0.5`), `lr_scheduler_min_lr` (`1e-6`)
- `early_stopping_patience` (`40`), `early_stopping_min_delta` (`1e-4`)
- `check_val_every_n_epoch` (`1`), `devices` (`[0]`), `io_workers` (`min(8, cpus)`)

BDT constraint (enabled by a non-empty `bdt_features` list):

- `bdt_features`: VBS branches the BDT is trained on
- `use_bdt_as_constraint` (`true`): use the BDT score as the DNN constraint
- `bdt`: XGBoost hyperparameters (`n_estimators`, `max_depth`, `learning_rate`,
  `subsample`, `colsample_bytree`, `min_child_weight`, `gamma`, `reg_lambda`,
  `reg_alpha`, `tree_method`, `n_jobs`, `scale_pos_weight`,
  `early_stopping_rounds`)
