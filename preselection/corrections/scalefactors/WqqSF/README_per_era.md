# Per-era Wqq SF files read by RDF

`Wqq_SF_<year>.json` are generated from the provider's file
`Wqq_SF_correctionlib_etainclusive.json` (described in `README.md`) by
`make_wqq_sf_jsons.py`, and are what `src/wqq_sf.cpp` reads.

| File | Content |
|---|---|
| `Wqq_SF_2016preVFP.json`, `_2016postVFP`, `_2017`, `_2018` | that era's own Run 2 fit |
| `Wqq_SF_2024Prompt.json`, `Wqq_SF_2025.json` | the 2024 and 2025 fits (2025: 2025 data, Summer24 MC relabelled) |
| `Wqq_SF_Run2Comb.json` | the combined Run 2 fit (`r2comb`), shipped alongside; not read by RDF |
| `Wqq_effMC_correctionlib_etainclusive.json` | provider's candle MC efficiency table, the `--eff` input |

The file name is the RDF `year` string, so RDF needs no era mapping.

## Corrections in each file

`Wqq_SF`: inputs `WP` (real), `pt` (real), `eta` (real), `systematic` (`nom`/`up`/`down`); output: the cumulative SF of score > WP, with score = (Xcs + Xqq/3)/(Xcs + Xqq/3 + QCD).

- **WP:** measured at 0.40, 0.60, 0.78, 0.90, 0.96, 0.98. The SF is piecewise linear between them (nom, up and down each interpolated), the measured points are reproduced exactly, and a WP outside [0.40, 0.98] is an error.
- **pt:** bins [250, 300, 400, 800, 13000), clamped. Above 800 GeV: the 400–800 nominal, with its up/down distances multiplied by 3.
- **eta:** unused (inclusive in |eta| < 2.5); kept so an |eta|-binned version can be dropped in later.

`Wqq_eff_MC`: inputs `WP`, `pt`; output: the candle MC efficiency of true W jets (tp2) to pass score > WP, built the same way. It's present only when `make_wqq_sf_jsons.py --eff` was given an efficiency table. RDF needs it for the `flags` use case (`applywqqsf.yaml`); `cut` works without it.

The current files carry it, from `Wqq_effMC_correctionlib_etainclusive.json` (correction `particleNet_Wqq_effMC`): the pre-fit tp2 yield pass / (pass + fail) of the same campaigns as the SFs (Run 2: `particlenetmd_sf_15oct_noov_Wnlo_msdjec_wpgrid_etaincl`, 2024/2025: `particlenetmd_sf_15oct_noov_Wlo_msdjec_wpgrid_etaincl`; overlap-free MC, true W = `fj_1_genMatch` 1 or 2, |eta| < 2.5, candle selection; W->lnu + jets has no true hadronic W, so the NLO/LO W choice does not change it), written by the provider's `get_SF_summary.py --export_json` next to the SF file. It has no uncertainty. When the SFs are remeasured, regenerate both files from the same campaign.

## Regenerate

```bash
python3 make_wqq_sf_jsons.py                    # SFs only
python3 make_wqq_sf_jsons.py --eff Wqq_effMC_correctionlib_etainclusive.json   # plus Wqq_eff_MC (current files)
```

The efficiency table has the same layout as the SF file (era → WP → pt bins), without the systematic level.
