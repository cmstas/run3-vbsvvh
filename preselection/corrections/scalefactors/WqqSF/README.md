# GloParT-3 Wqq tagging scale factors

Data/MC scale factors for the GloParT-3 Wqq tagger, measured for true W jets (`FatJet_genMatch` 1 = prompt W, 2 = W from a top). They are cumulative: SF(WP) corrects the efficiency of score > WP, with score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3` (`FatJet_Wqq`).

## Files

| File | Role |
|---|---|
| `Wqq_SF_correctionlib_etainclusive.json` | Original output of the SF fit framework: the SFs for every era in one file. Source for the per-era files. |
| `Wqq_effMC_correctionlib_etainclusive.json` | Original output of the SF fit framework: the candle MC efficiency of true W jets for score > WP. Source for the per-era files. |
| `make_wqq_sf_jsons.py` | Turns the two files above into the per-era files. |
| `Wqq_SF_2016preVFP.json`, `_2016postVFP`, `_2017`, `_2018`, `_2024Prompt`, `_2025` | The per-era files the RDF preselection reads (`src/wqq_sf.cpp`), one per RDF `year`, each holding that era's own fit. |
| `Wqq_SF_Run2Comb.json` | The combined Run 2 fit (`r2comb`) in the same per-era format. Kept for reference; not read by RDF. |

The two provider files are kept unchanged as delivered. Analyses never read them directly: RDF uses the per-era files, which differ from the provider files in three ways (details in [`README_per_era.md`](README_per_era.md)):

- **WP** is a real number: any WP in [0.40, 0.98] is valid, linearly interpolated between the measured points. In the provider files it is a text category limited to the six measured values.
- **pT above 800 GeV** gets the 400–800 GeV value with 3× its uncertainty. In the provider files it is simply clamped to the 400–800 GeV bin.
- **The efficiency** is included in each era file as `Wqq_eff_MC`.

When new SFs arrive, replace the two provider files and regenerate:

```bash
python3 make_wqq_sf_jsons.py --eff Wqq_effMC_correctionlib_etainclusive.json
```

## Provider files

- `Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF` (named for historical reasons; the tagger is GloParT-3). Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- `Wqq_effMC_correctionlib_etainclusive.json`, correction `particleNet_Wqq_effMC`. Inputs: `era`, `WP`, `eta`, `pt`. Pre-fit tp2 yield pass / (pass + fail) from the same fits as the SFs; no uncertainty.
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, and `r2comb` (combined Run 2 fit). No 2022/2023 yet.
- WPs: `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`, all from the same fits.
- pT bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.

## Measurement

- **Method:** 1-lepton ttbar tag-and-probe, single-pass pass/fail fit of the soft-drop mass (`fj_1_sdmass`, 50–220 GeV). The normalisation of the "other" (non-ttbar) MC floats freely with a lnU range of kappa = 2. Same datacard model and fit code (8963b15) at every WP.
- **Candle:** the 1lep_WSF RDF channel, with the fat-jet gen matching (`fj_1_genMatch`, true W = 1 or 2) and the soft-drop mass from JEC'd subjets (`feat/msd-subjet-jec`); V+jets samples after PR #127.
- **MC without overlapping samples:** the candle merge drops the duplicates in run3-vbsvvh's bkg/1lep_1FJ base list (inclusive WW/ZZ, split WZ, inclusive ttW/ttZ). W→lν: Run 2 uses `WJetsToLNu_HT-*` (8 HT bins, no inclusive sample); **Run 3 (2024, 2025) uses the LO `WtoLNu-4Jets` HT × MLNu samples (12 bins), not the NLO `WtoLNu-2Jets` PTLNu**.
- **Sources:** `wqq_sf/ParticleNetSF-WqqSF-nooverlap-genmatch/single_pass_results/<campaign>/SFsummary_jsons/` (`Wqq_SF_correctionlib.json`, `Wqq_effMC_correctionlib.json`). Run 2 and `r2comb` from `particlenetmd_sf_15oct_noov_Wnlo_msdjec_wpgrid_etaincl`; 2024 and 2025 from `particlenetmd_sf_15oct_noov_Wlo_msdjec_wpgrid_etaincl` (joined with `wqq_sf/small_tools/splice_sf_jsons.py`).
- **Fit quality:** every fit of both campaigns checked (all eras, WPs, pT bins): no failed fit, no SF / JMS / JMR uncertainty at a range edge.
- **Statistical limits:** at WP 0.96 / 0.98 the single Run 2 years at 400–800 GeV fluctuate by more than the WP trend (e.g. 2016preVFP 0.803 / 0.908 / 0.723 at WP 0.90 / 0.96 / 0.98, uncertainties 0.10–0.21). `r2comb` and Run 3 fall smoothly with the WP.
