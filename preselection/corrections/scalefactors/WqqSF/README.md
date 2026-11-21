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
- pT bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside. The fit uses exactly these bins.

## Measurement

- **Method:** 1-lepton ttbar tag-and-probe, single-pass pass/fail fit of the soft-drop mass (`fj_1_sdmass`, 50–220 GeV). Same datacard model and fit code (wqq_sf f898ec4) at every WP.
- **"Other" normalisation, lnU kappa = 4:** the "other" (non-ttbar: W+jets, diboson, ...) MC has a free normalisation with a range of x1/4 to x4. With kappa = 2 it sits at the x2 boundary in almost every fit (20 of 21 at WP 0.78), so the result would depend on an arbitrary boundary. The fit needs the larger normalisation because the MC underpredicts the data at the low-mass edge of the fit window (50–65 GeV, before the fit): in the fail region by up to 14 %, and for jets passing the WP by up to 76 % at 400–800 GeV. The ttbar processes share one normalisation fixed by the bulk of the distribution, and "other" is 1–6 % of the MC, so it is the only process that can absorb this. At kappa = 4 "other" comes out at x1.4–x4.0 and is at the limit in 2 of the 126 fits (WP 0.60 2016postVFP 400–800, WP 0.40 2025 250–300 GeV). The JMS/JMR in `../JMSR` come from the X2p-mass fit on the same candle and MC, with kappa = 6 (see there).
- **Candle:** the 1lep_WSF RDF channel, with the fat-jet gen matching (`fj_1_genMatch`, true W = 1 or 2) and the soft-drop mass from JEC'd subjets (`feat/msd-subjet-jec`); V+jets samples after PR #127.
- **MC without overlapping samples:** the candle merge drops the duplicates in run3-vbsvvh's bkg/1lep_1FJ base list (inclusive WW/ZZ, split WZ, inclusive ttW/ttZ). W→lν: Run 2 uses `WJetsToLNu_HT-*` (8 HT bins, no inclusive sample); **Run 3 (2024, 2025) uses the LO `WtoLNu-4Jets` HT × MLNu samples (12 bins), not the NLO `WtoLNu-2Jets` PTLNu**.
- **Sources:** `wqq_sf/ParticleNetSF-WqqSF-ptbinfix/single_pass_results/<campaign>/SFsummary_jsons/` (`Wqq_SF_correctionlib.json`, `Wqq_effMC_correctionlib.json`). Run 2 and `r2comb` from `particlenetmd_sf_14nov_ptfix_noov_Wnlo_msdjec_wpgrid_etaincl_otherlnU4`; 2024 and 2025 from `particlenetmd_sf_14nov_ptfix_noov_Wlo_msdjec_wpgrid_etaincl_otherlnU4` (joined with `wqq_sf/small_tools/splice_sf_jsons.py`).
- **Fit quality:** every fit of both campaigns checked (all eras, WPs, pT bins): no failed fit, no SF uncertainty at a range edge.
- **Statistical limits:** at WP 0.96 / 0.98 the single Run 2 years at 400–800 GeV fluctuate by more than the WP trend (e.g. 2016preVFP 0.803 / 0.912 / 0.726 at WP 0.90 / 0.96 / 0.98, uncertainties 0.10–0.21). `r2comb` and Run 3 fall smoothly with the WP.

## Changes in this version (2026-11-15)

Two changes against the files previously in this directory (from the 15oct campaigns, also on `main`):

1. **Exact pT bins.** The fit's pT projection used to include the 16.7 GeV slice above each upper bin edge (250–316.7 / 300–416.7 / 400–816.7 GeV). Fixed in wqq_sf 7737f2d. This moves mainly 250–300 GeV, which lost the 300–316.7 GeV slice (about 1/4 of its events).
2. **"Other" normalisation kappa 2 -> 4** (see Measurement).

Against the kappa-2 set with exact pT bins (`sf-corrlib-json` 6f87a3399), i.e. the effect of change 2 alone: mean |shift| 0.013. The SFs move down at 250–300 GeV (mean shift −0.017, mean |shift| 0.019) and little at 300–400 / 400–800 GeV (mean |shift| 0.010 each). The largest shift is −0.076 (2016postVFP WP 0.98 250–300, 0.344 -> 0.268). 5 of 126 points move by more than 1 sigma, all at 250–300 GeV (2025 WP 0.40 / 0.60 −2.1 sigma, 2025 WP 0.78 −1.5, 2018 WP 0.60 −1.3, 2018 WP 0.78 −1.1). At WP 0.78 no SF moves by more than 0.047, e.g. 2018 0.791/0.783/0.846 -> 0.746/0.775/0.845, 2024 0.740/0.784/0.760 -> 0.723/0.769/0.755, 2025 0.793/0.831/0.803 -> 0.747/0.816/0.801 (250–300 / 300–400 / 400–800).

Against the previous files here (both changes): 250–300 GeV mean shift −0.028 (14 of 42 points beyond 1 sigma), 300–400 / 400–800 GeV mean |shift| 0.009 / 0.008 (none beyond 1 sigma). The largest shift is −0.102 (2017 WP 0.96 250–300, 0.567 -> 0.465), and at WP 0.78 nothing moves by more than 0.047. `Wqq_effMC` changes only through the pT bins (pre-fit, independent of kappa).
