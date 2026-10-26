# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WP (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`, all from the same fits.
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside. Since 2026-10-26 the fit uses exactly these bins:
  before (up to 1027838ce) the fit's pt projection included the 16.7 GeV slice above each upper edge, i.e. 250-316.7,
  300-416.7, 400-816.7 GeV (fixed in wqq_sf 7737f2d).
- Soft-drop mass in the fit (`fj_1_sdmass`, 50-220 GeV); 1-lep ttbar tag-and-probe (single-pass pass/fail). The free
  normalisation of the "other" (non-ttbar) MC has a lnU range of kappa = 2. The JMS/JMR in `../JMSR` come from the
  GloParT-3 X2p-mass fit on the same candle and MC (as in d07a64def).
- Candle: 1lep_WSF producer with the new fat-jet gen matching (`fj_1_genMatch`, true W = 1/2) and the soft-drop mass
  from JEC'd subjets (`feat/msd-subjet-jec`), V+jets samples after PR #127.
- `Wqq_effMC_correctionlib_etainclusive.json`, correction `particleNet_Wqq_effMC` (inputs `era`, `WP`, `eta`, `pt`): candle MC
  efficiency of true W jets, pre-fit tp2 pass / (pass + fail), same campaigns as the SFs; no uncertainty. Per-era files for RDF: `README_per_era.md`.
- MC without overlapping samples: the candle merge drops the duplicates of run3-vbsvvh's bkg/1lep_1FJ base list
  (inclusive WW/ZZ, split WZ, inclusive ttW/ttZ dropped). W->lnu: Run 2 = `WJetsToLNu_HT-*` (8 HT bins, no inclusive
  sample); **Run 3 (2024, 2025) = LO `WtoLNu-4Jets` HT x MLNu (12 bins), not the NLO `WtoLNu-2Jets` PTLNu**.
- Sources (`wqq_sf/ParticleNetSF-WqqSF-ptbinfix/single_pass_results/<campaign>/SFsummary_jsons/Wqq_SF_correctionlib.json`
  and `Wqq_effMC_correctionlib.json`, all six WPs): `2016preVFP`, `2016postVFP`, `2017`, `2018`, `r2comb` from
  `particlenetmd_sf_25oct_ptfix_noov_Wnlo_msdjec_wpgrid_etaincl`; `2024`, `2025` from `particlenetmd_sf_25oct_ptfix_noov_Wlo_msdjec_wpgrid_etaincl`
  (joined with `wqq_sf/small_tools/splice_sf_jsons.py`). Every fit of both campaigns checked (all eras, WPs, pt bins):
  no failed fit, no SF / JMS / JMR uncertainty at a range edge.
- Same datacard model and fit code (7737f2d = 8963b15 + the pt-bin fix) at every WP. Statistically limited: at WP 0.96 / 0.98 the single Run 2
  years at 400-800 GeV fluctuate by more than the WP trend (e.g. 2016preVFP 0.803 / 0.908 / 0.723 at WP 0.90 / 0.96 /
  0.98, uncertainties 0.10-0.21); r2comb and Run 3 fall smoothly with the WP.
- vs 1027838ce (same campaigns and settings, old pt projection): only 250-300 GeV moves noticeably (it lost the
  300-316.7 GeV slice, ~1/4 of its events; its uncertainty grows by x1.23 on average): mean |shift| 0.029, max 0.104
  (2017 WP 0.96, 0.567 -> 0.463); 300-400 / 400-800 mean |shift| 0.007 / 0.004, all within 0.3 sigma. WP 0.78:
  <= 0.033 everywhere, e.g. 2018 0.758/0.785/0.844 -> 0.791/0.783/0.846, 2024 0.747/0.785/0.755 -> 0.740/0.784/0.760,
  2025 0.787/0.823/0.802 -> 0.793/0.831/0.803.
