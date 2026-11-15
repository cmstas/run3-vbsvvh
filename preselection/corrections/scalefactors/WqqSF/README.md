# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WP (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`, all from the same fits.
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside. Since 2026-10-26 the fit uses exactly these bins:
  before (up to 1027838ce) the fit's pt projection included the 16.7 GeV slice above each upper edge, i.e. 250-316.7,
  300-416.7, 400-816.7 GeV (fixed in wqq_sf 7737f2d).
- Soft-drop mass in the fit (`fj_1_sdmass`, 50-220 GeV); 1-lep ttbar tag-and-probe (single-pass pass/fail). The JMS/JMR
  in `../JMSR` come from the GloParT-3 X2p-mass fit on the same candle and MC.
- "Other" normalisation: the "other" (non-ttbar: W+jets, diboson, ...) MC has a free normalisation with a lnU range of
  kappa = 4 (x1/4 to x4; `WSF_OTHER_LNU=4` in the fit code). With the default kappa = 2 it sits at the x2 boundary in
  almost every fit (20 of 21 at WP 0.78), so the result would depend on an arbitrary boundary. The fit needs the larger
  normalisation because the MC underpredicts the data at the low-mass edge of the fit window (50-65 GeV, before the fit):
  in the fail region by up to 14 %, and for jets passing the WP by up to 76 % at 400-800 GeV. The ttbar processes share one
  normalisation and are fixed by the bulk of the distribution, and "other" is 1-6 % of the MC, so it is the only
  process that can absorb this. At kappa = 4 "other" comes out at x1.4-x4.0 and is at the limit in 2 of the 126 fits
  (WP 0.60 2016postVFP 400-800, WP 0.40 2025 250-300 GeV).
- Candle: 1lep_WSF producer with the new fat-jet gen matching (`fj_1_genMatch`, true W = 1/2) and the soft-drop mass
  from JEC'd subjets (`feat/msd-subjet-jec`), V+jets samples after PR #127.
- `Wqq_effMC_correctionlib_etainclusive.json`, correction `particleNet_Wqq_effMC` (inputs `era`, `WP`, `eta`, `pt`): candle MC
  efficiency of true W jets, pre-fit tp2 pass / (pass + fail), same campaigns as the SFs; no uncertainty. Per-era files for RDF: `README_per_era.md`.
- MC without overlapping samples: the candle merge drops the duplicates of run3-vbsvvh's bkg/1lep_1FJ base list
  (inclusive WW/ZZ, split WZ, inclusive ttW/ttZ dropped). W->lnu: Run 2 = `WJetsToLNu_HT-*` (8 HT bins, no inclusive
  sample); **Run 3 (2024, 2025) = LO `WtoLNu-4Jets` HT x MLNu (12 bins), not the NLO `WtoLNu-2Jets` PTLNu**.
- Sources (`wqq_sf/ParticleNetSF-WqqSF-ptbinfix/single_pass_results/<campaign>/SFsummary_jsons/Wqq_SF_correctionlib.json`
  and `Wqq_effMC_correctionlib.json`, all six WPs): `2016preVFP`, `2016postVFP`, `2017`, `2018`, `r2comb` from
  `particlenetmd_sf_14nov_ptfix_noov_Wnlo_msdjec_wpgrid_etaincl_otherlnU4`; `2024`, `2025` from
  `particlenetmd_sf_14nov_ptfix_noov_Wlo_msdjec_wpgrid_etaincl_otherlnU4` (joined with `wqq_sf/small_tools/splice_sf_jsons.py`).
  Every fit of both campaigns checked (all eras, WPs, pt bins): no failed fit, no SF uncertainty at a range edge.
- Same datacard model and fit code (wqq_sf f898ec4: 7737f2d pt-bin fix + `WSF_OTHER_LNU`) at every WP. Statistically limited: at WP 0.96 / 0.98 the single Run 2
  years at 400-800 GeV fluctuate by more than the WP trend (e.g. 2016preVFP 0.803 / 0.912 / 0.726 at WP 0.90 / 0.96 /
  0.98, uncertainties 0.10-0.21); r2comb and Run 3 fall smoothly with the WP.
- vs 6f87a3399 (same candle, MC and fit code, "other" kappa = 2): the SFs move down at 250-300 GeV (mean shift -0.017, mean |shift|
  0.019), little at 300-400 / 400-800 GeV (mean |shift| 0.010 each). Largest: -0.076 (2016postVFP WP 0.98 250-300,
  0.344 -> 0.268); 5 of 126 points move by more than 1 sigma, all at 250-300 GeV (2025 WP 0.40 / 0.60 -2.1 sigma, 2025 WP
  0.78 -1.5, 2018 WP 0.60 -1.3, 2018 WP 0.78 -1.1). WP 0.78: <= 0.047, e.g. 2018 0.791/0.783/0.846 -> 0.746/0.775/0.845,
  2024 0.740/0.784/0.760 -> 0.723/0.769/0.755, 2025 0.793/0.831/0.803 -> 0.747/0.816/0.801. The MC efficiencies
  (`Wqq_effMC`) are unchanged (pre-fit).
