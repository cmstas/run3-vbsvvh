# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WP (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`, all from the same fits (the 2026-10-21 version
  8af6b6b8b carried `0.78` only; its WP 0.78 values are unchanged here).
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.
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
- Sources (`wqq_sf/ParticleNetSF-WqqSF-nooverlap-genmatch/single_pass_results/<campaign>/SFsummary_jsons/Wqq_SF_correctionlib.json`
  and `Wqq_effMC_correctionlib.json`, all six WPs): `2016preVFP`, `2016postVFP`, `2017`, `2018`, `r2comb` from
  `particlenetmd_sf_15oct_noov_Wnlo_msdjec_wpgrid_etaincl`; `2024`, `2025` from `particlenetmd_sf_15oct_noov_Wlo_msdjec_wpgrid_etaincl`
  (joined with `wqq_sf/small_tools/splice_sf_jsons.py`). Every fit of both campaigns checked (all eras, WPs, pt bins):
  no failed fit, no SF / JMS / JMR uncertainty at a range edge.
- Same datacard model and fit code (8963b15) at every WP. Statistically limited: at WP 0.96 / 0.98 the single Run 2
  years at 400-800 GeV fluctuate by more than the WP trend (e.g. 2016preVFP 0.803 / 0.908 / 0.723 at WP 0.90 / 0.96 /
  0.98, uncertainties 0.10-0.21); r2comb and Run 3 fall smoothly with the WP.
- vs 116c67844 (X2p mass, kappa 6, NLO Run 3 W, WP 0.78): Run 2 and its WP 0.78 values are identical to
  d07a64def (same campaign); 2024 / 2025 differ from d07a64def by <= 0.02 (LO instead of NLO W). E.g. 2018:
  0.779/0.810/0.877 -> 0.758/0.785/0.844, 2024: 0.730/0.794/0.883 -> 0.747/0.785/0.755, 2025: 0.743/0.806/0.908 -> 0.787/0.823/0.802.
