# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WPs (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`.
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.
- Soft-drop mass in the fit, 1-lep ttbar tag-and-probe (single-pass pass/fail).
- Candle: 1lep_WSF producer with the new fat-jet gen matching (`fj_1_genMatch`, true W = 1/2) and the soft-drop mass
  from JEC'd subjets (`feat/msd-subjet-jec`), V+jets samples after PR #127.
- Source: `wqq_sf/ParticleNetSF-WqqSF-genmatch/single_pass_results/particlenetmd_sf_14oct_msdjec_wp0p78_etaincl/SFsummary_jsons/Wqq_SF_correctionlib.json`
  (fits 2026-10-14/15, all 126 era x WP x pt bins).
