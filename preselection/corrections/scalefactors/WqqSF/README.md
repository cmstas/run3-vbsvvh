# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WPs (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`.
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.
- Soft-drop mass in the fit, 1-lep ttbar tag-and-probe (single-pass pass/fail).
- Source: `wqq_sf/ParticleNetSF-WqqSF/single_pass_results/particlenetmd_sf_29sep_wpgrid_etapt/SFsummary_jsons/` (2026-09-30).
