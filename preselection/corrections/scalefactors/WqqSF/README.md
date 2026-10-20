# GloParT-3 Wqq tagging SF (correctionlib)

The correction is named `particleNet_Wqq_SF` for historical reasons; the tagger is GloParT-3.

`Wqq_SF_correctionlib_etainclusive.json`, correction `particleNet_Wqq_SF`.

- Inputs: `era` (string), `WP` (string), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real), `systematic` (`nom`/`up`/`down`, absolute SF values).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb` (Run 2 combined fit). No 2022/2023 yet.
- WP (score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD) from `globalParT3`): `0.78` only in this version (the
  previous version, d07a64def, had `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98` from a soft-drop-mass fit).
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.
- GloParT-3 X2p mass in the fit (`FatJet_globalParT3_massCorrX2p * FatJet_rawmass`), the same fit that measures the
  JMS/JMR in `../JMSR`; 1-lep ttbar tag-and-probe (single-pass pass/fail). The free normalisation of the "other"
  (non-ttbar) MC has a lnU range of kappa = 6 (was 2; the kappa scan showed "other" needs x2.9-x4.4).
- Candle: 1lep_WSF producer with the new fat-jet gen matching (`fj_1_genMatch`, true W = 1/2) and the soft-drop mass
  from JEC'd subjets (`feat/msd-subjet-jec`), V+jets samples after PR #127.
- `Wqq_effMC_correctionlib_etainclusive.json`, correction `particleNet_Wqq_effMC` (inputs `era`, `WP`, `eta`, `pt`): candle MC
  efficiency of true W jets, pre-fit tp2 pass / (pass + fail), same campaign; no uncertainty. Per-era files for RDF: `README_per_era.md`.
- MC without overlapping samples (2026-10-16): the candle merge drops the duplicates of run3-vbsvvh's bkg/1lep_1FJ base
  list (Run 3 W->lnu: NLO `WtoLNu-2Jets` PTLNu only, no LO `WtoLNu-4Jets`; inclusive WW/ZZ, split WZ, inclusive ttW/ttZ
  dropped). The earlier version of these files double counted Run 3 W+jets and the dibosons.
- Source: `wqq_sf/ParticleNetSF-WqqSF-nooverlap-genmatch/single_pass_results/particlenetmd_sf_16oct_noov_Wnlo_JMSR_x2p_msdjec_wp0p78_otherlnU6/SFsummary_jsons/Wqq_SF_correctionlib.json`
  (fits 2026-10-16/19, all 21 era x pt bins checked: no failed fit, no range-edge uncertainty). vs the previous
  version (soft-drop mass, kappa 2) at WP 0.78: 400-800 GeV up by 0.03-0.13 in every era, 250-400 GeV within
  -0.09..+0.04 (e.g. 2018: 0.758/0.785/0.844 -> 0.779/0.810/0.877, 2024: 0.740/0.779/0.764 -> 0.730/0.794/0.883).
