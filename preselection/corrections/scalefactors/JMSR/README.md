# JMS / JMR from the Wqq tag-and-probe fit (correctionlib)

Both files: inputs `era` (string), `WP` (string, `0.40`, `0.60`, `0.78`, `0.90`, `0.96`, `0.98`, same GloParT-3 Wqq
score as `../WqqSF`), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real). Eras `2016preVFP`, `2016postVFP`,
`2017`, `2018`, `2024`, `2025`, plus `r2comb`; no 2022/2023 yet. pt bins [250, 300), [300, 400), [400, 800) GeV,
clamped outside. Central values only, no uncertainty.

Mass: GloParT-3 X2p-corrected mass, `FatJet_globalParT3_massCorrX2p * FatJet_rawmass` (= `FatJet_massGloParT3` /
`fj_1_massGloParT3` of the RDF producer).

## `JMSR_rdf_correctionlib_etainclusive.json` (what to apply to MC)

Corrections `JMS_scale` (s), `JMR_factor` (f), `JMR_sigmaRel` (sigma_rel):
- scale: m' = s * m;
- resolution (hybrid): gen-matched jets m' = m_gen + f * (m - m_gen); unmatched jets m' = m * (1 + sqrt(f^2 - 1) *
  sigma_rel * N(0,1)) for f > 1.
Derived from the fitted nuisances on the fit's own true-W templates: s and f are the ratios of the Gaussian core mean /
width of the template morphed to the fitted theta vs the nominal (core bins from the nominal, +-2 sigma);
sigma_rel = core width / core mean of the nominal.

## `JMSR_correctionlib_etainclusive.json` (raw fit nuisances)

Corrections `JMS` and `JMR` = post-fit `tp2jms` / `tp2jmr`:
- JMS: theta = +1 corresponds to a mass scale of 1.05 (theta = -1 to 0.95).
- JMR: theta = +1 corresponds to an extra Gaussian smearing with sigma = 10%; the Down template is the
  symmetrized mirror (2*nominal - up), so negative theta is constrained.

Source: `wqq_sf/ParticleNetSF-WqqSF-genmatch/single_pass_results/particlenetmd_sf_15oct_JMSR_x2p_msdjec_wpgrid/SFsummary_jsons/`
(`JMSR_rdf_correctionlib.json`, `JMSR_correctionlib.json`; fits 2026-10-15, fit code 57651cb), same candle as `../WqqSF`
(new gen matching + soft-drop JEC).
