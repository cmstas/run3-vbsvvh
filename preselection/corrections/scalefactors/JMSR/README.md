# JMS / JMR from the Wqq tag-and-probe fit (correctionlib)

`JMSR_correctionlib_etainclusive.json`, corrections `JMS` and `JMR`.

- Inputs: `era` (string), `WP` (string, only `0.78`, same GloParT-3 Wqq score as `../WqqSF`), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real).
- Eras: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`, plus `r2comb`. No 2022/2023 yet.
- pt bins: [250, 300), [300, 400), [400, 800) GeV, clamped outside.
- Values are the post-fit nuisance parameters `tp2jms` / `tp2jmr` (central value only, no uncertainty):
  - JMS: theta = +1 corresponds to a mass scale of 1.05 (theta = -1 to 0.95).
  - JMR: theta = +1 corresponds to an extra Gaussian smearing with sigma = 10%; the Down template is the
    symmetrized mirror (2*nominal - up), so negative theta is constrained.
- Mass: GloParT-3 Generic-corrected mass (`globalParT3_massCorrGeneric * mass * (1 - rawFactor)`).
- Source: `wqq_sf/ParticleNetSF-WqqSF/single_pass_results/particlenetmd_sf_02oct_JMSR_generic_jmrsym/`
  (fits 2026-10-02). The |eta|-binned version is `SFsummary_jsons/JMSR_correctionlib.json` there.
