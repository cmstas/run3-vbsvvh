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
## How the fitted JMS/JMR (theta) are converted to s, f and sigma_rel

The fit does not measure s and f directly. It fits two shape nuisances of the true-W (tp2) mass template,
theta_JMS = `tp2jms` and theta_JMR = `tp2jmr` (the values in `JMSR_correctionlib_etainclusive.json`), built in the
Wqq SF code (`make2DTemplates.C`) as:

- JMS: `jmsUp` = template of m * 1.05, `jmsDown` = template of m * 0.95;
- JMR: `jmrUp` = template of m * G, G ~ N(1, 0.10) (seeded per lumi block and event), `jmrDown` = 2 * nominal - `jmrUp`
  (mirror, since a smear can only widen; negative bins clamped, normalisation kept at 2 N_nom - N_up).

combine morphs the template linearly in each bin (`shapeU`, symmetric up/down):

    T(theta) = T_nom + theta * (T_up - T_nom)        for theta >= 0
    T(theta) = T_nom + |theta| * (T_down - T_nom)    for theta < 0

So theta = +1 means "5 % scale" or "10 % extra smear", but a fractional theta is a bin-by-bin mix of two templates,
not a scale of 1 + 0.05 theta or a smear of 0.10 theta. For example the 10 % smear adds ~100-128 GeV^2 to the variance
of an 80 GeV peak (not (0.1 * 80)^2 = 64), and the mirrored `jmrDown` narrows the core less than `jmrUp` widens it.
The conversion therefore measures what the fitted template actually does to the peak, per era / WP / pt bin, in the
collector (`wqq_sf` `get_SF_summary.py --export_json`, `compute_jmsr_rdf`):

1. Take the fit's tp2 templates (pass + fail summed, 5-GeV bins, 50-220 GeV): T_nom and the four variations.
2. Find the Gaussian core of T_nom: fit a Gaussian (a parabola to log(counts), weighted by sqrt(counts)), keep the
   bins within mean +- 2 sigma, refit, iterate. This gives mu_0, sigma_0 and a fixed set of core bins (~7-8 bins).
3. Build T_JMS = T(theta_JMS) from `jmsUp`/`jmsDown` and T_JMR = T(theta_JMR) from `jmrUp`/`jmrDown` with the formula
   above (each with only its own nuisance moved), and fit the same Gaussian **on the same core bins** as the nominal:
   mean mu_JMS of T_JMS, width sigma_JMR of T_JMR.
4. Then

       s         = mu_JMS    / mu_0       (JMS_scale)
       f         = sigma_JMR / sigma_0    (JMR_factor)
       sigma_rel = sigma_0   / mu_0       (JMR_sigmaRel, MC relative mass resolution of true W jets)

These are exactly the s, f, sigma_rel used in the equations above: s scales the mass so that the MC peak position
moves as the fitted template does; f widens (f > 1) or narrows (f < 1) the MC peak by the same factor as the fitted
template, and sigma_rel sets the size of the random smear for jets without a gen match (sqrt(f^2 - 1) * sigma_rel is
the extra relative smear that turns a width sigma_rel into f * sigma_rel).

(Reusing the nominal's core bins matters: re-choosing the window for each template can move it by one 5-GeV bin and
change f by several percent between nearly identical templates.)

Examples at WP 0.78:

| era, pt [GeV] | theta_JMS | theta_JMR | s | f | sigma_rel |
|---|---|---|---|---|---|
| 2018, 300-400 | -0.200 | +0.001 | 0.991 | 1.000 | 0.101 |
| 2018, 400-800 | -0.040 | +0.228 | 0.998 | 1.076 | 0.098 |
| 2024, 300-400 | -0.211 | +0.159 | 0.990 | 1.043 | 0.110 |
| 2017, 250-300 | -0.306 | -0.124 | 0.986 | 0.970 | 0.108 |

## `JMSR_correctionlib_etainclusive.json` (raw fit nuisances)

Corrections `JMS` and `JMR` = post-fit `tp2jms` / `tp2jmr`:
- JMS: theta = +1 corresponds to a mass scale of 1.05 (theta = -1 to 0.95).
- JMR: theta = +1 corresponds to an extra Gaussian smearing with sigma = 10%; the Down template is the
  symmetrized mirror (2*nominal - up), so negative theta is constrained.

Source: `wqq_sf/ParticleNetSF-WqqSF-genmatch/single_pass_results/particlenetmd_sf_15oct_JMSR_x2p_msdjec_wpgrid/SFsummary_jsons/`
(`JMSR_rdf_correctionlib.json`, `JMSR_correctionlib.json`; fits 2026-10-15, fit code 57651cb), same candle as `../WqqSF`
(new gen matching + soft-drop JEC).
