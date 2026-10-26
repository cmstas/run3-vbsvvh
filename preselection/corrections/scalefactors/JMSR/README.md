# JMS / JMR from the Wqq tag-and-probe fit (correctionlib)

Both files: inputs `era` (string), `WP` (string, `0.78` only in this version; the previous one had `0.40`, `0.60`,
`0.78`, `0.90`, `0.96`, `0.98`; same GloParT-3 Wqq score as `../WqqSF`), `eta` (real, unused: inclusive in |eta| < 2.5), `pt` (real). Eras `2016preVFP`, `2016postVFP`,
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

### Core bins: what they are and how they are found

The true-W template is not a Gaussian over the whole 50-220 GeV range: it has a Gaussian-like peak around the W mass
and a long high-mass tail. Example, nominal tp2 template (pass + fail) of 2018, WP 0.78, pt 300-400 GeV, per 5-GeV bin
(bin centre: events):

    52:48  57:70  62:124  67:221  72:502  77:1226  82:2301  87:2718  92:2102  97:1320  102:821
    107:549  112:384  117:291  122:219  127:163  ...  182:6  ...  217:1

A single Gaussian fitted to all 34 bins gives sigma = 22.5 GeV, because the tail drags it wide; that says nothing about
the peak resolution. The **core bins** are the bins around the peak where the shape is Gaussian, here 70-105 GeV
(7 bins), where the fit gives sigma = 8.86 GeV.

The window depends on mu and sigma (bins with centres in mu +- 2 sigma), but mu and sigma come from the fit inside the
window, so `_core_window` iterates:

1. start with mu = centre of the highest bin, sigma = 10 GeV (a guess);
2. take the bins with centres in (mu - 2 sigma, mu + 2 sigma) and fit a Gaussian to those bins only, giving a new
   mu, sigma;
3. repeat step 2 with the new mu, sigma (10 times); once the window stops changing, the fit stops changing.

The Gaussian fit is a parabola fitted to log(counts), since the log of a Gaussian is a parabola
(ln N = c0 + c1 x + c2 x^2, so mu = -c1 / (2 c2), sigma = sqrt(-1 / (2 c2))), weighted by sqrt(counts). For the bin
above it settles at once, because the starting guess is already close:

    start:  mu = 87.5,  sigma = 10 (guess)
    iter 1: window 70-105 GeV (7 bins) -> mu = 88.09, sigma = 8.86
    iter 2: window 70-105 GeV (7 bins) -> mu = 88.09, sigma = 8.86   (unchanged from here on)

so mu_0 = 88.09 GeV, sigma_0 = 8.86 GeV, sigma_rel = 8.86 / 88.09 = 0.101. In bins where the first guess is further
off, the window moves for a step or two before it settles.

The window found on the nominal (here the 7 bins 70-105 GeV) is then **kept fixed**: T_JMS and T_JMR are fitted on
exactly these bins, s = mu_JMS / 88.09, f = sigma_JMR / 8.86. Re-finding the window for each morphed template could
make it gain or lose one 5-GeV bin at the edge and change f by several percent between nearly identical templates;
with a fixed window, s and f reflect only the change in shape.

Examples at WP 0.78:

| era, pt [GeV] | theta_JMS | theta_JMR | s | f | sigma_rel |
|---|---|---|---|---|---|
| 2018, 300-400 | -0.196 | +0.009 | 0.991 | 1.003 | 0.101 |
| 2018, 400-800 | -0.024 | +0.212 | 0.999 | 1.070 | 0.098 |
| 2024, 300-400 | -0.265 | +0.205 | 0.987 | 1.057 | 0.110 |
| 2017, 250-300 | -0.302 | -0.143 | 0.986 | 0.966 | 0.108 |

### Conversion code

Verbatim from `wqq_sf` `ParticleNetSF-WqqSF/get_SF_summary.py` (commit 57651cb), run by the collector with
`--export_json` for `mass_to_use = "globalpart"` campaigns. `m["templates"]` are the fit's pass and fail template ROOT
files of one era / WP / pt bin (histograms `tp2`, `tp2_tp2jmsUp`, `tp2_tp2jmsDown`, `tp2_tp2jmrUp`, `tp2_tp2jmrDown`),
`m["jms"]` / `m["jmr"]` the fitted `tp2jms` / `tp2jmr`.

```python
def _template_hist(paths, name):
    """Sum of histogram `name` over the pass and fail template files -> (contents, bin centres) as numpy arrays."""
    import numpy as np
    vals, centres = None, None
    for path in paths:
        f = ROOT.TFile.Open(path, "READ")
        h = f.Get(name) if f and not f.IsZombie() else None
        if not h:
            if f:
                f.Close()
            return None, None
        n = h.GetNbinsX()
        v = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
        c = np.array([h.GetBinCenter(i) for i in range(1, n + 1)])
        vals = v if vals is None else vals + v
        centres = c
        f.Close()
    return vals, centres

def _gauss_fit(vals, centres, sel):
    """Gaussian (mean, width) from a parabola fitted to log(counts) over the bins in sel, weighted by sqrt(counts).
    (None, None) if fewer than 3 bins or the parabola does not open downwards."""
    import numpy as np
    sel = sel & (vals > 0)
    if sel.sum() < 3:
        return None, None
    c2, c1, _c0 = np.polyfit(centres[sel], np.log(vals[sel]), 2, w=np.sqrt(vals[sel]))
    if c2 >= 0:
        return None, None
    return -c1 / (2.0 * c2), (-1.0 / (2.0 * c2)) ** 0.5

def _core_window(vals, centres, nsig=2.0, iters=10):
    """Bins of the Gaussian core of a peak: iterate the fit within mean +- nsig*width. Returns (sel, mean, width),
    sel None if the fit fails. Morphed templates are then fitted on the SAME bins: re-iterating the window on each
    template makes it jump by a bin (5-GeV bins, ~7 bins in the core) and the width ratio by several percent."""
    import numpy as np
    mu, sig = centres[int(np.argmax(vals))], 10.0
    sel = None
    for _ in range(iters):
        sel = (centres > mu - nsig * sig) & (centres < mu + nsig * sig)
        mu, sig = _gauss_fit(vals, centres, sel)
        if mu is None:
            return None, None, None
    return sel, mu, sig

def _morph(nom, up, down, theta):
    """The template combine fits for a shapeU nuisance at value theta: linear vertical interpolation, which is
    exact for symmetric up/down (down = 2*nom - up, as built for jmr) and beyond |theta| = 1. Clipped at 0."""
    import numpy as np
    shifted = nom + theta * (up - nom) if theta >= 0 else nom + (-theta) * (down - nom)
    return np.clip(shifted, 0.0, None)

def compute_jmsr_rdf(data):
    """Turns the fitted tp2jms / tp2jmr nuisance values into what RDF applies to a jet mass, measured on the fit's
    own true-W (tp2) templates (pass + fail summed), all fitted on the nominal's core bins (+-2 sigma):
        JMS_scale    = core mean of the JMS-morphed template / core mean of the nominal   (m' = s * m)
        JMR_factor   = core width of the JMR-morphed template / core width of the nominal (f, hybrid smearing)
        JMR_sigmaRel = core width / core mean of the nominal (MC relative resolution, unmatched-jet branch)
    Stored as data[key]['jms_scale'|'jmr_factor'|'jmr_sigma_rel']."""
    for key, m in data.items():
        if m.get("jms") is None or m.get("jmr") is None or not m.get("templates"):
            continue
        nom, c = _template_hist(m["templates"], "tp2")
        hists = {n: _template_hist(m["templates"], f"tp2_tp2{n}")[0] for n in ("jmsUp", "jmsDown", "jmrUp", "jmrDown")}
        if nom is None or any(h is None for h in hists.values()):
            continue
        win, mu0, s0 = _core_window(nom, c)
        if win is None:
            print(f"Warning: JMS/JMR core fit failed for {key}, no RDF values", file=sys.stderr)
            continue
        mu_s, _ = _gauss_fit(_morph(nom, hists["jmsUp"], hists["jmsDown"], m["jms"]), c, win)
        _, s_r = _gauss_fit(_morph(nom, hists["jmrUp"], hists["jmrDown"], m["jmr"]), c, win)
        if None in (mu_s, s_r):
            print(f"Warning: JMS/JMR core fit failed for {key}, no RDF values", file=sys.stderr)
            continue
        m["jms_scale"], m["jmr_factor"], m["jmr_sigma_rel"] = mu_s / mu0, s_r / s0, s0 / mu0
```

## `JMSR_correctionlib_etainclusive.json` (raw fit nuisances)

Corrections `JMS` and `JMR` = post-fit `tp2jms` / `tp2jmr`:
- JMS: theta = +1 corresponds to a mass scale of 1.05 (theta = -1 to 0.95).
- JMR: theta = +1 corresponds to an extra Gaussian smearing with sigma = 10%; the Down template is the
  symmetrized mirror (2*nominal - up), so negative theta is constrained.

Sources (`wqq_sf/ParticleNetSF-WqqSF-ptbinfix/single_pass_results/<campaign>/SFsummary_jsons/`
`JMSR_rdf_correctionlib.json`, `JMSR_correctionlib.json`, WP 0.78 kept; joined with `wqq_sf/small_tools/splice_sf_jsons.py`):
`2016preVFP`, `2016postVFP`, `2017`, `2018`, `r2comb` from `particlenetmd_sf_25oct_ptfix_noov_Wnlo_JMSR_x2p_msdjec_wpgrid`,
`2024`, `2025` from `particlenetmd_sf_25oct_ptfix_noov_Wlo_JMSR_x2p_msdjec_wpgrid` (LO `WtoLNu-4Jets` Run 3
W+jets instead of NLO). Both: X2p mass fit, free "other" normalisation with lnU kappa = 2, MC without overlapping samples,
same candle as `../WqqSF` (new gen matching + soft-drop JEC; the Wqq SFs there come from the soft-drop-mass fit on the
same MC). Every fit checked: no failed fit, deployed theta never at its range edge, theta +- 1 sigma inside [-1, 1].
Since 2026-10-26 the fit's pt bins are exactly [250, 300), [300, 400), [400, 800) GeV (before, up to 1027838ce, each bin
included the 16.7 GeV slice above its upper edge; fixed in wqq_sf 7737f2d). vs 1027838ce: s within 0.003, f within 0.023,
sigma_rel within 0.001 (largest changes at 250-300 GeV).
