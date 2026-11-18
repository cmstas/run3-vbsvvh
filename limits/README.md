# Expected limits for the 2L OS OF 1FJ and 2L OS SF 2FJ channels

Blind expected 95% CL upper limits on the signal strength `r`, computed from the datacards in `datacards/`, which a collaborator made.
`r` scales the `TotalSig_*` process of each card, i.e. the summed VBS VVH signal at that benchmark.

## Contents

| Path | What it is |
|---|---|
| `datacards/<benchmark>/R{2,3}_WWH_2L.txt` | 2L OS OF + 1 FJ cards (WWH), Run 2 and Run 3 |
| `datacards/<benchmark>/R{2,3}_WZH_ZZH_2L.txt` | 2L OS SF + 2 FJ cards (WZH and ZZH), Run 2 and Run 3 |
| `run_limits.sh` | Computes all limits (see below) |

The benchmarks are `vbs-hvv-c2v-1p0-c3-1p0` (SM), `vbs-hvv-c2v-1p5-c3-1p0` (c2V = 1.5) and `vbs-hvv-c2v-1p0-c3-10p0` (κλ = 10).

## How to reproduce

Requirements: `apptainer`, and `/cvmfs/unpacked.cern.ch` mounted. No CMSSW area is needed.

```bash
cd /users/jkil/run3-vbsvvh/limits
./run_limits.sh            # or: ./run_limits.sh /some/scratch/dir
```

What the script does:

1. Copies `datacards/` into a work directory (default `./work`). Combine writes its ROOT outputs there, so `datacards/` stays untouched. The work directory can be deleted afterwards.
2. Re-runs itself inside the `combine-standalone:latest` container from `/cvmfs/unpacked.cern.ch` (Combine v11.1.0 when these numbers were made).
3. For each benchmark and channel, merges the Run 2 and Run 3 cards with `combineCards.py R2=... R3=...` into `R23_<channel>_2L.txt`.
4. Runs `combine -M AsymptoticLimits <card> --run blind` on the Run 2, Run 3 and merged cards. `--run blind` builds the expected limits from the a-priori background and never reads the data, so the placeholder observation of 0 in region A plays no role.
5. Writes the 2.5%, 16%, 50%, 84% and 97.5% expected limits to `<workdir>/limits.txt`.

It runs in a few minutes. The `latest` tag of the container can change, so record the Combine version printed in the log if the numbers ever differ.

## Reference numbers (2026-10-08)

Median expected limit on `r`, with the 68% band in brackets.

| Channel | Benchmark | Run 2 | Run 3 | Run 2 + Run 3 |
|---|---|---|---|---|
| OS OF 1FJ (WWH) | SM | 538 [359, 827] | 450 [304, 685] | 321 [220, 478] |
| | c2V = 1.5 | 51.8 [33.6, 81.7] | 37.0 [24.0, 58.7] | 27.1 [18.1, 41.9] |
| | κλ = 10 | 240 [161, 367] | 117 [79.3, 178] | 99.8 [68.3, 149] |
| OS SF 2FJ (WZH/ZZH) | SM | 2131 [1532, 2989] | 900 [549, 1481] | 799 [501, 1267] |
| | c2V = 1.5 | 121 [85.8, 175] | 46.3 [26.3, 82.8] | 41.8 [24.8, 69.4] |
| | κλ = 10 | 536 [391, 750] | 277 [166, 455] | 237 [149, 371] |

## Known problems with these cards

These numbers come from the cards exactly as written. They are not yet suitable for the AN, for the following reasons.

1. **No ABCD background model.**
   The cards contain no `rateParam` linking region A to B·C/D. The background in every region is the simulated rate, constrained only by an MC-statistics `lnN`, so these are counting limits on simulation and not the data-driven ABCD estimate described in the AN.
   Regions B, C and D also enter as independent counting channels that carry signal, so their signal contamination is not handled as it would be in an ABCD fit.

2. **The regions do not match the AN signal-region tables.**
   For example, region A of the Run 3 1FJ card at c2V = 1.5 has 0.109 signal and 1.854 background events. The signal region in AN Table `2losof1fj:tab:sr_run3` (ABCDNet ≥ 0.87, VBS ≥ 0.71) has 1.689 and 72.27.
   The cards also reuse a single cut set: the 1FJ SM and κλ = 10 cards have identical backgrounds and observations, and the 2FJ backgrounds are identical for all three benchmarks. The AN instead optimizes the cuts per benchmark.
   Ask the collaborator which selection and cuts were used.

3. **The MC-statistics uncertainties look wrong.**
   Region-A background rate and its `mcstat` lnN, as read from the cards:

   | Card | Background in A | lnN |
   |---|---|---|
   | R2 WWH (SM, κλ = 10) | 3.341 | 1.0836 |
   | R2 WWH (c2V = 1.5) | 1.644 | 1.0843 |
   | R3 WWH (SM, κλ = 10) | 3.728 | 1.3344 |
   | R3 WWH (c2V = 1.5) | 1.854 | 1.0293 |
   | R2 WZH_ZZH (all) | 2.714 | 3.2809 |
   | R3 WZH_ZZH (all) | 0.101 | 1.0079 |

   A yield of 0.101 events cannot have a 0.8% statistical uncertainty, and 2.714 events with 228% is extreme.
   This pattern fits an uncertainty computed as variance/yield instead of √variance/yield. That happens when the squared variance is passed as `sumw2`, which is exactly what `VBSVVH/python/make_abcd.py` line 166 does (`h_proc.variances() ** 2`).
   Check how these cards were produced, and fix `make_abcd.py` before using it for datacards.

4. **Region A data is a placeholder.**
   The observation in A is 0 in every card. This is fine for blind expected limits, but the real count is needed for any unblinded or post-fit result.

## Next steps

- Build cards from the AN signal-region definitions (`abcd_scan.py`), one per benchmark and era.
- Model the background with the ABCD relation, e.g. free `rateParam`s for B, C and D and A expressed as B·C/D.
- Add the closure (non-closure) uncertainty and the systematic uncertainties.
