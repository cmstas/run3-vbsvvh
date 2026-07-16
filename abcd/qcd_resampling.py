import zlib

import numpy as np
import uproot
import awkward as ak


# Per-jet branches needed to rebuild the candidates (must exist in the ntuple).
PERJET_BRANCHES = [
    "fatjet_pt", "fatjet_eta", "fatjet_phi", "fatjet_mass", "fatjet_msoftdrop",
    "fatjet_globalParT3_mass", "fatjet_tau1", "fatjet_tau2",
    "fatjet_HvsQCD", "fatjet_VvsQCD",
]

MISSING = -999.0
VCAND_DR = 0.8  # dR of a V candidate jet from the H candidate (selections.cpp)


def _pt_label(lo, hi):
    hi_s = "Inf" if hi >= 1.0e5 else str(int(round(hi)))
    return f"{int(round(lo))}to{hi_s}"


def _cdf_from_counts(counts):
    """Normalised cumulative (length nbins, last entry == 1) or None if empty."""
    total = counts.sum()
    if total <= 0:
        return None
    return np.cumsum(counts) / total


def _invert_cdf(cdf, edges, u):
    """Vectorised inverse-CDF sampling with linear interpolation inside the bin.

    cdf   : (nbins,) monotonic, cdf[-1] == 1
    edges : (nbins+1,)
    u     : (N,) uniform draws in [0, 1)
    """
    k = np.searchsorted(cdf, u, side="right")
    k = np.clip(k, 0, len(cdf) - 1)
    cdf_lo = np.where(k > 0, cdf[np.clip(k - 1, 0, None)], 0.0)
    cdf_hi = cdf[k]
    denom = np.maximum(cdf_hi - cdf_lo, 1e-12)
    frac = np.clip((u - cdf_lo) / denom, 0.0, 1.0)
    return edges[k] + frac * (edges[k + 1] - edges[k])


class ScoreResampler:
    """Loads the data templates and resamples QCD fat-jet scores + candidates."""

    def __init__(self, template_path, qcd_shortname_prefix="QCD", seed=42):
        self.qcd_prefix = qcd_shortname_prefix
        self.base_seed = int(seed)

        with uproot.open(template_path) as f:
            h3 = f["h3_HvsQCD"]
            self.score_edges = np.asarray(h3.axis(0).edges(), dtype=np.float64)
            self.pt_edges = np.asarray(h3.axis(1).edges(), dtype=np.float64)
            self.eta_edges = np.asarray(h3.axis(2).edges(), dtype=np.float64)
            h_vals = h3.values()  # (nscore, npt, neta)

            self.n_pt = len(self.pt_edges) - 1
            self.n_eta = len(self.eta_edges) - 1

            # P(H | pT, |eta|) cumulative per (ipt, ieta)
            self.h_cdf = {}
            for ipt in range(self.n_pt):
                for ieta in range(self.n_eta):
                    self.h_cdf[(ipt, ieta)] = _cdf_from_counts(h_vals[:, ipt, ieta])
            # pT-integrated H marginal as a fallback for empty (pT,|eta|) cells
            self.h_cdf_pt = {
                ipt: _cdf_from_counts(h_vals[:, ipt, :].sum(axis=1))
                for ipt in range(self.n_pt)
            }

            # P(V | H, pT) from the eta-integrated joint TH2 per pT bin
            self.joint_hedges = None
            self.joint_vedges = None
            self.v_cdf = {}          # ipt -> (nHbins, nVbins) row-wise cumulative
            for ipt in range(self.n_pt):
                lab = _pt_label(self.pt_edges[ipt], self.pt_edges[ipt + 1])
                hj = f[f"joint_HV_pt{lab}"]
                if self.joint_hedges is None:
                    self.joint_hedges = np.asarray(hj.axis(0).edges(), dtype=np.float64)
                    self.joint_vedges = np.asarray(hj.axis(1).edges(), dtype=np.float64)
                jv = hj.values()  # (nHbins, nVbins)
                n_h = jv.shape[0]
                # fallback = V marginal of this pT bin (integrated over H)
                v_marg = _cdf_from_counts(jv.sum(axis=0))
                rows = np.zeros_like(jv, dtype=np.float64)
                for h in range(n_h):
                    row = _cdf_from_counts(jv[h])
                    rows[h] = row if row is not None else (
                        v_marg if v_marg is not None else np.linspace(0, 1, jv.shape[1]))
                self.v_cdf[ipt] = rows

    # ------------------------------------------------------------------ helpers
    def is_qcd(self, shortname):
        if shortname is None:
            return False
        if isinstance(shortname, bytes):
            shortname = shortname.decode("utf-8", "ignore")
        return str(shortname).startswith(self.qcd_prefix)

    def _bin_index(self, values, edges, n):
        # 0-based bin, clamped into [0, n-1] (overflow/underflow -> edge bins)
        idx = np.searchsorted(edges, values, side="right") - 1
        return np.clip(idx, 0, n - 1)

    def _sample_flat(self, flat_pt, flat_abseta, rng):
        """Draw (H*, V*) for a flat array of jets."""
        n = len(flat_pt)
        H = np.empty(n, dtype=np.float64)
        V = np.empty(n, dtype=np.float64)
        if n == 0:
            return H, V

        ipt = self._bin_index(flat_pt, self.pt_edges, self.n_pt)
        ieta = self._bin_index(flat_abseta, self.eta_edges, self.n_eta)
        u_h = rng.random(n)
        u_v = rng.random(n)

        # --- H* from P(H | pT, |eta|) (loop over the few (pT,|eta|) cells) ---
        for a in range(self.n_pt):
            for b in range(self.n_eta):
                sel = (ipt == a) & (ieta == b)
                if not sel.any():
                    continue
                cdf = self.h_cdf[(a, b)]
                if cdf is None:
                    cdf = self.h_cdf_pt[a]
                if cdf is None:  # utterly empty template bin -> leave scores at 0
                    H[sel] = 0.0
                    continue
                H[sel] = _invert_cdf(cdf, self.score_edges, u_h[sel])

        # --- V* from P(V | H*, pT) using the eta-integrated joint (loop over pT) ---
        jhbin = self._bin_index(H, self.joint_hedges, len(self.joint_hedges) - 1)
        for a in range(self.n_pt):
            sel = ipt == a
            if not sel.any():
                continue
            rows = self.v_cdf[a][jhbin[sel]]              # (Nsel, nVbins)
            u = u_v[sel]
            k = (rows < u[:, None]).sum(axis=1)
            k = np.clip(k, 0, rows.shape[1] - 1)
            cdf_lo = np.where(k > 0, rows[np.arange(len(k)), np.clip(k - 1, 0, None)], 0.0)
            cdf_hi = rows[np.arange(len(k)), k]
            denom = np.maximum(cdf_hi - cdf_lo, 1e-12)
            frac = np.clip((u - cdf_lo) / denom, 0.0, 1.0)
            V[sel] = self.joint_vedges[k] + frac * (self.joint_vedges[k + 1] - self.joint_vedges[k])
        return H, V

    # --------------------------------------------------------------- public API
    def resample_candidates(self, perjet, rng):
        """Given a dict of jagged per-jet arrays, return corrected candidate scalars.

        `perjet` maps each name in PERJET_BRANCHES to an awkward jagged array
        (one entry per event). Returns a dict of flat numpy arrays keyed by the
        boosted_{h,v1,v2}_candidate_* names to overwrite.
        """
        pt = ak.values_astype(perjet["fatjet_pt"], np.float64)
        eta = ak.values_astype(perjet["fatjet_eta"], np.float64)
        phi = ak.values_astype(perjet["fatjet_phi"], np.float64)
        mass = ak.values_astype(perjet["fatjet_mass"], np.float64)
        msd = ak.values_astype(perjet["fatjet_msoftdrop"], np.float64)
        pmass = ak.values_astype(perjet["fatjet_globalParT3_mass"], np.float64)
        tau1 = ak.values_astype(perjet["fatjet_tau1"], np.float64)
        tau2 = ak.values_astype(perjet["fatjet_tau2"], np.float64)
        tau21 = ak.where(tau1 > 0, tau2 / ak.where(tau1 > 0, tau1, 1.0), MISSING)

        counts = ak.num(pt, axis=1)
        flat_pt = ak.to_numpy(ak.flatten(pt))
        flat_abseta = np.abs(ak.to_numpy(ak.flatten(eta)))

        flat_H, flat_V = self._sample_flat(flat_pt, flat_abseta, rng)
        Hs = ak.unflatten(flat_H, counts)
        Vs = ak.unflatten(flat_V, counts)

        # --- H candidate = argmax H* ---
        h_loc = ak.argmax(Hs, axis=1, keepdims=True)
        h_eta = ak.firsts(eta[h_loc])
        h_phi = ak.firsts(phi[h_loc])

        # --- V candidates: dR >= 0.8 from H candidate, top-2 by V* ---
        deta = eta - h_eta
        dphi = (phi - h_phi + np.pi) % (2.0 * np.pi) - np.pi
        dR = np.sqrt(deta * deta + dphi * dphi)
        Vmask = ak.where(dR >= VCAND_DR, Vs, MISSING)
        order = ak.argsort(Vmask, axis=1, ascending=False)
        # pad so [:,0] and [:,1] always exist (events have >=3 jets, but be safe)
        order = ak.pad_none(order, 2, axis=1)
        v1_loc = order[:, 0:1]
        v2_loc = order[:, 1:2]

        def gather(arr, loc):
            return ak.to_numpy(ak.fill_none(ak.firsts(arr[loc]), MISSING))

        out = {}

        def fill(prefix, loc, score_src):
            score = ak.to_numpy(ak.fill_none(ak.firsts(score_src[loc]), MISSING))
            found = score > 0
            def kine(arr):
                return np.where(found, gather(arr, loc), MISSING)
            out[f"{prefix}_score"] = score.astype(np.float32)
            out[f"{prefix}_pt"] = kine(pt).astype(np.float32)
            out[f"{prefix}_eta"] = kine(eta).astype(np.float32)
            out[f"{prefix}_phi"] = kine(phi).astype(np.float32)
            out[f"{prefix}_mass"] = kine(mass).astype(np.float32)
            out[f"{prefix}_msoftdrop"] = kine(msd).astype(np.float32)
            out[f"{prefix}_part_mass"] = kine(pmass).astype(np.float32)
            out[f"{prefix}_tau21"] = kine(tau21).astype(np.float32)

        fill("boosted_h_candidate", h_loc, Hs)
        fill("boosted_v1_candidate", v1_loc, Vmask)
        fill("boosted_v2_candidate", v2_loc, Vmask)
        return out

    def rng_for(self, path):
        """A reproducible, per-file Generator (safe for threaded loading).

        Uses a stable hash (crc32) of the path so the same file resamples
        identically across processes -- e.g. training and a later --infer run.
        """
        h = (zlib.crc32(str(path).encode("utf-8")) ^ self.base_seed) & 0xFFFFFFFF
        return np.random.default_rng(h)
