import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


conv = load("btag_converter_test", "bEff-convert-to-correction.py")
plots = load("btag_plots_test", "plot-btag-eff-families.py")
global_plots = load("btag_global_test", "plot-btag-eff-global.py")


def counts(den=10., tight=2., loose=5., lt=3., untagged=5.):
    values = {}
    for flavor in conv.FLAVORS:
        values[flavor, "den"] = np.array([[den]])
        values[flavor, "T"] = np.array([[tight]])
        values[flavor, "L"] = np.array([[loose]])
        values[flavor, "LT"] = np.array([[lt]])
        values[flavor, "N"] = np.array([[untagged]])
    return values


class BTagEfficiencyTests(unittest.TestCase):
    def test_sample_lookup_prefers_family_then_exact(self):
        sample = "VBSWZH_c2v1p0"
        self.assertEqual(conv.select_efficiency_sample_key(sample, {"VBS_VVH", sample}), "VBS_VVH")
        self.assertEqual(conv.select_efficiency_sample_key(sample, {sample}), sample)
        with self.assertRaises(KeyError):
            conv.select_efficiency_sample_key(sample, set())

    def test_unweighted_category_uncertainty(self):
        uncertainty = conv.mcstat_efficiency_uncertainty(
            np.array([[2.]]), np.array([[10.]]), np.array([[2.]]), np.array([[10.]]))
        self.assertAlmostEqual(uncertainty[0, 0] ** 2, .2 * .8 / 10.)

    def test_weighted_category_uncertainty(self):
        uncertainty = conv.mcstat_efficiency_uncertainty(
            np.array([[3.]]), np.array([[5.]]), np.array([[5.]]), np.array([[13.]]))
        self.assertAlmostEqual(uncertainty[0, 0] ** 2, ((1 - 2 * .6) * 5 + .6 ** 2 * 13) / 25.)
        invalid = conv.mcstat_efficiency_uncertainty(
            np.array([[3.]]), np.array([[5.]]), np.array([[100.]]), np.array([[1.]]))
        self.assertTrue(np.isnan(invalid[0, 0]))

    def test_exclusive_category_identity_and_mask(self):
        valid = counts()
        for flavor in conv.FLAVORS:
            self.assertTrue(np.allclose(valid[flavor, "T"] + valid[flavor, "LT"] + valid[flavor, "N"],
                                        valid[flavor, "den"]))
            self.assertTrue(np.allclose(valid[flavor, "T"] + valid[flavor, "LT"], valid[flavor, "L"]))
        bad = counts(loose=4.)
        self.assertTrue(conv.invalid_count_bins(bad)["b"][0, 0])
        result = plots.efficiencies(conv, bad, {key: value.copy() for key, value in bad.items()})
        self.assertFalse(result[2]["b", "T"][0, 0])
        self.assertNotIn("L", conv.EXCLUSIVE_CATEGORIES)

    def test_manifest_completeness(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({"jobs": {
                "sample_0": {"sample": "sample", "job_idx": 0},
                "sample_1": {"sample": "sample", "job_idx": 1},
            }}))
            sample_dir = root / "outputs" / "sample"
            sample_dir.mkdir(parents=True)
            (sample_dir / "output_0.root").touch()
            with self.assertRaises(ValueError):
                conv.validate_job_manifest(root / "outputs", manifest)
            (sample_dir / "output_1.root").touch()
            self.assertEqual(conv.validate_job_manifest(root / "outputs", manifest)["sample"]["discovered_jobs"], 2)
            (sample_dir / "output_01.root").touch()
            with self.assertRaises(ValueError):
                conv.validate_job_manifest(root / "outputs", manifest)

    def test_sparse_summary_is_unavailable(self):
        values = {(flavor, wp): np.array([[.2]]) for flavor in conv.FLAVORS for wp in conv.WORKING_POINTS}
        errors = {(flavor, wp): np.array([[.1]]) for flavor in conv.FLAVORS for wp in conv.WORKING_POINTS}
        valid = {(flavor, wp): np.array([[flavor == "b" and wp == "T"]])
                 for flavor in conv.FLAVORS for wp in conv.WORKING_POINTS}
        pathological = {flavor: np.array([[False]]) for flavor in conv.FLAVORS}
        result = values, errors, valid, pathological
        with tempfile.TemporaryDirectory() as tmp:
            args = type("Args", (), {"plot_dir": Path(tmp)})()
            plots.summary_json({"a": None, "b": None}, {"a": result, "b": result}, args)
            summary = json.loads((Path(tmp) / "compatibility_summary.json").read_text())
            entry = summary["a__vs__b/all_flavours_all_categories"]
            self.assertFalse(entry["available"])
            self.assertIsNone(entry["fraction_gt2"])

    def test_overlapping_channels_rejected(self):
        with self.assertRaises(ValueError):
            global_plots.reject_overlapping_channels(["0lep_1FJ", "0lep_1FJ_met"])


if __name__ == "__main__":
    unittest.main()
