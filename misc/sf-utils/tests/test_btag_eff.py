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
families = load("btag_families_test", "btag_eff_families.py")


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
    def test_yaml_controls_strict_final_groups(self):
        sample = "VBSWZH_c2v1p0"
        self.assertEqual(families.sample_family(sample), "VBS_VVH")
        self.assertEqual(families.final_sample_family(sample), "VBS_VVH")
        self.assertEqual(families.final_channel("1lep_1FJ"), "leptonic")
        self.assertNotIn("0lep_1FJ_met", families.retained_source_channels())
        self.assertIn("0lep_1FJ_met", families.excluded_source_channels())
        with self.assertRaises(ValueError):
            families.final_channel("all_events")

    def test_yaml_rejects_incomplete_or_duplicate_final_groups(self):
        config = families.load_config()
        bad = json.loads(json.dumps(config))
        bad["final_merges"]["samples"]["VBS_VVH"] = []
        with self.assertRaises(ValueError):
            families.validate_config(bad)
        bad = json.loads(json.dumps(config))
        bad["final_merges"]["channels"]["0lep_0FJ"].append("1lep_1FJ")
        with self.assertRaises(ValueError):
            families.validate_config(bad)

    def test_runtime_lookup_uses_final_keys_without_exact_fallback(self):
        sample = "WplusH-Wto2Q-Hto2B_Par-M-125"
        self.assertEqual(families.final_runtime_keys("2024Prompt", "1lep_1FJ", sample),
                         ("btag_2024Prompt_leptonic", "VH"))
        with self.assertRaises(ValueError):
            families.final_runtime_keys("2024Prompt", "1lep_1FJ", "not_a_configured_sample")

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

    def test_final_input_discovery_requires_all_channels_and_manifests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError):
                conv.discover_final_inputs(root)
            for channel in families.retained_source_channels():
                directory = root / channel
                directory.mkdir()
                (directory / "manifest.json").write_text(json.dumps({"jobs": {}}))
            (root / "0lep_1FJ_met").mkdir()
            inputs, manifests, ignored = conv.discover_final_inputs(root)
            self.assertEqual(set(inputs), set(families.retained_source_channels()))
            self.assertEqual(set(manifests), set(inputs))
            self.assertEqual(ignored, ["0lep_1FJ_met"])
            (root / "1lep_1FJ" / "manifest.json").unlink()
            with self.assertRaises(ValueError):
                conv.discover_final_inputs(root)

    def test_failed_final_build_does_not_touch_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "btag_eff.json"
            output.write_text("original payload\n")
            previous_argv = sys.argv
            try:
                sys.argv = ["bEff-convert-to-correction.py", "--final", "--input-root", str(root),
                            "--year", "2024Prompt", "--output", str(output)]
                with self.assertRaises(ValueError):
                    conv.main()
            finally:
                sys.argv = previous_argv
            self.assertEqual(output.read_text(), "original payload\n")

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

    def test_global_diagnostics_accept_only_final_input_root_cli(self):
        parser = global_plots.parse_args
        previous_argv = sys.argv
        try:
            sys.argv = ["plot-btag-eff-global.py", "--input-root", "/tmp/inputs", "--year", "2024Prompt",
                        "--mode", "families", "--final"]
            args = parser()
        finally:
            sys.argv = previous_argv
        self.assertEqual(args.input_root, Path("/tmp/inputs"))

    def test_matrix_plots_use_module_categories(self):
        raw = counts()
        result = plots.efficiencies(conv, raw, {key: value.copy() for key, value in raw.items()})
        with tempfile.TemporaryDirectory() as tmp:
            args = type("Args", (), {"plot_dir": Path(tmp)})()
            plots.matrix_plots({"a": None}, {"a": result}, args, 1., 13.)
            self.assertTrue((Path(tmp) / "compatibility_b_T.pdf").exists())


if __name__ == "__main__":
    unittest.main()
