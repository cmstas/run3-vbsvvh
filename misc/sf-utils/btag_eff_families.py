"""Shared b-tag efficiency grouping configuration."""

from pathlib import Path

import yaml


DEFAULT_CONFIG = (Path(__file__).resolve().parents[2] / "preselection" / "corrections" /
                  "scalefactors" / "btagging" / "btag_eff_families.yaml")


def load_config(path=DEFAULT_CONFIG):
    config = yaml.safe_load(Path(path).read_text())
    if not isinstance(config, dict) or "preliminary_families" not in config or "final_merges" not in config:
        raise ValueError(f"Invalid b-tag family configuration: {path}")
    return config


def sample_family(sample, config_path=DEFAULT_CONFIG):
    """Return the first preliminary family match; YAML order resolves overlaps."""
    for family, needles in load_config(config_path)["preliminary_families"].items():
        if any(needle in sample for needle in needles):
            return family
    raise ValueError(f"No preliminary b-tag efficiency family is configured for sample {sample!r}")


def final_group(kind, preliminary_name, config_path=DEFAULT_CONFIG):
    groups = load_config(config_path)["final_merges"][kind]
    matches = [name for name, members in groups.items() if preliminary_name in members]
    if len(matches) != 1:
        raise ValueError(f"{preliminary_name!r} must occur in exactly one final {kind} group; found {matches}")
    return matches[0]


def final_sample_family(sample, config_path=DEFAULT_CONFIG):
    return final_group("samples", sample_family(sample, config_path), config_path)


def final_channel(channel, config_path=DEFAULT_CONFIG):
    return final_group("channels", channel, config_path)
