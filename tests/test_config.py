import pytest

from xenium_showcase.config import validate_config


def valid_config():
    return {"project_name": "test", "dataset_kind": "test", "input_zip": "input.zip",
            "output_dir": "results/test", "random_seed": 17, "markers": {"A": ["G1"]}}


def test_config_validation_accepts_minimum_contract():
    assert validate_config(valid_config())["random_seed"] == 17


def test_config_validation_rejects_invalid_contamination():
    cfg = valid_config()
    cfg["qc"] = {"anomaly": {"contamination": .8}}
    with pytest.raises(ValueError, match="contamination"):
        validate_config(cfg)

