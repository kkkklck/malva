import pytest
import yaml

from malva.index import _parse_bam_tag_name, load_flavor


class TestParseBamTagName:
    def test_extracts_the_tag(self):
        assert _parse_bam_tag_name("CB:{cell}") == "CB"
        assert _parse_bam_tag_name("CR:{cell}") == "CR"

    def test_defaults_to_cb(self):
        assert _parse_bam_tag_name("not a tag") == "CB"


class TestLoadFlavor:
    def test_from_flavor_config(self, tmp_path):
        config = tmp_path / "flavors.yaml"
        config.write_text(
            yaml.safe_dump({"barcode_flavors": {"myflavor": {"barcode_position": [0, 16]}}})
        )
        assert load_flavor("myflavor", str(config)) == {"barcode_position": [0, 16]}

    def test_unknown_flavor_raises(self, tmp_path):
        config = tmp_path / "flavors.yaml"
        config.write_text(yaml.safe_dump({"barcode_flavors": {}}))
        with pytest.raises(ValueError):
            load_flavor("missing", str(config))

    def test_from_custom_yaml(self, tmp_path):
        custom = tmp_path / "custom.yaml"
        custom.write_text(yaml.safe_dump({"barcode_position": [2, 27]}))
        assert load_flavor(str(custom), "unused") == {"barcode_position": [2, 27]}

    def test_missing_custom_yaml_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_flavor(str(tmp_path / "missing.yaml"), "unused")
