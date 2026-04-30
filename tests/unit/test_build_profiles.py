from __future__ import annotations

from pathlib import Path

from pyhmmer import plan7

from scripts.build_profiles import build_profiles


def test_build_profiles_writes_loadable_hmms(tmp_path: Path) -> None:
    manifest = build_profiles(
        config=Path("configs"),
        out_dir=tmp_path / "profiles",
        manifest_out=tmp_path / "profiles.yaml",
        source="seed",
    )

    assert manifest.height >= 1
    for profile in manifest["profile"].to_list():
        with plan7.HMMFile(profile) as hmm_file:
            hmm = hmm_file.read()
        assert hmm is not None
    assert (tmp_path / "profiles.yaml").exists()
    assert (tmp_path / "profiles.csv").exists()
    assert manifest["sha256"].str.len_chars().min() == 64
    assert manifest["nseq"].min() == 1
