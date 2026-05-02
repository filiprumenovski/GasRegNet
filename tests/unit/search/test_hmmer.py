from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from gasregnet.errors import MissingInputError
from gasregnet.search.hmmer import hmmsearch


class _FakeContext:
    def __init__(self, values: list[object]) -> None:
        self._values = values

    def __enter__(self) -> list[object]:
        return self._values

    def __exit__(self, *_: object) -> None:
        return None


class _FakeHit:
    name = b"target1"
    evalue = 1e-20
    score = 55.0
    bias = 0.1
    included = True


class _FakeTopHits:
    query_name = b"profile1"

    def __iter__(self) -> Any:
        return iter([_FakeHit()])


def test_hmmsearch_parses_pyhmmer_hits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile = tmp_path / "profile.hmm"
    sequences = tmp_path / "seqs.faa"
    profile.write_text("HMMER3/f\n", encoding="utf-8")
    sequences.write_text(">target1\nM\n", encoding="utf-8")

    monkeypatch.setattr(
        "gasregnet.search.hmmer.plan7.HMMFile",
        lambda _: _FakeContext([object()]),
    )
    monkeypatch.setattr(
        "gasregnet.search.hmmer.easel.SequenceFile",
        lambda *_, **__: _FakeContext([object()]),
    )
    monkeypatch.setattr(
        "gasregnet.search.hmmer.hmmer.hmmsearch",
        lambda *_args, **_kwargs: [_FakeTopHits()],
    )

    hits = hmmsearch(profile, sequences)

    assert hits.height == 1
    assert hits["query_id"].item() == "profile1"
    assert hits["target_id"].item() == "target1"
    assert hits["evalue"].item() == pytest.approx(1e-20)
    assert hits["included"].item() is True


def test_hmmsearch_streams_sequence_file_to_pyhmmer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile = tmp_path / "profile.hmm"
    sequences = tmp_path / "seqs.faa"
    profile.write_text("HMMER3/f\n", encoding="utf-8")
    sequences.write_text(">target1\nM\n", encoding="utf-8")
    streamed_sequences = object()
    seen_targets: list[object] = []

    monkeypatch.setattr(
        "gasregnet.search.hmmer.plan7.HMMFile",
        lambda _: _FakeContext([object()]),
    )
    monkeypatch.setattr(
        "gasregnet.search.hmmer.easel.SequenceFile",
        lambda *_, **__: _FakeContext(streamed_sequences),
    )

    def fake_hmmsearch(
        _queries: object,
        targets: object,
        **_kwargs: object,
    ) -> list[_FakeTopHits]:
        seen_targets.append(targets)
        return [_FakeTopHits()]

    monkeypatch.setattr("gasregnet.search.hmmer.hmmer.hmmsearch", fake_hmmsearch)

    hits = hmmsearch(profile, sequences)

    assert hits.height == 1
    assert seen_targets == [streamed_sequences]


def test_hmmsearch_batches_profiles_without_materializing_all(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile = tmp_path / "profile.hmm"
    sequences = tmp_path / "seqs.faa"
    profile.write_text("HMMER3/f\n", encoding="utf-8")
    sequences.write_text(">target1\nM\n", encoding="utf-8")
    seen_query_batch_sizes: list[int] = []
    sequence_file_opens = 0

    monkeypatch.setattr(
        "gasregnet.search.hmmer.plan7.HMMFile",
        lambda _: _FakeContext([object(), object(), object(), object(), object()]),
    )

    def fake_sequence_file(*_args: object, **_kwargs: object) -> _FakeContext:
        nonlocal sequence_file_opens
        sequence_file_opens += 1
        return _FakeContext(object())

    monkeypatch.setattr("gasregnet.search.hmmer.easel.SequenceFile", fake_sequence_file)

    def fake_hmmsearch(
        queries: object,
        _targets: object,
        **_kwargs: object,
    ) -> list[_FakeTopHits]:
        seen_query_batch_sizes.append(len(queries))  # type: ignore[arg-type]
        return [_FakeTopHits()]

    monkeypatch.setattr("gasregnet.search.hmmer.hmmer.hmmsearch", fake_hmmsearch)

    hits = hmmsearch(profile, sequences, query_batch_size=2)

    assert hits.height == 3
    assert seen_query_batch_sizes == [2, 2, 1]
    assert sequence_file_opens == 3


def test_hmmsearch_rejects_invalid_query_batch_size(tmp_path: Path) -> None:
    profile = tmp_path / "profile.hmm"
    sequences = tmp_path / "seqs.faa"
    profile.write_text("HMMER3/f\n", encoding="utf-8")
    sequences.write_text(">target1\nM\n", encoding="utf-8")

    with pytest.raises(ValueError, match="query_batch_size"):
        hmmsearch(profile, sequences, query_batch_size=0)


def test_hmmsearch_rejects_missing_profile(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError, match="HMM profile"):
        hmmsearch(tmp_path / "missing.hmm", tmp_path / "seqs.faa")


def test_hmmsearch_rejects_missing_sequences(tmp_path: Path) -> None:
    profile = tmp_path / "profile.hmm"
    profile.write_text("HMMER3/f\n", encoding="utf-8")

    with pytest.raises(MissingInputError, match="sequence FASTA"):
        hmmsearch(profile, tmp_path / "missing.faa")
