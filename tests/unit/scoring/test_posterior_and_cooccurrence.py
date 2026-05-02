from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.config import load_config
from gasregnet.scoring.candidates import score_candidates
from gasregnet.scoring.cooccurrence import assign_phylogenetic_profile_scores
from gasregnet.scoring.loci import score_loci
from gasregnet.scoring.posterior import assign_operon_regulation_score_bands
from gasregnet.simulation.synthetic_truth import simulate_synthetic_truth_corpus
from tests.unit.test_schemas import candidates_frame, loci_frame


def test_phylogenetic_profile_cooccurrence_scores_mutual_information() -> None:
    config = load_config(Path("configs"))
    loci = pl.concat(
        [
            loci_frame().with_columns(
                pl.lit("CO_1").alias("locus_id"),
                pl.lit("Genome A").alias("organism"),
            ),
            loci_frame().with_columns(
                pl.lit("CO_2").alias("locus_id"),
                pl.lit("Genome B").alias("organism"),
            ),
            loci_frame().with_columns(
                pl.lit("CN_1").alias("locus_id"),
                pl.lit("CN").alias("analyte"),
                pl.lit("Genome C").alias("organism"),
            ),
            loci_frame().with_columns(
                pl.lit("CN_2").alias("locus_id"),
                pl.lit("CN").alias("analyte"),
                pl.lit("Genome D").alias("organism"),
            ),
        ],
    )
    candidates = pl.concat(
        [
            candidates_frame().with_columns(
                pl.lit("CO_1").alias("locus_id"),
                pl.lit("Genome A").alias("organism"),
            ),
            candidates_frame().with_columns(
                pl.lit("CO_2").alias("locus_id"),
                pl.lit("Genome B").alias("organism"),
            ),
        ],
    )

    scored = assign_phylogenetic_profile_scores(
        candidates,
        loci,
        scoring=config.scoring,
    )

    assert scored["phylogenetic_profile_score"].min() == 1.0
    assert scored["candidate_score"].min() > candidates["candidate_score"].min()


def test_operon_score_band_adds_ordered_94_percent_bands() -> None:
    candidates = pl.concat(
        [
            candidates_frame().with_columns(pl.lit(20.0).alias("candidate_score")),
            candidates_frame().with_columns(
                pl.lit("cand2").alias("candidate_id"),
                pl.lit(5.0).alias("candidate_score"),
            ),
        ],
    )

    scored = assign_operon_regulation_score_bands(candidates)
    top = scored.sort("regulation_logit_score", descending=True).row(0, named=True)
    bottom = scored.sort("regulation_logit_score").row(0, named=True)

    assert top["candidate_id"] == "cand1"
    assert top["regulation_logit_score"] > bottom["regulation_logit_score"]
    assert (
        top["score_band_low"]
        < top["regulation_logit_score"]
        < top["score_band_high"]
    )
    assert top["score_band_model"] == "uncalibrated_sigmoid_beta_score_band_94"


def test_synthetic_truth_corpus_recovers_planted_regulators_above_decoys() -> None:
    config = load_config(Path("configs"))
    synthetic = simulate_synthetic_truth_corpus(
        n_genomes=12,
        positive_fraction=0.5,
        annotation_noise=0.0,
        seed=7,
    )

    loci = score_loci(synthetic.loci, config.scoring)
    candidates = score_candidates(loci, synthetic.genes, config.scoring)
    candidates = assign_phylogenetic_profile_scores(
        candidates,
        loci,
        scoring=config.scoring,
    )
    candidates = assign_operon_regulation_score_bands(candidates)
    truth = synthetic.ground_truth.filter(pl.col("is_true_regulator"))
    recovered = candidates.join(
        truth.select(["gene_accession", "is_true_regulator"]),
        on="gene_accession",
        how="left",
    )

    assert recovered.height == truth.height
    assert recovered["is_true_regulator"].all()
    assert recovered["regulation_logit_score"].min() > 0.4
