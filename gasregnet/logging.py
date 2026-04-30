"""Structured logging helpers."""

from __future__ import annotations

import logging

import structlog


def configure_logging(*, verbose: bool = False) -> None:
    """Configure structlog for command-line runs."""

    logging.basicConfig(
        force=True,
        format="%(message)s",
        level=logging.DEBUG if verbose else logging.INFO,
    )
    structlog.configure(
        processors=[
            structlog.processors.add_log_level,
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.JSONRenderer(),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(
            logging.DEBUG if verbose else logging.INFO,
        ),
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=True,
    )
