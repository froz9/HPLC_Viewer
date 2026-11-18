"""Logging configuration for HPLC Viewer."""

import logging
from pathlib import Path
from typing import Optional


def get_logger(name: str = "hplc_viewer", log_dir: Optional[Path] = None) -> logging.Logger:
    """Configure and return an application logger without duplicate handlers."""
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    try:
        if log_dir is None:
            log_dir = Path(__file__).resolve().parents[2] / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_dir / "hplc_viewer.log", encoding="utf-8")
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except OSError:
        # In read-only environments, fall back to console logging only.
        pass

    logger.propagate = False
    return logger
