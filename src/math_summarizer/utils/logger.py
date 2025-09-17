"""Logging utilities for the Math Summarizer."""

import logging
import sys
from pathlib import Path
from typing import Optional

from ..config import config


def setup_logger(name: str, log_file: Optional[str] = None) -> logging.Logger:
    """Set up a logger with console and file handlers."""
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, config.log_level.upper()))
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler
    if log_file or config.log_file:
        file_path = Path(log_file or config.log_file)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger
