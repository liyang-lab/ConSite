from __future__ import annotations
import hashlib
from pathlib import Path
from dataclasses import dataclass


def md5_text(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


@dataclass
class Hit:
    family: str  # Pfam accession, e.g., PF00018
    name: str    # Pfam name/description if available
    ali_start: int
    ali_end: int
    evalue: float
    score: float