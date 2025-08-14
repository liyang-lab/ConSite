from __future__ import annotations
from typing import Tuple
import numpy as np
from Bio import AlignIO


def read_stockholm(path) -> Tuple[np.ndarray, list]:
    aln = AlignIO.read(str(path), "stockholm")
    seq_ids = [rec.id for rec in aln]
    arr = np.array([list(str(rec.seq).upper()) for rec in aln], dtype='<U1')
    return arr, seq_ids