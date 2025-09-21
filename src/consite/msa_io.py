from __future__ import annotations
from typing import Tuple, Dict
import numpy as np
from Bio import AlignIO


def read_stockholm(path) -> Tuple[np.ndarray, list]:
    aln = AlignIO.read(str(path), "stockholm")
    seq_ids = [rec.id for rec in aln]
    arr = np.array([list(str(rec.seq).upper()) for rec in aln], dtype='<U1')
    return arr, seq_ids

def read_stockholm_with_meta(path) -> Tuple[np.ndarray, list, Dict[str, dict]]:
    aln = AlignIO.read(str(path), "stockholm")
    seq_ids = [rec.id for rec in aln]
    arr = np.array([list(str(rec.seq).upper()) for rec in aln], dtype='<U1')

    meta: Dict[str, dict] = {sid: {} for sid in seq_ids}
    # pass 2: parse GS lines for OS (organism), AC, etc.
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#=GS "):
                _, sid, key, *rest = line.strip().split(maxsplit=3)
                if sid in meta and rest:
                    val = rest[-1]
                    if key == "OS":   # organism species
                        meta[sid]["species"] = val
                    elif key == "AC":
                        meta[sid]["acc"] = val
    return arr, seq_ids, meta