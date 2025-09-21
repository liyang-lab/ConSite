from __future__ import annotations
from typing import Tuple, Dict
import numpy as np
from Bio import AlignIO


def read_stockholm(path) -> Tuple[np.ndarray, list]:
    aln = AlignIO.read(str(path), "stockholm")
    seq_ids = [rec.id for rec in aln]
    arr = np.array([list(str(rec.seq).upper()) for rec in aln], dtype='<U1')
    return arr, seq_ids

def read_stockholm_with_meta(path) -> Tuple[np.ndarray, list, Dict[str, dict], np.ndarray]:
    aln = AlignIO.read(str(path), "stockholm")
    seq_ids = [rec.id for rec in aln]
    arr = np.array([list(str(rec.seq).upper()) for rec in aln], dtype='<U1')

    meta: Dict[str, dict] = {sid: {} for sid in seq_ids}
    rf_mask = None

    # pass 2: parse GS lines for OS (organism), AC, etc. and GC RF line
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
            elif line.startswith("#=GC RF "):
                # Parse reference (match) annotation
                rf_line = line.strip().split(maxsplit=2)
                if len(rf_line) >= 3:
                    rf_annotation = rf_line[2]
                    # RF annotation: 'x' or uppercase = match column, '.' or lowercase = insert
                    rf_mask = np.array([c in 'x' or c.isupper() for c in rf_annotation], dtype=bool)

    # If no RF found, assume all columns are match columns
    if rf_mask is None:
        rf_mask = np.ones(arr.shape[1], dtype=bool)

    return arr, seq_ids, meta, rf_mask