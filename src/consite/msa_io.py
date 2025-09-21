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
    L = arr.shape[1]

    meta: Dict[str, dict] = {sid: {} for sid in seq_ids}
    rf_parts: list[str] = []

    # pass 2: parse GS lines (OS/AC etc.) AND accumulate GC RF across wrapped lines
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#=GS "):
                _, sid, key, *rest = line.strip().split(maxsplit=3)
                if sid in meta and rest:
                    val = rest[-1]
                    if key == "OS":
                        meta[sid]["species"] = val
                    elif key == "AC":
                        meta[sid]["acc"] = val
            elif line.startswith("#=GC RF"):
                # Collect all RF segments; they can be wrapped on multiple lines
                parts = line.strip().split(maxsplit=2)
                if len(parts) >= 3:
                    rf_parts.append(parts[2])

    if rf_parts:
        rf_str = "".join(rf_parts).replace(" ", "")
        # In RF, match columns are typically uppercase letters or 'x'; inserts are '.' (or lowercase)
        rf_mask = np.array([ (c != '.') for c in rf_str ], dtype=bool)
        # Make sure length matches the alignment width
        if rf_mask.size < L:
            # pad with False (treat missing tail as inserts)
            pad = np.zeros(L - rf_mask.size, dtype=bool)
            rf_mask = np.concatenate([rf_mask, pad])
        elif rf_mask.size > L:
            rf_mask = rf_mask[:L]
    else:
        rf_mask = np.ones(L, dtype=bool)  # fallback: treat all as match columns

    return arr, seq_ids, meta, rf_mask