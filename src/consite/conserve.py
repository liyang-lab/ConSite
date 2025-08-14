from __future__ import annotations
import numpy as np
from typing import Dict, Tuple

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA_SET = set(AA20)
BG = np.array([1/20.0]*20)  # uniform background; could swap for empirical


def column_counts(msa_col: str) -> np.ndarray:
    """Count AA in one MSA column (ignore gaps/unknown). Returns length-20 array.
    """
    counts = np.zeros(20, dtype=float)
    for ch in msa_col:
        if ch in AA_SET:
            counts[AA20.index(ch)] += 1
    return counts


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    m = 0.5*(p+q)
    def kld(a,b):
        mask = (a > 0)
        return np.sum(a[mask] * (np.log2(a[mask]) - np.log2(b[mask])))
    return 0.5*kld(p,m) + 0.5*kld(q,m)


def scores_from_msa(msa: np.ndarray) -> Dict[str, np.ndarray]:
    """msa: array of shape (Nseq, L) with single-letter codes ('-','A',...).
    Returns dict with entropy, jsd, consensus_freq per column.
    """
    L = msa.shape[1]
    ent = np.zeros(L)
    jsd = np.zeros(L)
    cons = np.zeros(L)
    for j in range(L):
        cnt = column_counts(''.join(msa[:, j]))
        total = cnt.sum()
        if total == 0:
            continue
        p = cnt / total
        ent[j] = -(p[p>0] * np.log2(p[p>0])).sum()
        jsd[j] = js_divergence(p, BG)
        cons[j] = p.max()
    # normalize entropy to [0,1] using max log2(20)
    ent_norm = ent / np.log2(20)
    jsd_norm = jsd / max(jsd.max(), 1e-9)
    return {"entropy": ent_norm, "jsd": jsd_norm, "consensus": cons}