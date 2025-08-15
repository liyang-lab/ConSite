from __future__ import annotations
from pathlib import Path
from typing import List
from .utils import Hit

def parse_domtbl(domtblout: Path, topn: int = 3) -> List[Hit]:
    """
    Parse HMMER --domtblout produced by `hmmsearch <HMMs> <seqs>`.
    Columns (HMMER v3.x): tname tacc tlen qname qacc qlen ... i-Evalue score ... hmm_from hmm_to ali_from ali_to ...
    We want the Pfam accession on the query (qacc) and the aligned coords on the target (ali_from/ali_to).
    """
    hits: List[Hit] = []
    with domtblout.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split()
            # target sequence (the protein) is columns 0/1
            # query model (Pfam) is columns 3/4
            qname = cols[3]          # e.g. Zn_ribbon_Top1
            qacc  = cols[4]          # e.g. PF03921.21 or "-"
            # alignment coords on the target sequence (1-based, inclusive)
            ali_start = int(cols[17])
            ali_end   = int(cols[18])
            # keep i-Evalue/score for ranking
            evalue = float(cols[12]) # i-Evalue
            score  = float(cols[13]) # score

            # Prefer the accession (PFxxxxx); fall back to the model name if missing
            family = qacc.split(".")[0] if qacc != "-" else qname
            hits.append(Hit(family, cols[0], ali_start, ali_end, evalue, score))

    hits.sort(key=lambda h: (h.evalue, -h.score))
    return hits[:topn]
