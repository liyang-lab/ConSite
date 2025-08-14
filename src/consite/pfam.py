from __future__ import annotations
from pathlib import Path
from typing import Optional
import re

_AC_LINE = re.compile(r'^#=GF\s+AC\s+(PF\d+)(?:\.\d+)?\s*$')

def extract_seed_for_accession(pfam_seed: Path, accession: str, out_sto: Path) -> Optional[Path]:
    """
    Pull a single family block from Pfam-A.seed by accession (PFxxxxx),
    writing the full Stockholm block to `out_sto`. Returns `out_sto` if found,
    else None.

    The Pfam-A.seed file is a concatenation of Stockholm blocks, each ending
    with a line:  //
    """
    want = accession.split('.')[0]  # tolerate version numbers in source/target
    in_block = False
    curr_acc: Optional[str] = None
    block: list[str] = []

    with pfam_seed.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            # Each block starts with the Stockholm header
            if line.startswith("# STOCKHOLM 1.0"):
                in_block = True
                curr_acc = None
                block = [line]
                continue

            if not in_block:
                continue

            block.append(line)

            if line.startswith("#=GF AC"):
                m = _AC_LINE.match(line)
                if m:
                    curr_acc = m.group(1)

            if line.strip() == "//":
                # End of the current block: write if it’s the one we want
                if curr_acc == want:
                    out_sto.write_text("".join(block), encoding="utf-8")
                    return out_sto
                # otherwise reset for the next block
                in_block = False
                curr_acc = None
                block = []

    return None
