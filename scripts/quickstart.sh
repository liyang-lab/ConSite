#!/usr/bin/env bash
set -euo pipefail

python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .


# Run the demo
consite \
  --fasta examples/P05362.fasta \
  --pfam-hmm pfam_db/Pfam-A.hmm \
  --pfam-seed pfam_db/Pfam-A.seed \
  --out results \
  --id P05362
