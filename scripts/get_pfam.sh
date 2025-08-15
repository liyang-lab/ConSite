#!/usr/bin/env bash
set -euo pipefail

# Where to place the local database
PFAM_DIR="${1:-pfam_db}"
mkdir -p "$PFAM_DIR"

echo "[1/3] Downloading Pfam-A.hmm.gz and Pfam-A.seed.gz (current release)…"
# Choose one mirror; EBI often hosts current release at:
#   https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/
# Pin a release if you want deterministic behavior.
curl -L -o "$PFAM_DIR/Pfam-A.hmm.gz"  https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
curl -L -o "$PFAM_DIR/Pfam-A.seed.gz" https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz

echo "[2/3] Uncompressing…"
gunzip -f "$PFAM_DIR/Pfam-A.hmm.gz"
gunzip -f "$PFAM_DIR/Pfam-A.seed.gz"

echo "[3/3] hmmpress Pfam-A.hmm…"
hmmpress "$PFAM_DIR/Pfam-A.hmm"

echo "Done. Files in: $PFAM_DIR"
