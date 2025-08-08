# ConSite

**ConSite** is a bioinformatics tool that takes a protein FASTA sequence as input, identifies conserved domains via NCBI's Conserved Domain Database (CDD) or local RPS-BLAST/HMMER searches, detects conserved sites within aligned regions, and outputs both static and interactive visualizations.

## Features
- **FASTA input → conserved domain search** (remote CDD API or local database)
- **Automatic domain alignment** using returned alignment blocks
- **Per-position conservation scoring** (entropy, Jensen–Shannon, consensus frequency)
- **Conserved site detection** with adjustable thresholds
- **Publication-quality visualization**:
  - Linear domain maps with highlighted conserved sites
  - Per-position conservation profiles
  - Interactive alignment viewers
- **Command-line interface (CLI)** and optional **web interface**
- **Reproducible outputs** (JSON, TSV, PNG, HTML report)

## Installation

### From source
```bash
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite
pip install -e .

### Roadmap
MVP CLI — FASTA input → CDD search → conserved-site calling → visualization
Streamlit web interface
Local RPS-BLAST/HMMER integration
AI-based conserved-site prediction (optional)

```bash
consite run \
  --fasta examples/P05362.fasta \
  --remote-cdd \
  --email yangli.evor@email.com \
  --out results/

```bash
consite run \
  --fasta myprotein.fasta \
  --rpsblast-db /path/to/Cdd_LE/ \
  --out results/

```bash
consite run \
  --fasta myprotein.fasta \
  --remote-cdd \
  --top-percent 10 \
  --out results/

Joey Wagner Yang Li. ConSite: conserved-domain alignment and conserved-site visualization from protein FASTA.
