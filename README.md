# ConSite

**ConSite** is a bioinformatics tool that takes a protein FASTA sequence as input, identifies conserved domains via local Pfam/HMMER searches, detects conserved sites within aligned regions, and outputs both structured data and publication-quality visualizations.

## Features

- **FASTA input → conserved domain search** using local Pfam database and HMMER
- **Automatic domain alignment** using Pfam SEED alignments
- **Per-position conservation scoring** (entropy, Jensen–Shannon divergence, consensus frequency)
- **Conserved site detection** with adjustable thresholds
- **Publication-quality visualization**:
  - Linear domain maps with highlighted conserved sites
  - Per-domain alignment panels with legible sequence display
  - Hollow red circles marking conserved positions
- **Command-line interface (CLI)** with comprehensive logging
- **Reproducible outputs** (JSON, TSV, PNG, Stockholm alignments)

## Installation

### Prerequisites

- Python 3.9 or higher
- HMMER 3.x installed and available in PATH
- Pfam database files (see Quick Start below)

### From Source

```bash
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite
pip install -e .
```

## Quick Start

### Option 1: Automatic Setup (Recommended)

We provide helper scripts to automate the entire setup process:

```bash
# Clone and enter the repository
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite

# Run the complete quickstart script
chmod +x scripts/quickstart.sh
./scripts/quickstart.sh
```

This script will:
1. Create a Python virtual environment
2. Install ConSite in development mode
3. Download and set up the Pfam database automatically
4. Run the demo with the included PO5362 example

### Option 2: Manual Setup

If you prefer to set up things manually or already have some components:

#### 1. Install ConSite

```bash
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e .
```

#### 2. Download Pfam Database

```bash
# Create directory for Pfam files
mkdir -p pfam_db

# Download Pfam-A HMM library
curl -L -o pfam_db/Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip pfam_db/Pfam-A.hmm.gz

# Download Pfam-A SEED alignments
curl -L -o pfam_db/Pfam-A.seed.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz
gunzip pfam_db/Pfam-A.seed.gz

# Press the HMM library for HMMER
hmmpress pfam_db/Pfam-A.hmm
```

### 3. Run ConSite

```bash
# Basic run with example protein
consite \
  --fasta examples/P05362.fasta \
  --pfam-hmm Pfam-A.hmm \
  --pfam-seed Pfam-A.seed \
  --out results \
  --id P05362

# With custom parameters
consite \
  --fasta myprotein.fasta \
  --pfam-hmm Pfam-A.hmm \
  --pfam-seed Pfam-A.seed \
  --out results \
  --topn 5 \
  --cpu 8 \
  --jsd-top-percent 15 \
  --log results/run.log
```

## Output Files

Each run produces a results folder containing:

- **`query.fasta`** - Input sequence used for analysis
- **`hits.json`** - Structured domain hit information
- **`scores.tsv`** - Per-position conservation scores and flags
- **`domain_map.png`** - Full sequence domain visualization
- **`*_panel.png`** - Individual domain alignment panels
- **`*_aligned.sto`** - Stockholm format alignments
- **`hmmsearch.domtblout`** - Raw HMMER domain table output

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--fasta` | Input protein FASTA file | Required |
| `--pfam-hmm` | Path to Pfam-A.hmm (pressed) | Required |
| `--pfam-seed` | Path to Pfam-A.seed | Required |
| `--out` | Output directory | Required |
| `--id` | Custom run ID (default: FASTA header) | Auto-detected |
| `--topn` | Number of top domains to analyze | 2 |
| `--cpu` | Number of CPU cores for HMMER | 4 |
| `--jsd-top-percent` | Top % of positions called conserved | 10.0 |
| `--log` | Log file for external tool output | `results/<id>/run.log` |
| `--quiet` | Suppress console output | False |
| `--keep` | Preserve existing output folder | False |

## How It Works

1. **Domain Detection**: Uses `hmmsearch` against Pfam-A.hmm to identify conserved domains
2. **SEED Extraction**: Pulls the corresponding Pfam SEED alignment for each hit
3. **HMM Building**: Creates a per-family HMM from the SEED using `hmmbuild`
4. **Sequence Alignment**: Aligns the query protein to the domain HMM using `hmmalign`
5. **Conservation Scoring**: Computes JSD and entropy scores for each position
6. **Visualization**: Generates domain maps and alignment panels

## Example Output

For the included ICAM1 example (P05362), ConSite identifies:
- **PF03921** (positions 25-115): Ig-like domain
- **PF21146** (positions 219-308): Ig-like domain

The tool produces clean visualizations showing domain boundaries and conserved sites as hollow red circles.

## Helper Scripts

We provide several helper scripts to automate common tasks:

### `scripts/quickstart.sh`
Complete one-command setup that:
- Creates Python virtual environment
- Installs ConSite
- Downloads and sets up Pfam database
- Runs the demo example

### `scripts/get_pfam.sh`
Downloads and prepares the Pfam database:
- Downloads current Pfam-A.hmm and Pfam-A.seed
- Uncompresses the files
- Runs `hmmpress` to prepare HMMER files

Usage:
```bash
./scripts/get_pfam.sh [directory]
# Default: downloads to ./pfam_db
# Custom: ./scripts/get_pfam.sh /path/to/custom/location
```

## Development

```
ConSite/
├── src/consite/          # Main package source
│   ├── cli.py            # Command-line interface
│   ├── hmmer_local.py    # HMMER tool wrappers
│   ├── parse_domtbl.py   # HMMER output parsing
│   ├── pfam.py           # Pfam SEED extraction
│   ├── msa_io.py         # Multiple sequence alignment I/O
│   ├── conserve.py       # Conservation scoring algorithms
│   ├── viz.py            # Visualization functions
│   └── utils.py          # Shared utilities
├── examples/              # Example input files
├── pfam_db/              # Pfam database files (not included)
└── results/               # Output directory
```

### Dependencies

- **biopython** ≥ 1.81 - Sequence and alignment handling
- **numpy** ≥ 1.24 - Numerical computations
- **matplotlib** ≥ 3.7 - Visualization generation
- **pandas** ≥ 2.0 - Data manipulation
- **scipy** ≥ 1.10 - Scientific computing
- **seaborn** ≥ 0.12 - Statistical visualization

## Citation

If you use ConSite in your research, please cite:

```
Joey Wagner, Yang Li. ConSite: conserved-domain alignment and conserved-site visualization from protein FASTA.
```

## License

This project is licensed under the terms specified in the LICENSE file.



