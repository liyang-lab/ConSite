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

- Python 3.10 or higher
- HMMER 3.x installed and available in PATH
- Pfam database files (see Quick Start below)

#### Installing HMMER

**macOS (Homebrew):**
```bash
brew install hmmer
```

**Linux (APT):**
```bash
sudo apt-get update
sudo apt-get install hmmer
```

**Windows (conda):**
```bash
conda install -c conda-forge hmmer
```

**Verify installation:**
```bash
hmmsearch --version
```

### From Source (1)

```bash
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite
```

## Quick Start

### Option 1: Automatic Setup (Recommended) (2)

We provide helper scripts to automate the setup process:

```bash
# Make scripts executable
chmod +x scripts/*.sh

# Download and set up Pfam database
./scripts/get_pfam.sh

# Run the demo
./scripts/quickstart.sh
```

**Note:** The scripts have different purposes:
- **`get_pfam.sh`**: Downloads and prepares the Pfam database files
- **`quickstart.sh`**: Sets up the Python environment and runs the demo

### Option 2: Manual Setup (2)

If you prefer to set up things manually or already have some components:

#### 1. Install ConSite

```bash
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
  --pfam-hmm pfam_db/Pfam-A.hmm \
  --pfam-seed pfam_db/Pfam-A.seed \
  --out results \
  --id P05362

# With custom parameters
consite \
  --fasta myprotein.fasta \
  --pfam-hmm pfam_db/Pfam-A.hmm \
  --pfam-seed pfam_db/Pfam-A.seed \
  --out results \
  --topn 5 \
  --cpu 8 \
  --jsd-top-percent 15 \
  --log results/run.log
```

## Output Files (4)

Each run produces a results folder containing:

- **`query.fasta`** - Input sequence used for analysis
- **`hits.json`** - Structured domain hit information
- **`scores.tsv`** - Per-position conservation scores and flags (columns: pos, in_domain, jsd, entropy, is_conserved)
- **`domain_map.png`** - Full sequence domain visualization with conserved sites marked
- **`*_panel.png`** - Individual domain alignment panels showing sequence and conserved positions
- **`*_aligned.sto`** - Stockholm format alignments of query to domain HMMs
- **`hmmsearch.domtblout`** - Raw HMMER domain table output
- **`run.log`** - Complete log of all external tool executions

## Command Line Options (5)

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

## How It Works (6)

1. **Domain Detection**: Uses `hmmsearch` against Pfam-A.hmm to identify conserved domains
2. **SEED Extraction**: Pulls the corresponding Pfam SEED alignment for each hit
3. **HMM Building**: Creates a per-family HMM from the SEED using `hmmbuild`
4. **Sequence Alignment**: Aligns the query protein to the domain HMM using `hmmalign`
5. **Conservation Scoring**: Computes JSD and entropy scores for each position
6. **Visualization**: Generates domain maps and alignment panels

### Current Status

**What works:**
- Complete domain detection and alignment pipeline
- Publication-quality visualizations
- Structured data outputs (JSON, TSV)
- Robust HMMER integration with logging

**Known limitations:**
- Conservation scoring currently uses query-only alignment (JSD/entropy are placeholders)
- Remote CDD mode is not yet implemented
- Conservation thresholds are relative (top X%) rather than absolute

**Next steps planned:**
- Integrate Pfam SEED sequences for real conservation scoring
- Add absolute conservation thresholds
- Implement remote CDD mode

## Example Output (7)

For the included ICAM1 example (P05362), ConSite identifies:
- **PF03921** (positions 25-115): Ig-like domain
- **PF21146** (positions 219-308): Ig-like domain

The tool produces clean visualizations showing domain boundaries and conserved sites as hollow red circles.

### Expected Results

After running the demo, you should see:
- A `results/P05362/` folder containing all outputs
- Two domain panels: `1_PF03921_panel.png` and `2_PF21146_panel.png`
- A full sequence map: `domain_map.png`
- Structured data: `hits.json` and `scores.tsv`

## Troubleshooting (8)

### Common Issues

**"command not found: hmmsearch"**
- Install HMMER using the instructions in Prerequisites above
- Ensure it's in your PATH: `which hmmsearch`

**"No such file or directory: pfam_db/Pfam-A.hmm"**
- Run the helper script: `./scripts/get_pfam.sh`
- Or manually download Pfam files as shown in Manual Setup

**"Permission denied" when running scripts**
- Make scripts executable: `chmod +x scripts/*.sh`

**Large log files**
- Use `--quiet` to suppress verbose output
- Check `--log` path is writable

## Development (9)

### Current Project Structure

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
├── scripts/               # Helper automation scripts
├── pfam_db/              # Pfam database files (not included in repo)
└── results/               # Output directory (not included in repo)
```

### Important Notes for Collaborators

- **Large files are excluded**: The `pfam_db/` and `results/` directories are in `.gitignore`
- **Pfam database**: You'll need to download this using the helper scripts
- **Virtual environment**: The `.venv/` folder is also excluded - you'll create your own
- **Example data**: The `examples/P05362.fasta` file is included for testing

### Dependencies

- **biopython** ≥ 1.81 - Sequence and alignment handling
- **numpy** ≥ 2.0 - Numerical computations
- **matplotlib** ≥ 3.7 - Visualization generation
- **pandas** ≥ 2.0 - Data manipulation
- **scipy** ≥ 1.16 - Scientific computing
- **Python** ≥ 3.10 - Required Python version

## Citation (10)

If you use ConSite in your research, please cite:

```
Joey Wagner, Yang Li. ConSite: conserved-domain alignment and conserved-site visualization from protein FASTA.
```

## License

This project is licensed under the terms specified in the LICENSE file.



