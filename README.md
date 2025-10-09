# ConSite

ConSite takes a protein FASTA sequence, finds conserved domains via local Pfam/HMMER, aligns your sequence to each domain, scores per-position conservation, and renders publication-quality figures and structured outputs.

## Features

- **FASTA → Pfam domain search** (local HMMER)
- **Automatic alignment to each hit** using the family's Pfam SEED HMM
- **Per-position conservation**: entropy, Jensen–Shannon divergence (JSD), consensus frequency, coverage
- **Conserved site calling** (top-X% by JSD, per domain)
- **Publication-quality visualization**
  - Linear domain map with hollow-red conserved sites
  - **MSA gradient panels** (from Pfam SEED): brightness-floored grayscale so letters remain legible; supports species labels, optional query row at top, gap glyphs (dash/dot/none), coverage masking/weighting
  - Per-domain alignment panels for the query with optional conservation background scale
- **Reproducible outputs** (JSON, TSV, PNG, Stockholm) and clear CLI logging

## Installation

### Option A: PyPI (recommended)

```bash
python -m pip install consite
```

### Prerequisites

- Python ≥ 3.10
- HMMER 3.x in your `PATH`
- Pfam database files (see Quick Start)

#### Install HMMER:
- **macOS (Homebrew)**: `brew install hmmer`
- **Linux (APT)**: `sudo apt-get update && sudo apt-get install hmmer`
- **Windows (conda)**: `conda install -c conda-forge hmmer`

Verify: `hmmsearch --version`

## Quick Start

### 1) Get the Pfam database

Use the helper script:

```bash
chmod +x scripts/*.sh
./scripts/get_pfam.sh
```

This downloads `Pfam-A.hmm` and `Pfam-A.seed`, uncompresses, and runs `hmmpress`.

(Manual alternative is in "Manual Setup" below.)

### 2) Run ConSite (demo)

```bash
consite \
  --fasta examples/P05362.fasta \
  --pfam-hmm pfam_db/Pfam-A.hmm \
  --pfam-seed pfam_db/Pfam-A.seed \
  --out results \
  --id P05362
```

### Manual Setup (from source)

```bash
git clone https://github.com/yangli-evo/ConSite.git
cd ConSite
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -e .

# (then run ./scripts/get_pfam.sh, or manually download + hmmpress)
```

### Manual Pfam download:

```bash
mkdir -p pfam_db
curl -L -o pfam_db/Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
curl -L -o pfam_db/Pfam-A.seed.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz
gunzip -f pfam_db/Pfam-A.hmm.gz pfam_db/Pfam-A.seed.gz
hmmpress pfam_db/Pfam-A.hmm
```

## Usage

### Basic:

```bash
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

### MSA panel tuned example (SEED gradient = JSD, labels with species+id, include query row, safe brightness):

```bash
consite \
  --fasta examples/GS2.fasta \
  --pfam-hmm pfam_db/Pfam-A.hmm \
  --pfam-seed pfam_db/Pfam-A.seed \
  --out results \
  --id GS2 \
  --msa-panel-nseq 8 \
  --msa-panel-metric jsd \
  --msa-labels species+id \
  --msa-include-query \
  --msa-min-brightness 0.28 \
  --panel-min-brightness 0.22
```

## Outputs

For each run you'll get a folder `results/<id>/` containing:

- **`query.fasta`** – input sequence
- **`hits.json`** – Pfam hits (family, coords, scores)
- **`scores.tsv`** – per-position tracks (columns: `pos`, `in_domain`, `jsd`, `entropy`, `is_conserved`)
- **`domain_map.png`** – linear domain map with conserved sites
- **`*_panel.png`** – per-domain query panels with conserved sites (hollow red)
- **`*_msa.png`** – Pfam SEED MSA panels with grayscale conservation gradient, labels, and optional query row
- **`*_sim.png`** – pairwise % identity heatmap among panel sequences (RF-masked columns)
- **`*_sim.tsv`** – pairwise % identity matrix (TSV)
- **`*_aligned.sto`** – Stockholm alignment of query to each family HMM
- **`hmmsearch.domtblout`** – raw HMMER domain table
- **`run.log`** – external tool logs

## Command-line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--fasta` | input protein FASTA | **Required** |
| `--pfam-hmm` | Path to `Pfam-A.hmm` (pressed) | **Required** |
| `--pfam-seed` | Path to `Pfam-A.seed` | **Required** |
| `--out` | Output directory | **Required** |
| `--id` | Custom run ID (subfolder name) | FASTA header |
| `--topn` | Number of top Pfam hits to analyze | 2 |
| `--cpu` | HMMER threads | 4 |
| `--jsd-top-percent` | Top % (by JSD) called conserved within each domain | 10.0 |
| `--log` | Append external tool output to this file | `results/<id>/run.log` |
| `--quiet` | Suppress console output | False |
| `--keep` | Keep existing output folder (don't overwrite) | False |

### MSA panel (SEED) controls

| Option | Description | Default |
|--------|-------------|---------|
| `--msa-panel-nseq` | Rows to display from the SEED alignment | 8 |
| `--msa-panel-metric` | Gradient metric: `entropy` → uses 1-entropy, or `jsd` | `entropy` |
| `--msa-labels` | Row labels: `id`, `species`, `species+id` | `species+id` |
| `--msa-include-query` | Prepend query row to the MSA panel | off |
| `--msa-min-brightness` | Floor for background brightness (keeps letters legible) | 0.26 |
| `--msa-min-coverage` | Mask columns below this coverage in panels (0–1) | 0.30 |
| `--cons-weight-coverage` | Weight conservation by coverage (alpha) | 1.0 |
| `--mask-inserts` | Use RF to mask insert columns | True |
| `--gap-glyph` | Gap rendering in MSA: `dash`, `dot`, or `none` | `dash` |
| `--gap-cell-brightness` | Brightness used in gap cells | 0.90 |
| `--write-sim-matrix` / `--no-write-sim-matrix` | Write pairwise % identity matrices for MSA panels | True |

### Per-domain panel (query) controls

| Option | Description | Default |
|--------|-------------|---------|
| `--panel-min-brightness` | Brightness floor for query panels | 0.18 |

## How It Works

1. **Domain detection:** `hmmsearch` against Pfam-A HMM library (GA thresholds).
2. **SEED extraction:** the matched family's block is pulled from `Pfam-A.seed`.
3. **Model building:** `hmmbuild` produces a family HMM.
4. **Query alignment:** `hmmalign` aligns your sequence to the family model.
5. **Scoring:**
   - MSA panels use SEED-based per-column scores (JSD or 1-entropy), with optional coverage masking/weighting.
   - Conserved sites are called on the query alignment per domain (top-X% by JSD).
6. **Visualization:** domain map, query panels, and SEED MSA gradient panels.

## Example (ICAM1)

For `examples/P05362.fasta`, ConSite typically finds:

- **PF03921** (~25–115): Ig-like domain
- **PF21146** (~219–308): Ig-like domain

You'll see two `*_panel.png`, two `*_msa.png`, a `domain_map.png`, plus `hits.json` and `scores.tsv`.

## Troubleshooting

- **`command not found: hmmsearch`**
  Install HMMER and ensure it's on your `PATH` ( `which hmmsearch` ).

- **`No such file or directory: pfam_db/Pfam-A.hmm`**
  Run `./scripts/get_pfam.sh` or download manually and `hmmpress`.

- **Large/verbose logs**
  Use `--quiet` and/or inspect `run.log`.

## Development

```bash
ConSite/
├── src/consite/
│   ├── cli.py            # CLI
│   ├── hmmer_local.py    # HMMER wrappers
│   ├── parse_domtbl.py   # HMMER output parsing
│   ├── pfam.py           # Pfam SEED extraction
│   ├── msa_io.py         # Stockholm I/O (+RF)
│   ├── conserve.py       # Scoring (entropy/JSD/coverage)
│   └── viz.py            # Plots (domain map, panels, MSA gradient)
├── scripts/              # Helper scripts (get_pfam, quickstart)
├── examples/             # Example FASTA files
└── pfam_db/              # Pfam files (user-provided)
```

**Dependencies:** biopython ≥1.81, numpy ≥2.0, pandas ≥2.0, scipy ≥1.16, matplotlib ≥3.7, Python ≥3.10.

## Notes & Roadmap

- Remote CDD mode is stubbed (local Pfam/HMMER path is fully supported).
- Current conserved-site calls are **relative** (top-X%); absolute thresholds are planned.

## Citation

**Joey Wagner, Yang Li.** ConSite: conserved-domain alignment and conserved-site visualization from protein FASTA.

## License

MIT (see `LICENSE`).