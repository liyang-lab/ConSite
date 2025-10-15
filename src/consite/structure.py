"""
Structure analysis module for ConSite.

Handles ColabFold prediction, Foldseek searches, and conservation mapping to PDB B-factors.
"""
from __future__ import annotations
import subprocess
import shutil
import csv
from pathlib import Path
from typing import Optional
import numpy as np


def run_colabfold(
    fasta: Path,
    out_pdb: Path,
    log: Optional[Path] = None,
    quiet: bool = False,
    colabfold_args: str = "--amber --templates off"
) -> bool:
    """
    Run ColabFold to predict protein structure.

    Args:
        fasta: Input FASTA file
        out_pdb: Output path for the predicted PDB
        log: Optional log file path
        quiet: Suppress stdout/stderr
        colabfold_args: Additional arguments to pass to colabfold_batch

    Returns:
        True if successful, False otherwise
    """
    colabfold_bin = shutil.which("colabfold_batch")
    if not colabfold_bin:
        raise RuntimeError(
            "colabfold_batch not found on PATH. "
            "Install ColabFold (pip install colabfold) or provide --pdb instead."
        )

    # Create temporary directory for ColabFold output
    tmp_dir = out_pdb.parent / "colabfold_tmp"
    tmp_dir.mkdir(exist_ok=True, parents=True)

    try:
        # Build command
        cmd = [
            colabfold_bin,
            str(fasta),
            str(tmp_dir),
        ]
        # Add extra arguments
        if colabfold_args:
            cmd.extend(colabfold_args.split())

        if not quiet:
            print(f"[INFO] Running ColabFold: {' '.join(cmd)}")

        # Run ColabFold
        stdout_dest = subprocess.DEVNULL if quiet else None
        stderr_dest = subprocess.DEVNULL if quiet else None

        if log:
            with log.open("a") as logf:
                logf.write(f"\n=== ColabFold run ===\n")
                logf.write(f"Command: {' '.join(cmd)}\n\n")
                result = subprocess.run(
                    cmd,
                    stdout=logf if quiet else stdout_dest,
                    stderr=logf if quiet else stderr_dest,
                    check=False
                )
        else:
            result = subprocess.run(
                cmd,
                stdout=stdout_dest,
                stderr=stderr_dest,
                check=False
            )

        if result.returncode != 0:
            if not quiet:
                print(f"[ERROR] ColabFold failed with return code {result.returncode}")
            return False

        # Find the top-ranked model (usually *_rank_001_*.pdb or *_relaxed_rank_001_*.pdb)
        pdb_candidates = list(tmp_dir.glob("*_rank_001_*.pdb")) + list(tmp_dir.glob("*_rank_1_*.pdb"))
        if not pdb_candidates:
            pdb_candidates = list(tmp_dir.glob("*.pdb"))

        if not pdb_candidates:
            if not quiet:
                print("[ERROR] No PDB file found in ColabFold output")
            return False

        # Use the first candidate (should be the top-ranked)
        top_pdb = sorted(pdb_candidates)[0]
        shutil.copy(top_pdb, out_pdb)

        if not quiet:
            print(f"[OK] Copied {top_pdb.name} → {out_pdb}")

        return True

    finally:
        # Clean up temporary directory
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)


def run_foldseek(
    query_pdb: Path,
    db: Path,
    out_tsv: Path,
    tmpdir: Path,
    log: Optional[Path] = None,
    quiet: bool = False,
    topk: int = 10,
    evalue: float = 1e-5,
    sensitivity: float = 9.5
) -> bool:
    """
    Run Foldseek to search for structurally similar proteins.

    Args:
        query_pdb: Input PDB structure
        db: Path to Foldseek database
        out_tsv: Output TSV file for results
        tmpdir: Temporary directory for Foldseek
        log: Optional log file path
        quiet: Suppress stdout/stderr
        topk: Number of top hits to return
        evalue: E-value cutoff
        sensitivity: Sensitivity parameter for Foldseek

    Returns:
        True if successful, False otherwise
    """
    foldseek_bin = shutil.which("foldseek")
    if not foldseek_bin:
        raise RuntimeError(
            "foldseek not found on PATH. "
            "Install Foldseek (https://github.com/steineggerlab/foldseek) or use --no-foldseek."
        )

    # Ensure tmpdir exists
    tmpdir.mkdir(exist_ok=True, parents=True)

    # Build command
    cmd = [
        foldseek_bin,
        "easy-search",
        str(query_pdb),
        str(db),
        str(out_tsv),
        str(tmpdir),
        "--format-output", "target,tdesc,evalue,tm,alntmscore,rmsd,qlen,tlen",
        "-e", str(evalue),
        "-s", str(sensitivity),
        "--max-seqs", str(topk)
    ]

    if not quiet:
        print(f"[INFO] Running Foldseek: {' '.join(cmd)}")

    stdout_dest = subprocess.DEVNULL if quiet else None
    stderr_dest = subprocess.DEVNULL if quiet else None

    if log:
        with log.open("a") as logf:
            logf.write(f"\n=== Foldseek search ===\n")
            logf.write(f"Command: {' '.join(cmd)}\n\n")
            result = subprocess.run(
                cmd,
                stdout=logf if quiet else stdout_dest,
                stderr=logf if quiet else stderr_dest,
                check=False
            )
    else:
        result = subprocess.run(
            cmd,
            stdout=stdout_dest,
            stderr=stderr_dest,
            check=False
        )

    if result.returncode != 0:
        if not quiet:
            print(f"[ERROR] Foldseek failed with return code {result.returncode}")
        return False

    if not quiet:
        print(f"[OK] Foldseek results written to {out_tsv}")

    return True


def write_bfactor_from_scores(
    pdb_in: Path,
    scores_tsv: Path,
    pdb_out: Path,
    track: str = "jsd"
) -> bool:
    """
    Write conservation scores to PDB B-factor column.

    Args:
        pdb_in: Input PDB file
        scores_tsv: Path to scores.tsv (contains pos, jsd, entropy columns)
        pdb_out: Output PDB file with B-factors set to conservation scores
        track: Which score to use ('jsd' or 'entropy')

    Returns:
        True if successful, False otherwise
    """
    if not pdb_in.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb_in}")
    if not scores_tsv.exists():
        raise FileNotFoundError(f"Scores TSV not found: {scores_tsv}")

    # Read scores into a dictionary {position: score}
    scores_dict = {}
    with scores_tsv.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos = int(row["pos"])
            if track == "jsd":
                score = float(row["jsd"])
            elif track == "entropy":
                score = float(row["entropy"])
            else:
                raise ValueError(f"Unknown track: {track}. Use 'jsd' or 'entropy'.")
            # Scale to 0-100 for B-factor range
            scores_dict[pos] = score * 100.0

    # Read PDB and update B-factors
    pdb_lines = []
    with pdb_in.open() as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # PDB format: residue number is columns 23-26 (1-indexed)
                try:
                    res_num = int(line[22:26].strip())
                    # B-factor is columns 61-66
                    new_bfactor = scores_dict.get(res_num, 0.0)
                    # Rewrite the line with new B-factor
                    new_line = line[:60] + f"{new_bfactor:6.2f}" + line[66:]
                    pdb_lines.append(new_line)
                except (ValueError, IndexError):
                    # If parsing fails, keep original line
                    pdb_lines.append(line)
            else:
                pdb_lines.append(line)

    # Write output PDB
    with pdb_out.open("w") as f:
        f.writelines(pdb_lines)

    return True


def render_static_pngs(
    pdb_in: Path,
    out_front: Path,
    out_back: Path,
    color_mode: str = "domain",
    domain_colors: Optional[dict] = None,
    log: Optional[Path] = None,
    quiet: bool = False
) -> bool:
    """
    Render static PNG images of the structure (front and back views).

    NOTE: This function requires either PyMOL or py3Dmol for rendering.
    Implementation will be added in a follow-up based on available dependencies.

    Args:
        pdb_in: Input PDB file
        out_front: Output path for front view PNG
        out_back: Output path for back view PNG
        color_mode: 'domain' for domain coloring, 'conservation' for B-factor gradient
        domain_colors: Optional dictionary mapping domain ranges to colors
        log: Optional log file path
        quiet: Suppress stdout/stderr

    Returns:
        True if successful, False otherwise
    """
    # TODO: Implement rendering using PyMOL or py3Dmol
    # This is a placeholder that will be implemented in the planning phase

    if not quiet:
        print("[WARN] Static PNG rendering not yet implemented. See structure.py:render_static_pngs")
        print(f"[WARN] Placeholder files will be created at {out_front} and {out_back}")

    # Create placeholder files for now
    out_front.write_text("# Placeholder for front view PNG\n")
    out_back.write_text("# Placeholder for back view PNG\n")

    return False  # Return False to indicate not implemented


def get_foldseek_hits_summary(foldseek_tsv: Path, max_hits: int = 10) -> list[dict]:
    """
    Read Foldseek results and return a summary of top hits.

    Args:
        foldseek_tsv: Path to Foldseek output TSV
        max_hits: Maximum number of hits to return

    Returns:
        List of dictionaries with hit information
    """
    if not foldseek_tsv.exists():
        return []

    hits = []
    with foldseek_tsv.open() as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "target_id", "target_desc", "evalue", "tm", "alntmscore", "rmsd", "qlen", "tlen"
        ])
        for i, row in enumerate(reader):
            if i >= max_hits:
                break
            hits.append({
                "target_id": row["target_id"],
                "target_desc": row["target_desc"],
                "evalue": row["evalue"],
                "tm": row["tm"],
                "alntmscore": row["alntmscore"],
                "rmsd": row["rmsd"],
                "qlen": row["qlen"],
                "tlen": row["tlen"]
            })
    return hits
