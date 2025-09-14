from __future__ import annotations
import json
import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional, List
from Bio import SeqIO
import numpy as np
import re
import shutil

from .utils import ensure_dir, Hit
from .hmmer_local import run_hmmsearch, run_hmmbuild, run_hmmalign
from .parse_domtbl import parse_domtbl
from .pfam import extract_seed_for_accession
from .msa_io import read_stockholm
from .conserve import scores_from_msa
from .viz import plot_domain_map, plot_alignment_panel, plot_msa_with_gradient  # <-- added import

def _ensure_hmmer_or_exit():
    missing = [t for t in ("hmmsearch","hmmbuild","hmmalign") if shutil.which(t) is None]
    if missing:
        tools = ", ".join(missing)
        raise SystemExit(f"[ERROR] Required HMMER tool(s) not found on PATH: {tools}. "
                         "Install HMMER 3.x (brew/apt/conda) and try again.")


def _write_scores_tsv(
    seq_len: int,
    jsd_global: np.ndarray,
    entropy_global: np.ndarray,
    conserved: set[int],
    hits: list[Hit],
    out_tsv: Path,
) -> None:
    """Write per-position scores and indicators to TSV."""
    in_domain = np.zeros(seq_len, dtype=bool)
    for h in hits:
        a, b = max(1, h.ali_start), min(seq_len, h.ali_end)
        if a <= b:
            in_domain[a - 1 : b] = True

    with out_tsv.open("w") as f:
        f.write("pos\tin_domain\tjsd\tentropy\tis_conserved\n")
        for pos in range(1, seq_len + 1):
            f.write(
                f"{pos}\t{int(in_domain[pos-1])}\t"
                f"{float(jsd_global[pos-1]):.6g}\t"
                f"{float(entropy_global[pos-1]):.6g}\t"
                f"{int(pos in conserved)}\n"
            )


def run_pipeline(
    fasta: Path,
    outdir: Path,
    pfam_hmm: Optional[Path] = None,
    pfam_seed: Optional[Path] = None,
    *,
    remote_cdd: bool = False,
    email: Optional[str] = None,
    topn: int = 2,
    cpu: int = 4,
    jsd_top_percent: float = 10.0,
    log: Optional[Path] = None,
    quiet: bool = False,
    run_id: Optional[str] = None,
    keep: bool = False,
    msa_panel_nseq: int = 8,                 # <-- added
    msa_panel_metric: str = "entropy",       # <-- added ("entropy" means use 1-entropy for shading)
) -> None:
    """Run either local Pfam/HMMER pipeline or (placeholder) remote CDD mode."""
    ensure_dir(outdir)

    # Decide run id (default: FASTA header’s first token) and sanitize
    records_peek = list(SeqIO.parse(str(fasta), "fasta"))
    if not records_peek:
        raise SystemExit("No sequences in FASTA")
    default_id = records_peek[0].id.split()[0]
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", run_id or default_id)

    # Make the final output dir: <out>/<safe_id>
    outdir = (outdir / safe_id).resolve()
    if outdir.exists() and not keep:
        shutil.rmtree(outdir)
    ensure_dir(outdir)

    # Reset FASTA iterator (we consumed it for the peek above)
    records = [records_peek[0]]
    if len(records_peek) > 1 and not quiet:
        print("[WARN] Multiple sequences provided; using the first.")
    rec = records[0]
    seq_str = str(rec.seq)                 # keep full sequence string
    seq_len = len(seq_str)
    seq_fa = outdir / "query.fasta"
    SeqIO.write([rec], str(seq_fa), "fasta")

    # Remote CDD (placeholder so CLI runs without Pfam files)
    if remote_cdd:
        if not quiet:
            print("[INFO] Remote CDD mode selected.")
            if email:
                print(f"[INFO] Using email: {email}")
        (outdir / "hits.json").write_text("[]")
        plot_domain_map(seq_len, [], [], outdir / "domain_map.png")
        if not quiet:
            print("[INFO] Remote CDD results not implemented in this build; "
                  "provide --pfam-hmm/--pfam-seed for local Pfam/HMMER mode.")
        return

    # Local Pfam/HMMER
    if pfam_hmm is None or pfam_seed is None:
        raise SystemExit("Local mode requires --pfam-hmm and --pfam-seed (or use --remote-cdd).")

    domtbl = outdir / "hmmsearch.domtblout"
    log_path = log or (outdir / "run.log")
    run_hmmsearch(seq_fa, pfam_hmm, domtbl, cpu=cpu, log_path=log_path, quiet=quiet)

    hits: List[Hit] = parse_domtbl(domtbl, topn=topn)
    (outdir / "hits.json").write_text(json.dumps([h.__dict__ for h in hits], indent=2))

    # If no hits, still draw an empty domain map and exit gracefully
    if not hits:
        if not quiet:
            print("[INFO] hmmsearch returned no reportable domains.")
        plot_domain_map(seq_len, [], [], outdir / "domain_map.png")
        return

    # accumulate per-position tracks across domains
    jsd_global = np.zeros(seq_len, dtype=float)
    entropy_global = np.zeros(seq_len, dtype=float)

    conserved_positions_global: List[int] = []
    total = len(hits)

    # Per-hit: SEED → hmmbuild → hmmalign → score → call conserved sites
    for i, h in enumerate(hits, 1):
        if not quiet:
            print(f"[{i}/{total}] {h.family}  ali:{h.ali_start}-{h.ali_end}", end="\r", flush=True)

        with TemporaryDirectory() as td_str:
            td = Path(td_str)

            # Extract SEED for this Pfam accession
            seed_path = td / f"{h.family}.seed.sto"
            seed_ok = extract_seed_for_accession(pfam_seed, h.family, seed_path) is not None
            if not seed_ok:
                if not quiet:
                    print(f"\n[WARN] SEED not found for {h.family}; skipping.")
                continue

            # ---- NEW: MSA gradient panel from SEED ----
            seed_msa, seed_ids = read_stockholm(seed_path)
            n_show = max(1, min(len(seed_ids), int(msa_panel_nseq)))
            idx = np.linspace(0, len(seed_ids) - 1, n_show, dtype=int)  # simple spread
            msa_sub = seed_msa[idx, :]
            names_sub = [seed_ids[k] for k in idx]

            seed_scores = scores_from_msa(seed_msa)
            if msa_panel_metric == "jsd":
                col_metric = seed_scores["jsd"]                    # 0..1
                title_metric = "bg=JSD"
            else:
                col_metric = 1.0 - seed_scores["entropy"]          # 1 - entropy in 0..1
                title_metric = "1 - entropy"

            msa_png = outdir / f"{i}_{h.family}_msa.png"
            plot_msa_with_gradient(
                msa_sub, names_sub, msa_png,
                title=f"{h.family}  ({title_metric})",
                metric_values=col_metric
            )
            # ---- END NEW ----

            # Build a temporary HMM from the SEED
            fam_hmm = td / f"{h.family}.hmm"
            run_hmmbuild(seed_path, fam_hmm, log_path=log_path, quiet=quiet)

            # Align the query to that model
            q_fa = td / "query.fa"
            SeqIO.write([rec], str(q_fa), "fasta")
            sto = outdir / f"{i}_{h.family}_aligned.sto"
            run_hmmalign(fam_hmm, q_fa, sto, log_path=log_path, quiet=quiet)

            # Score conservation (JSD/entropy) and call top X% within domain span
            msa, _ = read_stockholm(sto)
            scores = scores_from_msa(msa)

            jsd = scores["jsd"]
            dom_range = np.arange(max(1, h.ali_start), min(seq_len, h.ali_end) + 1)
            if dom_range.size > 0:
                vals = jsd[dom_range - 1]
                k = max(1, int(len(vals) * (jsd_top_percent / 100.0)))
                thr = np.partition(vals, -k)[-k]
                conserved_local = dom_range[vals >= thr].tolist()
                conserved_positions_global.extend(conserved_local)

                # Update global tracks inside this domain span
                entropy = scores.get("entropy", np.zeros(seq_len))
                jsd_global[dom_range - 1] = np.maximum(jsd_global[dom_range - 1], jsd[dom_range - 1])
                entropy_global[dom_range - 1] = np.maximum(entropy_global[dom_range - 1], entropy[dom_range - 1])

                # Render a per-domain panel PNG (correct signature)
                panel_png = outdir / f"{i}_{h.family}_panel.png"
                plot_alignment_panel(
                    seq=seq_str,
                    hit=h,
                    conserved=set(conserved_local),
                    out_png=panel_png,
                    cons_values=scores["jsd"],  # <— per-position conservation
                    cons_clip=(5,95), cons_gamma=0.7, cons_smooth=3, cons_show_scale=True
                )

    if not quiet and total:
        print()  # newline after the progress line

    # Write per-position scores/flags
    _write_scores_tsv(
        seq_len=seq_len,
        jsd_global=jsd_global,
        entropy_global=entropy_global,
        conserved=set(conserved_positions_global),
        hits=hits,
        out_tsv=outdir / "scores.tsv",
    )

    plot_domain_map(seq_len, hits, conserved_positions_global, outdir / "domain_map.png")
    if len(conserved_positions_global) == 0 and not quiet:
        print("[INFO] No conserved positions were called (check JSD cutoff or alignment).")


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ConSite CLI")
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)

    # Remote CDD (optional)
    p.add_argument("--remote-cdd", action="store_true",
                   help="Use remote NCBI CD-Search instead of local Pfam/HMMER.")
    p.add_argument("--email", default=None,
                   help="Contact email for remote CDD submissions (recommended).")

    # Local Pfam/HMMER (optional; required if not using --remote-cdd)
    p.add_argument("--pfam-hmm", type=Path, default=None,
                   help="Path to Pfam-A.hmm (pressed).")
    p.add_argument("--pfam-seed", type=Path, default=None,
                   help="Path to Pfam-A.seed (Stockholm).")

    p.add_argument("--topn", type=int, default=2)
    p.add_argument("--cpu", type=int, default=4)
    p.add_argument("--jsd-top-percent", type=float, default=10.0)

    # NEW: MSA gradient panel controls
    p.add_argument("--msa-panel-nseq", type=int, default=8,
                   help="Rows to show in the MSA gradient panel (from SEED).")
    p.add_argument("--msa-panel-metric", choices=["entropy", "jsd"], default="entropy",
                   help="Column metric for gradient: 1-entropy (default) or JSD.")

    # Logging / verbosity
    p.add_argument("--log", type=Path, default=None, help="Append external tool logs here.")
    p.add_argument("--quiet", action="store_true", help="Suppress tool stdout/stderr.")
    p.add_argument("--id", dest="run_id", default=None,
                   help="Subfolder name under --out (default: FASTA record id)")
    p.add_argument("--keep", action="store_true",
                   help="Do not delete an existing output folder (default: overwrite)")
    return p


def main():
    args = build_argparser().parse_args()

    # allow either remote CDD OR local Pfam/HMMER
    if not args.remote_cdd and not (args.pfam_hmm and args.pfam_seed):
        raise SystemExit(
            "Either use --remote-cdd (remote CDD mode) OR provide both "
            "--pfam-hmm and --pfam-seed for local Pfam/HMMER mode."
        )
    
    if not args.remote_cdd:
        _ensure_hmmer_or_exit()

    run_pipeline(
        fasta=args.fasta,
        outdir=args.out,
        pfam_hmm=args.pfam_hmm,
        pfam_seed=args.pfam_seed,
        remote_cdd=args.remote_cdd,
        email=args.email,
        topn=args.topn,
        cpu=args.cpu,
        jsd_top_percent=args.jsd_top_percent,
        log=args.log,
        quiet=args.quiet,
        run_id=args.run_id,
        keep=args.keep,
        msa_panel_nseq=args.msa_panel_nseq,             # <-- pass through
        msa_panel_metric=args.msa_panel_metric,         # <-- pass through
    )


if __name__ == "__main__":
    main()
