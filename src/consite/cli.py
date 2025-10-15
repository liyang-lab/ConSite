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
import html
from typing import Tuple


from .utils import ensure_dir, Hit
from .hmmer_local import run_hmmsearch, run_hmmbuild, run_hmmalign
from .parse_domtbl import parse_domtbl
from .pfam import extract_seed_for_accession
from .msa_io import read_stockholm, read_stockholm_with_meta
from .conserve import scores_from_msa
from .viz import plot_domain_map, plot_alignment_panel, plot_msa_with_gradient, plot_similarity_matrix
from .structure import (
    run_colabfold,
    run_foldseek,
    write_bfactor_from_scores,
    render_static_pngs,
)

def _split_id_range(s: str) -> Tuple[str, str]:
    """Split "A0A...._HUMAN/24-115" -> ("A0A...._HUMAN", "24-115")"""
    m = re.match(r"^(.+?)/(\d+-\d+)$", s)
    return (m.group(1), m.group(2)) if m else (s, "")

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


def _generate_html_report(run_dir: Path, quiet: bool = False) -> None:
    """Generate HTML report by calling the consite_make_report module."""
    import csv
    from .structure import get_foldseek_hits_summary

    domain_map = run_dir / "domain_map.png"
    hits_json = run_dir / "hits.json"
    scores_tsv = run_dir / "scores.tsv"
    query_fa = run_dir / "query.fasta"
    domtbl = run_dir / "hmmsearch.domtblout"

    def read_first_fasta_header(p: Path) -> tuple[str, int]:
        hdr, length = "?", 0
        if not p.exists():
            return hdr, length
        with p.open() as f:
            for line in f:
                if line.startswith(">"):
                    hdr = line.strip()[1:]
                else:
                    length += len(line.strip())
        return hdr, length

    def read_hits(p: Path):
        if not p.exists():
            return []
        return json.loads(p.read_text())

    def sniff_domains(rd: Path):
        panels = sorted(rd.glob("*_panel.png"))
        domains = []
        for panel in panels:
            stem = panel.stem
            parts = stem.split("_")
            if len(parts) >= 2 and parts[0].isdigit():
                idx = int(parts[0])
                pf = parts[1]
                msa = panel.with_name(f"{idx}_{pf}_msa.png")
                sto = panel.with_name(f"{idx}_{pf}_aligned.sto")
                domains.append({
                    "idx": idx,
                    "pf": pf,
                    "panel": panel.name,
                    "msa": msa.name if msa.exists() else None,
                    "sto": sto.name if sto.exists() else None,
                })
        domains.sort(key=lambda d: d["idx"])
        return domains

    def small_table_from_scores(scores_path: Path, max_rows=200):
        rows = []
        if not scores_path.exists():
            return rows
        with scores_path.open() as f:
            r = csv.DictReader(f, delimiter="\t")
            for i, row in enumerate(r):
                if i >= max_rows:
                    break
                rows.append({
                    "pos": row["pos"],
                    "in_domain": row["in_domain"],
                    "jsd": row["jsd"],
                    "entropy": row["entropy"],
                    "is_conserved": row["is_conserved"],
                })
        return rows

    def sniff_structure(rd: Path):
        """Check for structure-related files and return metadata."""
        struct_dir = rd / "structure"
        if not struct_dir.exists():
            return None

        structure_data = {
            "has_structure": False,
            "model_pdb": None,
            "consurf_pdb": None,
            "foldseek_tsv": None,
            "domain_front_png": None,
            "domain_back_png": None,
            "cons_front_png": None,
            "cons_back_png": None,
        }

        for pdb in struct_dir.glob("*_model.pdb"):
            structure_data["model_pdb"] = f"structure/{pdb.name}"
            structure_data["has_structure"] = True
            break
        for pdb in struct_dir.glob("*_model_consurf.pdb"):
            structure_data["consurf_pdb"] = f"structure/{pdb.name}"
            break

        foldseek = struct_dir / "foldseek.tsv"
        if foldseek.exists():
            structure_data["foldseek_tsv"] = f"structure/{foldseek.name}"

        for png in struct_dir.glob("model_domain_front.png"):
            structure_data["domain_front_png"] = f"structure/{png.name}"
        for png in struct_dir.glob("model_domain_back.png"):
            structure_data["domain_back_png"] = f"structure/{png.name}"
        for png in struct_dir.glob("model_cons_front.png"):
            structure_data["cons_front_png"] = f"structure/{png.name}"
        for png in struct_dir.glob("model_cons_back.png"):
            structure_data["cons_back_png"] = f"structure/{png.name}"

        return structure_data if structure_data["has_structure"] else None

    query_hdr, query_len = read_first_fasta_header(query_fa)
    hits = read_hits(hits_json)
    domains = sniff_domains(run_dir)
    scores_preview = small_table_from_scores(scores_tsv)
    structure = sniff_structure(run_dir)

    # Read Foldseek hits if available
    foldseek_hits = []
    if structure and structure.get("foldseek_tsv"):
        foldseek_path = run_dir / structure["foldseek_tsv"]
        foldseek_hits = get_foldseek_hits_summary(foldseek_path)

    # Build HTML
    title = f"ConSite report — {run_dir.name}"
    css = """
    body { font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 24px; color:#111; }
    header { margin-bottom: 20px; }
    h1 { font-size: 1.6rem; margin: 0 0 4px 0; }
    .sub { color:#666; font-size: 0.95rem; }
    .section { margin: 26px 0; }
    .grid { display:grid; gap:14px; }
    .two { grid-template-columns: 1fr 1fr; }
    img { max-width: 100%; height:auto; border:1px solid #e5e7eb; border-radius:10px; }
    .card { border:1px solid #e5e7eb; border-radius:12px; padding:14px; background:#fff; box-shadow:0 1px 2px rgba(0,0,0,0.03); }
    .muted { color:#6b7280; }
    .kvs { display:grid; grid-template-columns: max-content 1fr; gap:6px 12px; }
    .k { color:#6b7280; }
    details summary { cursor:pointer; font-weight:600; }
    table { border-collapse: collapse; width:100%; font-size: 0.9rem; }
    th, td { border-bottom:1px solid #eee; padding:6px 8px; text-align:left; }
    code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, "Liberation Mono", monospace; }
    footer { margin-top: 28px; color:#6b7280; font-size:0.85rem; }
    .pf { font-weight:600; }
    """

    def esc(s):
        return html.escape(str(s))

    # Hits table rows
    hits_rows = ""
    for h in hits:
        hits_rows += f"<tr><td class='pf'>{esc(h.get('family',''))}</td><td>{esc(h.get('ali_start',''))}-{esc(h.get('ali_end',''))}</td><td>{esc(h.get('evalue',''))}</td><td>{esc(h.get('score',''))}</td></tr>"

    # Domains blocks
    dom_blocks = []
    for d in domains:
        items = []
        items.append(
            f"<div class='card'><div class='muted'>Per-domain panel</div><img src='{esc(d['panel'])}' alt='{esc(d['panel'])}'></div>"
        )
        if d["msa"]:
            items.append(
                f"<div class='card'><div class='muted'>SEED MSA panel</div><img src='{esc(d['msa'])}' alt='{esc(d['msa'])}'></div>"
            )
        sto_link = (
            f"<a href='{esc(d['sto'])}' download>{esc(d['sto'])}</a>"
            if d["sto"]
            else ""
        )
        dom_blocks.append(f"""
        <section class='section'>
          <h3>Domain {d['idx']}: <span class='pf'>{esc(d['pf'])}</span></h3>
          <div class='grid two'>
            {''.join(items)}
          </div>
          <div class='muted' style="margin-top:8px;">{sto_link}</div>
        </section>
        """)

    # Structure section
    structure_section = ""
    if structure:
        struct_items = []

        # Interactive viewer (Mol*)
        if structure.get("consurf_pdb"):
            viewer_html = f"""
            <div class='card'>
              <div class='muted'>Interactive 3D viewer (Mol*)</div>
              <div id="molstar-viewer" style="width:100%; height:500px; position:relative;"></div>
              <div class='muted' style="margin-top:8px;">
                Download: <a href="{esc(structure['consurf_pdb'])}" download>model_consurf.pdb</a>
                {' &middot; <a href="' + esc(structure['model_pdb']) + '" download>model.pdb</a>' if structure.get('model_pdb') else ''}
              </div>
            </div>
            """
            struct_items.append(viewer_html)

        # Static renders (if available)
        render_items = []
        if structure.get("domain_front_png"):
            render_items.append(f"""
            <div class='card'>
              <div class='muted'>Domain-colored structure (front)</div>
              <img src='{esc(structure["domain_front_png"])}' alt='Structure front view'>
            </div>
            """)
        if structure.get("cons_front_png"):
            render_items.append(f"""
            <div class='card'>
              <div class='muted'>Conservation-colored structure</div>
              <img src='{esc(structure["cons_front_png"])}' alt='Conservation front view'>
            </div>
            """)

        # Foldseek hits table
        if foldseek_hits:
            foldseek_rows = ""
            for hit in foldseek_hits[:10]:
                target_link = hit['target_id']
                if len(hit['target_id']) == 4 or hit['target_id'].startswith('AF-'):
                    pdb_url = f"https://www.rcsb.org/structure/{hit['target_id'][:4]}"
                    target_link = f"<a href='{pdb_url}' target='_blank'>{esc(hit['target_id'])}</a>"

                foldseek_rows += f"""<tr>
                  <td>{target_link}</td>
                  <td>{esc(hit.get('target_desc', ''))[:60]}</td>
                  <td>{esc(hit.get('evalue', ''))}</td>
                  <td>{esc(hit.get('tm', ''))}</td>
                  <td>{esc(hit.get('rmsd', ''))}</td>
                </tr>"""

            foldseek_table = f"""
            <div class='card'>
              <div class='muted'>Foldseek structural similarity hits</div>
              <div style="max-height:320px; overflow:auto; margin-top:8px;">
                <table>
                  <thead><tr><th>Target</th><th>Description</th><th>E-value</th><th>TM-score</th><th>RMSD</th></tr></thead>
                  <tbody>{foldseek_rows}</tbody>
                </table>
              </div>
              <div class='muted' style="margin-top:8px;">
                Download: <a href="{esc(structure['foldseek_tsv'])}" download>foldseek.tsv</a>
              </div>
            </div>
            """
            struct_items.append(foldseek_table)

        # Build the structure section
        structure_section = f"""
        <section class='section'>
          <h2>Protein Structure</h2>
          <div class='grid {"two" if len(render_items) > 0 else ""}'>
            {''.join(struct_items)}
          </div>
          {('<div class="grid two" style="margin-top:14px;">' + ''.join(render_items) + '</div>') if render_items else ''}
        </section>
        """

    # Scores preview
    score_rows = ""
    for r in scores_preview[:100]:
        score_rows += f"<tr><td>{r['pos']}</td><td>{r['in_domain']}</td><td>{r['jsd']}</td><td>{r['entropy']}</td><td>{r['is_conserved']}</td></tr>"

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{esc(title)}</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>{css}</style>
</head>
<body>
<header>
  <h1>{esc(title)}</h1>
  <div class="sub">{esc(run_dir)}</div>
</header>

<section class="section card">
  <div class="kvs">
    <div class="k">Query</div><div>{esc(query_hdr)}</div>
    <div class="k">Length</div><div>{query_len}</div>
  </div>
</section>

<section class="section">
  <h2>Overview</h2>
  <div class="grid two">
    <div class="card">
      <div class="muted">Domain map</div>
      <img src="{esc(domain_map.name) if domain_map.exists() else ''}" alt="domain_map">
    </div>
    <div class="card">
      <div class="muted">Hits</div>
      <table>
        <thead><tr><th>Pfam</th><th>Aligned range</th><th>i-Evalue</th><th>Score</th></tr></thead>
        <tbody>{hits_rows if hits_rows else "<tr><td colspan='4' class='muted'>No hits</td></tr>"}</tbody>
      </table>
      <div style="margin-top:8px" class="muted">
        Downloads:
        {'<a href="'+esc(hits_json.name)+'" download>hits.json</a>' if hits_json.exists() else ''}
        &nbsp;&middot;&nbsp;
        {'<a href="'+esc(domtbl.name)+'" download>hmmsearch.domtblout</a>' if domtbl.exists() else ''}
        &nbsp;&middot;&nbsp;
        {'<a href="'+esc(scores_tsv.name)+'" download>scores.tsv</a>' if scores_tsv.exists() else ''}
        &nbsp;&middot;&nbsp;
        {'<a href="'+esc(query_fa.name)+'" download>query.fasta</a>' if query_fa.exists() else ''}
      </div>
    </div>
  </div>
</section>

{structure_section}

{"".join(dom_blocks)}

<section class="section card">
  <details>
    <summary>Scores preview (first 100 rows)</summary>
    <div style="margin-top:10px; max-height:320px; overflow:auto;">
      <table>
        <thead><tr><th>pos</th><th>in_domain</th><th>jsd</th><th>entropy</th><th>is_conserved</th></tr></thead>
        <tbody>{score_rows if score_rows else "<tr><td colspan='5' class='muted'>scores.tsv missing</td></tr>"}</tbody>
      </table>
    </div>
  </details>
</section>

<footer>Generated by ConSite static report builder.</footer>

{'<script type="text/javascript" src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"></script>' if structure and structure.get('consurf_pdb') else ''}
{'<link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css" />' if structure and structure.get('consurf_pdb') else ''}

<script>
// Initialize Mol* viewer if structure is available
(function() {{
  const viewerElem = document.getElementById('molstar-viewer');
  if (!viewerElem) return;

  molstar.Viewer.create('molstar-viewer', {{
    layoutIsExpanded: false,
    layoutShowControls: true,
    layoutShowRemoteState: false,
    layoutShowSequence: true,
    layoutShowLog: false,
    layoutShowLeftPanel: true,
    viewportShowExpand: true,
    viewportShowSelectionMode: true,
    viewportShowAnimation: false,
  }}).then(viewer => {{
    viewer.loadPdb('{structure.get("consurf_pdb", "") if structure else ""}');
  }});
}})();
</script>
</body>
</html>
"""
    out_html = run_dir / "report.html"
    out_html.write_text(html_doc, encoding="utf-8")
    if not quiet:
        print(f"[OK] Wrote {out_html}")


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
    msa_panel_nseq: int = 8,
    msa_panel_metric: str = "entropy",
    msa_labels: str = "species+id",
    msa_include_query: bool = False,
    msa_min_brightness: float = 0.25,
    panel_min_brightness: float = 0.18,
    panel_bg: str = "jsd",
    msa_min_coverage: float = 0.3,
    mask_inserts: bool = True,
    gap_glyph: str = "dash",
    gap_cell_brightness: float = 0.9,
    cons_weight_coverage: float = 1.0,
    write_sim_matrix: bool = True,
    html_report: bool = False,
    # Structure analysis parameters
    predict_structure: bool = False,
    pdb: Optional[Path] = None,
    colabfold_args: str = "--amber --templates off",
    foldseek_db: Optional[Path] = None,
    foldseek_topk: int = 10,
    no_foldseek: bool = False,
    cons_to_bfactor: str = "jsd",
    no_structure_renders: bool = False,
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
            seed_msa, seed_ids, seed_meta, rf_mask = read_stockholm_with_meta(seed_path)
            n_show = max(1, min(len(seed_ids), int(msa_panel_nseq)))
            idx = np.linspace(0, len(seed_ids) - 1, n_show, dtype=int)  # simple spread
            msa_sub = seed_msa[idx, :]
            names_sub = []
            for k in idx:
                raw = seed_ids[k]
                if msa_labels == "id":
                    names_sub.append(raw)
                elif msa_labels == "species":
                    base, rng = _split_id_range(raw)
                    species = seed_meta.get(base, {}).get("species") or base.split("_")[-1]
                    label = f"{species}/{rng}" if rng else species
                    names_sub.append(label)
                else:  # species+id (default)
                    base, rng = _split_id_range(raw)
                    species = seed_meta.get(base, {}).get("species") or base.split("_")[-1]
                    label = f"{species}/{rng} ({base})" if rng else f"{species} ({base})"
                    names_sub.append(label)

            seed_scores = scores_from_msa(seed_msa)
            if msa_panel_metric == "jsd":
                col_metric = seed_scores["jsd"]
                title_metric = "bg=JSD"
            else:
                col_metric = 1.0 - seed_scores["entropy"]
                title_metric = "1 - entropy"

            # Weight or mask by coverage
            col_metric = col_metric * (seed_scores["coverage"] ** cons_weight_coverage)

            # Optionally zero out very low-coverage columns outright:
            col_metric[seed_scores["coverage"] < msa_min_coverage] = np.nan

            # Mask insert columns (RF annotation) if enabled
            if mask_inserts:
                col_metric[~rf_mask] = np.nan


            # Build a temporary HMM from the SEED
            fam_hmm = td / f"{h.family}.hmm"
            run_hmmbuild(seed_path, fam_hmm, log_path=log_path, quiet=quiet)

            # Align the query to that model
            q_fa = td / "query.fa"
            SeqIO.write([rec], str(q_fa), "fasta")
            sto = outdir / f"{i}_{h.family}_aligned.sto"
            run_hmmalign(fam_hmm, q_fa, sto, log_path=log_path, quiet=quiet)

            # Read the query alignment with RF
            q_msa, q_ids, q_meta, q_rf = read_stockholm_with_meta(sto)

            # Keep only match columns in each alignment
            # (lengths differ because of inserts, but the number of match states should agree)
            seed_match_mask = rf_mask               # from the SEED block you already parsed
            q_match_mask    = q_rf                  # from hmmalign output

            # Helper: extract match-only array and the model-state index (1..M) for each kept column
            def _keep_match_with_index(msa_arr: np.ndarray, rf: np.ndarray):
                assert msa_arr.shape[1] == rf.size, "RF mask length must equal alignment width"
                kept_idx = np.where(rf)[0]                # column indices to keep
                kept_msa = msa_arr[:, kept_idx]           # match-only
                # model-state ordinal for each match column: 1,2,3,... (used to intersect alignments)
                model_pos = np.cumsum(rf)[rf]             # shape = (#match,)
                return kept_msa, kept_idx, model_pos

            seed_match_msa, seed_kept_idx, seed_pos = _keep_match_with_index(msa_sub, seed_match_mask)
            q_match_msa,    q_kept_idx,    q_pos    = _keep_match_with_index(q_msa,    q_match_mask)

            # Intersect by model-state index to guarantee identical column count & order
            common_pos = np.intersect1d(seed_pos, q_pos)
            seed_common_mask = np.isin(seed_pos, common_pos)
            q_common_mask    = np.isin(q_pos,    common_pos)

            seed_match_msa  = seed_match_msa[:, seed_common_mask]
            q_row_match     = q_match_msa[0:1, q_common_mask]

            # Trim the background metric to the same seed match columns (then to the common set)
            col_metric_match = col_metric[seed_match_mask]
            col_metric_match = col_metric_match[seed_common_mask]

            # Build the final panel MSA
            if msa_include_query:
                final_msa   = np.vstack([q_row_match, seed_match_msa])
                final_names = [f"QUERY: {rec.id}"] + names_sub
            else:
                final_msa   = seed_match_msa
                final_names = names_sub

            # Render with the trimmed metric so columns align
            msa_png = outdir / f"{i}_{h.family}_msa.png"
            plot_msa_with_gradient(
                final_msa, final_names, msa_png,
                title=f"{h.family}  ({title_metric})",
                metric_values=col_metric_match,
                min_brightness=msa_min_brightness,
                gap_glyph=gap_glyph,
                gap_cell_brightness=gap_cell_brightness
            )

            # Similarity matrix (pairwise % identity) for the same rows/cols
            if write_sim_matrix:
                import pandas as pd
                sim_png = outdir / f"{i}_{h.family}_sim.png"
                M = plot_similarity_matrix(final_msa, final_names, sim_png)

                # Also write TSV
                sim_tsv = outdir / f"{i}_{h.family}_sim.tsv"
                pd.DataFrame(M, index=final_names, columns=final_names).to_csv(sim_tsv, sep="\t", float_format="%.2f")

            # Score conservation (JSD/entropy) and call top X% within domain span
            msa, _ = read_stockholm(sto)
            scores = scores_from_msa(msa)

            jsd = scores["jsd"]
            dom_range = np.arange(max(1, h.ali_start), min(seq_len, h.ali_end) + 1)
            if dom_range.size > 0:
                # For conserved site calling, we need to map the domain range to the RF positions
                # Since we don't have RF from hmmalign output, we use the domain positions as-is
                # but could filter by the SEED RF mask if we had a mapping
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

                # For per-domain panel, we use JSD from query alignment
                # but could mask low coverage regions or disable background entirely
                domain_cons_values = scores["jsd"] if panel_bg == "jsd" else None

                plot_alignment_panel(
                    seq=seq_str,
                    hit=h,
                    conserved=set(conserved_local),
                    out_png=panel_png,
                    cons_values=domain_cons_values,  # <— per-position conservation or None
                    cons_clip=(5,95), cons_gamma=0.7, cons_smooth=3, cons_show_scale=True,
                    cons_min_brightness=panel_min_brightness
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

    # Structure analysis pipeline
    if predict_structure or pdb is not None:
        if not quiet:
            print("[INFO] Running structure analysis pipeline...")

        # Create structure directory
        struct_dir = outdir / "structure"
        ensure_dir(struct_dir)

        # Step 1: Get or predict PDB
        model_pdb = struct_dir / f"{safe_id}_model.pdb"
        if pdb is not None:
            # Use provided PDB
            if not pdb.exists():
                if not quiet:
                    print(f"[ERROR] Provided PDB file not found: {pdb}")
            else:
                shutil.copy(pdb, model_pdb)
                if not quiet:
                    print(f"[OK] Using provided PDB: {pdb}")
        elif predict_structure:
            # Run ColabFold
            if not quiet:
                print("[INFO] Running ColabFold structure prediction...")
            success = run_colabfold(
                fasta=seq_fa,
                out_pdb=model_pdb,
                log=log_path,
                quiet=quiet,
                colabfold_args=colabfold_args
            )
            if not success:
                if not quiet:
                    print("[ERROR] ColabFold prediction failed; skipping structure analysis.")
                model_pdb = None

        # Step 2: Map conservation to B-factors
        if model_pdb and model_pdb.exists():
            consurf_pdb = struct_dir / f"{safe_id}_model_consurf.pdb"
            scores_file = outdir / "scores.tsv"
            if scores_file.exists():
                try:
                    write_bfactor_from_scores(
                        pdb_in=model_pdb,
                        scores_tsv=scores_file,
                        pdb_out=consurf_pdb,
                        track=cons_to_bfactor
                    )
                    if not quiet:
                        print(f"[OK] Conservation scores mapped to B-factors: {consurf_pdb}")
                except Exception as e:
                    if not quiet:
                        print(f"[ERROR] Failed to map conservation to B-factors: {e}")

        # Step 3: Run Foldseek if database is provided
        if model_pdb and model_pdb.exists() and foldseek_db and not no_foldseek:
            if not quiet:
                print("[INFO] Running Foldseek structural search...")
            foldseek_tsv = struct_dir / "foldseek.tsv"
            foldseek_tmp = struct_dir / "tmp"
            try:
                success = run_foldseek(
                    query_pdb=model_pdb,
                    db=foldseek_db,
                    out_tsv=foldseek_tsv,
                    tmpdir=foldseek_tmp,
                    log=log_path,
                    quiet=quiet,
                    topk=foldseek_topk
                )
                if success and not quiet:
                    print(f"[OK] Foldseek results: {foldseek_tsv}")
            except Exception as e:
                if not quiet:
                    print(f"[ERROR] Foldseek search failed: {e}")

        # Step 4: Render static PNGs (placeholder for now)
        if model_pdb and model_pdb.exists() and not no_structure_renders:
            if not quiet:
                print("[INFO] Rendering static structure images...")
            front_png = struct_dir / "model_domain_front.png"
            back_png = struct_dir / "model_domain_back.png"
            try:
                render_static_pngs(
                    pdb_in=model_pdb,
                    out_front=front_png,
                    out_back=back_png,
                    color_mode="domain",
                    log=log_path,
                    quiet=quiet
                )
            except Exception as e:
                if not quiet:
                    print(f"[WARN] Static rendering not yet implemented: {e}")

    # Generate HTML report if requested
    if html_report:
        if not quiet:
            print("[INFO] Generating HTML report...")
        _generate_html_report(outdir, quiet=quiet)


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
    p.add_argument("--msa-labels", choices=["id", "species", "species+id"],
                   default="species+id",
                   help="How to label MSA rows.")
    p.add_argument("--msa-include-query", action="store_true",
                   help="Add the query as the first row in the MSA panel.")
    p.add_argument("--msa-min-brightness", type=float, default=0.25,
                   help="Floor for background brightness in MSA panels (0..1).")
    p.add_argument("--panel-min-brightness", type=float, default=0.18,
                   help="Floor for background brightness in per-domain panels (0..1).")
    p.add_argument("--panel-bg", choices=["none","jsd"], default="jsd",
                   help="Background conservation values for per-domain panels.")

    # NEW: Coverage and insert handling
    p.add_argument("--msa-min-coverage", type=float, default=0.3,
                   help="Mask columns below this coverage in MSA panels (0..1).")
    p.add_argument("--mask-inserts", action="store_true", default=True,
                   help="Use RF to mask inserts everywhere (default: True).")
    p.add_argument("--gap-glyph", choices=["dash", "dot", "none"], default="dash",
                   help="Glyph to show for gaps in MSA: dash (–), dot (·), or none.")
    p.add_argument("--gap-cell-brightness", type=float, default=0.9,
                   help="Brightness for gap cells (0..1, higher = brighter).")
    p.add_argument("--cons-weight-coverage", type=float, default=1.0,
                   help="Alpha for coverage weighting in conservation (1.0 = linear).")

    # Similarity matrix
    p.add_argument("--write-sim-matrix", action="store_true", default=True,
                   help="Write pairwise % identity matrices for MSA panels (default: True).")
    p.add_argument("--no-write-sim-matrix", dest="write_sim_matrix", action="store_false",
                   help="Do not write similarity matrices.")

    # HTML report
    p.add_argument("--html-report", action="store_true",
                   help="Generate a static HTML report after pipeline completes.")

    # Structure prediction and analysis
    p.add_argument("--predict-structure", action="store_true",
                   help="Run ColabFold to predict protein structure.")
    p.add_argument("--pdb", type=Path, default=None,
                   help="Use existing PDB file (skip structure prediction).")
    p.add_argument("--colabfold-args", type=str, default="--amber --templates off",
                   help="Additional arguments to pass to colabfold_batch.")

    p.add_argument("--foldseek-db", type=Path, default=None,
                   help="Path to Foldseek database for structural similarity search.")
    p.add_argument("--foldseek-topk", type=int, default=10,
                   help="Number of top Foldseek hits to report.")
    p.add_argument("--no-foldseek", action="store_true",
                   help="Skip Foldseek search even if PDB is present.")

    p.add_argument("--cons-to-bfactor", choices=["jsd", "entropy"], default="jsd",
                   help="Map conservation scores to PDB B-factors (jsd or entropy).")
    p.add_argument("--no-structure-renders", action="store_true",
                   help="Skip generating static PNG renders of structure.")

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
        msa_panel_nseq=args.msa_panel_nseq,
        msa_panel_metric=args.msa_panel_metric,
        msa_labels=args.msa_labels,
        msa_include_query=args.msa_include_query,
        msa_min_brightness=args.msa_min_brightness,
        panel_min_brightness=args.panel_min_brightness,
        panel_bg=args.panel_bg,
        msa_min_coverage=args.msa_min_coverage,
        mask_inserts=args.mask_inserts,
        gap_glyph=args.gap_glyph,
        gap_cell_brightness=args.gap_cell_brightness,
        cons_weight_coverage=args.cons_weight_coverage,
        write_sim_matrix=args.write_sim_matrix,
        html_report=args.html_report,
        # Structure analysis parameters
        predict_structure=args.predict_structure,
        pdb=args.pdb,
        colabfold_args=args.colabfold_args,
        foldseek_db=args.foldseek_db,
        foldseek_topk=args.foldseek_topk,
        no_foldseek=args.no_foldseek,
        cons_to_bfactor=args.cons_to_bfactor,
        no_structure_renders=args.no_structure_renders,
    )


if __name__ == "__main__":
    main()
