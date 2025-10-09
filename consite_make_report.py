#!/usr/bin/env python3
import json, pathlib, html, csv
from pathlib import Path

def read_first_fasta_header(p: Path) -> tuple[str,int]:
    hdr, length = "?", 0
    if not p.exists(): return hdr, length
    with p.open() as f:
        for line in f:
            if line.startswith(">"):
                hdr = line.strip()[1:]
            else:
                length += len(line.strip())
    return hdr, length

def read_hits(p: Path):
    if not p.exists(): return []
    return json.loads(p.read_text())

def sniff_domains(run_dir: Path):
    panels = sorted(run_dir.glob("*_panel.png"))
    domains = []
    for panel in panels:
        # expect "<idx>_<PFxxxxx>_panel.png"
        stem = panel.stem
        parts = stem.split("_")
        if len(parts) >= 2 and parts[0].isdigit():
            idx = int(parts[0])
            pf = parts[1]
            msa = panel.with_name(f"{idx}_{pf}_msa.png")
            sto = panel.with_name(f"{idx}_{pf}_aligned.sto")
            sim_png = panel.with_name(f"{idx}_{pf}_sim.png")
            sim_tsv = panel.with_name(f"{idx}_{pf}_sim.tsv")
            domains.append({
                "idx": idx,
                "pf": pf,
                "panel": panel.name,
                "msa": msa.name if msa.exists() else None,
                "sto": sto.name if sto.exists() else None,
                "sim_png": sim_png.name if sim_png.exists() else None,
                "sim_tsv": sim_tsv.name if sim_tsv.exists() else None,
            })
    domains.sort(key=lambda d: d["idx"])
    return domains

def small_table_from_scores(scores_tsv: Path, max_rows=200):
    rows = []
    if not scores_tsv.exists(): return rows
    with scores_tsv.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(r):
            if i >= max_rows: break
            rows.append({
                "pos": row["pos"],
                "in_domain": row["in_domain"],
                "jsd": row["jsd"],
                "entropy": row["entropy"],
                "is_conserved": row["is_conserved"],
            })
    return rows

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Make a static HTML report for a ConSite run folder.")
    ap.add_argument("run_dir", type=Path, help="Path to results/<ID>/ (folder containing domain_map, panels, etc.)")
    ap.add_argument("-o", "--out", type=Path, default=None, help="Output HTML path (default: report.html in run_dir)")
    args = ap.parse_args()

    run_dir = args.run_dir.resolve()
    out_html = args.out or (run_dir / "report.html")

    domain_map = (run_dir / "domain_map.png")
    hits_json  = (run_dir / "hits.json")
    scores_tsv = (run_dir / "scores.tsv")
    query_fa   = (run_dir / "query.fasta")
    domtbl     = (run_dir / "hmmsearch.domtblout")

    query_hdr, query_len = read_first_fasta_header(query_fa)
    hits = read_hits(hits_json)
    domains = sniff_domains(run_dir)
    scores_preview = small_table_from_scores(scores_tsv)

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
    /* lightbox */
    #viewer {
      position: fixed; inset: 0; background: rgba(0,0,0,.6);
      display: none; align-items: center; justify-content: center; z-index: 9999;
    }
    #viewer.open { display: flex; }
    #viewer .frame {
      background: #fff; border-radius: 12px; max-width: 96vw; max-height: 90vh;
      box-shadow: 0 10px 30px rgba(0,0,0,.25); display:flex; flex-direction:column;
    }
    #viewer header { display:flex; align-items:center; justify-content:space-between;
      padding:10px 14px; border-bottom:1px solid #eee; }
    #viewer header .ttl { font-weight:600; }
    #viewer button.close { border:0; background:#eef; padding:6px 10px; border-radius:8px; cursor:pointer; }
    #viewer .scroller {
      overflow: auto; /* enables horizontal + vertical scroll */
      padding: 10px;
    }
    #viewer img { height: auto; max-height: 75vh; /* natural width -> horizontal scroll */
      display:block; }
    """

    def esc(s): return html.escape(str(s))

    # Hits table rows
    hits_rows = ""
    for h in hits:
        hits_rows += f"<tr><td class='pf'>{esc(h.get('family',''))}</td><td>{esc(h.get('ali_start',''))}-{esc(h.get('ali_end',''))}</td><td>{esc(h.get('evalue',''))}</td><td>{esc(h.get('score',''))}</td></tr>"

    # Domains blocks
    dom_blocks = []
    for d in domains:
        items = []
        items.append(f"""<div class='card'><div class='muted'>Per-domain panel</div>
          <a class="zoom" href="{esc(d['panel'])}" data-title="Per-domain panel — {esc(d['pf'])}">
            <img src='{esc(d['panel'])}' alt='{esc(d['panel'])}'>
          </a>
        </div>""")
        if d["msa"]:
            items.append(f"""<div class='card'><div class='muted'>SEED MSA panel</div>
          <a class="zoom" href="{esc(d['msa'])}" data-title="SEED MSA panel — {esc(d['pf'])}">
            <img src='{esc(d['msa'])}' alt='{esc(d['msa'])}'>
          </a>
        </div>""")

        # Add similarity matrix if available
        if d["sim_png"]:
            sim_tsv_link = f"<a href='{esc(d['sim_tsv'])}' download>matrix.tsv</a>" if d["sim_tsv"] else ""
            items.append(f"""<div class='card'>
              <div class='muted'>Similarity matrix</div>
              <a class="zoom" href="{esc(d['sim_png'])}" data-title="Similarity matrix — {esc(d['pf'])}">
                <img src='{esc(d['sim_png'])}' alt='{esc(d['sim_png'])}'>
              </a>
              <div class='muted' style="margin-top:8px;">{sim_tsv_link}</div>
            </div>""")

        sto_link = f"<a href='{esc(d['sto'])}' download>{esc(d['sto'])}</a>" if d["sto"] else ""
        dom_blocks.append(f"""
        <section class='section'>
          <h3>Domain {d['idx']}: <span class='pf'>{esc(d['pf'])}</span></h3>
          <div class='grid two'>
            {''.join(items)}
          </div>
          <div class='muted' style="margin-top:8px;">{sto_link}</div>
        </section>
        """)

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
      <a class="zoom" href="{esc(domain_map.name) if domain_map.exists() else ''}" data-title="Domain map">
        <img src="{esc(domain_map.name) if domain_map.exists() else ''}" alt="domain_map">
      </a>
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

<div id="viewer" aria-modal="true" role="dialog">
  <div class="frame">
    <header>
      <div class="ttl" id="viewer-ttl"></div>
      <button class="close" id="viewer-close">Close ✕</button>
    </header>
    <div class="scroller">
      <img id="viewer-img" src="" alt="">
    </div>
  </div>
</div>

<script>
(function () {
  const viewer = document.getElementById('viewer');
  const vimg   = document.getElementById('viewer-img');
  const vttl   = document.getElementById('viewer-ttl');
  const close  = document.getElementById('viewer-close');

  function openViewer(src, title) {
    vimg.src = src;
    vimg.alt = title || '';
    vttl.textContent = title || '';
    viewer.classList.add('open');
    // ensure we start scrolled to the left each time
    viewer.querySelector('.scroller').scrollLeft = 0;
  }
  function dismiss() {
    viewer.classList.remove('open');
    vimg.src = '';
  }

  document.addEventListener('click', function (e) {
    const a = e.target.closest('a.zoom');
    if (!a) return;
    e.preventDefault();
    openViewer(a.getAttribute('href'), a.dataset.title || a.getAttribute('href'));
  });

  close.addEventListener('click', dismiss);
  viewer.addEventListener('click', (e) => { if (e.target === viewer) dismiss(); });
  document.addEventListener('keydown', (e) => { if (e.key === 'Escape') dismiss(); });
})();
</script>
</body>
</html>
"""
    out_html.write_text(html_doc, encoding="utf-8")
    print(f"[OK] Wrote {out_html}")

if __name__ == "__main__":
    main()
