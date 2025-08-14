# src/consite/hmmer_local.py
from __future__ import annotations
import subprocess as sp
from pathlib import Path
from typing import Optional, Tuple, Union, IO

# What subprocess.run accepts for stdout/stderr
LogTarget = Union[int, IO[bytes]]  # int covers sp.DEVNULL; IO[bytes] is an open file


def _open_log_sink(log_path: Optional[Path], quiet: bool) -> Tuple[Optional[LogTarget], Optional[IO[bytes]]]:
    """
    Decide where tool output should go.

    Returns (sink, to_close):
      - sink: what to pass to subprocess.run for stdout/stderr (or None to inherit).
      - to_close: the file object we need to close afterward (if any).
    """
    if log_path is not None:
        f = open(log_path, "ab")  # append binary; safe for subprocess fd redirection
        return f, f
    if quiet:
        return sp.DEVNULL, None
    # inherit parent's stdout/stderr (i.e., let it print to the console)
    return None, None


def _run(cmd: list[str], *, log_path: Optional[Path] = None, quiet: bool = False) -> None:
    """
    Run a command, sending output to a log, to /dev/null (quiet), or to the console.
    """
    sink, to_close = _open_log_sink(log_path, quiet)
    try:
        if sink is None:
            # inherit stdout/stderr
            sp.run(cmd, check=True)
        else:
            sp.run(cmd, check=True, stdout=sink, stderr=sink)
    finally:
        if to_close is not None:
            to_close.close()


def run_hmmsearch(
    seq_fa: Path,
    pfam_hmm: Path,
    domtblout: Path,
    *,
    cpu: int = 4,
    log_path: Optional[Path] = None,
    quiet: bool = False,
) -> None:
    """
    hmmsearch --cut_ga --noali --notextw --cpu N --domtblout <out> <Pfam-A.hmm> <query.fa>
    """
    cmd = [
        "hmmsearch",
        "--cut_ga",      # use Pfam GA thresholds
        "--noali",       # don't print alignments in text output
        "--notextw",     # no fixed line width
        "--cpu", str(cpu),
        "--domtblout", str(domtblout),
        str(pfam_hmm),
        str(seq_fa),
    ]
    _run(cmd, log_path=log_path, quiet=quiet)


def run_hmmbuild(
    seed_sto: Path,
    out_hmm: Path,
    *,
    log_path: Optional[Path] = None,
    quiet: bool = False,
) -> None:
    """
    hmmbuild <out.hmm> <seed.sto>
    """
    cmd = ["hmmbuild", str(out_hmm), str(seed_sto)]
    _run(cmd, log_path=log_path, quiet=quiet)


def run_hmmalign(
    hmm: Path,
    fasta: Path,
    out_sto: Path,
    *,
    log_path: Optional[Path] = None,
    quiet: bool = False,
) -> None:
    """
    hmmalign <model.hmm> <query.fa>  > aligned.sto
    We capture alignment (stdout) to out_sto, while sending stderr to the chosen sink.
    """
    cmd = ["hmmalign", str(hmm), str(fasta)]
    sink, to_close = _open_log_sink(log_path, quiet)
    try:
        with open(out_sto, "wb") as outfh:
            if sink is None:
                # inherit stdout would go to terminal; we override to file and let stderr inherit
                sp.run(cmd, check=True, stdout=outfh)
            else:
                sp.run(cmd, check=True, stdout=outfh, stderr=sink)
    finally:
        if to_close is not None:
            to_close.close()
