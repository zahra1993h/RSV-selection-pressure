"""
Selection Pressure Analysis Pipeline (Section 7.5) — v6 (RESUMABLE)

WHAT CHANGED FROM v5, AND WHY RSV-B WAS BROKEN
==============================================
Three genuine bugs in v5 explain the RSV-B failure:

  BUG 1 (crash) -- `--run` raised AttributeError immediately.
      run_all(..., dataset=args.dataset) was called, but `--dataset` was
      never registered with argparse. Fixed: --dataset and --gene now exist.

  BUG 2 (the "defective" results) -- the resume logic was unsafe:
          if out_json.exists() and not force:  skip
      HyPhy creates the output file early and fills it as it goes. If a job
      was killed (OOM, closed terminal, laptop sleep), a TRUNCATED .json was
      left behind. Every later run then SKIPPED it because the file existed,
      so the corruption became permanent and --parse read garbage.
      Fixed: jobs write to <name>.json.tmp and are renamed to <name>.json
      only after HyPhy exits 0 AND the JSON parses AND the number of rows
      matches the codon count of the alignment. A .json present here is a
      .json you can trust. Broken files are detected and rerun automatically,
      no --force needed.

  BUG 3 (CPU oversubscription -- this is why MEME felt so slow) --
      conda's HyPhy is the HYPHYMP build: a SINGLE job already grabs every
      core via OpenMP. v5 launched 3 such jobs concurrently, so on an 8-core
      machine 24 threads fought over 8 cores. Fixed: each job is given an
      explicit CPU=N share, sized so parallel x CPU <= total cores.

OTHER IMPROVEMENTS
  * --status : shows exactly how far each dataset/gene/method has got,
    with a progress bar, before you commit to another long run.
  * Global longest-first job queue across ALL datasets instead of one
    dataset at a time, so cores never sit idle waiting for a slow MEME.
  * Per-job log files (hyphy_results/<dataset>/<gene>.<METHOD>.log) instead
    of one shared errors.log that parallel threads were deleting from
    under each other.
  * MEME speed controls: --meme-rates (1 is ~2x faster than the default 2)
    and --mpi N to run MEME under hyphy-mpi, which splits sites across
    processes and is the largest single speed-up available.
  * A .lock file per job, so two terminals cannot start the same run twice.

All extraction, QC, frame-detection, dedup, tree-pruning and parsing logic
from v5 is UNCHANGED.

Inherited from previous versions:

  1. Auto-detects the correct reading frame per gene/dataset.
  2. Detects a UNIVERSALLY-SHARED stop codon position (gene-end signal)
     and trims the extraction window there for ALL sequences uniformly.
  3. Unconditional positional trim of the final codon -> guarantees
     identical sequence length across the alignment.
  4. Per-sequence QC drops genuinely defective sequences.
  5. Deduplicates identical sequences to unique haplotypes before
     writing the final FASTA/tree -- this is what let FUBAR/MEME finish;
     dramatically reduces tree size and optimization time for HyPhy.
  6. Clears MAPLE's internal "_MinorSeqsClade" node labels (invalid
     HyPhy identifiers due to '/' characters).
  7. On HyPhy failure, prints the full errors.log content.
  8. NEW: runs multiple independent HyPhy jobs CONCURRENTLY (via
     --parallel N) instead of strictly one at a time, using idle CPU
     cores since different dataset/gene/method jobs don't depend on
     each other at all.

DATASETS: RSV-A/B x pre/post-vaccine (4 total), "adopter countries" only.

REQUIREMENTS: HyPhy>=2.5 (WSL/Linux/macOS only), biopython, pandas, dendropy
    pip install biopython pandas dendropy

USAGE
    # 0. WHERE AM I?  (safe, reads only, takes seconds)
    python3 selection_pressure_pipeline_6.py --status
    python3 selection_pressure_pipeline_6.py --status --dataset rsvB_prevaccine

    # 1. extraction (only needed if the gene FASTAs are missing)
    python3 selection_pressure_pipeline_6.py --extract --dataset rsvB_prevaccine

    # 2. CONTINUE the interrupted work -- reruns only what is missing or
    #    corrupt, leaves finished jobs alone. Safe to Ctrl-C and relaunch.
    python3 selection_pressure_pipeline_6.py --run --dataset rsvB_prevaccine
    python3 selection_pressure_pipeline_6.py --run --dataset rsvB_postvaccine

    # 3. see what would run, without running it
    python3 selection_pressure_pipeline_6.py --run --dataset rsvB_prevaccine --dry-run

    # 4. MEME only, as fast as possible (needs hyphy-mpi + mpirun)
    python3 selection_pressure_pipeline_6.py --run --method meme --mpi 8

    # 5. parse once --status shows everything DONE
    python3 selection_pressure_pipeline_6.py --parse

DO NOT `rm -rf hyphy_results` any more. v6 resumes correctly, so deleting
finished RSV-A results just makes you run them again for nothing.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------
ROOT = Path("/mnt/c/Users/zahra/Desktop")

GFF_PATHS = {
    "A": ROOT / "adopter countries" / "clean" / "results_rsvA_postvaccine" / "nextclade.gff",
    "B": ROOT / "adopter countries" / "clean" / "results_rsvB_postvaccine" / "nextclade.gff",
}

FILTERED_DIR = ROOT / "adopter countries" / "clean" / "filtered"
MAPLE_DIR = ROOT / "adopter countries" / "maple file"

DATASETS = {}
for subtype in ["A", "B"]:
    for period in ["prevaccine", "postvaccine"]:
        prefix = f"rsv{subtype}_{period}"
        DATASETS[prefix] = {
            "fasta": FILTERED_DIR / f"{prefix}_FILTERED.fasta",
            "tree": MAPLE_DIR / f"{prefix}_MAPLE_nexusTree.tree",
            "subtype": subtype,
        }

GENES = ["F", "G"]

# Where the per-gene FASTAs, pruned trees and HyPhy JSONs live.
# ABSOLUTE on purpose: in v5 this was Path("hyphy_results"), a relative path,
# so results landed in whatever folder you happened to be standing in.
OUTPUT_DIR = ROOT / "RSV seq data" / "hyphy_results"

def autodiscover_output_dir(default: Path) -> Path:
    """If the configured results folder is not there, look for it.

    Handles the folder having been moved or renamed (hyphy_results /
    hyphy_result / HyPhy_results ...) without another round of guessing.
    """
    if default.exists():
        return default
    if not ROOT.exists():
        return default
    hits = []
    for d in ROOT.glob("*/"):
        if not d.is_dir():
            continue
        for sub in d.iterdir() if d.exists() else []:
            try:
                if sub.is_dir() and "hyphy" in sub.name.lower() \
                        and "result" in sub.name.lower():
                    hits.append(sub)
            except OSError:
                continue
    for sub in ROOT.iterdir() if ROOT.exists() else []:
        try:
            if sub.is_dir() and "hyphy" in sub.name.lower() \
                    and "result" in sub.name.lower():
                hits.append(sub)
        except OSError:
            continue
    if hits:
        # prefer the one that actually contains dataset subfolders
        scored = sorted(hits, key=lambda p: -sum(
            1 for lbl in DATASETS if (p / lbl).is_dir()))
        print(f"[info] configured results dir not found; using "
              f"auto-discovered {scored[0]}")
        return scored[0]
    return default


HYPHY_BIN = "hyphy"
HYPHY_MPI_BIN = "hyphy-mpi"

# --- MEME tuning ------------------------------------------------------
# MEME_BRANCHES: v5 ran MEME on internal branches only while FEL and FUBAR
#   used all branches. That is a defensible choice (terminal branches carry
#   most sequencing error), but it means your >=2/3 consensus is combining
#   tests asked of DIFFERENT branch sets. Either keep "Internal" and say so
#   explicitly in the methods section, or set this to "All" for consistency.
MEME_BRANCHES = "Internal"
# MEME_RATES: number of synonymous rate classes. HyPhy default is 2; 1 is
#   roughly twice as fast and is reasonable on shallow within-subtype data.
MEME_RATES = 2
# A job is treated as abandoned if its .lock has not been touched for this
# long, so a crashed run does not block the next one forever.
LOCK_STALE_SECONDS = 3600

FUBAR_POSTERIOR_THRESHOLD = 0.9
MEME_PVALUE_THRESHOLD = 0.1
FEL_PVALUE_THRESHOLD = 0.1
MIN_METHODS_AGREEING = 2

STOP_CODONS = {"TAA", "TAG", "TGA"}
OFFSET_CANDIDATES = [-2, -1, 0, 1, 2]
SAMPLE_SIZE_FOR_DETECTION = 500
SHARED_STOP_MAJORITY_THRESHOLD = 0.50
DEDUPLICATE_IDENTICAL_SEQUENCES = True


# ----------------------------------------------------------------------
# STEP 0 -- GFF coordinates
# ----------------------------------------------------------------------
def parse_gff_coords(gff_path: Path) -> dict:
    coords = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            start, end = int(parts[3]), int(parts[4])
            attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
            name = attrs.get("gene_name") or attrs.get("Name") or attrs.get("ID")
            if name:
                coords[name] = (start, end)
    return coords


_gff_cache = {}
def get_gff_coords(subtype: str) -> dict:
    if subtype not in _gff_cache:
        gff_path = GFF_PATHS[subtype]
        if not gff_path.exists():
            raise FileNotFoundError(f"GFF not found for subtype {subtype}: {gff_path}")
        _gff_cache[subtype] = parse_gff_coords(gff_path)
    return _gff_cache[subtype]


# ----------------------------------------------------------------------
# frame-offset detection
# ----------------------------------------------------------------------
def count_stops_for_offset(records, start, end, offset):
    total_stops, total_codons = 0, 0
    for rec_seq in records:
        a, b = start + offset, end + 1 + offset
        if a < 0 or b > len(rec_seq):
            continue
        sub = rec_seq[a:b].upper()
        sub = sub[: len(sub) - (len(sub) % 3)]
        codons = [sub[i:i + 3] for i in range(0, len(sub), 3)]
        for codon in codons[:-1]:
            if "-" in codon or "N" in codon:
                continue
            total_codons += 1
            if codon in STOP_CODONS:
                total_stops += 1
    return total_stops, total_codons


def detect_frame_offset(sample, start, end, gene):
    best_offset, best_rate = 0, None
    print(f"  [info] auto-detecting frame for {gene}:")
    for offset in OFFSET_CANDIDATES:
        stops, codons = count_stops_for_offset(sample, start, end, offset)
        rate = stops / codons if codons else 1.0
        print(f"    offset {offset:+d}: {stops} stops / {codons} codons (rate {rate:.4f})")
        if best_rate is None or rate < best_rate:
            best_rate, best_offset = rate, offset
    print(f"  [info] selected offset {best_offset:+d} for {gene} (stop rate {best_rate:.4f})")
    return best_offset


# ----------------------------------------------------------------------
# universally-shared stop codon detection (gene-end signal)
# ----------------------------------------------------------------------
def detect_shared_stop_codon_index(sample, start, end, offset):
    a, b = start + offset, end + 1 + offset
    position_hits = Counter()
    position_total = Counter()

    for rec_seq in sample:
        if a < 0 or b > len(rec_seq):
            continue
        sub = rec_seq[a:b].upper()
        sub = sub[: len(sub) - (len(sub) % 3)]
        codons = [sub[i:i + 3] for i in range(0, len(sub), 3)]
        for idx, codon in enumerate(codons):
            if "-" in codon or "N" in codon:
                continue
            position_total[idx] += 1
            if codon in STOP_CODONS:
                position_hits[idx] += 1

    n_codons = (b - a) // 3
    for idx in range(n_codons):
        total = position_total.get(idx, 0)
        hits = position_hits.get(idx, 0)
        if total == 0:
            continue
        frac = hits / total
        if frac >= SHARED_STOP_MAJORITY_THRESHOLD:
            print(f"  [info] found universally shared stop codon at codon index {idx} "
                  f"({hits}/{total} = {frac:.1%} of sequences) -- treating as true CDS end")
            return idx
    return None


# ----------------------------------------------------------------------
# per-sequence QC -- unconditional positional trim of last codon
# ----------------------------------------------------------------------
def clean_and_check_stop(nt_seq: str):
    seq = nt_seq.upper()
    seq = seq[: len(seq) - (len(seq) % 3)]
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]

    if codons:
        codons = codons[:-1]

    for codon in codons:
        if "-" in codon or "N" in codon:
            continue
        if codon in STOP_CODONS:
            return "".join(codons), True

    return "".join(codons), False


# ----------------------------------------------------------------------
# STEP 1 -- extract per-gene nucleotide fasta (frame + universal-trim + QC + dedup)
# ----------------------------------------------------------------------
def extract_gene_nt_fasta(aligned_fasta: Path, coords: dict, gene: str, out_path: Path,
                           min_nongap_frac: float = 0.5):
    if gene not in coords:
        print(f"  [WARN] gene '{gene}' not in GFF: available = {list(coords.keys())}")
        return None

    start, end = coords[gene]

    sample = []
    for i, rec in enumerate(SeqIO.parse(str(aligned_fasta), "fasta")):
        if i >= SAMPLE_SIZE_FOR_DETECTION:
            break
        sample.append(str(rec.seq))

    offset = detect_frame_offset(sample, start, end, gene)
    shared_stop_idx = detect_shared_stop_codon_index(sample, start, end, offset)

    a = start + offset
    if shared_stop_idx is not None:
        b = a + (shared_stop_idx + 1) * 3
    else:
        b = end + 1 + offset
        print(f"  [info] no universally shared stop found for {gene}; using full annotated window")

    expected_len = b - a
    expected_len -= expected_len % 3

    kept = []
    dropped_stop = 0
    dropped_gap = 0
    dropped_length = 0

    for rec in SeqIO.parse(str(aligned_fasta), "fasta"):
        full = str(rec.seq)
        if a < 0 or b > len(full):
            dropped_gap += 1
            continue
        sub = full[a:b]
        if len(sub) != expected_len:
            dropped_length += 1
            continue
        nongap_frac = sum(1 for c in sub if c not in "-N") / len(sub)
        if nongap_frac < min_nongap_frac:
            dropped_gap += 1
            continue

        cleaned_seq, bad_stop = clean_and_check_stop(sub)
        if bad_stop:
            dropped_stop += 1
            continue

        kept.append((rec.id, cleaned_seq))

    if not kept:
        print(f"  [WARN] no sequences passed QC for gene {gene}")
        return None

    lengths = {len(s) for _, s in kept}
    if len(lengths) > 1:
        print(f"  [ERROR] inconsistent lengths still present after fix: {lengths} -- aborting {gene}")
        return None
    seq_len = lengths.pop()

    n_before_dedup = len(kept)
    if DEDUPLICATE_IDENTICAL_SEQUENCES:
        seen = {}
        representatives = []
        for seq_id, seq in kept:
            if seq not in seen:
                seen[seq] = seq_id
                representatives.append((seq_id, seq))
        kept = representatives
        n_collapsed = n_before_dedup - len(kept)
    else:
        n_collapsed = 0

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for seq_id, seq in kept:
            f.write(f">{seq_id}\n{seq}\n")
    print(f"  [ok] {gene}: {len(kept)} unique seqs -> {out_path.name} (window {a}-{b}, offset {offset:+d}, "
          f"length {seq_len} nt); dropped {dropped_gap} (low coverage), "
          f"{dropped_length} (length mismatch), {dropped_stop} (internal stop codon), "
          f"{n_collapsed} (duplicate, kept 1 representative each)")
    return [k[0] for k in kept]


# ----------------------------------------------------------------------
# STEP 2 -- prune MAPLE tree via dendropy (with internal-label fix)
# ----------------------------------------------------------------------
def load_tree_any_format(tree_path: Path):
    import dendropy
    last_err = None
    for schema in ("nexus", "newick"):
        try:
            tree = dendropy.Tree.get(path=str(tree_path), schema=schema, preserve_underscores=True)
            print(f"  [info] loaded {tree_path.name} as {schema.upper()} format")
            return tree
        except Exception as e:
            last_err = e
            continue
    raise ValueError(f"Could not parse {tree_path} as Nexus or Newick: {last_err}")


def prune_tree(tree_path: Path, keep_ids: list, out_tree_path: Path):
    try:
        import dendropy  # noqa: F401
    except ImportError:
        print("  [ERROR] pip install dendropy")
        return None

    if not tree_path.exists():
        print(f"  [WARN] tree not found: {tree_path}")
        return None

    try:
        tree = load_tree_any_format(tree_path)
    except ValueError as e:
        print(f"  [ERROR] {e}")
        return None

    taxon_labels = {t.label for t in tree.taxon_namespace}
    keep_set = set(keep_ids) & taxon_labels
    missing = set(keep_ids) - taxon_labels
    if missing:
        print(f"  [WARN] {len(missing)} gene-fasta IDs not found in tree e.g. {list(missing)[:3]}")

    if len(keep_set) < 4:
        print(f"  [ERROR] only {len(keep_set)} matching taxa -- aborting")
        return None

    tree.retain_taxa_with_labels(list(keep_set))

    n_cleared = 0
    for node in tree.preorder_node_iter():
        if not node.is_leaf() and node.label is not None:
            node.label = None
            n_cleared += 1
    if n_cleared:
        print(f"  [info] cleared {n_cleared} internal node labels (MAPLE clade names, not needed by HyPhy)")

    out_tree_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(path=str(out_tree_path), schema="newick", suppress_rooting=True,
               unquoted_underscores=True, suppress_internal_node_labels=True)
    print(f"  [ok] pruned tree -> {out_tree_path.name} ({len(keep_set)} taxa)")
    return keep_set


def filter_fasta_to_ids(fasta_path: Path, keep_ids: set, out_path: Path) -> int:
    kept = [r for r in SeqIO.parse(str(fasta_path), "fasta") if r.id in keep_ids]
    with open(out_path, "w") as f:
        for r in kept:
            f.write(f">{r.id}\n{str(r.seq)}\n")
    return len(kept)


def extract_all(dataset=None, gene=None):
    labels = [dataset] if dataset else list(DATASETS.keys())
    genes = [gene] if gene else list(GENES)
    for label in labels:
        d = DATASETS[label]
        print(f"\n=== {label} ===")
        if not d["fasta"].exists():
            print(f"  [WARN] missing whole-genome FILTERED fasta: {d['fasta']}")
            continue

        coords = get_gff_coords(d["subtype"])
        for gene in genes:
            raw_fasta = OUTPUT_DIR / label / f"{gene}.raw.fasta"
            seq_ids = extract_gene_nt_fasta(d["fasta"], coords, gene, raw_fasta)
            if not seq_ids:
                continue

            pruned_tree = OUTPUT_DIR / label / f"{gene}.tree.nwk"
            kept_ids = prune_tree(d["tree"], seq_ids, pruned_tree)
            if not kept_ids:
                continue

            final_fasta = OUTPUT_DIR / label / f"{gene}.fasta"
            n = filter_fasta_to_ids(raw_fasta, kept_ids, final_fasta)
            print(f"  [ok] final {gene}: {n} sequences ready for HyPhy")


# ----------------------------------------------------------------------
# STEP 3 -- run HyPhy (with full errors.log reporting on failure)
# ----------------------------------------------------------------------
def count_codons(aln: Path) -> int:
    """Alignment length / 3. All sequences are guaranteed equal length by
    the extraction step, so reading the first record is enough."""
    if not aln.exists():
        return 0
    length, started = 0, False
    with open(aln) as f:
        for line in f:
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif started:
                length += len(line.strip())
    return length // 3


def inspect_result(out_json: Path, n_codons: int):
    """Classify an existing HyPhy .json.

    Returns (status, n_rows, note) where status is one of
    DONE / MISSING / CORRUPT / PARTIAL.
    This is the fix for BUG 2: v5 only asked 'does the file exist'.
    """
    if not out_json.exists():
        return ("MISSING", 0, "")
    if out_json.stat().st_size == 0:
        return ("CORRUPT", 0, "empty file")
    try:
        with open(out_json, errors="replace") as f:
            d = json.load(f)
    except json.JSONDecodeError as e:
        return ("CORRUPT", 0, f"truncated JSON at line {e.lineno}")
    except OSError as e:
        return ("CORRUPT", 0, str(e))

    try:
        rows = len(d["MLE"]["content"]["0"])
    except (KeyError, TypeError):
        return ("CORRUPT", 0, "no MLE/content/0 table")
    if rows == 0:
        return ("CORRUPT", 0, "empty MLE table")
    if n_codons and rows < n_codons:
        return ("PARTIAL", rows, f"{rows}/{n_codons} sites")
    return ("DONE", rows, "")


def lock_path_for(out_json: Path) -> Path:
    return out_json.with_suffix(out_json.suffix + ".lock")


def lock_is_live(lock: Path) -> bool:
    if not lock.exists():
        return False
    if time.time() - lock.stat().st_mtime > LOCK_STALE_SECONDS:
        return False
    try:
        pid = int(lock.read_text().split()[0])
    except (ValueError, IndexError, OSError):
        return True
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def build_hyphy_cmd(method, aln, tree, out_tmp, extra_args, cpu, mpi):
    if mpi and mpi > 1 and method == "meme":
        base = ["mpirun", "-np", str(mpi), HYPHY_MPI_BIN]
    else:
        base = [HYPHY_BIN]
        if cpu:
            # BUG 3 fix: cap this job's OpenMP threads so concurrent jobs
            # do not oversubscribe the machine.
            base.append(f"CPU={cpu}")
    cmd = base + [method,
                  "--alignment", str(aln),
                  "--tree", str(tree),
                  "--output", str(out_tmp)]
    return cmd + list(extra_args or [])


def run_hyphy_method(method, aln, tree, out_json, extra_args=None, force=False,
                     cpu=0, mpi=0, dry_run=False, log_prefix=""):
    n_codons = count_codons(aln)
    status, rows, note = inspect_result(out_json, n_codons)

    if status == "DONE" and not force:
        print(f"  [skip] {out_json.name} -- already complete ({rows} sites)")
        return True
    if status in ("CORRUPT", "PARTIAL"):
        # This is the case v5 silently accepted forever.
        print(f"  [redo] {out_json.name} is {status} ({note}) -- rerunning")
        try:
            out_json.unlink()
        except OSError:
            pass

    if not aln.exists() or not tree.exists():
        print(f"  [WARN] missing input for {method}: {aln} / {tree}")
        return False

    lock = lock_path_for(out_json)
    if lock_is_live(lock):
        print(f"  [busy] {out_json.name} is already running in another "
              f"process -- skipping")
        return False

    out_tmp = out_json.with_suffix(out_json.suffix + ".tmp")
    logf = out_json.with_suffix(".log")
    cmd = build_hyphy_cmd(method, aln, tree, out_tmp, extra_args, cpu, mpi)

    if dry_run:
        print(f"  [dry ] {' '.join(cmd)}")
        return True

    print(f"  [run ] {' '.join(cmd)}")
    lock.write_text(f"{os.getpid()} {datetime.now().isoformat()}\n")
    t0 = time.time()
    ok = False
    try:
        with open(logf, "a") as lf:
            lf.write(f"\n===== {datetime.now().isoformat()} =====\n")
            lf.write(" ".join(cmd) + "\n")
            lf.flush()
            proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)
            wait = 1.0
            while proc.poll() is None:
                time.sleep(wait)
                wait = min(wait * 1.5, 30.0)
                try:
                    os.utime(lock, None)      # keep the lock fresh
                except OSError:
                    pass
            rc = proc.returncode

        elapsed = time.time() - t0
        if rc != 0:
            tail = ""
            if logf.exists():
                tail = logf.read_text(errors="replace")[-2000:]
            print(f"  [FAIL] {log_prefix}{method} rc={rc} after "
                  f"{elapsed/60:.1f} min -- last of {logf.name}:\n{tail}")
        else:
            st, rws, nt = inspect_result(out_tmp, n_codons)
            if st == "DONE":
                # atomic: the final .json only ever appears complete
                os.replace(out_tmp, out_json)
                ok = True
                print(f"  [ok  ] {out_json.name} -- {rws} sites in "
                      f"{elapsed/60:.1f} min")
            else:
                print(f"  [FAIL] {log_prefix}{method} produced {st} output "
                      f"({nt}); left as {out_tmp.name}")
    except FileNotFoundError:
        print(f"  [FAIL] cannot find '{cmd[0]}'. Is HyPhy on your PATH? "
              f"Try `which hyphy`, or activate your conda env.")
    except KeyboardInterrupt:
        print("  [stop] interrupted -- partial output kept as .tmp, "
              "rerun --run to redo this job cleanly")
        raise
    finally:
        if lock.exists():
            try:
                lock.unlink()
            except OSError:
                pass

    with open(OUTPUT_DIR / "progress.tsv", "a") as f:
        f.write(f"{datetime.now().isoformat()}\t{out_json.parent.name}\t"
                f"{out_json.stem}\t{method}\t{'ok' if ok else 'fail'}\t"
                f"{int(time.time() - t0)}\n")
    return ok


def collect_jobs(dataset=None, gene=None, method=None, meme_rates=None,
                 meme_branches=None):
    """Build the full job list across every dataset/gene/method.

    v5 processed one dataset at a time, so at the end of each dataset the
    machine sat with 2 idle workers waiting for a single slow MEME. One
    global queue keeps every core busy until the very last job.
    """
    labels = [dataset] if dataset else list(DATASETS.keys())
    genes = [gene] if gene else list(GENES)
    if method:
        methods = [method] if isinstance(method, str) else list(method)
    else:
        methods = ["meme", "fel", "fubar"]

    rates = MEME_RATES if meme_rates is None else meme_rates
    branches = meme_branches or MEME_BRANCHES

    jobs = []
    for label in labels:
        gd = OUTPUT_DIR / label
        for g in genes:
            aln, tree = gd / f"{g}.fasta", gd / f"{g}.tree.nwk"
            spec = {
                "fubar": (gd / f"{g}.FUBAR.json", None),
                "fel":   (gd / f"{g}.FEL.json", None),
                "meme":  (gd / f"{g}.MEME.json",
                          ["--branches", branches, "--rates", str(rates)]),
            }
            for m in methods:
                out_json, extra = spec[m]
                jobs.append({"method": m, "aln": aln, "tree": tree,
                             "out": out_json, "extra": extra,
                             "label": label, "gene": g})
    return jobs


def run_all(force=False, parallel=None, dataset=None, gene=None, method=None,
            mpi=0, meme_rates=None, meme_branches=None, dry_run=False):
    if dataset and dataset not in DATASETS:
        print(f"[ERROR] unknown dataset: {dataset}. "
              f"Valid options: {list(DATASETS.keys())}")
        return

    n_cpu = os.cpu_count() or 2
    if mpi and mpi > 1:
        parallel = 1                      # MPI already uses the whole box
        print(f"[info] MEME will run under mpirun -np {mpi}; "
              f"one job at a time")
    elif parallel is None:
        parallel = min(3, max(1, n_cpu - 1))

    # BUG 3 fix: divide the cores between workers instead of letting every
    # HYPHYMP job try to grab all of them.
    cpu_per_job = max(1, n_cpu // parallel) if parallel > 1 else 0
    if cpu_per_job:
        print(f"[info] {n_cpu} CPUs detected -> {parallel} concurrent jobs "
              f"x CPU={cpu_per_job} each")

    jobs = collect_jobs(dataset, gene, method, meme_rates, meme_branches)

    todo, done, blocked = [], [], []
    for j in jobs:
        st, rows, note = inspect_result(j["out"], count_codons(j["aln"]))
        if st == "DONE" and not force:
            done.append(j)
        elif lock_is_live(lock_path_for(j["out"])):
            blocked.append(j)
        elif not j["aln"].exists() or not j["tree"].exists():
            blocked.append(j)
            print(f"[skip] {j['label']}/{j['gene']}/{j['method']}: inputs "
                  f"missing -- run --extract --dataset {j['label']} first")
        else:
            j["_status"] = st
            j["_note"] = note
            todo.append(j)

    print(f"[info] {len(done)} already complete, {len(todo)} to run, "
          f"{len(blocked)} skipped")
    if not todo:
        print("[info] nothing to do. Run --status to see the full picture.")
        return

    if not dry_run:
        need = HYPHY_MPI_BIN if (mpi and mpi > 1) else HYPHY_BIN
        if shutil.which(need) is None and not Path(need).exists():
            print(f"[ERROR] '{need}' not found on PATH. Activate your conda "
                  f"env, or edit HYPHY_BIN at the top of this file.")
            return
        if mpi and mpi > 1 and shutil.which("mpirun") is None:
            print("[ERROR] --mpi requested but mpirun is not installed. "
                  "Try: conda install -c bioconda hyphy  (MPI build), "
                  "or drop --mpi.")
            return

    # longest first: MEME dominates the wall clock, so start those before
    # the cheap FEL/FUBAR jobs rather than after.
    rank = {"meme": 0, "fel": 1, "fubar": 2}
    todo.sort(key=lambda j: (rank.get(j["method"], 9), j["label"], j["gene"]))

    for j in todo:
        tag = "rerun " + j["_status"] if j["_status"] != "MISSING" else "new"
        print(f"       - {j['label']}/{j['gene']}/{j['method'].upper()} ({tag})")

    t_start = time.time()
    results = []

    def _run(j):
        print()
        print(f"=== {j['label']} / {j['gene']} / {j['method']} (starting) ===")
        ok = run_hyphy_method(j["method"], j["aln"], j["tree"], j["out"],
                              extra_args=j["extra"], force=force,
                              cpu=cpu_per_job, mpi=mpi, dry_run=dry_run,
                              log_prefix=f"{j['label']}/{j['gene']} ")
        print(f"=== {j['label']} / {j['gene']} / {j['method']} -> "
              f"{'OK' if ok else 'FAILED'} ===")
        return j, ok

    if parallel <= 1:
        for j in todo:
            results.append(_run(j))
    else:
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = [executor.submit(_run, j) for j in todo]
            for future in as_completed(futures):
                results.append(future.result())

    n_ok = sum(1 for _, ok in results if ok)
    print()
    print(f"########## {n_ok}/{len(results)} jobs succeeded in "
          f"{(time.time() - t_start)/60:.1f} min ##########")
    failed = [j for j, ok in results if not ok]
    if failed:
        print("Failed jobs (rerun --run to retry; logs are next to the JSONs):")
        for j in failed:
            print(f"  - {j['label']}/{j['gene']}/{j['method']}  "
                  f"see {j['out'].with_suffix('.log')}")
    print("Run --status for the current picture.")


def _show_dir(path: Path, label: str, limit=40, pattern="*"):
    print(f"\n--- {label}")
    print(f"    {path}")
    if not path.exists():
        print("    [MISSING] this directory does not exist")
        return []
    try:
        items = sorted(path.glob(pattern))
    except OSError as e:
        print(f"    [ERROR] {e}")
        return []
    if not items:
        print("    (empty)")
        return []
    for p in items[:limit]:
        kind = "dir " if p.is_dir() else "file"
        size = "" if p.is_dir() else f"  {p.stat().st_size/1e6:>8.2f} MB"
        print(f"    [{kind}] {p.name}{size}")
    if len(items) > limit:
        print(f"    ... and {len(items) - limit} more")
    return items


def doctor():
    """Print exactly what is on disk, so paths can be fixed from facts
    rather than from guesses about how the folders are named."""
    print("=" * 66)
    print("DOCTOR -- checking every path this pipeline depends on")
    print("=" * 66)

    print(f"\nROOT        {ROOT}   {'OK' if ROOT.exists() else '[MISSING]'}")
    print(f"OUTPUT_DIR  {OUTPUT_DIR}   "
          f"{'OK' if OUTPUT_DIR.exists() else '[MISSING]'}")

    if ROOT.exists():
        print("\n--- folders directly under ROOT (Desktop)")
        for p in sorted(ROOT.iterdir()):
            try:
                if p.is_dir():
                    print(f"    [dir ] {p.name}")
            except OSError:
                pass

    _show_dir(FILTERED_DIR, "FILTERED_DIR (whole-genome input FASTAs)",
              pattern="*.fasta")
    _show_dir(MAPLE_DIR, "MAPLE_DIR (input trees)")

    print("\n--- expected input files, one row per dataset")
    for label, d in DATASETS.items():
        f_ok = "OK" if d["fasta"].exists() else "MISSING"
        t_ok = "OK" if d["tree"].exists() else "MISSING"
        print(f"    {label:<20} fasta={f_ok:<8} tree={t_ok}")
        if f_ok == "MISSING":
            print(f"        looked for: {d['fasta']}")
        if t_ok == "MISSING":
            print(f"        looked for: {d['tree']}")

    _show_dir(OUTPUT_DIR, "OUTPUT_DIR (results root)")
    for label in DATASETS:
        _show_dir(OUTPUT_DIR / label, f"OUTPUT_DIR / {label}")

    print("\n--- what the pipeline expects inside each dataset folder")
    print("    F.fasta  F.tree.nwk  F.FUBAR.json  F.MEME.json  F.FEL.json")
    print("    G.fasta  G.tree.nwk  G.FUBAR.json  G.MEME.json  G.FEL.json")
    print("\n    If the files above exist under different names, tell me the")
    print("    real names and I will adjust the script rather than you")
    print("    renaming hundreds of files.")

    print("\n--- tools")
    for tool in ("hyphy", "hyphy-mpi", "mpirun"):
        w = shutil.which(tool)
        print(f"    {tool:<10} {w or '[not on PATH]'}")
    print(f"    CPUs       {os.cpu_count()}")
    print("\n" + "=" * 66)


# ----------------------------------------------------------------------
# STEP 3b -- STATUS: how far have I actually got?
# ----------------------------------------------------------------------
def status_all(dataset=None, gene=None):
    labels = [dataset] if dataset else list(DATASETS.keys())
    if dataset and dataset not in DATASETS:
        print(f"[ERROR] unknown dataset: {dataset}. "
              f"Valid options: {list(DATASETS.keys())}")
        return
    genes = [gene] if gene else list(GENES)

    print()
    print("=== EXTRACTION (step 1) ===")
    print(f"{'dataset':<22}{'gene':<6}{'seqs':>8}{'codons':>9}  tree")
    ready = {}
    for label in labels:
        gd = OUTPUT_DIR / label
        for g in genes:
            aln, tree = gd / f"{g}.fasta", gd / f"{g}.tree.nwk"
            if aln.exists():
                n_seq = sum(1 for line in open(aln) if line.startswith(">"))
                n_cod = count_codons(aln)
            else:
                n_seq, n_cod = 0, 0
            ok_tree = "yes" if tree.exists() else "MISSING"
            shown = n_seq if aln.exists() else "-"
            print(f"{label:<22}{g:<6}{str(shown):>8}{str(n_cod or '-'):>9}"
                  f"  {ok_tree}")
            ready[(label, g)] = (aln.exists() and tree.exists(), n_cod)

    print()
    print("=== HyPhy JOBS (step 2) ===")
    print(f"{'dataset':<22}{'gene':<6}{'method':<7}{'status':<9}{'sites':>7}  note")
    counts = Counter()
    todo = []
    for label in labels:
        gd = OUTPUT_DIR / label
        for g in genes:
            inputs_ok, n_cod = ready[(label, g)]
            for m, fn in (("FUBAR", "FUBAR"), ("MEME", "MEME"), ("FEL", "FEL")):
                out_json = gd / f"{g}.{fn}.json"
                if not inputs_ok:
                    st, rows, note = "NO-INPUT", 0, "run --extract first"
                else:
                    st, rows, note = inspect_result(out_json, n_cod)
                    if st != "DONE" and lock_is_live(lock_path_for(out_json)):
                        st, note = "RUNNING", "another process has it"
                    elif st == "MISSING":
                        tmp = out_json.with_suffix(out_json.suffix + ".tmp")
                        if tmp.exists():
                            note = (f"interrupted run left "
                                    f"{tmp.stat().st_size/1e6:.1f} MB .tmp")
                counts[st] += 1
                if st in ("MISSING", "CORRUPT", "PARTIAL"):
                    todo.append((label, g, m))
                print(f"{label:<22}{g:<6}{m:<7}{st:<9}"
                      f"{str(rows or ''):>7}  {note}")

    total = sum(counts.values())
    done = counts["DONE"]
    print()
    print("=== SUMMARY ===")
    for k in ("DONE", "RUNNING", "PARTIAL", "CORRUPT", "MISSING", "NO-INPUT"):
        if counts[k]:
            print(f"  {k:<9} {counts[k]}")
    bar = 34
    filled = int(bar * done / total) if total else 0
    pct = 100.0 * done / total if total else 0.0
    print()
    print(f"  [{'#' * filled}{'.' * (bar - filled)}] {done}/{total} jobs "
          f"({pct:.0f}%)")

    if todo:
        print()
        print("  still to run:")
        for label, g, m in todo:
            print(f"    - {label} / {g} / {m}")
        ds = f" --dataset {dataset}" if dataset else ""
        print()
        print(f"  next:  python3 {Path(sys.argv[0]).name} --run{ds}")
    else:
        print()
        print(f"  Nothing left to run. Next: python3 "
              f"{Path(sys.argv[0]).name} --parse")

    # timing hints from the ledger, so you can predict the remaining wait
    ledger = OUTPUT_DIR / "progress.tsv"
    if ledger.exists():
        times = {}
        for line in ledger.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 6 and parts[4] == "ok":
                times.setdefault(parts[3], []).append(int(parts[5]))
        if times:
            print()
            print("  observed runtimes so far (median):")
            for m, vals in sorted(times.items()):
                vals.sort()
                med = vals[len(vals) // 2]
                print(f"    {m.upper():<7} {med/60:.1f} min  "
                      f"(n={len(vals)})")

    if MEME_BRANCHES != "All":
        print()
        print(f"  NOTE: MEME is configured with --branches {MEME_BRANCHES} "
              f"while FEL and FUBAR use all branches.")
        print("        Keep it if intended, but state it in your methods "
              "section, since the >=2/3")
        print("        consensus then combines tests over different branch "
              "sets.")


# ----------------------------------------------------------------------
# STEP 4 -- parse + >=2/3 consensus
# ----------------------------------------------------------------------
def _mle_df(path):
    """Read a HyPhy MLE table into a DataFrame.

    HyPhy does not always list every data column in the 'headers' block.
    FUBAR in particular appends MCMC convergence diagnostics (PSRF, N_eff)
    to each row without declaring them, so the row width can exceed the
    header count and pandas raises
        'N columns passed, passed data had M columns'.
    Reconcile the two rather than assuming they match.
    """
    d = json.load(open(path))
    rows = d["MLE"]["content"]["0"]
    h = [x[0] for x in d["MLE"]["headers"]]
    ncol = max(len(r) for r in rows) if rows else len(h)

    if len(h) < ncol:
        # name the known FUBAR diagnostics, then fall back to generic names
        extras = ["PSRF", "N_eff"]
        add = []
        for i in range(ncol - len(h)):
            add.append(extras[i] if i < len(extras) else f"extra_{i}")
        h = h + add
    elif len(h) > ncol:
        h = h[:ncol]

    rows = [list(r) + [None] * (ncol - len(r)) for r in rows]
    df = pd.DataFrame(rows, columns=h)
    df.insert(0, "site", range(1, len(df) + 1))
    return df


def _col(df, *candidates, contains=None):
    """Find a column by exact name, then by substring, so small changes in
    HyPhy's column labelling between versions do not break parsing."""
    for c in candidates:
        if c in df.columns:
            return c
    if contains:
        for c in df.columns:
            if contains.lower() in str(c).lower():
                return c
    raise KeyError(f"none of {candidates} (or '*{contains}*') in {list(df.columns)}")


def parse_fubar(p):
    df = _mle_df(p)
    pcol = _col(df, "Prob[alpha<beta]", contains="alpha<beta")
    a = _col(df, "alpha", contains="alpha")
    b = _col(df, "beta", contains="beta")
    df["FUBAR_flag"] = df[pcol] > FUBAR_POSTERIOR_THRESHOLD
    out = df[["site", a, b, pcol, "FUBAR_flag"]]
    return out.rename(columns={a: "alpha", b: "beta",
                               pcol: "Prob[alpha<beta]"})


def parse_meme(p):
    df = _mle_df(p)
    pcol = _col(df, contains="p-value")
    df["MEME_flag"] = df[pcol] < MEME_PVALUE_THRESHOLD
    return df[["site", pcol, "MEME_flag"]].rename(columns={pcol: "MEME_pvalue"})


def parse_fel(p):
    df = _mle_df(p)
    pcol = _col(df, contains="p-value")
    a = _col(df, "alpha", contains="alpha")
    b = _col(df, "beta", contains="beta")
    df["FEL_flag"] = (df[pcol] < FEL_PVALUE_THRESHOLD) & (df[b] > df[a])
    out = df[["site", a, b, pcol, "FEL_flag"]]
    return out.rename(columns={a: "alpha", b: "beta", pcol: "FEL_pvalue"})


def combine_one(label, gene):
    gd = OUTPUT_DIR / label
    paths = [gd / f"{gene}.FUBAR.json", gd / f"{gene}.MEME.json", gd / f"{gene}.FEL.json"]
    # v6: validate before parsing. v5 would happily read a truncated file
    # and either crash with a JSONDecodeError or silently under-report.
    n_cod = count_codons(gd / f"{gene}.fasta")
    bad = []
    for p in paths:
        st, rows, note = inspect_result(p, n_cod)
        if st != "DONE":
            bad.append(f"{p.name} [{st}{': ' + note if note else ''}]")
    if bad:
        print(f"  [WARN] {label}/{gene} not ready: {', '.join(bad)}")
        print(f"         fix with: --run --dataset {label} --gene {gene}")
        return pd.DataFrame()
    f, m, e = parse_fubar(paths[0]), parse_meme(paths[1]), parse_fel(paths[2])
    merged = f.merge(m, on="site").merge(e, on="site", suffixes=("_fubar", "_fel"))
    merged["n_methods_flagging"] = (merged["FUBAR_flag"].astype(int) +
                                     merged["MEME_flag"].astype(int) +
                                     merged["FEL_flag"].astype(int))
    merged["positively_selected"] = merged["n_methods_flagging"] >= MIN_METHODS_AGREEING
    merged.insert(0, "gene", gene)
    merged.insert(0, "dataset", label)
    return merged


def parse_all():
    dfs = [combine_one(l, g) for l in DATASETS for g in GENES]
    dfs = [d for d in dfs if not d.empty]
    if not dfs:
        print("Nothing parsed.")
        return
    full = pd.concat(dfs, ignore_index=True)
    full.to_csv(OUTPUT_DIR / "all_sites_selection_results.csv", index=False)
    summary = (full[full["positively_selected"]].groupby(["dataset", "gene"]).size()
               .reset_index(name="n_positively_selected_sites"))
    summary.to_csv(OUTPUT_DIR / "summary_positive_sites_by_dataset_gene.csv", index=False)
    detail = full[full["positively_selected"]][
        ["dataset", "gene", "site", "n_methods_flagging", "Prob[alpha<beta]", "MEME_pvalue", "FEL_pvalue"]
    ].sort_values(["dataset", "gene", "site"])
    detail.to_csv(OUTPUT_DIR / "positively_selected_sites_detail.csv", index=False)
    print(f"\n{len(full)} site-records; {full['positively_selected'].sum()} flagged (>=2/3 methods).")
    print(f"Outputs -> {OUTPUT_DIR.resolve()}")


# ----------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="RSV selection-pressure pipeline (v6, resumable). "
                    "Start with --status.")
    ap.add_argument("--status", action="store_true",
                    help="report how far each dataset/gene/method has got")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--parse", action="store_true")
    ap.add_argument("--force", action="store_true",
                    help="redo jobs even if their JSON is already valid")
    ap.add_argument("--dataset", default=None,
                    help=f"restrict to one dataset: {list(DATASETS.keys())}")
    ap.add_argument("--gene", default=None, choices=GENES,
                    help="restrict to one gene")
    ap.add_argument("--method", default=None, nargs="+",
                    choices=["fel", "fubar", "meme"],
                    help="restrict to one or more methods, e.g. "
                         "--method fel fubar")
    ap.add_argument("--parallel", type=int, default=None,
                    help="concurrent HyPhy jobs (default: min(3, CPUs-1)); "
                         "cores are split between them")
    ap.add_argument("--mpi", type=int, default=0,
                    help="run MEME under 'mpirun -np N hyphy-mpi'; biggest "
                         "single speed-up if you have the MPI build")
    ap.add_argument("--meme-rates", type=int, default=None,
                    help=f"MEME synonymous rate classes (default "
                         f"{MEME_RATES}); 1 is about twice as fast")
    ap.add_argument("--meme-branches", default=None,
                    choices=["All", "Internal", "Leaves"],
                    help=f"branch set for MEME (default {MEME_BRANCHES})")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the HyPhy commands without running them")
    ap.add_argument("--doctor", action="store_true",
                    help="print every path this pipeline uses and what is "
                         "actually on disk -- run this when --status says "
                         "MISSING but you know the files exist")
    ap.add_argument("--output-dir", default=None,
                    help=f"where the per-gene FASTAs and HyPhy JSONs live "
                         f"(default: {OUTPUT_DIR})")
    ap.add_argument("--clear-locks", action="store_true",
                    help="delete leftover .lock files (use only if you are "
                         "sure no HyPhy job is actually running)")
    args = ap.parse_args()

    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir).expanduser()
    else:
        OUTPUT_DIR = autodiscover_output_dir(OUTPUT_DIR)
    print(f"[info] results directory: {OUTPUT_DIR}")
    if not OUTPUT_DIR.exists():
        print("[warn] that directory does not exist. Run --doctor to see "
              "what is actually on disk.")

    if args.clear_locks:
        n = 0
        for lk in OUTPUT_DIR.rglob("*.lock"):
            lk.unlink()
            n += 1
        print(f"[info] removed {n} lock file(s)")

    # Only create it when we are going to write. Creating it during --status
    # or --doctor would mask a wrong path with a new empty folder.
    if args.extract or args.run:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if not (args.status or args.extract or args.run or args.parse
            or args.doctor):
        print("Use --doctor / --status / --extract / --run / --parse. "
              "See docstring at top of file.")
        sys.exit(0)
    if args.doctor:
        doctor()
    if args.status:
        status_all(dataset=args.dataset, gene=args.gene)
    if args.extract:
        if args.dataset and args.dataset not in DATASETS:
            print(f"[ERROR] unknown dataset: {args.dataset}. "
                  f"Valid: {list(DATASETS.keys())}")
            sys.exit(1)
        extract_all(dataset=args.dataset, gene=args.gene)
    if args.run:
        run_all(force=args.force, parallel=args.parallel,
                dataset=args.dataset, gene=args.gene, method=args.method,
                mpi=args.mpi, meme_rates=args.meme_rates,
                meme_branches=args.meme_branches, dry_run=args.dry_run)
    if args.parse:
        parse_all()
