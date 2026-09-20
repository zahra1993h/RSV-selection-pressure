#!/usr/bin/env python3
"""
map_G_positions.py -- map your HyPhy G-gene alignment columns onto the
reference numbering used everywhere else in the manuscript.

RUN (plain python, from the Desktop, inside WSL):

    python3 map_G_positions.py

WHY THIS IS NOT JUST AN OFFSET
    F needed a single offset (+2 for RSV-A, +6 for RSV-B) because F is
    the same length in every strain. G is not. Contemporary RSV-A is
    genotype ON1, carrying a 72-nucleotide (24 aa) duplication in the
    C-terminal hypervariable region; contemporary RSV-B is genotype BA,
    carrying a 60-nucleotide (20 aa) duplication. A duplication is an
    INSERTION, so residues before it and after it need different
    corrections. No single number can be right for the whole gene.

    This script therefore does a proper global alignment (Needleman-
    Wunsch) between your consensus and the reference consensus, and maps
    every selected site individually through that alignment.

WHICH REFERENCE
    It prefers your own Nextclade CDS translation, because that is what
    Tables 11-16 are numbered against. Keeping Tables 17/18 on the same
    numbering is what allows the residue-209 convergence to be stated.

    Put the Nextclade G translation next to this script as:
        nextclade.cds_translation.G.fasta            (used for both), or
        nextclade.cds_translation.G.rsvA.fasta  and
        nextclade.cds_translation.G.rsvB.fasta       (subtype-specific)

    If neither is present the script says so and stops, rather than
    silently falling back to a different reference and producing numbers
    that disagree with your other tables.
"""

import os
from collections import Counter

BASE = "RSV seq data/hyphy_results"
DATASETS = {"A": ("rsvA_prevaccine", "rsvA_postvaccine"),
            "B": ("rsvB_prevaccine", "rsvB_postvaccine")}

# G-gene sites from all_sites_selection_results.csv (alignment columns)
G_SITES = {
    "A": {"shared": [69, 212, 245, 253, 265, 271, 272, 282, 295, 296, 306,
                     308, 312],
          "pre":    [19, 95, 101, 159, 208, 221, 240, 259, 311],
          "post":   [92, 113, 119, 127, 135, 176, 214, 223, 231, 260, 297,
                     317]},
    "B": {"shared": [82, 93, 109, 146, 149, 212, 220, 244, 246, 251, 260,
                     280, 286, 289],
          "pre":    [1, 14, 247, 263, 270, 276, 298],
          "post":   [65, 95, 125, 188, 197, 205, 216, 231, 258, 268, 305]},
}

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, c1 in enumerate(_b):
    for j, c2 in enumerate(_b):
        for k, c3 in enumerate(_b):
            CODON[c1 + c2 + c3] = _aa[i * 16 + j * 4 + k]


def read_fasta(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur)); cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def consensus_from_nt(seqs):
    ncod = len(seqs[0]) // 3
    out = []
    for c in range(ncod):
        col = [CODON.get(s[c*3:c*3+3].upper().replace("U", "T")) for s in seqs]
        col = [x for x in col if x]
        out.append(Counter(col).most_common(1)[0][0] if col else "X")
    return "".join(out)


def consensus_from_aa(seqs):
    n = max(len(s) for s in seqs)
    out = []
    for i in range(n):
        col = [s[i] for s in seqs if i < len(s) and s[i] not in "-*X"]
        out.append(Counter(col).most_common(1)[0][0] if col else "X")
    return "".join(out)


def needleman_wunsch(a, b, match=2, mismatch=-1, gap=-6):
    """Global alignment. Returns list of (i, j) index pairs, -1 for gaps."""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = score[i-1][0] + gap
    for j in range(1, m + 1):
        score[0][j] = score[0][j-1] + gap
    for i in range(1, n + 1):
        ai = a[i-1]
        row, prev = score[i], score[i-1]
        for j in range(1, m + 1):
            d = prev[j-1] + (match if ai == b[j-1] else mismatch)
            row[j] = d if d >= prev[j] + gap else prev[j] + gap
            if row[j] < row[j-1] + gap:
                row[j] = row[j-1] + gap
    pairs, i, j = [], n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + \
                (match if a[i-1] == b[j-1] else mismatch):
            pairs.append((i-1, j-1)); i -= 1; j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap:
            pairs.append((i-1, -1)); i -= 1
        else:
            pairs.append((-1, j-1)); j -= 1
    return pairs[::-1]


def find_reference(subtype):
    for name in (f"nextclade.cds_translation.G.rsv{subtype}.fasta",
                 f"nextclade.cds_translation.G.fasta"):
        if os.path.exists(name):
            return name
    return None


print("=" * 72)
print("MAPPING G-GENE ALIGNMENT COLUMNS TO REFERENCE NUMBERING")
print("=" * 72)

missing = []
for subtype, (pre_ds, _post_ds) in DATASETS.items():
    print(f"\n{'='*72}\nRSV-{subtype}")
    ref_path = find_reference(subtype)
    if not ref_path:
        print("  No Nextclade G translation found. Looked for:")
        print(f"    nextclade.cds_translation.G.rsv{subtype}.fasta")
        print( "    nextclade.cds_translation.G.fasta")
        print("  Export it from your Nextclade run and put it beside this")
        print("  script. Not falling back to another reference, because the")
        print("  numbering must match Tables 11-16.")
        missing.append(subtype)
        continue

    qpath = os.path.join(BASE, pre_ds, "G.fasta")
    if not os.path.exists(qpath):
        print(f"  missing {qpath}")
        continue

    q = consensus_from_nt(read_fasta(qpath))
    r = consensus_from_aa(read_fasta(ref_path))
    print(f"  your G consensus : {len(q)} aa   ({qpath})")
    print(f"  reference        : {len(r)} aa   ({ref_path})")

    pairs = needleman_wunsch(q, r)
    qi2ri = {i: j for i, j in pairs if i >= 0}
    ident = sum(1 for i, j in pairs if i >= 0 and j >= 0 and q[i] == r[j])
    aligned = sum(1 for i, j in pairs if i >= 0 and j >= 0)
    print(f"  aligned {aligned} positions, {100*ident/max(aligned,1):.1f}% identity")

    # describe the offset structure
    offs = sorted({(j - i) for i, j in pairs if i >= 0 and j >= 0})
    print(f"  distinct offsets across the gene: {offs}")
    if len(offs) == 1:
        print(f"  -> a single offset of {offs[0]:+d} applies to the whole gene")
    else:
        print("  -> offset CHANGES along the gene (duplication/indel present);")
        print("     each site must be mapped individually, as below")

    print(f"\n  {'period':<9}{'your col':>9}  ->  {'reference':<10}{'residue'}")
    print("  " + "-" * 46)
    for period, sites in G_SITES[subtype].items():
        for s in sites:
            j = qi2ri.get(s - 1, -1)
            if j < 0:
                print(f"  {period:<9}{s:>9}  ->  {'no match':<10}"
                      f"(insertion relative to reference)")
            else:
                print(f"  {period:<9}{s:>9}  ->  {j+1:<10}{r[j]}")

print("\n" + "=" * 72)
if missing:
    print("INCOMPLETE: no reference found for " + ", ".join(missing))
    print("Tables 17/18 cannot be finalised until G is mapped.")
else:
    print("Send this whole output back and the corrected Tables 17 and 18")
    print("can be generated with F and G both in reference numbering.")
print("=" * 72)
