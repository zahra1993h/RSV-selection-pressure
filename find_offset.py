#!/usr/bin/env python3
"""
find_offset.py -- determine the exact offset between YOUR alignment
columns and canonical F0 residue numbering, using the UniProt reference
protein rather than a PDB entry.

RUN (plain python, from the Desktop, inside WSL - needs internet):

    python3 find_offset.py

WHY NOT USE THE PDB
    check_numbering.py compared your alignment to 4MMV / 5UDE. That test
    is unreliable, because both are engineered constructs:

      4MMV (DS-Cav1) is the F ectodomain TRUNCATED AT RESIDUE 513 and
      fused to a fibritin "foldon" trimerisation tag. Anything numbered
      above ~513 in that file is fibritin, NOT RSV F. So the RSV-A sites
      541 and 552 were being compared against a bacteriophage tag, which
      is why no offset could ever satisfy them.

      Both constructs also carry engineered point mutations (S155C,
      S190F, V207L, S290C), so individual residues can disagree with a
      wild-type alignment even when the numbering is perfectly correct.

    The canonical F0 precursor sequence is the right reference. This
    script downloads it and slides your consensus against it to find the
    offset that maximises identity across the WHOLE protein, not just at
    a handful of sites.

REFERENCES USED
    P03420  Fusion glycoprotein F0, human RSV strain A2       (574 aa)
    Q6V2E7  Fusion glycoprotein F0, human RSV strain B 9320   (574 aa)

OUTPUT
    For each subtype: the best offset, the identity it achieves, and a
    corrected list of your selected sites in true F0 numbering.
"""

import os
import sys
import urllib.request
from collections import Counter

BASE = "RSV seq data/hyphy_results"
REF = {"A": ("P03420", "rsvA_prevaccine", "rsvA_postvaccine"),
       "B": ("Q6V2E7", "rsvB_prevaccine", "rsvB_postvaccine")}
SITES = {"A": {"pre": [115, 513, 552], "post": [118, 541]},
         "B": {"shared": [6, 112, 203], "pre": [62, 119, 523]}}

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


def consensus_protein(seqs):
    """Majority amino acid at every codon column."""
    ncod = len(seqs[0]) // 3
    out = []
    for c in range(ncod):
        col = []
        for s in seqs:
            cdn = s[c * 3:c * 3 + 3].upper().replace("U", "T")
            aa = CODON.get(cdn)
            if aa:
                col.append(aa)
        out.append(Counter(col).most_common(1)[0][0] if col else "X")
    return "".join(out)


def fetch_uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            txt = r.read().decode()
    except Exception as e:
        print(f"  could not download {acc}: {e}")
        print(f"  open {url} in a browser, save it, and place it next to")
        print(f"  this script as {acc}.fasta")
        p = f"{acc}.fasta"
        if os.path.exists(p):
            txt = open(p).read()
            print(f"  using local {p}")
        else:
            return None
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))


def best_offset(query, ref, lo=-30, hi=30):
    """Slide query along ref; return sorted (identity, offset, n_compared)."""
    res = []
    for off in range(lo, hi + 1):
        same = tot = 0
        for i, aa in enumerate(query):
            j = i + off                      # 0-based index into ref
            if 0 <= j < len(ref) and aa != "X":
                tot += 1
                if aa == ref[j]:
                    same += 1
        if tot >= 50:
            res.append((same / tot, off, tot, same))
    res.sort(reverse=True)
    return res


print("=" * 70)
print("FINDING TRUE F0 NUMBERING OFFSET (reference: UniProt, not PDB)")
print("=" * 70)

summary = {}
for subtype, (acc, pre_ds, post_ds) in REF.items():
    print(f"\nRSV-{subtype}   reference {acc}")
    ref = fetch_uniprot(acc)
    if not ref:
        continue
    print(f"          reference length {len(ref)} aa")

    path = os.path.join(BASE, pre_ds, "F.fasta")
    if not os.path.exists(path):
        print(f"  missing {path} -- run this from the Desktop")
        continue
    seqs = read_fasta(path)
    q = consensus_protein(seqs)
    print(f"          your alignment {len(seqs)} seqs, {len(q)} codons")

    ranked = best_offset(q, ref)
    if not ranked:
        print("  could not compare")
        continue
    ident, off, tot, same = ranked[0]
    print(f"\n  best alignment: offset {off:+d}  ->  "
          f"{100*ident:.1f}% identity ({same}/{tot} residues)")
    for i2, o2, t2, s2 in ranked[1:3]:
        print(f"  next best    : offset {o2:+d}  ->  {100*i2:.1f}%")

    # column N of your alignment == F0 residue N + off
    if ident < 0.80:
        print("\n  WARNING: identity below 80%. Either the wrong reference")
        print("  strain, or the extraction is not a clean F reading frame.")
        print("  Resolve this before trusting any position.")
    elif off == 0:
        print("\n  -> your columns ALREADY equal F0 numbering. No change needed.")
    else:
        print(f"\n  -> F0 residue = your column {off:+d}")
    summary[subtype] = (off, ident)

    print("\n  corrected site numbers:")
    for period, sites in SITES[subtype].items():
        fixed = [f"{s}->{s+off}" for s in sites]
        print(f"    {period:<8}{', '.join(fixed)}")

print("\n" + "=" * 70)
print("WHAT TO DO WITH THIS")
print("=" * 70)
for st, (off, ident) in summary.items():
    if off == 0:
        print(f"RSV-{st}: no change needed.")
    else:
        print(f"RSV-{st}: add {off:+d} to every position before quoting it "
              f"against\n        antigenic sites, structures or the literature.")
print("\nThe HyPhy statistics themselves are unaffected -- only the residue")
print("LABELS change. Tables 17/18, both figures and every antigenic-site")
print("assignment in the Results text must be renumbered consistently.")
print("Re-check the Site O / Site V assignments after renumbering: a shift")
print("can move a site into or out of an epitope boundary.")
