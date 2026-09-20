#!/usr/bin/env python3
"""
rsv_F_domain_map.py -- linear domain map of the RSV F0 precursor showing
EVERY positively selected codon site, including the ones that cannot
appear on a prefusion structure.

RUN
    python3 rsv_F_domain_map.py

Writes rsv_F_domain_map.png / .pdf (300 dpi, publication ready).

WHY THIS FIGURE EXISTS
    Four of five RSV-A sites and three of six RSV-B sites fall in the
    signal peptide, p27, transmembrane region or cytoplasmic tail. None
    of those are present in a soluble prefusion ectodomain, so a PyMOL
    figure alone would silently drop most of your F results. This panel
    shows all of them against the F0 precursor architecture, and makes
    the interpretive point for you: sites outside the ectodomain cannot
    be antibody-escape signals, because antibodies never see them.

    Pair it with rsv_F_pymol.py as panels A and B of one figure.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

# ------------------------------------------------------------------ data
F0_LEN = 574

DOMAINS = [
    (1, 26, "SP", "#9e9e9e"),
    (27, 109, "F2", "#7fb3d5"),
    (110, 136, "p27", "#e8a33d"),
    (137, 529, "F1 ectodomain", "#82c09a"),
    (530, 550, "TM", "#b9846f"),
    (551, 574, "CT", "#9e9e9e"),
]

ANTIGENIC = [
    (62, 69, "Ø"), (196, 212, "Ø"),
    (155, 206, "V"), (305, 310, "V"),
    (254, 277, "II"),
    (422, 471, "IV"),
]

# Positions are REFERENCE F0 numbering. Offsets applied to the raw HyPhy
# alignment columns: RSV-A +2, RSV-B +6 (determined by global alignment of
# the consensus translation against the reference, 97% identity, single
# offset across the whole gene). Do not edit these back to column numbers.
SITES = {
    "RSV-A": {"shared": [], "pre": [117, 515, 554], "post": [120, 543]},
    "RSV-B": {"shared": [12, 118, 209], "pre": [68, 125, 529], "post": []},
}

COLORS = {"shared": "#7b3fa0", "pre": "#1f6fb4", "post": "#d13b2e"}
LABELS = {"shared": "Both periods", "pre": "Pre-rollout only",
          "post": "Post-rollout only"}

ECTO = [(27, 109), (137, 529)]


def in_ecto(r):
    return any(lo <= r <= hi for lo, hi in ECTO)


# ------------------------------------------------------------------ plot
fig, axes = plt.subplots(2, 1, figsize=(11.5, 6.8), sharex=True)
fig.subplots_adjust(hspace=0.55, top=0.86, bottom=0.17, left=0.07, right=0.985)

for ax, (subtype, groups) in zip(axes, SITES.items()):
    # antigenic site shading behind everything
    for lo, hi, name in ANTIGENIC:
        ax.add_patch(Rectangle((lo, -0.62), hi - lo, 1.5,
                               facecolor="#ffe9a8", edgecolor="none",
                               alpha=0.55, zorder=0))

    # domain bar
    for lo, hi, name, colr in DOMAINS:
        ax.add_patch(Rectangle((lo, -0.20), hi - lo, 0.40,
                               facecolor=colr, edgecolor="black",
                               linewidth=0.7, zorder=2))
        if hi - lo > 28:
            ax.text((lo + hi) / 2, 0.0, name, ha="center", va="center",
                    fontsize=8.5, zorder=3,
                    color="white" if name in ("F2", "F1 ectodomain") else "black")

    # furin cleavage markers
    for pos in (109, 136):
        ax.plot([pos, pos], [-0.20, 0.20], color="black", lw=1.4,
                linestyle=(0, (2, 1.4)), zorder=4)
    ax.text(122.5, 0.46, "furin", ha="center", fontsize=7, style="italic")

    # selected sites -- stagger vertically where positions are close
    # together, otherwise labels such as 115/118 and 541/552 collide
    placed = []
    allsites = sorted((r, g) for g, rs in groups.items() for r in rs)
    for r, group in allsites:
        tier = 0
        while any(abs(r - pr) < 34 and pt == tier for pr, pt in placed):
            tier += 1
        placed.append((r, tier))
        top = 0.62 + tier * 0.58
        solid = in_ecto(r)
        ax.plot([r, r], [0.20, top], color=COLORS[group], lw=1.2, zorder=5)
        ax.plot(r, top + 0.08, marker="o", markersize=7.5,
                markerfacecolor=COLORS[group] if solid else "white",
                markeredgecolor=COLORS[group], markeredgewidth=1.6,
                zorder=6)
        ax.text(r, top + 0.28, str(r), ha="center", va="bottom",
                fontsize=7.6, color=COLORS[group], zorder=6)

    # call out the residue where selection and entropy analyses converge
    if subtype == "RSV-B":
        ax.annotate("209\nSite Ø\n(Q209R)", xy=(209, 1.05), xytext=(272, 1.78),
                    fontsize=7.4, ha="center", va="center", color="#7b3fa0",
                    arrowprops=dict(arrowstyle="-", lw=0.8, color="#7b3fa0"))

    n_tot = sum(len(v) for v in groups.values())
    n_map = sum(1 for v in groups.values() for r in v if in_ecto(r))
    ax.set_title(f"{subtype} F  —  {n_tot} positively selected sites, "
                 f"{n_map} mappable onto a prefusion ectodomain structure",
                 fontsize=9.5, loc="left", pad=6)

    ax.set_xlim(-8, F0_LEN + 8)
    ax.set_ylim(-0.75, 2.30)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_position(("data", -0.62))
    ax.tick_params(axis="x", labelsize=8)

axes[1].set_xlabel("F0 precursor residue number", fontsize=9)
axes[1].set_xticks([1, 100, 200, 300, 400, 500, 574])

handles = [
    Line2D([], [], marker="o", ls="none", markersize=7.5,
           markerfacecolor=COLORS[k], markeredgecolor=COLORS[k],
           label=LABELS[k]) for k in ("shared", "pre", "post")
] + [
    Line2D([], [], marker="o", ls="none", markersize=7.5,
           markerfacecolor="white", markeredgecolor="black",
           markeredgewidth=1.6,
           label="Open symbol: outside mature ectodomain\n(not antibody-accessible)"),
    Patch(facecolor="#ffe9a8", edgecolor="none",
          label="Antigenic sites Ø, II, IV, V"),
]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
           fontsize=8, bbox_to_anchor=(0.5, -0.005))

fig.suptitle("Positively selected codon sites (≥2 of 3 HyPhy methods) "
             "on the RSV F0 precursor", fontsize=11, y=0.975)

fig.savefig("rsv_F_domain_map.png", dpi=300, bbox_inches="tight")
fig.savefig("rsv_F_domain_map.pdf", bbox_inches="tight")
print("wrote rsv_F_domain_map.png and rsv_F_domain_map.pdf")

for subtype, groups in SITES.items():
    out = []
    for g, rs in groups.items():
        for r in rs:
            out.append(f"{r}{'' if in_ecto(r) else '*'}")
    print(f"{subtype}: " + ", ".join(sorted(out, key=lambda s: int(s.rstrip('*')))))
print("* = outside the mature ectodomain (SP, p27, TM or cytoplasmic tail)")
