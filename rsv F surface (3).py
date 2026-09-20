#!/usr/bin/env python3
"""
rsv_F_surface.py -- RSV prefusion F trimer as a SOLID SURFACE with
positively selected codon sites highlighted by rollout period.

=======================================================================
HOW TO RUN
=======================================================================
INTERACTIVE (opens the PyMOL window, keeps it open so you can rotate,
zoom and adjust before saving) -- this is what you want:

    pymol rsv_F_surface.py -- A          # RSV-A, structure 4MMV
    pymol rsv_F_surface.py -- B          # RSV-B, structure 5UDE

Note: NO -c and NO -q. Those flags mean "headless, no GUI" and are why
the previous run only wrote a PNG file.

Once the window is open you can keep typing commands, e.g.
    set transparency, 0.3        # see interior
    show cartoon, F*             # cartoon under the surface
    hide surface, F*             # switch back to cartoon only
    png myfigure.png, width=2400, height=2000, dpi=300, ray=1

HEADLESS (batch, writes PNGs only):
    pymol -cq rsv_F_surface.py -- A png

=======================================================================
WHAT WAS WRONG WITH THE PREVIOUS IMAGE
=======================================================================
The trimer is built by C3 crystal symmetry, so PyMOL loads the assembly
as 3 STATES of one object. The old script merged those states into a
single object, which put three chains with identical chain IDs on top of
each other. PyMOL then tried to draw one continuous backbone through all
three and produced the shattered ribbon you saw. Here the states are
split into three separate objects (F_1, F_2, F_3) and grouped, so each
protomer keeps its own connectivity and the surface closes properly.

=======================================================================
STRUCTURES (verified on RCSB)
    4MMV  prefusion RSV F DS-Cav1, strain A2       -> RSV-A
    5UDE  RSV F B9320 DS-Cav1, strain 9320 (B)     -> RSV-B
          5UDE carries S155C, S190F, V207L, S290C. V207L is 4 residues
          from site 203 and S190F is inside Site O, so describe it as a
          "prefusion-stabilised DS-Cav1 variant", not "RSV F".

NUMBERING: resolved. These are reference F0 positions. The coverage
report still prints the residue identity found at each site; it should
now agree with the reference residue, allowing for the engineered
mutations listed above.
=======================================================================
"""

import sys
from pymol import cmd

# REFERENCE F0 numbering (NOT raw HyPhy alignment columns).
# Offsets applied: RSV-A +2, RSV-B +6, determined by global alignment of
# each consensus translation against the reference (97% identity, one
# offset across the whole gene). Matches corrected Tables 17 and 18.
SITES = {
    "A": {"shared": [], "pre": [117, 515, 554], "post": [120, 543]},
    "B": {"shared": [12, 118, 209], "pre": [68, 125, 529], "post": []},
}
DEFAULT_PDB = {"A": "4MMV", "B": "5UDE"}

NON_ECTO = [(1, 26, "signal peptide (cleaved)"),
            (110, 136, "p27 (excised by furin)"),
            (530, 550, "transmembrane"),
            (551, 574, "cytoplasmic tail")]

# Order matters: later entries overwrite earlier ones where ranges
# overlap. Sites O (196-212) and V (155-206) share residues 196-206, so
# Site O is painted LAST and displayed intact; state this in the legend.
# Colours are deliberately pale so that the saturated site colours below
# remain distinguishable - do not use a blue here, it collides with the
# pre-rollout marine.
ANTIGENIC = [
    ("site_V",  "155-206+305-310", "paleyellow"),
    ("site_II", "254-277",         "lightpink"),
    ("site_IV", "422-471",         "palegreen"),
    ("site_O",  "62-69+196-212",   "lightorange"),
]
COLORS = {"shared": "purple", "pre": "marine", "post": "red"}


def non_ecto(r):
    for lo, hi, name in NON_ECTO:
        if lo <= r <= hi:
            return name
    return None


args = [a for a in sys.argv[1:] if a != "--"]
subtype = (args[0].upper() if args else "B")
if subtype not in SITES:
    sys.exit("first argument must be A or B")
headless = len(args) > 1 and args[1].lower() == "png"
pdb = DEFAULT_PDB[subtype]

print("=" * 64)
print(f"RSV-{subtype} F   structure {pdb}")
print("=" * 64)

cmd.reinitialize()

# ---- load the biological assembly as separate protomers -------------
cmd.set("assembly", "1")
cmd.fetch(pdb, "asm", async_=0)
n = cmd.count_states("asm")
print(f"assembly states (protomers): {n}")
if n > 1:
    cmd.split_states("asm", prefix="F_")
    cmd.delete("asm")
else:
    cmd.set_name("asm", "F_1")
cmd.set("assembly", "")
cmd.group("trimer", "F_*")
cmd.remove("not polymer")
cmd.remove("solvent")

prot = "F_*"
# split_states zero-pads (F_0001, not F_1) -- take the real name
FIRST = cmd.get_object_list("F_*")[0]
print(f"protomers loaded: {cmd.get_object_list('F_*')}  (first = {FIRST})")
print(f"CA atoms total  : {cmd.count_atoms(prot + ' and name CA')}")

# ---- solid surface --------------------------------------------------
cmd.hide("everything")
cmd.show("surface", prot)
cmd.color("grey70", prot)
cmd.bg_color("white")
cmd.set("surface_quality", 1)
cmd.set("transparency", 0.0)          # fully solid
cmd.set("surface_smooth_edges", 1)
cmd.set("two_sided_lighting", 1)
cmd.set("ambient", 0.18)
cmd.set("specular", 0.15)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("ray_opaque_background", 1)

# ---- antigenic sites as coloured surface patches --------------------
for name, resi, colr in ANTIGENIC:
    cmd.select(name, f"{prot} and resi {resi}")
    cmd.color(colr, name)

# ---- selected sites -------------------------------------------------
report = []
for group, resis in SITES[subtype].items():
    keep = []
    for r in resis:
        bad = non_ecto(r)
        if bad:
            report.append((r, group, "EXCLUDED", bad))
            continue
        if cmd.count_atoms(f"{prot} and resi {r} and name CA") == 0:
            report.append((r, group, "MISSING", "unresolved in this structure"))
        else:
            aa = []
            cmd.iterate(f"{FIRST} and resi {r} and name CA",
                        "aa.append(oneletter)", space={"aa": aa})
            report.append((r, group, "PRESENT",
                           f"modelled as {aa[0] if aa else '?'}{r}"))
            keep.append(r)
    sel = f"sel_{group}"
    if keep:
        cmd.select(sel, f"{prot} and resi " + "+".join(map(str, keep)))
        cmd.color(COLORS[group], sel)        # colour the surface patch
        # spheres on top so small sites stay visible against the surface
        cmd.create(f"spheres_{group}", sel + " and name CA")
        cmd.show("spheres", f"spheres_{group}")
        cmd.color(COLORS[group], f"spheres_{group}")
    else:
        cmd.select(sel, "none")
cmd.set("sphere_scale", 1.1, "spheres_*")

# ---- labels, one copy per residue -----------------------------------
cmd.select("lab", f"(sel_shared or sel_pre or sel_post) and {FIRST} and name CA")
cmd.label("lab", '"%s%s" % (oneletter, resi)')
cmd.set("label_size", 20)
cmd.set("label_color", "black")
cmd.set("label_outline_color", "white")
cmd.set("label_position", (0, 3.0, 0))
cmd.deselect()

cmd.orient(prot)
cmd.turn("x", -12)
cmd.zoom(prot, 3)

# ---- report ---------------------------------------------------------
print("\n" + "-" * 64)
print(f"{'site':<8}{'period':<10}{'status':<10}note")
for r, g, st, note in sorted(report):
    print(f"F{r:<7}{g:<10}{st:<10}{note}")
n_ok = sum(1 for _, _, s, _ in report if s == "PRESENT")
print(f"\n{n_ok} of {len(report)} RSV-{subtype} F sites shown on the surface.")
print("Excluded sites are not in the mature ectodomain -- show them on the")
print("linear F0 domain map instead, and say so in the legend.")
print("Verify the one-letter codes above against your alignment consensus.")

if headless:
    for i in range(3):
        cmd.png(f"rsv{subtype}_F_surface_view{i+1}.png",
                width=2400, height=2000, dpi=300, ray=1)
        cmd.turn("y", 120)
    print(f"\nwrote rsv{subtype}_F_surface_view1-3.png")
else:
    print("\nPyMOL window is open. Rotate with the mouse. Useful commands:")
    print("   set transparency, 0.35        # ONLY for exploring - a")
    print("                                 # semi-transparent surface")
    print("                                 # is not a publication figure")
    print("   show cartoon, F_*             # ribbon under the surface")
    print("   png fig.png, width=2400, height=2000, dpi=300, ray=1")
cmd.save(f"rsv{subtype}_F_surface.pse")
print(f"session saved: rsv{subtype}_F_surface.pse")
