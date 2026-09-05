# fig1_structures.py -- Figure 1(b): the two structures at ONE SHARED SCALE
#
# Run headless from analysis/fig1/ (next to the two .pdb files):
#     pymol -cq fig1_structures.py
# then combine the panels:
#     python stitch_fig1_structures.py
#
# Shared-scale guarantee: both molecules are centered on the origin, a single
# camera view is locked with get_view/set_view, and each is rendered at that
# identical view AND identical pixel size -- so Angstroms-per-pixel is provably
# the same in both panels. The collapsed 135-mer (Rg 24.0 A) must therefore look
# LARGER than the confident 258-mer (Rg 22.2 A); that size contrast is the figure.
from pymol import cmd

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)
cmd.set("cartoon_transparency", 0)

cmd.load("fig1_ubiquitin_repeat.pdb", "confident")   # 258 res, pLDDT 85.87 (compact)
cmd.load("fig1_eerhv_repeat.pdb",     "collapsed")    # 135 res, pLDDT 28.85 (extended)
cmd.hide("everything")
cmd.show("cartoon")

# color by per-residue pLDDT (output_to_pdb wrote it into the B-factor column),
# SAME fixed range on both. blue = confident (high pLDDT), red = low-confidence.
for obj in ("confident", "collapsed"):
    cmd.spectrum("b", "red_white_blue", obj, minimum=20, maximum=90)

# center BOTH molecules on the origin so one shared camera frames each centrally
for obj in ("confident", "collapsed"):
    c = cmd.centerofmass(obj)
    cmd.translate([-c[0], -c[1], -c[2]], obj, camera=0)

# fix the shared SCALE on the LARGER (collapsed) molecule, then lock the view
cmd.orient("collapsed")
cmd.zoom("collapsed", 10)          # 10 A padding
view = cmd.get_view()              # the camera distance stored here == the shared scale

W = H = 900
cmd.disable("collapsed"); cmd.enable("confident")
cmd.set_view(view); cmd.ray(W, H); cmd.png("panel_confident.png", dpi=300)

cmd.disable("confident"); cmd.enable("collapsed")
cmd.set_view(view); cmd.ray(W, H); cmd.png("panel_collapsed.png", dpi=300)

print("wrote panel_confident.png and panel_collapsed.png")
