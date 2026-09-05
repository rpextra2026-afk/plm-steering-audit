# stitch_fig1_structures.py -- combine the two equal-scale panels into one PDF.
# Run after fig1_structures.py, from the same folder.
#
# The assert is the shared-scale tripwire: the two PyMOL panels were rendered at
# an identical view and identical pixel size, so they MUST be the same size here.
# If they are not, the scale drifted -- stop and do not use the output.
from PIL import Image

a = Image.open("panel_confident.png").convert("RGB")   # left  = confident 258-mer
b = Image.open("panel_collapsed.png").convert("RGB")    # right = collapsed 135-mer
assert a.size == b.size, "panels differ in size -> scale is NOT shared, stop"

canvas = Image.new("RGB", (a.width + b.width, a.height), "white")
canvas.paste(a, (0, 0))
canvas.paste(b, (a.width, 0))
canvas.save("fig1_structures.pdf", "PDF", resolution=300.0)
print("wrote fig1_structures.pdf")

# To stack the two molecules vertically instead (a taller, narrower panel), use:
#     canvas = Image.new("RGB", (a.width, a.height + b.height), "white")
#     canvas.paste(a, (0, 0)); canvas.paste(b, (0, a.height))
