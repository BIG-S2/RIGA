#!/usr/bin/env python3
"""supplementary_figure10_partial_residuals.py (2026-10-14; panels e and f moved to Supplementary Fig. 10 at Prof. Zhu's suggestion) -- compose panels a-d of main
Figure 1 (revision/figures/fig1_v3/{a..f}.png from figure1_overview_panels.R, each drawn at its final
size at 300 dpi) WITHOUT rescaling, so the type on the page is the type in the panels (6-8 pt).
Layout (7.2 in wide): band 1 = a+b (full width); band 2 = c | d. Writes revision/manuscript/figures/data_v4.png and the
supplementary composite of e | f as revision/manuscript/figures/supp_fig_S10_partial_residuals.png."""
from PIL import Image, ImageDraw, ImageFont
import pathlib
D = pathlib.Path("${ADSP_ROOT}/revision/figures/fig1_v3")
OUT = pathlib.Path("${ADSP_ROOT}/revision/manuscript/figures/data_v4.png")
OUT_SUPP = pathlib.Path("${ADSP_ROOT}/revision/manuscript/figures/supp_fig_S10_partial_residuals.png")
DPI = 300; W = int(7.2 * DPI); LAB = 42; GAP = 12
font = ImageFont.truetype("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans-Bold.ttf", 40)   # ~9.5 pt on the page

def white_bg(path):
    im = Image.open(path).convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    return Image.alpha_composite(bg, im).convert("RGB")
def fit_w(img, w):
    r = w / img.width
    return img.resize((w, int(img.height * r)), Image.LANCZOS)

P = {k: white_bg(D / f"{k}.png") for k in "acdef"}     # a.png holds the a+b block with its own tags
ab = P["a"]; d = P["d"]; e = P["e"]; f = P["f"]
if ab.width > W: ab = fit_w(ab, W)
c = fit_w(P["c"], int(2.7 * DPI))                        # the BIC curve is a pre-rendered image (2.7 in wide)
band1_h = ab.height
band2_h = LAB + max(c.height, d.height)
band3_h = LAB + max(e.height, f.height)
H = band1_h + GAP + band2_h + GAP
canvas = Image.new("RGB", (W, H), "white"); draw = ImageDraw.Draw(canvas)
def put(img, x, y, label):
    draw.text((x + 6, y), label, fill="black", font=font)
    canvas.paste(img, (x, y + LAB))
canvas.paste(ab, (0, 0))
y2 = band1_h + GAP
put(c, 0, y2, "c")
put(d, W - d.width, y2, "d")
canvas.save(OUT, dpi=(DPI, DPI))
supp = Image.new("RGB", (W, band3_h + GAP), "white"); draw = ImageDraw.Draw(supp)
draw.text((6, 0), "a", fill="black", font=font); supp.paste(e, (0, LAB)); draw.text((W - f.width + 6, 0), "b", fill="black", font=font); supp.paste(f, (W - f.width, LAB))
supp.save(OUT_SUPP, dpi=(DPI, DPI)); print("wrote", OUT_SUPP, supp.size)
print("wrote", OUT, canvas.size, "=", round(canvas.width / DPI, 2), "x", round(canvas.height / DPI, 2), "in")
