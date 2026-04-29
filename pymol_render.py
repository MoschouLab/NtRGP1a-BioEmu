"""
NtRGP1a – Figure 4: PyMOL 3.1.0 publication render.

Renders WT and C8S average structures colored by B-factor (RMSF in Å),
blue = rigid → white → red = flexible, ray-traced cartoon ribbon.

Steps:
  1. Parse RMSF range from both PDB B-factor columns (shared scale).
  2. Run a headless PyMOL script (pymol -cq) to produce individual PNGs.
  3. Compose the two panels into a side-by-side publication figure.

Output:
  figures/WT_fig4_pymol_render.png       – ray-traced WT
  figures/C8S_fig4_pymol_render.png      – ray-traced C8S
  figures/fig4_publication_render.png    – composite with shared colorbar
"""

import os
import subprocess
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.image as mpimg

WORK_DIR  = os.path.abspath(".")
PYMOL_BIN = "/home/moschou/miniconda3/envs/bioemu/bin/pymol"
VARIANTS  = {
    "WT":  os.path.join(WORK_DIR, "structures/WT_average_rmsf_allatom.pdb"),
    "C8S": os.path.join(WORK_DIR, "structures/C8S_average_rmsf_allatom.pdb"),
}
OUT_DIR = os.path.join(WORK_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# B-factor parser – used to set a shared color scale across both structures
# ---------------------------------------------------------------------------

def _parse_bfactors(pdb_path: str) -> np.ndarray:
    bfac = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                bfac.append(float(line[60:66]))
    return np.array(bfac)


# ---------------------------------------------------------------------------
# Step 1 – determine global RMSF range
# ---------------------------------------------------------------------------

all_bfac   = np.concatenate([_parse_bfactors(p) for p in VARIANTS.values()])
global_min = float(all_bfac.min())
global_max = float(all_bfac.max())
print(f"RMSF range (global): {global_min:.2f} – {global_max:.2f} Å")


# ---------------------------------------------------------------------------
# Step 2 – headless PyMOL render of individual structures
# ---------------------------------------------------------------------------

_pymol_script = f"""\
from pymol import cmd
import os

variants = {repr(dict(VARIANTS))}
out_dir   = {repr(OUT_DIR)}
b_min     = {global_min:.4f}
b_max     = {global_max:.4f}

# Quality settings (applied once; cmd.delete("all") preserves them)
cmd.set("antialias",            2)
cmd.set("ray_opaque_background", 1)
cmd.set("ray_shadows",           0)
cmd.set("ray_trace_mode",        1)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_smooth_loops",  1)
cmd.bg_color("white")

for label, pdb in variants.items():
    cmd.delete("all")
    cmd.load(pdb, label)
    cmd.hide("everything", "all")
    cmd.show("cartoon", "all")
    cmd.hide("lines",    "all")
    cmd.set("cartoon_gap_cutoff",  0)
    cmd.set("cartoon_loop_radius", 0.4)
    cmd.set("cartoon_sampling",    14)
    cmd.spectrum(expression="b", palette="blue_white_red",
                 selection="all", minimum=b_min, maximum=b_max)
    cmd.orient(label)
    out = os.path.join(out_dir, f"{{label}}_fig4_pymol_render.png")
    cmd.ray(2400, 1800)
    cmd.png(out, dpi=300)
    print(f"  Saved: {{out}}", flush=True)
"""

with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as tf:
    tf.write(_pymol_script)
    _tmp = tf.name

print("Running PyMOL headless...")
result = subprocess.run(
    [PYMOL_BIN, "-cq", _tmp],
    capture_output=True, text=True, cwd=WORK_DIR
)
os.unlink(_tmp)

if result.stdout:
    print(result.stdout.rstrip())
if result.returncode != 0:
    print("--- PyMOL stderr ---")
    print(result.stderr[-3000:])
    raise RuntimeError(f"PyMOL exited with code {result.returncode}")


# ---------------------------------------------------------------------------
# Step 3 – composite publication figure (matplotlib, shared colorbar)
# ---------------------------------------------------------------------------

wt_img  = mpimg.imread(os.path.join(OUT_DIR, "WT_fig4_pymol_render.png"))
c8s_img = mpimg.imread(os.path.join(OUT_DIR, "C8S_fig4_pymol_render.png"))

fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="white",
                         gridspec_kw={"wspace": 0.04})
fig.suptitle(
    "NtRGP1a  –  Average structure colored by RMSF (Å)\n"
    "blue = rigid  ·  red = flexible",
    fontsize=13, fontweight="bold"
)

for ax, img, label in zip(axes, [wt_img, c8s_img], ["WT", "C8S"]):
    ax.imshow(img)
    ax.set_title(label, fontsize=15, fontweight="bold", pad=6)
    ax.axis("off")

# bwr matches PyMOL's blue_white_red palette
norm = mcolors.Normalize(vmin=global_min, vmax=global_max)
sm   = plt.cm.ScalarMappable(cmap=plt.cm.bwr, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes.tolist(), orientation="vertical",
                    shrink=0.65, pad=0.03, aspect=28)
cbar.set_label("RMSF (Å)", fontsize=11)
cbar.ax.tick_params(labelsize=9)

out_pub = os.path.join(OUT_DIR, "fig4_publication_render.png")
fig.savefig(out_pub, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"  Saved: {out_pub}")
print("Done.")
