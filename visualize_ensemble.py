"""
NtRGP1a BioEmu conformational ensemble visualisation.
Produces 3 figures per variant (WT and C8S):
  Fig 1 – average Cα structure (3-D ribbon)
  Fig 2 – ensemble overlay / spaghetti plot (blue→red by frame)
  Fig 3 – average structure colored by per-residue RMSF (main deliverable)
"""

import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

os.makedirs("figures", exist_ok=True)
os.makedirs("structures", exist_ok=True)

VARIANTS = {
    "WT":  ("output_2000/topology.pdb",     "output_2000/samples.xtc"),
    "C8S": ("output_2000_C8S/topology.pdb", "output_2000_C8S/samples.xtc"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _patch_bfactors(pdb_path: str, bfactors: np.ndarray) -> None:
    """Overwrite the B-factor column (cols 61-66) of a PDB with given values."""
    lines = open(pdb_path).readlines()
    out = []
    atom_idx = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and atom_idx < len(bfactors):
            bf = float(bfactors[atom_idx])
            line = line[:60] + f"{bf:6.2f}" + line[66:]
            atom_idx += 1
        out.append(line)
    with open(pdb_path, "w") as f:
        f.writelines(out)


def _residue_ids(traj_ca: md.Trajectory) -> np.ndarray:
    return np.array([a.residue.resSeq for a in traj_ca.topology.atoms])


# ---------------------------------------------------------------------------
# Per-variant pipeline
# ---------------------------------------------------------------------------

def process_variant(label: str, top_path: str, traj_path: str) -> None:
    print(f"\n=== {label} ===")
    traj = md.load(traj_path, top=top_path)
    print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms")

    ca_idx = traj.topology.select("name CA")
    traj_ca = traj.atom_slice(ca_idx)

    # Align every frame to frame 0 (Cα RMSD)
    traj_ca.superpose(traj_ca, frame=0)

    res_ids = _residue_ids(traj_ca)
    n_res = traj_ca.n_atoms
    n_frames = traj_ca.n_frames

    # Mean Cα coordinates (nm → Å for display; keep nm for MDTraj ops)
    mean_xyz_nm = traj_ca.xyz.mean(axis=0)          # (n_res, 3) nm
    mean_xyz_A  = mean_xyz_nm * 10                   # Å

    # ------------------------------------------------------------------
    # Figure 1 — Average structure
    # ------------------------------------------------------------------
    avg_traj = md.Trajectory(mean_xyz_nm[np.newaxis], traj_ca.topology)
    avg_pdb = f"structures/{label}_average_structure.pdb"
    avg_traj.save_pdb(avg_pdb)
    print(f"  Saved: {avg_pdb}")

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    x, y, z = mean_xyz_A[:, 0], mean_xyz_A[:, 1], mean_xyz_A[:, 2]
    ax.plot(x, y, z, color="steelblue", lw=1.8, zorder=2)
    ax.scatter(x, y, z, c="steelblue", s=12, alpha=0.85, zorder=3)
    ax.set_title(f"{label}  –  Average Cα Structure", fontsize=12, fontweight="bold")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.tick_params(labelsize=7)
    plt.tight_layout()
    fig1_path = f"figures/{label}_fig1_average_structure.png"
    fig.savefig(fig1_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {fig1_path}")

    # ------------------------------------------------------------------
    # Figure 2 — Ensemble overlay (spaghetti plot)
    # ------------------------------------------------------------------
    cmap2 = plt.cm.coolwarm
    norm2 = mcolors.Normalize(vmin=0, vmax=n_frames - 1)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    for i in range(n_frames):
        c = cmap2(norm2(i))
        xi = traj_ca.xyz[i, :, 0] * 10
        yi = traj_ca.xyz[i, :, 1] * 10
        zi = traj_ca.xyz[i, :, 2] * 10
        ax.plot(xi, yi, zi, color=c, lw=0.25, alpha=0.12)

    sm = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.1)
    cbar.set_label("Frame index", fontsize=9)
    cbar.set_ticks([0, n_frames // 2, n_frames - 1])

    ax.set_title(f"{label}  –  Cα Ensemble Overlay  ({n_frames} conformations)",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.tick_params(labelsize=7)
    plt.tight_layout()
    fig2_path = f"figures/{label}_fig2_ensemble_overlay.png"
    fig.savefig(fig2_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {fig2_path}")

    # ------------------------------------------------------------------
    # Figure 3 — Average structure colored by per-residue RMSF
    # ------------------------------------------------------------------
    mean_traj = md.Trajectory(mean_xyz_nm[np.newaxis], traj_ca.topology)
    rmsf_nm = md.rmsf(traj_ca, mean_traj, frame=0,
                      atom_indices=np.arange(n_res))   # (n_res,) in nm
    rmsf_A = rmsf_nm * 10                               # Å

    # Save RMSF-annotated Cα PDB (RMSF stored in B-factor column)
    rmsf_pdb = f"structures/{label}_average_rmsf.pdb"
    avg_traj.save_pdb(rmsf_pdb)
    _patch_bfactors(rmsf_pdb, rmsf_A)
    print(f"  Saved: {rmsf_pdb}")

    # All-atom average structure for PyMOL cartoon rendering
    # Align full trajectory with the same Cα-RMSD superposition, then
    # take the mean over all heavy atoms and broadcast per-residue RMSF
    # onto every atom of that residue via resSeq lookup.
    heavy_idx = traj.topology.select("not element H")
    traj.superpose(traj, frame=0, atom_indices=ca_idx)
    traj_heavy   = traj.atom_slice(heavy_idx)
    mean_heavy_nm = traj_heavy.xyz.mean(axis=0)
    avg_heavy_traj = md.Trajectory(mean_heavy_nm[np.newaxis], traj_heavy.topology)

    res_to_rmsf  = {a.residue.resSeq: rmsf_A[i]
                    for i, a in enumerate(traj_ca.topology.atoms)}
    rmsf_allatom = np.array([res_to_rmsf.get(a.residue.resSeq, 0.0)
                              for a in traj_heavy.topology.atoms])

    allatom_pdb = f"structures/{label}_average_rmsf_allatom.pdb"
    avg_heavy_traj.save_pdb(allatom_pdb)
    _patch_bfactors(allatom_pdb, rmsf_allatom)
    print(f"  Saved: {allatom_pdb}")

    # 3-D structure colored by RMSF + RMSF profile sub-panel
    cmap3 = plt.cm.coolwarm
    norm3 = mcolors.Normalize(vmin=rmsf_A.min(), vmax=rmsf_A.max())

    fig = plt.figure(figsize=(12, 5))
    ax3d = fig.add_subplot(121, projection="3d")

    # Draw backbone segments colored by RMSF of the N-terminal residue of each bond
    for i in range(n_res - 1):
        seg_x = mean_xyz_A[i:i+2, 0]
        seg_y = mean_xyz_A[i:i+2, 1]
        seg_z = mean_xyz_A[i:i+2, 2]
        ax3d.plot(seg_x, seg_y, seg_z, color=cmap3(norm3(rmsf_A[i])), lw=2.5)

    sc = ax3d.scatter(
        mean_xyz_A[:, 0], mean_xyz_A[:, 1], mean_xyz_A[:, 2],
        c=rmsf_A, cmap=cmap3, norm=norm3, s=25, zorder=5
    )
    cbar3 = fig.colorbar(sc, ax=ax3d, shrink=0.5, pad=0.1)
    cbar3.set_label("RMSF (Å)", fontsize=9)
    ax3d.set_title(f"{label}  –  Average Structure\nColored by RMSF",
                   fontsize=11, fontweight="bold")
    ax3d.set_xlabel("x (Å)"); ax3d.set_ylabel("y (Å)"); ax3d.set_zlabel("z (Å)")
    ax3d.tick_params(labelsize=7)

    # RMSF profile (residue index on x-axis)
    ax2d = fig.add_subplot(122)
    ax2d.fill_between(res_ids, rmsf_A, alpha=0.3,
                      color=plt.cm.coolwarm(0.85))
    ax2d.plot(res_ids, rmsf_A, color=plt.cm.coolwarm(0.85), lw=1.2)
    ax2d.set_xlabel("Residue number", fontsize=10)
    ax2d.set_ylabel("RMSF (Å)", fontsize=10)
    ax2d.set_title(f"{label}  –  Per-residue RMSF profile", fontsize=11,
                   fontweight="bold")
    ax2d.set_xlim(res_ids[0], res_ids[-1])
    ax2d.tick_params(labelsize=8)
    ax2d.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig3_path = f"figures/{label}_fig3_rmsf_colored.png"
    fig.savefig(fig3_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {fig3_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    for label, (top, traj) in VARIANTS.items():
        process_variant(label, top, traj)
    print("\nDone.  Figures → ./figures/    Structures → ./structures/")
