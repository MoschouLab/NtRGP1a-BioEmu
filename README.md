# NtRGP1a BioEmu Conformational Analysis (WT & C8S)

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19881948-blue.svg)](https://doi.org/10.5281/zenodo.19881948)
[![BioEmu](https://img.shields.io/badge/BioEmu-v1.1-orange.svg)](https://github.com/microsoft/bioemu)

Structural ensemble analysis of **NtRGP1a** (*Nicotiana tabacum* Glycine-rich
RNA-Binding Protein 1a, UniProt [D6PZY5](https://www.uniprot.org/uniprot/D6PZY5),
156 aa) and its C8S point mutant, generated with
[BioEmu](https://github.com/microsoft/bioemu) and visualised with MDTraj and PyMOL.

Part of the **PLANTEX ERC Consolidator Grant** project —
[Moschou Lab](https://www.imbb.forth.gr/en/research/labs/moschou-lab),
IMBB-FORTH, Heraklion, Greece.

---

## Repository

https://github.com/MoschouLab/NtRGP1a-BioEmu

---

## Biological context

NtRGP1a is a plant glycine-rich RNA-binding protein with an N-terminal RRM
(RNA Recognition Motif) and a C-terminal glycine-rich disordered (GRD) tail.
The C8S point mutation exchanges the single cysteine for serine, abolishing
potential redox-sensing behaviour. BioEmu was used to sample the conformational
landscape of both variants (826 conformations for WT, 837 for C8S); this
repository contains the downstream ensemble analysis and publication figures.

---

## Requirements

| Software | Version | Install |
|---|---|---|
| conda | ≥ 23 | [miniconda](https://docs.conda.io/en/latest/miniconda.html) |
| MDTraj | 1.11.0 | `conda install -c conda-forge mdtraj` |
| PyMOL Open-Source | 3.1.0 | `conda install -c conda-forge pymol-open-source` |
| matplotlib | ≥ 3.7 | included in `bioemu` env |
| numpy | ≥ 1.24 | included in `bioemu` env |

All Python dependencies are available in the `bioemu` conda environment.

---

## Raw trajectory data

Raw BioEmu output trajectories are deposited on Zenodo:

> **DOI: [10.5281/zenodo.19881948](https://doi.org/10.5281/zenodo.19881948)**

| File | Variant | Conformations |
|---|---|---|
| `samples.xtc` | WT | 826 |
| `samplesC8S.xtc` | C8S | 837 |
| `topology.pdb` | shared topology | — |

Download and place the trajectories as follows before running the scripts:

```
output_2000/
    topology.pdb
    samples.xtc
output_2000_C8S/
    topology.pdb
    samples.xtc
```

---

## Reproducing the figures

Activate the conda environment first:

```bash
conda activate bioemu
```

**Step 1 — Ensemble analysis (Figures 1–3 + intermediate PDBs)**

```bash
python visualize_ensemble.py
```

**Step 2 — PyMOL publication render (Figure 4)**

```bash
python pymol_render.py
```

Both scripts are self-contained and write all outputs automatically.

---

## Output files

### `figures/`

| File | Description |
|---|---|
| `WT_fig1_average_structure.png` | Mean Cα structure of the WT ensemble (3-D ribbon) |
| `WT_fig2_ensemble_overlay.png` | All 826 WT conformations superimposed, blue→red by frame index |
| `WT_fig3_rmsf_colored.png` | Mean WT structure colored by per-residue RMSF + RMSF profile |
| `WT_fig4_pymol_render.png` | Ray-traced PyMOL render of WT, cartoon colored blue→red by RMSF |
| `C8S_fig1_average_structure.png` | Mean Cα structure of the C8S ensemble |
| `C8S_fig2_ensemble_overlay.png` | All 837 C8S conformations superimposed |
| `C8S_fig3_rmsf_colored.png` | Mean C8S structure colored by per-residue RMSF + RMSF profile |
| `C8S_fig4_pymol_render.png` | Ray-traced PyMOL render of C8S colored by RMSF |
| `fig4_publication_render.png` | WT and C8S side-by-side composite with shared RMSF colorbar |

All figures are saved at 300 dpi.

### `structures/`

| File | Description |
|---|---|
| `WT_average_structure.pdb` | Mean Cα coordinates of the WT ensemble |
| `WT_average_rmsf.pdb` | Mean Cα structure; B-factor = per-residue RMSF (Å) |
| `WT_average_rmsf_allatom.pdb` | Mean heavy-atom structure; B-factor = per-residue RMSF (Å) |
| `C8S_average_structure.pdb` | Mean Cα coordinates of the C8S ensemble |
| `C8S_average_rmsf.pdb` | Mean Cα structure; B-factor = per-residue RMSF (Å) |
| `C8S_average_rmsf_allatom.pdb` | Mean heavy-atom structure; B-factor = per-residue RMSF (Å) |

The `*_allatom.pdb` files are used by `pymol_render.py` for cartoon rendering;
the Cα-only PDBs are used by `visualize_ensemble.py` for the matplotlib figures
and can be opened in PyMOL/ChimeraX directly (`spectrum b` to reproduce the
RMSF coloring).

---

## Methods summary

1. Trajectories loaded with MDTraj 1.11.0.
2. All frames aligned to frame 0 by Cα RMSD (`mdtraj.Trajectory.superpose`).
3. Per-residue RMSF computed against the mean structure (`mdtraj.rmsf`).
4. RMSF values stored in the B-factor column of the average-structure PDBs.
5. PyMOL 3.1.0 headless (`pymol -cq`) used for ray-traced cartoon renders
   with a shared blue–white–red color scale across WT and C8S.

---

## Citation

If you use this analysis or the deposited trajectories, please cite:

> *Citation pending* — manuscript in preparation, Moschou Lab, IMBB-FORTH.

BioEmu:
> Wayment-Steele H.K. et al. (2024). *Predicting multiple conformations via
> sequence clustering and AlphaFold2*. **Nature** 625, 832–839.
> https://doi.org/10.1038/s41586-023-06832-9

---

## License

MIT License — see [LICENSE](LICENSE) for details.
