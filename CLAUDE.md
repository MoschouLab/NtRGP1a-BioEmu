# NtRGP1a BioEmu Conformational Analysis

## Project context
Visualisation of BioEmu conformational ensembles for NtRGP1a 
(Nicotiana tabacum Glycine-rich RNA-binding protein, 156 aa, UniProt D6PZY5)
and its C8S point mutant. Part of the PLANTEX/ERC project at IMBB-FORTH.

## Environment
- conda env: bioemu
- Working directory: /home/moschou/Documents/bioemu_analysis/NtRGP1a
- Python packages available: mdtraj, mdanalysis, matplotlib, numpy, pymol (3.1.0 Open-Source)

## Input data
- output_2000/topology.pdb + samples.xtc  → WT  (826 conformations)
- output_2000_C8S/topology.pdb + samples.xtc → C8S (837 conformations)

## Goal: generate 3 figures for each variant (WT and C8S)

### Figure 1 — Average structure
- Compute mean Cα coordinates across all conformations
- Save as average_structure.pdb
- Render as image (use matplotlib 3D or py3Dmol)

### Figure 2 — Ensemble overlay (spaghetti plot)
- Align all conformations to the first frame (Cα RMSD alignment)
- Plot all Cα traces superimposed
- Use transparency and color gradient (blue→red by frame index)

### Figure 3 — Average structure colored by per-residue RMSF
- Calculate per-residue RMSF across the ensemble
- Map RMSF values onto the average structure
- Color scale: blue (rigid) → red (flexible)
- This is the main deliverable figure

### Figure 4 — PyMOL publication render
- Load *_average_rmsf.pdb for both WT and C8S
- Color by B-factor (blue=rigid, red=flexible) using spectrum
- Ray-trace and export as PNG 300 dpi
- Use PyMOL headless mode (pymol -cq script.py)
- PyMOL 3.1.0 Open-Source is installed in the bioemu conda env — no fallback needed

## Output
All figures saved to ./figures/ as high-resolution PNG (300 dpi)
All intermediate PDBs saved to ./structures/

## Code style
- Single self-contained Python scripts
- Use MDTraj as primary library (already validated in this project)
- Always activate conda env bioemu before running