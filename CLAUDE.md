# CLAUDE.md — invariants for the BsCM cluster-model project

**Read this file first, every session.** This is a standalone repository (a cluster
model of chorismate mutase in ORCA 6.x). It is *not* related to OrcaStudio. Full task
spec lives in `protocol.md` (Ukrainian). Structural ground truth lives in `notes.md`
and `scripts/active_site.json`.

## Goal
Reproduce the methodology of **Burschowsky et al., FEBS Open Bio 7 (2017)** in ORCA 6.x
(the original used Gaussian09). Compute the activation free energy ΔG‡ for three systems:
`water` (uncatalysed), `WT` (wild-type, Arg90), `Arg90Cit` (neutral citrulline mutant).
Key deliverables:
- `ΔΔG‡_cat = ΔG‡(water) − ΔG‡(WT) > 0` (enzyme lowers the barrier)
- `ΔΔG‡(Cit − WT) = ΔG‡(Arg90Cit) − ΔG‡(WT) > 0` (removing one positive charge raises
  the barrier — the "money" result showing catalysis here is electrostatic).

## Charges — RE-VERIFY against the actual cluster EVERY time
| species | charge |
|---|---|
| chorismate / prephenate (ligand) | −2 |
| Arg7 | +1 |
| Arg63 | +1 |
| Glu78 | −1 |
| Cit90 (citrulline) | 0 |
| Arg90 (WT) | +1 |
| Thr74, Cys75, Tyr108, H₂O | 0 |

Totals: **Cit90 cluster = −1**, **WT cluster = 0**. **Multiplicity = 1 everywhere**
(all closed-shell). Confirmed against the 3ZP7 reference active site (see `notes.md`).

## Freezing (Burschowsky §3.3)
Capping H atoms **and the atom they attach to** are fixed at crystallographic
coordinates. The freeze list is a set of **0-based atom indices** and MUST be
**identical across R / P / TS / WT / Cit**.

## Atom ordering
Identical atom ordering across R / P / WT / Cit is **mandatory** — otherwise NEB and
the Cartesian constraints break. Never reorder atoms between states.

## ORCA conventions
- Invoke ORCA by its **full absolute path**: `/opt/orca/orca` (verified present).
- **One isolated directory per calculation** under `calc/` (create later, one per run).
- Use **B3LYP/G** (the Gaussian variant, VWN3 correlation) for reproducibility vs the
  paper. Reference method budget: B3LYP/6-31G(d), ~100–170 atoms.
- Logs: stream/tail only, never load whole `.out` files.
- Completion criterion: process `.exit_code == 0` **AND** `.out` contains
  `ORCA TERMINATED NORMALLY`.

## Structure ground truth (PDB 3ZP7, established session 1)
- 3ZP7 = Arg90Cit BsCM + chorismate + prephenate, 1.70 Å, space group P1, 6 protein
  chains A–F (two homotrimers in the AU). Ligand codes: **ISJ = chorismate (substrate)**,
  **PRE = prephenate (product)**, **CIR = citrulline at residue 90**, HOH = water.
- Reference active site: **ISJ chain A / resnum 1119** (substrate copy).
- The active site is at a **subunit interface**: home chain A donates Arg63, Thr74,
  Cys75; partner chain B donates Arg7, Glu78, Cit90, Tyr108. Keep this A/B pairing
  consistent when extracting the cluster.
- Reactive ether oxygen of chorismate is **O11** (bonded to ring C1 and enolpyruvyl
  C12; the C–O bond that breaks). Carboxylates: O8/O9 (on C7) and O14/O15 (on C13).
  Hydroxyl: O10 (on C6). Structural water: **HOH A/2119** (bridges both ligand
  carboxylates — do NOT discard).

## Built cluster (session 2) — Cit reactant
- `structures/clusters/cit_react.xyz` — 114 atoms, charge **−1**, mult 1 (422 e⁻).
- Atom order is **LAW**: `[ISJ 0–23][H₂O 24–26][Arg7 27–42][Arg63 43–58][Thr74 59–67]
  [Cys75 68–73][Glu78 74–83][Tyr108 84–99][Cit90 100–113]`. Full map in
  `scripts/cit_react_order.json`. Cit90 is the trailing block so WT swaps only it.
- Freeze list: `scripts/constraints_wt_cit.txt` (0-based; cap-H + anchor C per trimmed
  residue; substrate + water free). 14 frozen atoms.
- **Reactive atoms** (for NEB/TS later): bond breaking **C1–O11** (ether), bond forming
  **C3–C16** (Claisen chair TS).
- Rebuild: `.venv/bin/python scripts/build_cluster.py`. Product = edit reactant geometry
  (same atoms/order); WT = regenerate the Cit90 block as ethylguanidinium (+1).

## Languages
- Scripts and commit messages: **English**.
- Chemistry notes (`notes.md`): **Ukrainian**.

## Environment
- Python venv at `.venv/` (gemmi installed). Run scripts with `.venv/bin/python`.
- Working dir / repo root is this directory.

## Discipline
Do not cut, cap, or run ORCA until the structural selection is reviewed and approved.
This session (1) was reconnaissance only.
