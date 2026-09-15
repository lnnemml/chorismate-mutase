#!/usr/bin/env python3
"""
build_cluster.py -- Build ONE cluster: Cit-variant, REACTANT (chorismate substrate).

Trimming + capping + freeze-list only. Does NOT build product / WT, does NOT run
xTB/ORCA. Source of truth: scripts/active_site.json (recon session).

Selection (all from 3ZP7, reference active site):
  substrate : ISJ chain A/1119 (chorismate dianion, -2)   -- kept whole, H from CCD ideal
  water     : HOH chain A/2119 (structural, bridges both ligand carboxylates) -- whole
  Arg7  (B) -> ethylguanidinium (+1)
  Arg63 (A) -> ethylguanidinium (+1)
  Thr74 (A) -> ethanol-type alcohol (0)
  Cys75 (A) -> methanethiol (0)
  Glu78 (B) -> propionate (-1)
  Tyr108(B) -> p-cresol (0)
  Cit90 (B) -> neutral ethyl-urea (0)

Net charge = -1, closed shell (mult 1).

Capping: each residue is cut at ONE C-C bond; the removed heavy atom is replaced by a
capping H placed along the original bond vector at 1.09 A from the anchor carbon.
Freeze: the capping H AND its anchor carbon (Burschowsky Sec. 3.3). Substrate + water
are fully free.

ATOM ORDER (this file defines the law reused one-to-one by product & WT):
  [ISJ substrate] [HOH water] [Arg7] [Arg63] [Thr74] [Cys75] [Glu78] [Tyr108] [Cit90]
Residue 90 is placed LAST so the WT build swaps only that trailing block, leaving all
other indices (and their constraints) identical.

Outputs:
  structures/clusters/cit_react.xyz
  scripts/constraints_wt_cit.txt   (0-based frozen indices, one per line)
  scripts/cit_react_order.json     (global_idx -> residue / CCD atom name / element)
"""
import json
import os
import sys

import numpy as np
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CIF = os.path.join(REPO, "structures", "raw", "3zp7.cif")
ISJ_CCD = os.path.join(REPO, "structures", "raw", "ISJ.cif")
XYZ_OUT = os.path.join(REPO, "structures", "clusters", "cit_react.xyz")
CONS_OUT = os.path.join(HERE, "constraints_wt_cit.txt")
ORDER_OUT = os.path.join(HERE, "cit_react_order.json")

CH = 1.09  # A, C-H capping bond length

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def get_residue(model, chain, resname, resnum):
    for ch in model:
        if ch.name != chain:
            continue
        for res in ch:
            if res.name == resname and res.seqid.num == resnum:
                return res
    raise KeyError(f"{resname} {chain}/{resnum} not found")


def apos(res, name):
    for a in res:
        if a.name == name:
            return np.array([a.pos.x, a.pos.y, a.pos.z])
    raise KeyError(f"atom {name} not in {res.name}")


def cap_h(anchor_xyz, removed_xyz, length=CH):
    v = removed_xyz - anchor_xyz
    v = v / np.linalg.norm(v)
    return anchor_xyz + length * v


def elem_of(name):
    n = name.lstrip("0123456789")
    if n[:1] in ("C", "N", "O", "S", "H"):
        return n[:1]
    raise ValueError(f"cannot infer element from {name}")


# ---------------------------------------------------------------------------
# ISJ substrate: crystal heavy-atom coords + CCD topology template; H placed by
# ideal LOCAL geometry (RDKit AddHs). A rigid whole-molecule CCD->crystal fit is
# NOT used -- chorismate's flexible enolpyruvyl tail would distort transferred H.
# Dianion (-2): both carboxylates deprotonated (no HO9/HO15); ring O10-H kept.
# ---------------------------------------------------------------------------
ISJ_HEAVY = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "O8", "O9", "O10",
             "O11", "C12", "C13", "O14", "O15", "C16"]
ISJ_BONDS = [("C1", "C2", "S"), ("C1", "C6", "S"), ("C1", "O11", "S"),
             ("C2", "C3", "D"), ("C3", "C4", "S"), ("C3", "C7", "S"),
             ("C4", "C5", "D"), ("C5", "C6", "S"), ("C6", "O10", "S"),
             ("C7", "O8", "D"), ("C7", "O9", "S"),
             ("O11", "C12", "S"), ("C12", "C13", "S"), ("C12", "C16", "D"),
             ("C13", "O14", "D"), ("C13", "O15", "S")]
ISJ_CHARGES = {"O9": -1, "O15": -1}  # deprotonated carboxylates -> dianion


def build_isj(model):
    res = get_residue(model, "A", "ISJ", 1119)
    spec = {"keep": ISJ_HEAVY, "bonds": ISJ_BONDS, "charges": ISJ_CHARGES}
    atoms, qtot = build_group(res, spec, cap=None)
    if qtot != -2:
        sys.exit(f"ISJ charge {qtot} != -2")
    return atoms, res


# ---------------------------------------------------------------------------
# generic RDKit builder: heavy atoms (crystal coords) + CCD-template bonds +
# formal charges; optional cap H along the cut bond; then AddHs(addCoords).
# ---------------------------------------------------------------------------
def build_group(res, spec, cap):
    """cap: None (no capping, e.g. substrate) or dict {anchor, removed}.
    Returns (atoms, qtot) if cap is None,
    else (atoms, qtot, anchor_idx, cap_idx)."""
    rw = Chem.RWMol()
    idx = {}
    conf_xyz = {}
    for name in spec["keep"]:
        a = Chem.Atom(elem_of(name))
        a.SetProp("ccd", name)
        i = rw.AddAtom(a)
        idx[name] = i
        conf_xyz[i] = apos(res, name)
    cap_lbl = None
    if cap is not None:
        anc = cap["anchor"]
        cxyz = cap_h(apos(res, anc), apos(res, cap["removed"]))
        h = Chem.Atom("H")
        cap_lbl = "Hcap"
        h.SetProp("ccd", cap_lbl)
        hi = rw.AddAtom(h)
        conf_xyz[hi] = cxyz
        rw.AddBond(idx[anc], hi, Chem.BondType.SINGLE)
    bt = {"S": Chem.BondType.SINGLE, "D": Chem.BondType.DOUBLE,
          "A": Chem.BondType.AROMATIC}
    for a, b_, o in spec["bonds"]:
        rw.AddBond(idx[a], idx[b_], bt[o])
    for name, q in spec.get("charges", {}).items():
        rw.GetAtomWithIdx(idx[name]).SetFormalCharge(q)
    m = rw.GetMol()
    conf = Chem.Conformer(m.GetNumAtoms())
    for i in range(m.GetNumAtoms()):
        x, y, z = conf_xyz[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    m.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(m)
    m = Chem.AddHs(m, addCoords=True)
    conf = m.GetConformer()
    out = []
    for atom in m.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        xyz = np.array([p.x, p.y, p.z])
        if atom.HasProp("ccd"):
            label = atom.GetProp("ccd")
        else:  # AddHs-generated H: label by its heavy parent
            par = atom.GetNeighbors()[0]
            plab = par.GetProp("ccd") if par.HasProp("ccd") else par.GetSymbol()
            label = f"H@{plab}"
        out.append((atom.GetSymbol(), xyz, label))
    qtot = Chem.GetFormalCharge(m)
    if cap is None:
        return out, qtot
    anchor_idx = next(k for k, (_, _, l) in enumerate(out) if l == cap["anchor"])
    cap_idx = next(k for k, (_, _, l) in enumerate(out) if l == cap_lbl)
    return out, qtot, anchor_idx, cap_idx


# ---------------------------------------------------------------------------
# water: O from crystal, 2 H oriented toward the two ligand carboxylate O it bridges
# ---------------------------------------------------------------------------
def build_water(model, isj_res):
    w = get_residue(model, "A", "HOH", 2119)
    Ow = np.array([w[0].pos.x, w[0].pos.y, w[0].pos.z])
    acc = [apos(isj_res, "O15"), apos(isj_res, "O9")]  # its bridging partners
    atoms = [("O", Ow, "O")]
    for k, a in enumerate(acc):
        d = a - Ow
        d = d / np.linalg.norm(d)
        atoms.append(("H", Ow + 0.96 * d, f"H{k+1}"))
    return atoms


# ---------------------------------------------------------------------------
# fragment specifications
# ---------------------------------------------------------------------------
ARG = lambda: {
    "keep": ["CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "anchor": "CG", "removed": "CB",
    "bonds": [("CG", "CD", "S"), ("CD", "NE", "S"), ("NE", "CZ", "S"),
              ("CZ", "NH1", "D"), ("CZ", "NH2", "S")],
    "charges": {"NH1": +1},
}
GLU = {
    "keep": ["CB", "CG", "CD", "OE1", "OE2"],
    "anchor": "CB", "removed": "CA",
    "bonds": [("CB", "CG", "S"), ("CG", "CD", "S"),
              ("CD", "OE1", "D"), ("CD", "OE2", "S")],
    "charges": {"OE2": -1},
}
TYR = {
    "keep": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "anchor": "CB", "removed": "CA",
    # explicit Kekule benzene: CG=CD1, CE1=CZ, CE2=CD2
    "bonds": [("CB", "CG", "S"),
              ("CG", "CD1", "D"), ("CD1", "CE1", "S"), ("CE1", "CZ", "D"),
              ("CZ", "CE2", "S"), ("CE2", "CD2", "D"), ("CD2", "CG", "S"),
              ("CZ", "OH", "S")],
}
CYS = {
    "keep": ["CB", "SG"],
    "anchor": "CB", "removed": "CA",
    "bonds": [("CB", "SG", "S")],
}
THR = {
    "keep": ["CB", "OG1", "CG2"],
    "anchor": "CB", "removed": "CA",
    "bonds": [("CB", "OG1", "S"), ("CB", "CG2", "S")],
}
CIT = {
    "keep": ["C4", "C5", "N6", "C7", "O7", "N8"],
    "anchor": "C4", "removed": "C3",
    "bonds": [("C4", "C5", "S"), ("C5", "N6", "S"), ("N6", "C7", "S"),
              ("C7", "O7", "D"), ("C7", "N8", "S")],
}

# (label, chain, resname, resnum, spec, expected_charge)
FRAGMENTS = [
    ("Arg7",  "B", "ARG", 7,   ARG(), +1),
    ("Arg63", "A", "ARG", 63,  ARG(), +1),
    ("Thr74", "A", "THR", 74,  THR,    0),
    ("Cys75", "A", "CYS", 75,  CYS,    0),
    ("Glu78", "B", "GLU", 78,  GLU,   -1),
    ("Tyr108","B", "TYR", 108, TYR,    0),
    ("Cit90", "B", "CIR", 90,  CIT,    0),
]


def main():
    st = gemmi.read_structure(CIF)
    model = st[0]

    records = []   # (global_idx, residue_label, atom_label, element, xyz)
    frozen = []    # 0-based global indices
    charge_tally = []
    gidx = 0

    # ---- substrate ISJ (free) ----
    isj_atoms, isj_res = build_isj(model)
    for e, xyz, lab in isj_atoms:
        records.append((gidx, "ISJ_substrate", lab, e, xyz)); gidx += 1
    charge_tally.append(("ISJ_substrate(chorismate)", -2))

    # ---- water (free) ----
    for e, xyz, lab in build_water(model, isj_res):
        records.append((gidx, "HOH_water", lab, e, xyz)); gidx += 1
    charge_tally.append(("HOH_water", 0))

    # ---- protein fragments ----
    for label, chain, rn, num, spec, qexp in FRAGMENTS:
        res = get_residue(model, chain, rn, num)
        cap = {"anchor": spec["anchor"], "removed": spec["removed"]}
        out, qtot, a_i, c_i = build_group(res, spec, cap)
        if qtot != qexp:
            sys.exit(f"CHARGE MISMATCH {label}: got {qtot}, expected {qexp}")
        base = gidx
        for e, xyz, lab in out:
            records.append((gidx, label, lab, e, xyz)); gidx += 1
        frozen.append(base + a_i)   # anchor carbon
        frozen.append(base + c_i)   # capping H
        charge_tally.append((f"{label}({rn}{num} chain {chain})", qtot))

    # ---- validation ----
    natoms = len(records)
    net_charge = sum(q for _, q in charge_tally)
    nelec = sum(gemmi.Element(e).atomic_number for _, _, _, e, _ in records) - net_charge
    assert net_charge == -1, f"net charge {net_charge} != -1"
    assert nelec % 2 == 0, f"odd electron count {nelec}"

    # ---- write xyz ----
    os.makedirs(os.path.dirname(XYZ_OUT), exist_ok=True)
    with open(XYZ_OUT, "w") as fh:
        fh.write(f"{natoms}\n")
        fh.write("Cit-variant reactant (chorismate); charge=-1 mult=1; "
                 "atom order=LAW (see cit_react_order.json)\n")
        for _, _, _, e, xyz in records:
            fh.write(f"{e:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}\n")

    # ---- write constraints ----
    frozen.sort()
    with open(CONS_OUT, "w") as fh:
        fh.write("# 0-based frozen atom indices for cit_react.xyz (Burschowsky 3.3):\n")
        fh.write("# capping H + its anchor carbon, per trimmed residue. "
                 "Substrate + water fully free.\n")
        for i in frozen:
            fh.write(f"{i}\n")

    # ---- write order map ----
    order_map = [{"global_idx": g, "residue": rl, "atom": al, "element": e,
                  "frozen": g in set(frozen)}
                 for g, rl, al, e, _ in records]
    with open(ORDER_OUT, "w") as fh:
        json.dump({"natoms": natoms, "net_charge": net_charge,
                   "multiplicity": 1, "n_electrons": nelec,
                   "frozen_indices": frozen,
                   "reactive_atoms": {"bond_breaking": "C1-O11",
                                      "bond_forming": "C3-C16"},
                   "atoms": order_map}, fh, indent=2)

    # ---- report ----
    print(f"total atoms          : {natoms}")
    print(f"net charge           : {net_charge}  (expected -1)")
    print(f"electrons            : {nelec}  ({'even/closed-shell' if nelec%2==0 else 'ODD'})")
    print(f"multiplicity         : 1")
    print(f"frozen atoms         : {len(frozen)}  indices={frozen}")
    print("\ncharge tally:")
    for name, q in charge_tally:
        print(f"  {name:32s} {q:+d}")
    print(f"  {'-'*40}")
    print(f"  {'NET':32s} {net_charge:+d}")
    print("\nper-residue atom ranges:")
    ranges = {}
    for g, rl, _, _, _ in records:
        ranges.setdefault(rl, [g, g])
        ranges[rl][1] = g
    for rl, (a, b) in ranges.items():
        print(f"  {rl:20s} idx {a:3d}..{b:3d}  ({b-a+1} atoms)")
    print(f"\nwrote: {XYZ_OUT}")
    print(f"wrote: {CONS_OUT}")
    print(f"wrote: {ORDER_OUT}")


if __name__ == "__main__":
    main()
