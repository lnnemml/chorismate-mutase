#!/usr/bin/env python3
"""
recon.py -- Structural reconnaissance of PDB 3ZP7 (Arg90Cit BsCM + chorismate + prephenate).

Read-only. Establishes the structural "ground truth" for the cluster model:
  a) chains and residue counts
  b) all HETATM residues; ligand identification; residue-90 representation
  c) pick ONE reference active site around a SUBSTRATE (chorismate / ISJ) copy
  d) protein residues within 5 A of that ligand; localize the 7 catalytic residues,
     flag those donated by a neighbouring subunit; warn on any missing
  e) crystallographic waters near the ligand carboxylates / ether oxygen
  f) H-bond table: ligand carboxylate + ether O -> polar atoms of the 7 residues

Writes scripts/active_site.json (machine-readable selection for the next session).
Does NOT cut, cap, or run ORCA.

Ligand chemistry (PDB chemical component dictionary):
  ISJ = chorismate  (SUBSTRATE)   -- 16 heavy atoms
  PRE = prephenate  (PRODUCT)     -- 16 heavy atoms
  CIR = citrulline  (residue 90, isosteric neutral replacement of catalytic Arg90)
  HOH = water
"""
import json
import os
import sys
from collections import defaultdict

import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CIF = os.path.join(REPO, "structures", "raw", "3zp7.cif")
JSON_OUT = os.path.join(HERE, "active_site.json")

CONTACT_CUTOFF = 5.0   # A, residue-selection shell
HB_CUTOFF = 3.6        # A, heavy-atom donor/acceptor distance for H-bond table
WATER_CUTOFF = 3.5     # A, water-to-ligand-O contact

# The 7 catalytic residues expected around the substrate (Burschowsky et al. 2017).
CATALYTIC = [("ARG", 7), ("ARG", 63), ("THR", 74), ("CYS", 75),
             ("GLU", 78), ("CIR", 90), ("TYR", 108)]

# Polar / charged atoms per catalytic residue used for the H-bond table.
POLAR_ATOMS = {
    "ARG": ["NE", "NH1", "NH2"],           # guanidinium
    "THR": ["OG1"],                        # hydroxyl
    "CYS": ["SG"],                         # thiol
    "GLU": ["OE1", "OE2"],                 # carboxylate
    "CIR": ["N6", "O7", "N8"],             # ureido group (delta-N, C=O, terminal NH2)
    "TYR": ["OH"],                         # phenol hydroxyl
}

# ISJ (chorismate) oxygen roles, VERIFIED by bond connectivity in this structure
# (see scripts/recon.py output / notes.md):
#   ring carboxylate      : O8, O9   (both on C7)
#   enolpyruvyl carboxylate: O14, O15 (both on C13)
#   enol-ether bridging O : O11  (bonded to ring C1 AND enolpyruvyl C12) -- the
#                                 C-O bond that breaks in the Claisen rearrangement
#   secondary hydroxyl    : O10  (on ring C6)
ISJ_CARBOXYL = ["O8", "O9", "O14", "O15"]
ISJ_ETHER = ["O11"]
ISJ_HYDROXYL = ["O10"]
# carboxylate carbon -> its two oxygens (used for intra-ligand charge-bridge test)
ISJ_COO_GROUPS = {"C7": ["O8", "O9"], "C13": ["O14", "O15"]}
CRYO_CODES = {"GOL", "EDO", "PEG", "MPD", "SO4", "PO4", "CL", "NA", "K",
              "MG", "CA", "ZN", "ACT", "DMS", "FMT", "NO3"}


def banner(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


def dist(a, b):
    return a.pos.dist(b.pos)


STD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
          "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
          "TYR", "VAL", "MSE", "SEC", "PYL"}


def is_aa(res):
    """Standard/modified amino-acid test that works across gemmi versions."""
    if res.name in STD_AA:
        return True
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info) and info.is_amino_acid()


def main():
    if not os.path.exists(CIF):
        sys.exit(f"CIF not found: {CIF}")
    st = gemmi.read_structure(CIF)
    st.setup_entities()
    model = st[0]

    banner("3ZP7 reconnaissance")
    print(f"file            : {CIF}")
    print(f"structure name  : {st.name}")
    print(f"space group     : {st.spacegroup_hm}")
    print(f"resolution      : {st.resolution} A")
    print(f"models          : {len(st)}")

    # ---- (a) chains and residue counts -------------------------------------
    banner("(a) Chains and residue counts")
    for ch in model:
        poly = sum(1 for r in ch if is_aa(r))
        het = sum(1 for r in ch if r.het_flag == "H" and r.name != "HOH")
        wat = sum(1 for r in ch if r.name == "HOH")
        print(f"  chain {ch.name}: {len(ch):4d} residues "
              f"(amino acids={poly}, het(non-water)={het}, water={wat})")

    # ---- (b) HETATM inventory ----------------------------------------------
    banner("(b) HETATM residues (non-water)")
    lig_index = defaultdict(list)   # code -> list of (chain, seqid, natoms)
    for ch in model:
        for res in ch:
            if res.het_flag == "H" and res.name != "HOH":
                lig_index[res.name].append((ch.name, res.seqid.num, len(res)))
    names = {"ISJ": "chorismate (SUBSTRATE)",
             "PRE": "prephenate (PRODUCT)",
             "CIR": "citrulline (residue 90; neutral isostere of Arg90)"}
    for code in sorted(lig_index):
        role = names.get(code, "cryoprotectant/ion -> EXCLUDE"
                         if code in CRYO_CODES else "unknown ligand")
        insts = lig_index[code]
        print(f"  {code}  [{role}]  copies={len(insts)}")
        for cn, sn, na in insts:
            print(f"       chain {cn}  resnum {sn}  atoms {na}")
    nwater = sum(1 for ch in model for r in ch if r.name == "HOH")
    print(f"  HOH  [water]  copies={nwater}")

    # residue-90 representation
    banner("(b') Residue 90 representation")
    for ch in model:
        for res in ch:
            if res.seqid.num == 90 and res.name in ("CIR", "ARG"):
                print(f"  chain {ch.name}: residue 90 = {res.name} "
                      f"(het_flag={res.het_flag!r}), atoms: "
                      f"{' '.join(a.name for a in res)}")
                break
    print("  -> Residue 90 is modelled as HET residue CIR (citrulline), a distinct")
    print("     3-letter code, NOT a modified ARG. Neutral ureido side chain.")

    # ---- (c) choose reference active site ----------------------------------
    banner("(c) Reference active-site selection")
    # Prefer a chain that has a clean single ISJ (no alt-conf duplicates) and a
    # single PRE, choosing the lexicographically first such chain.
    isj_by_chain = defaultdict(list)
    for ch in model:
        for res in ch:
            if res.name == "ISJ":
                isj_by_chain[ch.name].append(res)
    clean = sorted(cn for cn, lst in isj_by_chain.items() if len(lst) == 1)
    if not clean:
        sys.exit("No chain with a single clean ISJ copy found.")
    ref_chain = clean[0]
    ref_isj = isj_by_chain[ref_chain][0]
    print(f"  chains with a single clean ISJ copy : {clean}")
    print(f"  CHOSEN reference ligand : ISJ (chorismate/substrate) "
          f"chain {ref_chain} resnum {ref_isj.seqid.num}")
    print("  Rationale: substrate copy (per task); clean single occupancy; "
          "first such chain alphabetically.")

    ref_atoms = list(ref_isj)

    # ---- (d) protein residues within 5 A -----------------------------------
    banner(f"(d) Protein residues within {CONTACT_CUTOFF} A of the reference ISJ")
    # min distance from each protein residue to any ISJ atom
    contacts = {}   # (chain, resname, seqid) -> min_dist
    for ch in model:
        for res in ch:
            if not is_aa(res) and res.name != "CIR":
                continue
            dmin = min((dist(a, la) for a in res for la in ref_atoms),
                       default=1e9)
            if dmin <= CONTACT_CUTOFF:
                contacts[(ch.name, res.name, res.seqid.num)] = dmin

    for key in sorted(contacts, key=lambda k: contacts[k]):
        cn, rn, sn = key
        tag = "  <-- other subunit" if cn != ref_chain else ""
        print(f"  {rn:3s} {sn:<4d} chain {cn}   dmin={contacts[key]:.2f} A{tag}")

    banner("(d') Localization of the 7 catalytic residues")
    found = {}
    for (resn, seqn) in CATALYTIC:
        hits = [(cn, rn, sn, contacts[(cn, rn, sn)])
                for (cn, rn, sn) in contacts if rn == resn and sn == seqn]
        if hits:
            hits.sort(key=lambda x: x[3])
            cn, rn, sn, d = hits[0]
            same = (cn == ref_chain)
            found[f"{resn}{seqn}"] = {"resname": rn, "chain": cn, "resnum": sn,
                                      "min_dist": round(d, 2),
                                      "same_subunit": same}
            tag = "same subunit" if same else "OTHER SUBUNIT (cross-interface)"
            print(f"  {resn}{seqn:<4d}: chain {cn} resnum {sn}  "
                  f"dmin={d:.2f} A  [{tag}]")
        else:
            found[f"{resn}{seqn}"] = None
            print(f"  {resn}{seqn:<4d}: *** NOT within {CONTACT_CUTOFF} A "
                  f"of this ligand copy -- WARNING ***")
    missing = [k for k, v in found.items() if v is None]
    if missing:
        print(f"\n  WARNING: catalytic residues not in contact: {missing}")
    else:
        print("\n  All 7 catalytic residues located within the contact shell.")

    # ---- (e) crystallographic waters ---------------------------------------
    banner(f"(e) Crystallographic waters within {WATER_CUTOFF} A of ligand "
           f"carboxylate / ether O")
    key_o = {a.name: a for a in ref_atoms
             if a.name in ISJ_CARBOXYL + ISJ_ETHER + ISJ_HYDROXYL}
    water_cands = []
    for ch in model:
        for res in ch:
            if res.name != "HOH":
                continue
            ow = res[0]
            near = [(oname, dist(ow, oat))
                    for oname, oat in key_o.items()
                    if dist(ow, oat) <= WATER_CUTOFF]
            if near:
                near.sort(key=lambda x: x[1])
                water_cands.append((ch.name, res.seqid.num, near))
    if not water_cands:
        print("  (none within cutoff)")
    for cn, sn, near in sorted(water_cands, key=lambda x: x[2][0][1]):
        pretty = ", ".join(f"{on}={d:.2f}" for on, d in near)
        bridge = ""
        # a water contacting a ligand carboxylate AND being near Glu78 bridges
        # two negative charges -> flag as candidate "structural" water
        cnames = [on for on, _ in near]
        if any(o in ISJ_CARBOXYL for o in cnames):
            bridge = "  [near ligand carboxylate -> check charge-bridge]"
        print(f"  HOH {sn:<5d} chain {cn}: {pretty}{bridge}")

    # A "structural" water separates negative charges. Test each candidate for
    # simultaneous contact with (i) BOTH ligand carboxylate groups, and/or
    # (ii) a ligand carboxylate AND the Glu78 carboxylate.
    banner("(e') Structural-water test (water bridging two negative charges)")
    fv = found.get("GLU78")
    glu_o = []
    if fv:
        for ch in model:
            if ch.name == fv["chain"]:
                for res in ch:
                    if res.name == "GLU" and res.seqid.num == fv["resnum"]:
                        glu_o = [a for a in res if a.name in ("OE1", "OE2")]
    struct_waters = []
    for cn, sn, near in water_cands:
        wr = None
        for ch in model:
            if ch.name == cn:
                for res in ch:
                    if res.name == "HOH" and res.seqid.num == sn:
                        wr = res
        if wr is None:
            continue
        ow = wr[0]
        # which distinct ligand carboxylate GROUPS does it contact?
        grp_hits = {}
        for cc, oxys in ISJ_COO_GROUPS.items():
            dd = min((d for o, d in near if o in oxys), default=1e9)
            if dd <= WATER_CUTOFF:
                grp_hits[cc] = round(dd, 2)
        dglu = min((dist(ow, go) for go in glu_o), default=1e9)
        bridges_ligand_coos = len(grp_hits) >= 2
        bridges_glu = (len(grp_hits) >= 1 and dglu <= WATER_CUTOFF)
        if bridges_ligand_coos or bridges_glu:
            entry = {"chain": cn, "resnum": sn,
                     "ligand_COO_groups": grp_hits,
                     "d_glu78_carboxyl": (round(dglu, 2) if dglu < 1e8 else None),
                     "bridges_two_ligand_carboxylates": bridges_ligand_coos,
                     "bridges_ligand_and_glu78": bridges_glu}
            struct_waters.append(entry)
            msg = []
            if bridges_ligand_coos:
                msg.append("bridges BOTH ligand carboxylates "
                           f"({grp_hits})")
            if bridges_glu:
                msg.append(f"also near Glu78-COO ({dglu:.2f} A)")
            print(f"  HOH {sn} chain {cn}: " + "; ".join(msg))
            print("      -> separates two negative charges; KEEP (structural water).")
    if not struct_waters:
        print("  (no water simultaneously bridges two negative charges within cutoff)")

    # ---- (f) H-bond table --------------------------------------------------
    banner(f"(f) H-bond distances (<= {HB_CUTOFF} A): ligand O -> catalytic polar atoms")
    print(f"  {'ligand_O':10s} {'role':16s} {'residue':10s} {'atom':6s} {'d(A)':>6s}  subunit")
    print("  " + "-" * 66)
    hbonds = []
    lig_polar = [a for a in ref_atoms
                 if a.name in ISJ_CARBOXYL + ISJ_ETHER + ISJ_HYDROXYL]
    for key, info in found.items():
        if info is None:
            continue
        cn, rn, sn = info["chain"], info["resname"], info["resnum"]
        res = None
        for ch in model:
            if ch.name == cn:
                for r in ch:
                    if r.name == rn and r.seqid.num == sn:
                        res = r
        if res is None:
            continue
        for pa_name in POLAR_ATOMS.get(rn, []):
            pa = next((a for a in res if a.name == pa_name), None)
            if pa is None:
                continue
            for lo in lig_polar:
                d = dist(lo, pa)
                if d <= HB_CUTOFF:
                    role = ("carboxylate" if lo.name in ISJ_CARBOXYL
                            else "ether-O" if lo.name in ISJ_ETHER
                            else "hydroxyl")
                    sub = "same" if cn == ref_chain else "OTHER"
                    hbonds.append({"ligand_atom": lo.name, "role": role,
                                   "residue": f"{rn}{sn}", "res_atom": pa_name,
                                   "chain": cn, "dist": round(d, 2),
                                   "subunit": sub})
    for hb in sorted(hbonds, key=lambda h: h["dist"]):
        print(f"  {hb['ligand_atom']:10s} {hb['role']:16s} "
              f"{hb['residue']:10s} {hb['res_atom']:6s} {hb['dist']:6.2f}  "
              f"{hb['subunit']}")
    if not hbonds:
        print("  (no contacts within cutoff)")

    # ---- charge tally (informational) --------------------------------------
    banner("(g) Preliminary cluster charge tally (per actual residues found)")
    charge_map = {"ISJ": -2, "PRE": -2, "ARG": +1, "GLU": -1, "CIR": 0,
                  "THR": 0, "CYS": 0, "TYR": 0}
    total = charge_map["ISJ"]
    print(f"  ligand ISJ (chorismate)         : {charge_map['ISJ']:+d}")
    for (resn, seqn) in CATALYTIC:
        info = found.get(f"{resn}{seqn}")
        if info is None:
            print(f"  {resn}{seqn:<4d} : (missing)")
            continue
        q = charge_map.get(resn, 0)
        total += q
        note = " (WT would be ARG90 = +1)" if resn == "CIR" else ""
        print(f"  {resn}{seqn:<4d} chain {info['chain']:1s}          : {q:+d}{note}")
    print(f"  ----")
    print(f"  Cit90 cluster total (as-is)     : {total:+d}   (expected -1)")
    print(f"  WT cluster (swap Cit90->Arg90)  : {total + 1:+d}   (expected  0)")
    print("  Multiplicity = 1 (all closed-shell).")

    # ---- write JSON --------------------------------------------------------
    out = {
        "pdb": "3ZP7",
        "cif": os.path.relpath(CIF, REPO),
        "space_group": st.spacegroup_hm,
        "resolution_A": st.resolution,
        "reference_ligand": {
            "code": "ISJ", "name": "chorismate", "role": "substrate",
            "chain": ref_chain, "resnum": ref_isj.seqid.num,
            "carboxylate_O": ISJ_CARBOXYL, "ether_O": ISJ_ETHER,
            "hydroxyl_O": ISJ_HYDROXYL,
        },
        "product_ligand_note": "PRE = prephenate (product); not used as reference",
        "catalytic_residues": found,
        "contact_cutoff_A": CONTACT_CUTOFF,
        "structural_water_candidates": struct_waters,
        "all_water_contacts": [
            {"chain": cn, "resnum": sn,
             "contacts": [{"ligand_O": on, "dist": round(d, 2)} for on, d in near]}
            for cn, sn, near in sorted(water_cands, key=lambda x: x[2][0][1])
        ],
        "hbonds": sorted(hbonds, key=lambda h: h["dist"]),
        "charges": {
            "map": charge_map,
            "cit90_cluster_total": total,
            "wt_cluster_total": total + 1,
            "multiplicity": 1,
        },
        "notes": [
            "Read-only reconnaissance; no atoms cut/capped; ORCA not run.",
            "Atom ordering / freeze-list to be fixed identically across R/P/TS/WT/Cit "
            "in a later session.",
        ],
    }
    with open(JSON_OUT, "w") as fh:
        json.dump(out, fh, indent=2)
    banner("Wrote machine-readable selection")
    print(f"  {JSON_OUT}")


if __name__ == "__main__":
    main()
