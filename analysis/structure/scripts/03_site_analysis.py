"""Environment of residue 486 (R-type numbering) in a human R-type PK crystal structure.
Usage: python -I scripts/03_site_analysis.py PDBID RESNAME_EXPECTED [OUTPREFIX]
Author numbering in 2VGB/2VGI equals UniProt P30613 (R-type) numbering (checked against _struct_ref_seq and the supplied FASTA).
Outputs JSON to data/site_<PDBID>.json and a readable text report to stdout."""
import sys, os, json, tempfile, itertools
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, freesasa
from Bio import SeqIO
from Bio.PDB import PDBIO, Select, NeighborSearch
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pkcommon import *

pid = sys.argv[1]; expected = sys.argv[2]; POS = int(sys.argv[3]) if len(sys.argv) > 3 else 486
OFF = int(sys.argv[4]) if len(sys.argv) > 4 else 0   # R-type number = author number + OFF
m = load(pid)
rep = {"pdb": pid, "position_author": POS, "rtype_offset": OFF}

# 1. sequence check against the supplied FASTA (human only)
if OFF == 0 and pid in ("2VGB", "2VGI", "8XFD", "8TBS"):
    wt = str(next(SeqIO.parse("../sequence/supplied/PKD_protein_WT.fasta", "fasta")).seq)
    mut = str(next(SeqIO.parse("../sequence/supplied/PKD_protein_mutant.fasta", "fasta")).seq)
    seqchk = {}
    for ch in m:
        mism = []
        n = 0
        for r in protein_residues(ch):
            i = r.id[1]
            if 1 <= i <= len(wt):
                n += 1
                aa = index_to_one(three_to_index(r.resname)) if r.resname not in ("SEP", "TPO") else "X"
                if aa != wt[i - 1]:
                    mism.append(f"{wt[i-1]}{i}{aa}")
        seqchk[ch.id] = {"n_compared": n, "differences_vs_supplied_WT": mism}
    rep["sequence_check"] = seqchk
    rep["supplied_WT_486"] = wt[485]; rep["supplied_mut_486"] = mut[485]

# 2. secondary structure from the deposited mmCIF annotation
d = MMCIF2Dict(f"raw/{pid}.cif")
def L(k):
    v = d.get(k, []); return v if isinstance(v, list) else [v]
ss = []
for cid, typ, ch1, b, ch2, e in zip(L("_struct_conf.id"), L("_struct_conf.conf_type_id"), L("_struct_conf.beg_auth_asym_id"),
                                     L("_struct_conf.beg_auth_seq_id"), L("_struct_conf.end_auth_asym_id"), L("_struct_conf.end_auth_seq_id")):
    if int(b) <= POS <= int(e):
        ss.append(f"{ch1}:{typ} {cid} {b}-{e}")
for sid, rid, ch1, b, e in zip(L("_struct_sheet_range.sheet_id"), L("_struct_sheet_range.id"), L("_struct_sheet_range.beg_auth_asym_id"),
                                L("_struct_sheet_range.beg_auth_seq_id"), L("_struct_sheet_range.end_auth_seq_id")):
    if int(b) <= POS <= int(e):
        ss.append(f"{ch1}:STRAND sheet {sid} {rid} {b}-{e}")
rep["deposited_secondary_structure_covering_pos"] = ss
# neighbouring helices
near = []
for cid, typ, ch1, b, e in zip(L("_struct_conf.id"), L("_struct_conf.conf_type_id"), L("_struct_conf.beg_auth_asym_id"),
                               L("_struct_conf.beg_auth_seq_id"), L("_struct_conf.end_auth_seq_id")):
    if ch1 == "A" and (POS - 15 <= int(e) and int(b) <= POS + 15):
        near.append(f"{typ} {cid} {b}-{e}")
rep["deposited_ss_elements_near_pos_chainA"] = near

# 3. SASA (freesasa, Lee-Richards default, ProtOr radii; protein atoms only, HETATM and waters excluded)
class ChainSel(Select):
    def __init__(self, ids): self.ids = ids
    def accept_chain(self, c): return c.id in self.ids
    def accept_residue(self, r): return r.id[0] == " " or r.resname in ("SEP", "TPO")
    def accept_atom(self, a): return a.element != "H"
tmpd = tempfile.mkdtemp()
io = PDBIO(); io.set_structure(m)
chains = [c.id for c in m]
sasa = {}
def run_sasa(ids):
    f = os.path.join(tmpd, "".join(ids) + ".pdb"); io.save(f, ChainSel(ids))
    s = freesasa.Structure(f); res = freesasa.calc(s)
    return res.residueAreas()
tet = run_sasa(chains)
for c in chains:
    mono = run_sasa([c])
    ra_m = mono[c][str(POS)]; ra_t = tet[c][str(POS)]
    sasa[c] = {"resname": ra_m.residueType,
               "monomer_total_A2": round(ra_m.total, 1), "monomer_rel_total": round(ra_m.relativeTotal, 3),
               "monomer_rel_sidechain": round(ra_m.relativeSideChain, 3),
               "tetramer_total_A2": round(ra_t.total, 1), "tetramer_rel_total": round(ra_t.relativeTotal, 3),
               "tetramer_rel_sidechain": round(ra_t.relativeSideChain, 3),
               "buried_by_other_chains_A2": round(ra_m.total - ra_t.total, 1)}
rep["sasa"] = sasa

# 4. interactions of the side chain of residue POS in each chain
inter = {}
for ch in m:
    if POS not in ch: continue
    r = ch[POS]
    donors = {"ARG": ARG_DONORS, "TRP": TRP_DONORS, "GLN": {"NE2", "OE1"}}.get(r.resname, set())
    inter[ch.id] = {"resname": r.resname,
                    "polar_contacts_le_3.5A": polar_contacts(r, m, donors, 3.5),
                    "polar_contacts_3.5_4.0A": [x for x in polar_contacts(r, m, donors, 4.0) if x[2] > 3.5],
                    "salt_bridges_le_4.0A": salt_bridges(r, m, 4.0) if r.resname == "ARG" else [],
                    "cation_pi_le_6.0A": cation_pi(r, m, 6.0) if r.resname == "ARG" else [],
                    "all_contacts_le_4.0A": all_contacts(r, m, 4.0),
                    "sidechain_contacts_le_4.0A": all_contacts(r, m, 4.0, sidechain_only=True),
                    "B_factors": {a.get_id(): round(a.get_bfactor(), 1) for a in r},
                    "occupancy_min": min(a.get_occupancy() for a in r),
                    "atoms": [a.get_id() for a in r]}
rep["interactions"] = inter

# 5. inter-chain interfaces: residues with any heavy atom within 4.0 A of another chain; classify by domain
atoms_by_chain = {c.id: [a for r in protein_residues(c) for a in r] for c in m}
pairs = {}
iface_res = {}
for c1, c2 in itertools.permutations(chains, 2):
    ns = NeighborSearch(atoms_by_chain[c2])
    res = set()
    for a in atoms_by_chain[c1]:
        if ns.search(a.coord, 4.0):
            res.add(a.get_parent().id[1])
    iface_res[(c1, c2)] = sorted(res)
for c1, c2 in itertools.combinations(chains, 2):
    rs = iface_res[(c1, c2)]
    from collections import Counter
    pairs[f"{c1}-{c2}"] = {"n_interface_residues_in_" + c1: len(rs), "domains": dict(Counter(domain_of(x, OFF) for x in rs))}
rep["chain_pair_interfaces"] = pairs
# label pair types
ptype = {}
for k, v in pairs.items():
    dom = v["domains"]
    if not dom: ptype[k] = "none"; continue
    top = max(dom, key=dom.get)
    ptype[k] = {"A": "A/A' (large interface, A domains)", "C": "C/C' (large interface, C domains)"}.get(top, f"minor ({top})") if v[list(v)[0]] >= 10 else f"minor contact ({dom})"
rep["chain_pair_type"] = ptype

# 6. distances from residue POS (any atom) to sites
dist = {}
for ch in m:
    if POS not in ch: continue
    r = ch[POS]
    dd = {}
    # contacts to other chains
    for c2 in chains:
        if c2 == ch.id: continue
        best = (1e9, None, None)
        for a in r:
            for b in atoms_by_chain[c2]:
                x = a - b
                if x < best[0]: best = (x, a, b)
        dd[f"min_to_chain_{c2}"] = (round(float(best[0]), 2), fmt_atom(best[1]), fmt_atom(best[2]))
        # nearest interface residue of own chain facing c2
        rs = iface_res[(ch.id, c2)]
        if rs:
            bi = (1e9, None)
            for n in rs:
                if n == POS: bi = (0.0, n); break
                x = min_dist(r, ch[n])[0]
                if x < bi[0]: bi = (x, n)
            dd[f"min_to_interface_residue_facing_{c2}"] = (round(float(bi[0]), 2), f"{ch[bi[1]].resname}{bi[1]}")
    for lig in ("PGA", "MN", "K", "FBP", "OXL", "ATP", "MG", "PYR", "WV2", "AFP"):
        L_ = hetero(m, lig)
        if not L_: continue
        own = [x for x in L_ if x.get_parent().id == ch.id]
        allb = sorted(((min_dist(r, x)[0], x) for x in L_), key=lambda t: t[0])
        if own:
            o = sorted(((min_dist(r, x)[0], x) for x in own), key=lambda t: t[0])[0]
            dd[f"min_to_{lig}_same_chain"] = (round(float(o[0]), 2), f"{o[1].get_parent().id}:{lig}{o[1].id[1]}")
        dd[f"min_to_{lig}_any_chain"] = (round(float(allb[0][0]), 2), f"{allb[0][1].get_parent().id}:{lig}{allb[0][1].id[1]}")
        if "CA" in r:
            dd[f"CA_to_nearest_{lig}_atom_same_chain"] = round(float(min(r['CA'] - b for x in own for b in x)), 2) if own else None
    # catalytic / substrate-binding residues listed by UniProt P30613 (R-type numbering) and by Valentini 2002
    for n in (116, 118, 120, 156, 157, 313, 315, 338, 339, 371):
        nn = n - OFF
        if nn in ch:
            dd[f"min_to_{ch[nn].resname}{nn}"] = round(float(min_dist(r, ch[nn])[0]), 2)
    inter_dom = {}
    dist[ch.id] = dd
rep["distances"] = dist

# 7. intra-chain domain contacts of residue POS (4.0 A), by domain
dom_contacts = {}
for ch in m:
    if POS not in ch: continue
    c = all_contacts(ch[POS], m, 4.0)
    out = {}
    for k, v in c.items():
        cid, rn = k.split(":")
        num = int("".join(x for x in rn[3:] if x.isdigit() or x == "-") or 0)
        dom = domain_of(num, OFF) if cid == ch.id and rn[:3] not in ("HOH",) else ("other chain" if cid != ch.id else "water")
        out.setdefault(dom, []).append(rn)
    dom_contacts[ch.id] = out
rep["contacts_by_domain_4A"] = dom_contacts
rep["domain_of_pos"] = domain_of(POS, OFF)

json.dump(rep, open(f"data/site_{pid}_{POS}.json", "w"), indent=1, default=str)
print(json.dumps(rep, indent=1, default=str))
