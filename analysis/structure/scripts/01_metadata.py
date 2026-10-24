"""Extract metadata from mmCIF files (method, resolution, R-free, entities, ligands, chains, observed residues)."""
import sys, json
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

def L(d,k):
    v=d.get(k,[]); return v if isinstance(v,list) else [v]

for pid in sys.argv[1:]:
    f=f"raw/{pid}.cif"; d=MMCIF2Dict(f)
    print("=====",pid)
    for k in ["_struct.title","_exptl.method","_refine.ls_d_res_high","_refine.ls_R_factor_R_work","_refine.ls_R_factor_R_free",
              "_pdbx_database_status.recvd_initial_deposition_date","_symmetry.space_group_name_H-M","_cell.length_a","_cell.length_b","_cell.length_c",
              "_citation.title","_citation.journal_abbrev","_citation.year","_citation.journal_volume","_citation.page_first","_citation.page_last",
              "_citation.pdbx_database_id_DOI","_citation.pdbx_database_id_PubMed"]:
        print(f"  {k}: {L(d,k)}")
    print("  revisions:", list(zip(L(d,'_pdbx_audit_revision_history.revision_date'),L(d,'_pdbx_audit_revision_history.major_revision'),L(d,'_pdbx_audit_revision_history.minor_revision'))))
    print("  authors(citation):", [a for a,c in zip(L(d,'_citation_author.name'),L(d,'_citation_author.citation_id')) if c=='primary'])
    print("  entities:")
    for i,t,desc,n,mw in zip(L(d,'_entity.id'),L(d,'_entity.type'),L(d,'_entity.pdbx_description'),L(d,'_entity.pdbx_number_of_molecules'),L(d,'_entity.formula_weight')):
        print(f"    {i} {t} {desc} n={n} mw={mw}")
    print("  entity_poly strand ids:", list(zip(L(d,'_entity_poly.entity_id'),L(d,'_entity_poly.pdbx_strand_id'))))
    print("  src organism:", L(d,'_entity_src_gen.pdbx_gene_src_scientific_name'), L(d,'_entity_src_gen.pdbx_gene_src_gene'), L(d,'_entity_src_gen.pdbx_host_org_scientific_name'))
    print("  struct_ref:", list(zip(L(d,'_struct_ref.db_name'),L(d,'_struct_ref.pdbx_db_accession'),L(d,'_struct_ref.db_code'))))
    print("  struct_ref_seq:", list(zip(L(d,'_struct_ref_seq.pdbx_strand_id'),L(d,'_struct_ref_seq.seq_align_beg'),L(d,'_struct_ref_seq.seq_align_end'),L(d,'_struct_ref_seq.db_align_beg'),L(d,'_struct_ref_seq.db_align_end'),L(d,'_struct_ref_seq.pdbx_auth_seq_align_beg'),L(d,'_struct_ref_seq.pdbx_auth_seq_align_end'))))
    print("  seq_dif:", list(zip(L(d,'_struct_ref_seq_dif.pdbx_pdb_strand_id'),L(d,'_struct_ref_seq_dif.mon_id'),L(d,'_struct_ref_seq_dif.pdbx_auth_seq_num'),L(d,'_struct_ref_seq_dif.db_mon_id'),L(d,'_struct_ref_seq_dif.pdbx_seq_db_seq_num'),L(d,'_struct_ref_seq_dif.details'))))
    print("  assemblies:", list(zip(L(d,'_pdbx_struct_assembly.id'),L(d,'_pdbx_struct_assembly.details'),L(d,'_pdbx_struct_assembly.oligomeric_details'))))
    print("  assembly_gen:", list(zip(L(d,'_pdbx_struct_assembly_gen.assembly_id'),L(d,'_pdbx_struct_assembly_gen.oper_expression'),L(d,'_pdbx_struct_assembly_gen.asym_id_list'))))
    print("  oper:", list(zip(L(d,'_pdbx_struct_oper_list.id'),L(d,'_pdbx_struct_oper_list.type'))))
    s=MMCIFParser(QUIET=True).get_structure(pid,f)
    m=s[0]
    for ch in m:
        aa=[r for r in ch if is_aa(r,standard=False) and r.id[0] in (' ','H_SEP','H_TPO')]
        aa=[r for r in ch if r.id[0]==' ' or r.resname in ('SEP','TPO')]
        het=[r for r in ch if r.id[0].startswith('H_') and r.resname not in ('SEP','TPO')]
        nums=[r.id[1] for r in aa]
        # gaps
        gaps=[]
        for a,b in zip(nums,nums[1:]):
            if b-a>1: gaps.append(f"{a+1}-{b-1}")
        from collections import Counter
        print(f"  chain {ch.id}: {len(aa)} aa residues, first {nums[0] if nums else None} last {nums[-1] if nums else None}, gaps {gaps}; het {dict(Counter(r.resname for r in het))}; waters {sum(1 for r in ch if r.id[0]=='W')}")
