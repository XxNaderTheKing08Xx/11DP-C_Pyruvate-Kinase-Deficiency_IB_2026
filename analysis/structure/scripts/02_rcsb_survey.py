"""Survey all RCSB entries mapped to UniProt P30613 (human PKLR) and P12928 (rat Pklr) via RCSB GraphQL API.
Output: data/rcsb_survey.tsv"""
import json, requests, sys
ids=[]
for acc in ["P30613","P12928"]:
    d=json.load(open(f"data/rcsb_search_{acc}.json"))
    ids+= [(acc,r['identifier']) for r in d['result_set']]
q='''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} exptl{method}
 rcsb_entry_info{resolution_combined deposited_polymer_monomer_count polymer_entity_count_protein}
 refine{ls_R_factor_R_free} rcsb_accession_info{deposit_date initial_release_date}
 rcsb_primary_citation{pdbx_database_id_PubMed pdbx_database_id_DOI year journal_abbrev title}
 polymer_entities{rcsb_polymer_entity{pdbx_mutation pdbx_description} entity_poly{pdbx_seq_one_letter_code_can} rcsb_polymer_entity_container_identifiers{reference_sequence_identifiers{database_accession}} rcsb_entity_source_organism{scientific_name}}
 nonpolymer_entities{nonpolymer_comp{chem_comp{id name}}}}}'''
r=requests.post("https://data.rcsb.org/graphql",json={"query":q,"variables":{"ids":[i for _,i in ids]}},timeout=120)
r.raise_for_status(); data=r.json()["data"]["entries"]
json.dump(data,open("data/rcsb_survey.json","w"),indent=1)
acc_of=dict((i,a) for a,i in ids)
with open("data/rcsb_survey.tsv","w") as f:
    f.write("pdb\tuniprot\ttitle\tmethod\tresolution_A\tRfree\tdeposited\tmutation\tligands\tPMID\tDOI\n")
    for e in sorted(data,key=lambda e:e['rcsb_accession_info']['deposit_date']):
        muts=";".join(str(p['rcsb_polymer_entity'].get('pdbx_mutation')) for p in e['polymer_entities'])
        ligs=";".join(f"{n['nonpolymer_comp']['chem_comp']['id']}({n['nonpolymer_comp']['chem_comp']['name']})" for n in (e.get('nonpolymer_entities') or []))
        res=e['rcsb_entry_info']['resolution_combined']; rf=(e.get('refine') or [{}])[0].get('ls_R_factor_R_free')
        c=e.get('rcsb_primary_citation') or {}
        f.write("\t".join(map(str,[e['rcsb_id'],acc_of[e['rcsb_id']],e['struct']['title'],e['exptl'][0]['method'],res[0] if res else None,rf,e['rcsb_accession_info']['deposit_date'][:10],muts,ligs,c.get('pdbx_database_id_PubMed'),c.get('pdbx_database_id_DOI')]))+"\n")
print(len(data))
