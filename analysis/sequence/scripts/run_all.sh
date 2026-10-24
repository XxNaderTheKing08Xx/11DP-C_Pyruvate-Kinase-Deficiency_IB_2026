#!/bin/sh
# Re-run the whole sequence analysis. Run from the analysis/sequence folder:
#   sh scripts/run_all.sh <bin dir with clustalo, needle and blastp/blastn>
# The bin dir is that of a conda environment with Clustal Omega 1.2.4, EMBOSS 6.6.0 and BLAST+ 2.17.0
# (see env/bioenv_tools.txt), e.g. created with
#   micromamba create -p ./bioenv -c conda-forge -c bioconda clustalo=1.2.4 emboss=6.6.0 blast=2.17.0
# Python: set PYTHON to an interpreter with the packages in requirements.txt (default: python3).
# Scripts 02 (ClinVar) and 03 (Ensembl gene-symbol lookup) query web services; all other inputs are in raw/.
# Script 10 writes Figure 1 to the repository folder figures/ (../../figures/).
set -e
BIN="$1"; PY="${PYTHON:-python3}"
$PY -I scripts/01_characterize.py
$PY -I scripts/02_clinvar.py
$PY -I scripts/03_orthologs.py
mkdir -p results/msa_local
for p in core_panel broad_panel; do
  "$BIN/clustalo" -i results/$p.fasta -t Protein --outfmt=clustal --resno --wrap=60 -o results/msa_local/$p.clustalo.aln --force --distmat-out=results/msa_local/$p.distmat --full --percent-id --guidetree-out=results/msa_local/$p.guidetree
  "$BIN/clustalo" -i results/$p.fasta -t Protein --outfmt=fasta -o results/msa_local/$p.clustalo.afa --force
done
$PY -I scripts/05_pairwise.py "$BIN"
$PY -I scripts/06_conservation.py
$PY -I scripts/07_compara_survey.py
$PY -I scripts/08_primate_refseq_check.py
$PY -I scripts/09_verify.py "$BIN"
$PY -I scripts/10_figure_table.py
