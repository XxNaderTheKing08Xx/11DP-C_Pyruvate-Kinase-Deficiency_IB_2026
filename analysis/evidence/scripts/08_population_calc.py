"""Simple population-genetic arithmetic for PKLR c.1456C>T (R486W), using gnomAD v4.1 counts.

Usage: python scripts/08_population_calc.py raw/gnomad_gnomad_r4_R486W_1-155291918-G-A.json > raw/population_calc_R486W.txt

All results are INFERENCES under stated assumptions (Hardy-Weinberg equilibrium, random mating,
no selection, gnomAD individuals as a sample of the general population); they are not measurements.
"""
import json
import math
import sys

v = json.load(open(sys.argv[1]))["data"]["variant"]
print("Assumptions: Hardy-Weinberg equilibrium; individuals = AN/2 (autosomal locus).\n")
for part in ("exome", "genome", "joint"):
    x = v[part]
    q = x["ac"] / x["an"]
    n = x["an"] / 2
    exp_hom = n * q * q
    carriers = 2 * q * (1 - q)
    print(f"{part}: AC={x['ac']} AN={x['an']} q={q:.5f} (={q*100:.3f}%); individuals~{n:.0f}; "
          f"expected homozygotes n*q^2={exp_hom:.2f}; observed={x['homozygote_count']}; "
          f"carrier frequency 2pq={carriers:.4f} (1 in {1/carriers:.0f}); q^2 per million={q*q*1e6:.1f}")

q = v["joint"]["ac"] / v["joint"]["an"]
print("\nComparison with disease prevalence (Secrest et al. 2020: diagnosed 3.2-8.5 per million; "
      "diagnosed + undiagnosed up to 51 per million) and with the share of R486W among mutated alleles "
      "in the PKD Natural History Study (Bianchi et al. 2020: 31/354 = 8.8%).")
share = 31 / 354
for prev in (3.2e-6, 8.5e-6, 51e-6):
    q_all = math.sqrt(prev)  # total frequency of pathogenic PKLR alleles if all affected individuals are biallelic
    q_r486w = share * q_all
    print(f"  prevalence {prev*1e6:>4.1f}/million -> total pathogenic allele freq {q_all:.5f}; "
          f"expected R486W freq if it contributes {share*100:.1f}% of alleles = {q_r486w:.5f}; "
          f"observed/expected = {q/q_r486w:.1f}")
print(f"\nExpected frequency of R486W homozygotes in the population (q^2): {q*q*1e6:.1f} per million "
      f"(compare with diagnosed PKD prevalence of all genotypes, 3.2-8.5 per million).")
