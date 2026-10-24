"""Parse ClinVar VCV XML (from 01_clinvar.py) into a readable summary.

Usage: python scripts/02_parse_clinvar.py raw/clinvar_efetch_vcv_R486W_1513.xml > raw/clinvar_summary_R486W_1513.txt
"""
import sys
import xml.etree.ElementTree as ET
from collections import Counter

root = ET.parse(sys.argv[1]).getroot()
va = root.find("VariationArchive")
print("VCV", va.get("Accession"), "version", va.get("Version"), "| name:", va.get("VariationName"))
print("NumberOfSubmissions", va.get("NumberOfSubmissions"), "NumberOfSubmitters", va.get("NumberOfSubmitters"),
      "DateLastUpdated", va.get("DateLastUpdated"), "DateCreated", va.get("DateCreated"))
cr = va.find("ClassifiedRecord")
for sl in cr.find("SimpleAllele").find("Location").findall("SequenceLocation"):
    print("Location", sl.get("Assembly"), sl.get("Accession"), sl.get("positionVCF"), sl.get("referenceAlleleVCF"), ">", sl.get("alternateAlleleVCF"))
print("SPDI", cr.find("SimpleAllele").findtext("CanonicalSPDI"))
for x in cr.find("SimpleAllele").iter("XRef"):
    print("XRef", x.get("DB"), x.get("ID"), x.get("Type") or "")
af = cr.find("SimpleAllele").find("AlleleFrequencyList")
if af is not None:
    for a in af:
        print("AlleleFrequency", a.attrib)
print("\n== Aggregate germline classification ==")
cl = cr.find("Classifications")
g = cl.find("GermlineClassification")
print("Attributes:", g.attrib)
print("Review status:", g.findtext("ReviewStatus"), "| Description:", g.findtext("Description"),
      "| Explanation:", g.findtext("Explanation"))
print("Aggregate citations (PubMed):", [c.findtext("ID") for c in g.findall("Citation") if c.find("ID") is not None and c.find("ID").get("Source") == "PubMed"])
print("\n== RCV records ==")
for rcv in cr.find("RCVList").findall("RCVAccession"):
    conds = [t.findtext("Name/ElementValue") for t in rcv.iter("ClassifiedCondition")]
    conds = [c.text for c in rcv.iter("ClassifiedCondition")]
    rc = rcv.find("RCVClassifications/GermlineClassification")
    if rc is None:
        rc = rcv.find("RCVClassifications/NoClassification")
    print(rcv.get("Accession"), "v" + rcv.get("Version"), "|", "; ".join(conds), "|",
          rc.findtext("Description"), "|", rc.findtext("ReviewStatus"),
          "| submissions:", rc.find("Description").get("SubmissionCount"))
print("\n== Submitted records (SCV) ==")
counts = Counter()
for ca in cr.find("ClinicalAssertionList").findall("ClinicalAssertion"):
    acc = ca.find("ClinVarAccession")
    cls = ca.find("Classification")
    desc = cls.findtext("GermlineClassification") if cls is not None else None
    counts[desc] += 1
    traits = [t.findtext("Name/ElementValue") for t in ca.iter("Trait")]
    method = ca.findtext("ObservedInList/ObservedIn/Method/MethodType")
    origin = ca.findtext("ObservedInList/ObservedIn/Sample/Origin")
    print(acc.get("Accession"), "v" + acc.get("Version"), "|", acc.get("SubmitterName"), "|",
          desc, "|", cls.findtext("ReviewStatus") if cls is not None else "", "|",
          cls.get("DateLastEvaluated") if cls is not None else "", "|", "; ".join([t for t in traits if t]),
          "|", method, "|", origin)
print("\nClassification counts across SCVs:", dict(counts))
