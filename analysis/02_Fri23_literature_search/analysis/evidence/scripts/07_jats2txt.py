"""Convert JATS full-text XML (Europe PMC / PMC) to plain text, keeping section titles and table text.
Usage: python scripts/07_jats2txt.py raw/fulltext_PMCxxxx.xml ... (writes .txt next to each)"""
import sys, re
import xml.etree.ElementTree as ET
for f in sys.argv[1:]:
    root = ET.parse(f).getroot()
    out = []
    for el in root.iter():
        tag = el.tag.split('}')[-1]
        if tag in ("title", "p", "caption", "tr", "label", "article-title", "abstract"):
            t = " ".join("".join(el.itertext()).split())
            if t:
                out.append(("\n## " if tag == "title" else "") + t)
    open(f[:-4] + ".txt", "w").write("\n".join(out))
    print(f, len(out))
