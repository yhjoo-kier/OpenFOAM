"""Check that all cite keys in main.tex exist in Zotero_YHJoo.bib."""
import re

BIB = r"C:\Vaults\Research\Zotero_YHJoo.bib"
TEX = r"C:\Users\YHJoo\Projects\OpenFOAM\papers\cfd_visual_qa\latex\main.tex"

with open(TEX, "r", encoding="utf-8") as f:
    tex = f.read()

cites = set()
for m in re.findall(r"\\cite\{([^}]+)\}", tex):
    for k in m.split(","):
        cites.add(k.strip())

with open(BIB, "r", encoding="utf-8") as f:
    bib = f.read()

bib_keys = set(re.findall(r"@\w+\{([^,]+),", bib))

found = cites & bib_keys
missing = cites - bib_keys

print(f"Cite keys in tex: {len(cites)}")
print(f"Found in .bib:    {len(found)}")
print(f"Missing from .bib: {len(missing)}")
if missing:
    for m in sorted(missing):
        print(f"  {m}")
else:
    print("All cite keys resolved!")
