"""Match our 33 temporary citekeys to actual Zotero BetterBibTeX citekeys."""
import re

BIB_PATH = r"C:\Vaults\Research\Zotero_YHJoo.bib"
TEX_PATH = r"C:\Users\YHJoo\Projects\OpenFOAM\papers\cfd_visual_qa\latex\main.tex"

# Our temp keys and identifying info (author fragment, year, title fragment or DOI)
MATCH_RULES = [
    ("calzolari2021deep", "calzolari", "2021", "deep learning", "10.1016/j.buildenv.2021.108315"),
    ("liu2019review", "ding", "2019", "cross ventilation", "10.1016/j.buildenv.2019.106394"),
    ("joo2026image2cfd", "joo", "2026", "indoor cfd", None),
    ("pandey2025openfoamgpt", "pandey", "2025", "openfoamgpt", "10.1063/5.0257555"),
    ("dong2025finetuning", "dong", "2025", "fine-tuning", "10.1016/j.taml.2025.100594"),
    ("chow2025physbench", "chow", "2025", "physbench", None),
    ("burgess2025microvqa", "burgess", "2025", "microvqa", "10.1109/CVPR52734.2025.01821"),
    ("doris2025designqa", "doris", "2025", "designqa", "10.1115/1.4067333"),
    ("picard2025concept", "picard", "2025", "concept to manufacturing", "10.1007/s10462-025-11290-y"),
    ("pandeyottley2025visliteracy", "pandey", "2025", "visualization literacy", None),
    ("chen2024climateIQA", "chen", "2024", "climateiq", None),
    ("ezemba2025simulation", "ezemba", "2025", "simulation", "10.1145/3722573.3727826"),
    ("yue2024mmmu", "yue", "2024", "mmmu", "10.1109/CVPR52733.2024.00913"),
    ("roberts2024", "roberts", "2024", "scifibench", None),
    ("lu2024mathvista", "lu", "2024", "mathvista", None),
    ("wang2024charxiv", "wang", "2024", "charxiv", None),
    ("alampara2025macbench", "alampara", "2025", "macbench", "10.1038/s43588-025-00836-3"),
    ("pramanick2024spiqa", "pramanick", "2024", "spiqa", None),
    ("li2024arxivcap", "li", "2024", "multimodal arxiv", "10.18653/v1/2024.acl-long.775"),
    ("somasekharan2025", "somasekharan", "2025", "cfdllmbench", None),
    ("wang2025taml", "wang", "2025", "evaluations.*computational fluid", "10.1016/j.taml.2025.100597"),
    ("zhang2025autoturb", "zhang", "2025", "autoturb", "10.1063/5.0247759"),
    ("mdpi2026vlmcfd", None, "2026", "mixed.?reality.*indoor", "10.3390/technologies14030157"),
    ("elrefaie2025drivaernet", "elrefaie", "2025", "drivaernet", None),
    ("dai2025gpt4o", "dai", "2025", "indoor air distribution", "10.1016/j.buildenv.2025.113647"),
    ("alanisruiz2025cdcgan", "alanis", "2025", "pollutant dispersion", "10.1016/j.buildenv.2025.112856"),
    ("kashefi2024", "kashefi", "2024", "misleading gallery", None),
    ("schulzebuschoff2025visual", "schulze", "2025", "visual cognition", "10.1038/s42256-024-00963-y"),
    ("jassim2024grasp", "jassim", "2024", "grasp", None),
    ("wang2025dynsuperCLEVR", "wang", "2025", "compositional 4d", None),
    ("ghaffari2024failure", "ghaffari", "2024", "failure cases", "10.1609/aaaiss.v3i1.31189"),
    ("liu2022dlflowviz", "liu", "2022", "flow visualization", "10.1186/s42774-022-00113-1"),
    ("banerjee2024physicsinformedcv", "banerjee", "2024", "physics-informed computer vision", "10.1145/3689037"),
]

# Parse bib file into entries
with open(BIB_PATH, "r", encoding="utf-8") as f:
    bib_text = f.read()

# Extract entries: citekey, full text block
entries = []
for m in re.finditer(r'@\w+\{([^,]+),\s*\n(.*?)(?=\n@|\Z)', bib_text, re.DOTALL):
    citekey = m.group(1).strip()
    block = m.group(2).lower()
    entries.append((citekey, block))

print(f"Loaded {len(entries)} bib entries\n")
print(f"{'Our Key':<40} {'Zotero Key':<50} {'Match Method'}")
print("-" * 120)

mapping = {}
not_found = []

for our_key, author_frag, year, title_frag, doi in MATCH_RULES:
    found = None
    method = ""

    # Try DOI match first (most reliable)
    if doi:
        doi_lower = doi.lower()
        for ck, block in entries:
            if doi_lower in block:
                found = ck
                method = "DOI"
                break

    # Try author + year + title fragment
    if not found and author_frag:
        for ck, block in entries:
            if (author_frag.lower() in block
                    and year in block
                    and re.search(title_frag.lower(), block)):
                found = ck
                method = "author+year+title"
                break

    # Try title fragment + year only
    if not found:
        for ck, block in entries:
            if year in block and re.search(title_frag.lower(), block):
                found = ck
                method = "year+title"
                break

    if found:
        mapping[our_key] = found
        print(f"{our_key:<40} {found:<50} {method}")
    else:
        not_found.append(our_key)
        print(f"{our_key:<40} {'*** NOT FOUND ***':<50}")

print(f"\n=== Summary ===")
print(f"Matched: {len(mapping)}/{len(MATCH_RULES)}")
print(f"Not found: {len(not_found)}")
if not_found:
    print(f"  Missing: {not_found}")

# Output sed-like replacement commands
if mapping:
    print(f"\n=== Replacement Map (our_key -> zotero_key) ===")
    for old, new in sorted(mapping.items()):
        if old != new:
            print(f"  {old}  ->  {new}")
