"""Fetch RIS entries via DOI content negotiation for all cited papers."""
import urllib.request
import time
import sys

# All papers with known DOIs
papers = [
    # Key, DOI
    ("calzolari2021deep", "10.1016/j.buildenv.2021.108315"),
    ("joo2026image2cfd", None),  # submitted, no DOI yet
    ("pandey2025openfoamgpt", "10.1063/5.0257555"),
    ("dong2025finetuning", "10.1016/j.taml.2025.100594"),
    ("chow2025physbench", None),  # ICLR, no DOI in Scopus
    ("burgess2025microvqa", "10.1109/CVPR52734.2025.01821"),
    ("doris2025designqa", "10.1115/1.4067333"),
    ("picard2025concept", "10.1007/s10462-025-11290-y"),
    ("pandeyottley2025visliteracy", None),  # EuroVis 2025, arXiv:2503.16632
    ("chen2024climateIQA", None),  # arXiv:2406.09838
    ("ezemba2025simulation", "10.1145/3722573.3727826"),
    ("yue2024mmmu", "10.1109/CVPR52733.2024.00913"),
    ("roberts2024", None),  # SciFIBench, arXiv:2405.08807
    ("lu2024mathvista", None),  # ICLR 2024
    ("wang2024charxiv", None),  # NeurIPS 2024
    ("alampara2025macbench", "10.1038/s43588-025-00836-3"),
    ("pramanick2024spiqa", None),  # NeurIPS 2024
    ("li2024arxivcap", "10.18653/v1/2024.acl-long.775"),
    ("somasekharan2025", None),  # arXiv:2509.20374
    ("wang2025taml", "10.1016/j.taml.2025.100597"),
    ("zhang2025autoturb", "10.1063/5.0247759"),
    ("mdpi2026vlmcfd", "10.3390/technologies14030157"),
    ("elrefaie2025drivaernet", None),  # NeurIPS 2025
    ("dai2025gpt4o", "10.1016/j.buildenv.2025.113647"),
    ("alanisruiz2025cdcgan", "10.1016/j.buildenv.2025.112856"),
    ("kashefi2024", None),  # arXiv:2405.15406
    ("schulzebuschoff2025visual", "10.1038/s42256-024-00963-y"),
    ("jassim2024grasp", None),  # IJCAI 2024
    ("wang2025dynsuperCLEVR", None),  # ICLR 2025
    ("ghaffari2024failure", "10.1609/aaaiss.v3i1.31189"),
    ("liu2022dlflowviz", "10.1186/s42774-022-00113-1"),
    ("banerjee2024physicsinformedcv", "10.1145/3689037"),
    ("liu2019review", "10.1016/j.buildenv.2019.106394"),
]

out_path = sys.argv[1] if len(sys.argv) > 1 else "all_references.ris"
success = 0
failed = []
manual = []

with open(out_path, "w", encoding="utf-8") as f:
    for key, doi in papers:
        if doi is None:
            manual.append(key)
            continue
        url = f"https://doi.org/{doi}"
        req = urllib.request.Request(url, headers={
            "Accept": "application/x-research-info-systems",
            "User-Agent": "PaperSearch/1.0 (mailto:yhjoo@kier.re.kr)"
        })
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                ris_text = resp.read().decode("utf-8", errors="replace")
                # Add custom note with our citekey
                ris_text = ris_text.rstrip()
                if not ris_text.endswith("ER  -"):
                    ris_text += "\nER  -"
                f.write(f"// Key: {key}\n")
                f.write(ris_text + "\n\n")
                success += 1
                print(f"  OK: {key} ({doi})")
        except Exception as e:
            failed.append((key, doi, str(e)))
            print(f"  FAIL: {key} ({doi}) - {e}")
        time.sleep(0.5)

print(f"\n=== Results ===")
print(f"DOI fetched: {success}")
print(f"DOI failed: {len(failed)}")
for k, d, e in failed:
    print(f"  {k}: {d} -> {e}")
print(f"No DOI (manual): {len(manual)}")
for k in manual:
    print(f"  {k}")
print(f"\nRIS saved to: {out_path}")
