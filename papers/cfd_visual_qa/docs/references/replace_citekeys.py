"""Replace temporary citekeys with Zotero BetterBibTeX citekeys in main.tex."""
import re

TEX_PATH = r"C:\Users\YHJoo\Projects\OpenFOAM\papers\cfd_visual_qa\latex\main.tex"

# Mapping: our temp key -> Zotero key
REPLACEMENTS = {
    "calzolari2021deep": "calzolari.liu2021DeepLearning",
    "liu2019review": "ding.lam2019DatadrivenModel",
    "joo2026image2cfd": "joo2026AutomatedIndoor",
    "pandey2025openfoamgpt": "pandey.etal2025OpenFOAMGPTRetrievalaugmented",
    "dong2025finetuning": "dong.etal2025FinetuningLarge",
    "chow2025physbench": "chow.etal2025PhysBenchBenchmarking",
    "burgess2025microvqa": "burgess.etal2025MicroVQAMultimodal",
    "doris2025designqa": "doris.etal2024DesignQAMultimodal",
    "picard2025concept": "picard.etal2025ConceptManufacturing",
    "pandeyottley2025visliteracy": "pandey.ottley2025BenchmarkingVisionLanguage",
    "chen2024climateIQA": "chen.etal2024ClimateIQAMultimodal",
    "ezemba2025simulation": "ezemba.etal2025SimulationVs",
    "yue2024mmmu": "yue.etal2024MMMUMassive",
    "roberts2024": "roberts.etal2024SciFIBenchBenchmarking",
    "lu2024mathvista": "lu.etal2024MathVistaEvaluating",
    "wang2024charxiv": "wang.etal2024CharXivCharting",
    "alampara2025macbench": "alampara.etal2025ProbingLimitations",
    "pramanick2024spiqa": "pramanick.etal2024SPIQADataset",
    "li2024arxivcap": "li.etal2024MultimodalArXiv",
    "somasekharan2025": "somasekharan.panigrahi2025CFDLLMBenchBenchmark",
    "wang2025taml": "wang.etal2025EvaluationsLarge",
    "zhang2025autoturb": "zhang.etal2025AutoTurbUsing",
    "mdpi2026vlmcfd": "futamura.fukuda2026QuantitativeEvaluation",
    "elrefaie2025drivaernet": "elrefaie.etal2025DrivaerNetLargeScale",
    "dai2025gpt4o": "dai.etal2025ApplicationLarge",
    "alanisruiz2025cdcgan": "alanisruiz.etal2025DeepConvolutional",
    "kashefi2024": "kashefi2024MisleadingGallery",
    "schulzebuschoff2025visual": "schulzebuschoff.etal2025VisualCognition",
    "jassim2024grasp": "jassim.etal2024GRASPNovel",
    "wang2025dynsuperCLEVR": "wang2025dynsuperCLEVR",  # keep as-is, not in Zotero
    "ghaffari2024failure": "ghaffari.krishnaswamy2024ExploringFailure",
    "liu2022dlflowviz": "liu.etal2022DeepLearning",
    "banerjee2024physicsinformedcv": "banerjee.etal2024PhysicsInformedComputer",
}

with open(TEX_PATH, "r", encoding="utf-8") as f:
    text = f.read()

count = 0
for old, new in REPLACEMENTS.items():
    if old == new:
        continue
    n = text.count(old)
    if n > 0:
        text = text.replace(old, new)
        count += n
        print(f"  {old} -> {new}  ({n} replacements)")

# Also remove the entire thebibliography block
bib_start = text.find("%% Temporary bibliography")
bib_end = text.find("\\end{thebibliography}")
if bib_start != -1 and bib_end != -1:
    bib_end = bib_end + len("\\end{thebibliography}")
    bib_block = text[bib_start:bib_end]
    bib_replacement = (
        "\\bibliographystyle{elsarticle-num}\n"
        "\\bibliography{Zotero_YHJoo}"
    )
    text = text[:bib_start] + bib_replacement + text[bib_end:]
    print(f"\n  Replaced thebibliography block with \\bibliography{{Zotero_YHJoo}}")

with open(TEX_PATH, "w", encoding="utf-8") as f:
    f.write(text)

print(f"\nDone: {count} citekey replacements in {TEX_PATH}")
