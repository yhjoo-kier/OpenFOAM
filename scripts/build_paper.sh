#!/bin/bash
# Build paper PDF from Markdown source
# Usage: bash scripts/build_paper.sh

set -e

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LATEX_DIR="$PROJECT_DIR/latex"
SCRIPTS_DIR="$PROJECT_DIR/scripts"
MD_SOURCE="$PROJECT_DIR/docs/paper_draft_v1.md"
TEX_FILE="$LATEX_DIR/paper.tex"
BIB_SOURCE="C:/Vaults/Research/Zotero_YHJoo.bib"
FIG_SOURCE="$PROJECT_DIR/results/paper_figures_phase2"

# 1. Create latex directory
mkdir -p "$LATEX_DIR/figures"

# 2. Convert Markdown -> LaTeX
echo "=== Converting Markdown to LaTeX ==="
# Use PaperSearch venv Python (Windows Store python3 is a stub)
PYTHON="C:/Users/YHJoo/Projects/PaperSearch/.venv/Scripts/python.exe"
"$PYTHON" "$SCRIPTS_DIR/md2latex.py" "$MD_SOURCE" "$TEX_FILE"

# 3. Copy figures (PDF from phase2, then PaperBanana PNGs override)
echo "=== Copying figures ==="
cp "$FIG_SOURCE"/fig_*.pdf "$LATEX_DIR/figures/" 2>/dev/null || echo "  (no PDF figures found)"

# PaperBanana-generated figures (PNG, override corresponding PDFs)
BANANA_FIGS="$PROJECT_DIR/results/paperbanana_selected"
if [ -d "$BANANA_FIGS" ]; then
    for png in "$BANANA_FIGS"/*.png; do
        base=$(basename "$png" .png)
        # Remove the PDF so LaTeX picks the PNG
        rm -f "$LATEX_DIR/figures/${base}.pdf"
        cp "$png" "$LATEX_DIR/figures/"
    done
    echo "  PaperBanana figures applied"
fi

# 4. Copy BibTeX
echo "=== Copying BibTeX ==="
cp "$BIB_SOURCE" "$LATEX_DIR/Zotero_YHJoo.bib"

# 5. Build PDF
echo "=== Building PDF (pdflatex -> bibtex -> pdflatex -> pdflatex) ==="
cd "$LATEX_DIR"
pdflatex -interaction=nonstopmode paper.tex > build_pass1.log 2>&1 || true
bibtex paper > build_bibtex.log 2>&1 || true
pdflatex -interaction=nonstopmode paper.tex > build_pass2.log 2>&1 || true
pdflatex -interaction=nonstopmode paper.tex > build_pass3.log 2>&1 || true

# 6. Check result
if [ -f paper.pdf ]; then
    echo "=== SUCCESS: $LATEX_DIR/paper.pdf ==="
    ls -lh paper.pdf
else
    echo "=== FAILED: No PDF generated ==="
    echo "Check logs: build_pass1.log, build_bibtex.log"
    exit 1
fi
