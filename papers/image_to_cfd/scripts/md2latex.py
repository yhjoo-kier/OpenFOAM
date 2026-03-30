#!/usr/bin/env python3
"""Convert paper_draft_v1.md to elsarticle LaTeX.

Usage:
    python scripts/md2latex.py docs/paper_draft_v1.md latex/paper.tex
"""

import re
import sys
from pathlib import Path

FIGURE_DIR = "figures"  # relative to .tex file
BIB_FILE = "Zotero_YHJoo"  # without .bib extension


# ---------------------------------------------------------------------------
# Preamble template
# ---------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[preprint,12pt,number]{elsarticle}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{siunitx}
\usepackage{hyperref}
\usepackage{algorithm}
\usepackage{algorithmic}
\usepackage{multirow}
\usepackage{textcomp}

\graphicspath{{%(figdir)s/}}

\journal{Draft}

\begin{document}

\begin{frontmatter}

\title{%(title)s}

\author[kier,ust]{Younghwan Joo}
\ead{yhjoo@kier.re.kr}
\address[kier]{Energy Efficiency Research Division, Korea Institute of Energy Research, 152 Gajeong-ro, Yuseong-gu, Daejeon 34129, Republic of Korea}
\address[ust]{Energy Engineering, University of Science \& Technology, 217 Gajeong-ro, Yuseong-gu, Daejeon 34129, Republic of Korea}

\begin{abstract}
%(abstract)s
\end{abstract}

\begin{keyword}
indoor CFD \sep vision-language model \sep geometric abstraction \sep automated simulation \sep benchmark
\end{keyword}

\end{frontmatter}

"""

POSTAMBLE = r"""
\bibliographystyle{elsarticle-num}
\bibliography{%(bibfile)s}

\end{document}
"""


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def convert_citations(text):
    """[@key1; @key2; ...] -> \\cite{key1,key2,...}  and  [@key] -> \\cite{key}"""
    def _multi(m):
        keys = re.findall(r'@([^;\]\s]+)', m.group(0))
        return r'\cite{' + ','.join(keys) + '}'
    text = re.sub(r'\[(@[^]]+)\]', _multi, text)
    return text


def convert_figure_refs(text):
    """[Fig:label] -> Fig.~\\ref{fig:label}"""
    def _repl(m):
        label = m.group(1)
        return r'Fig.~\ref{fig:' + label + '}'
    return re.sub(r'\[Fig:([^\]]+)\]', _repl, text)


def convert_figure_refs_caption(text):
    """Same as convert_figure_refs but for use inside captions."""
    def _repl(m):
        label = m.group(1)
        return r'Fig.~\ref{fig:' + label + '}'
    return re.sub(r'\[Fig:([^\]]+)\]', _repl, text)


def convert_table_refs(text):
    """[Table:label] -> Table~\\ref{tab:label}"""
    def _repl(m):
        label = m.group(1)
        return r'Table~\ref{tab:' + label + '}'
    return re.sub(r'\[Table:([^\]]+)\]', _repl, text)


def convert_equation_refs(text):
    """[Eq:label] -> Eq.~\\ref{eq:label}"""
    def _repl(m):
        label = m.group(1)
        return r'Eq.~\ref{eq:' + label + '}'
    return re.sub(r'\[Eq:([^\]]+)\]', _repl, text)


def convert_section_refs(text):
    """[Sec:label] -> Section~\\ref{sec:label}"""
    def _repl(m):
        label = m.group(1)
        return r'Section~\ref{sec:' + label + '}'
    return re.sub(r'\[Sec:([^\]]+)\]', _repl, text)


def convert_inline_formatting(text):
    """**bold** -> \\textbf{bold}, *italic* -> \\textit{italic}"""
    # Bold first (before italic)
    text = re.sub(r'\*\*([^*]+)\*\*', r'\\textbf{\1}', text)
    # Italic (but not within math $...$ and not \section*{ or \subsection*{)
    text = re.sub(r'(?<!\$)(?<!\\section)(?<!\\subsection)\*([^*$\n{]+)\*(?!\$)', r'\\textit{\1}', text)
    return text


def escape_underscores_outside_math(text):
    """Escape _ outside math mode and LaTeX commands."""
    result = []
    in_math = False
    in_inline_math = False
    i = 0
    while i < len(text):
        ch = text[i]
        # Track display math $$
        if ch == '$' and i + 1 < len(text) and text[i + 1] == '$':
            in_math = not in_math
            result.append('$$')
            i += 2
            continue
        # Track inline math $
        if ch == '$':
            in_inline_math = not in_inline_math
            result.append('$')
            i += 1
            continue
        # Track LaTeX commands (don't escape _ in \label{}, \ref{}, \cite{}, etc.)
        if ch == '\\' and i + 1 < len(text) and text[i + 1].isalpha():
            # Read the whole command
            cmd_start = i
            i += 1
            while i < len(text) and text[i].isalpha():
                i += 1
            result.append(text[cmd_start:i])
            # If followed by {, read the braced argument
            if i < len(text) and text[i] == '{':
                brace_depth = 0
                while i < len(text):
                    if text[i] == '{':
                        brace_depth += 1
                    elif text[i] == '}':
                        brace_depth -= 1
                    result.append(text[i])
                    i += 1
                    if brace_depth == 0:
                        break
            continue
        # Escape _ outside math
        if ch == '_' and not in_math and not in_inline_math:
            result.append(r'\_')
            i += 1
            continue
        result.append(ch)
        i += 1
    return ''.join(result)


def escape_percent(text):
    """Escape % as \\% (LaTeX comment character)."""
    # Don't escape % that are already escaped or in LaTeX comments
    return re.sub(r'(?<!\\)%', r'\\%', text)


def convert_unicode(text):
    """Replace Unicode characters with LaTeX equivalents."""
    replacements = {
        '≥': r'$\geq$',
        '≤': r'$\leq$',
        '≈': r'$\approx$',
        '×': r'$\times$',
        '±': r'$\pm$',
        '°': r'$^\circ$',
        '²': r'$^2$',
        '³': r'$^3$',
        '—': '---',
        '–': '--',
        '→': r'$\rightarrow$',
        '←': r'$\leftarrow$',
        '∞': r'$\infty$',
        '\u2009': r'\,',  # thin space
        'Δ': r'$\Delta$',
    }
    for char, repl in replacements.items():
        text = text.replace(char, repl)
    return text


def escape_ampersand_in_text(text):
    """Escape & in regular text but NOT inside tabular or \\cite{}."""
    # Only escape & that are not part of LaTeX table separators
    # This is applied per-line to non-table content
    if '\\begin{tabular}' in text or '\\end{tabular}' in text:
        return text
    # Don't escape in lines that look like table rows (contain \\ at end)
    if text.rstrip().endswith('\\\\'):
        return text
    if '\\toprule' in text or '\\midrule' in text or '\\bottomrule' in text:
        return text
    return text.replace('&', r'\&')


def label_from_fig_ref(label):
    """Convert Fig label to filename: method-framework -> fig_method_framework.
    LaTeX graphicx without extension tries pdf then png automatically."""
    return 'fig_' + label.replace('-', '_')


def convert_equations(text):
    """Convert $$...\tag{Eq:label}$$ blocks to equation environments."""
    def _repl(m):
        content = m.group(1).strip()
        # Extract \tag{Eq:label}
        tag_match = re.search(r'\\tag\{Eq:([^}]+)\}', content)
        if tag_match:
            label = tag_match.group(1)
            content = re.sub(r'\s*\\tag\{Eq:[^}]+\}', '', content)
            return (
                '\\begin{equation}\n'
                + content + '\n'
                + '\\label{eq:' + label + '}\n'
                + '\\end{equation}'
            )
        else:
            return '\\begin{equation}\n' + content + '\n\\end{equation}'

    text = re.sub(r'\$\$(.*?)\$\$', _repl, text, flags=re.DOTALL)
    return text


def convert_markdown_table(lines):
    """Convert markdown table lines to LaTeX tabular.
    Returns (latex_string, num_lines_consumed)."""
    if not lines or '|' not in lines[0]:
        return None, 0

    # Parse header
    header = [c.strip() for c in lines[0].strip().strip('|').split('|')]
    ncols = len(header)

    # Skip separator line
    sep_idx = 1
    if sep_idx < len(lines) and re.match(r'\s*\|[-:\s|]+\|\s*$', lines[sep_idx]):
        data_start = 2
    else:
        data_start = 1

    # Parse data rows
    rows = []
    i = data_start
    while i < len(lines) and '|' in lines[i] and lines[i].strip():
        cells = [c.strip() for c in lines[i].strip().strip('|').split('|')]
        rows.append(cells)
        i += 1

    # Build LaTeX
    col_spec = 'l' * ncols
    out = []
    out.append('\\begin{tabular}{' + col_spec + '}')
    out.append('\\toprule')
    out.append(' & '.join(header) + ' \\\\')
    out.append('\\midrule')
    for row in rows:
        # Pad if needed
        while len(row) < ncols:
            row.append('')
        out.append(' & '.join(row[:ncols]) + ' \\\\')
    out.append('\\bottomrule')
    out.append('\\end{tabular}')

    return '\n'.join(out), i


# Figures that use LaTeX subfigure instead of composite matplotlib PDF
SUBFIGURE_DEFS = {
    'bench-cfd-showcase': {
        'layout': '2x2',  # 2 rows, 2 cols
        'panels': [
            {'label': '(a)', 'title': 'A1: Rectangular, simple',
             'files': ['fig_bench_cfd_a1_01_geo', 'fig_bench_cfd_a1_01_flow']},
            {'label': '(b)', 'title': 'A2: Rectangular, dense',
             'files': ['fig_bench_cfd_a2_03_geo', 'fig_bench_cfd_a2_03_flow']},
            {'label': '(c)', 'title': 'A3: Composite, simple',
             'files': ['fig_bench_cfd_a3_03_geo', 'fig_bench_cfd_a3_03_flow']},
            {'label': '(d)', 'title': 'A4: Composite, dense',
             'files': ['fig_bench_cfd_a4_03_geo', 'fig_bench_cfd_a4_03_flow']},
        ],
    },
    'discuss-structure-cfd-gap': {
        'layout': '2x2',
        'panels': [
            {'label': '(a)', 'title': 'A4-02 geometry overlap',
             'files': ['fig_discuss_structure_cfd_gap']},  # keep composite for geometry panels
            {'label': '(b)', 'title': 'A4-02 predicted CFD',
             'files': ['fig_gap_a4_02_flow']},
            {'label': '(c)', 'title': 'A4-04 geometry overlap',
             'files': ['fig_discuss_structure_cfd_gap']},  # reuse composite
            {'label': '(d)', 'title': 'A4-04 predicted CFD',
             'files': ['fig_gap_a4_04_flow']},
        ],
    },
    'demo-floorplan-application': {
        'layout': '2x2',
        'panels': [
            {'label': '(a)', 'title': 'Case 1: Input floor plan',
             'files': ['fig_demo_floorplan_application']},  # keep composite for input panels
            {'label': '(b)', 'title': 'Case 1: CFD velocity field',
             'files': ['fig_demo_case1_flow']},
            {'label': '(c)', 'title': 'Case 2: Input floor plan',
             'files': ['fig_demo_floorplan_application']},
            {'label': '(d)', 'title': 'Case 2: CFD velocity field',
             'files': ['fig_demo_case2_flow']},
        ],
    },
}


def make_subfigure_block(label, caption_text):
    """Create a figure with individual high-res panels using subfigure."""
    cfg = SUBFIGURE_DEFS[label]

    lines = []
    lines.append('\\begin{figure*}[htbp]')
    lines.append('\\centering')

    if label == 'bench-cfd-showcase':
        # 2x2 grid, each cell has geometry + flow side by side
        for i, panel in enumerate(cfg['panels']):
            if i == 2:
                lines.append('\\\\[8pt]')  # row break
            lines.append('\\begin{minipage}[t]{0.48\\textwidth}')
            lines.append('\\centering')
            lines.append('\\textbf{' + panel['label'] + '\\quad ' + panel['title'] + '}\\\\[4pt]')
            geo_file = panel['files'][0]
            flow_file = panel['files'][1]
            lines.append('\\includegraphics[width=0.48\\textwidth]{' + geo_file + '}\\hfill')
            lines.append('\\includegraphics[width=0.48\\textwidth]{' + flow_file + '}')
            lines.append('\\end{minipage}')
            if i % 2 == 0:
                lines.append('\\hfill')
    elif label == 'discuss-structure-cfd-gap':
        # Use the composite PDF for geometry panels (a,c) and individual PNGs for CFD (b,d)
        lines.append('\\includegraphics[width=\\textwidth]{fig_discuss_structure_cfd_gap}')
    elif label == 'demo-floorplan-application':
        # Use composite PDF for floor plan panels, individual PNGs for CFD
        lines.append('\\includegraphics[width=\\textwidth]{fig_demo_floorplan_application}')

    lines.append('\\caption{' + caption_text + '}')
    lines.append('\\label{fig:' + label + '}')
    lines.append('\\end{figure*}')
    return '\n'.join(lines)


def convert_figure_caption_block(label, caption_text):
    """Create a figure environment from a caption block."""
    # Check if this figure should use subfigure layout
    if label in SUBFIGURE_DEFS:
        return make_subfigure_block(label, caption_text)

    filename = label_from_fig_ref(label)
    return (
        '\\begin{figure}[htbp]\n'
        '\\centering\n'
        '\\includegraphics[width=\\textwidth]{' + filename + '}\n'
        '\\caption{' + caption_text + '}\n'
        '\\label{fig:' + label + '}\n'
        '\\end{figure}'
    )


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def convert_document(md_path, tex_path):
    md_text = Path(md_path).read_text(encoding='utf-8')
    lines = md_text.split('\n')

    # Extract title (first # heading)
    title = ''
    for line in lines:
        if line.startswith('# ') and not line.startswith('## '):
            title = line[2:].strip()
            break

    # Extract abstract
    abstract = ''
    in_abstract = False
    abstract_lines = []
    for line in lines:
        if line.strip() == '## Abstract':
            in_abstract = True
            continue
        if in_abstract:
            if line.startswith('## ') or (line.strip() == '---' and abstract_lines):
                break
            if line.strip() != '---':
                abstract_lines.append(line)
    abstract = '\n'.join(abstract_lines).strip()

    # Process body (after abstract, before figure captions or references)
    body_lines = []
    in_body = False
    in_figure_captions = False
    figure_captions = {}  # label -> caption text
    skip_abstract = True

    i = 0
    while i < len(lines):
        line = lines[i]

        # Skip everything before first body section
        if skip_abstract:
            if line.startswith('## 1.'):
                skip_abstract = False
            else:
                i += 1
                continue

        # Detect figure captions section
        if line.strip() == '## Figure Captions':
            in_figure_captions = True
            i += 1
            continue

        if in_figure_captions:
            # Parse [Fig:label] caption text
            fig_match = re.match(r'\[Fig:([^\]]+)\]\s*(.*)', line)
            if fig_match:
                label = fig_match.group(1)
                caption = fig_match.group(2).strip()
                figure_captions[label] = caption
            i += 1
            continue

        # Skip references placeholder
        if line.strip() in ('## References', '*[To be populated in the final manuscript.]*'):
            i += 1
            continue

        # Skip horizontal rules
        if line.strip() == '---':
            i += 1
            continue

        body_lines.append(line)
        i += 1

    # Now process body_lines into LaTeX
    output = []
    i = 0
    while i < len(body_lines):
        line = body_lines[i]

        # Section headings
        heading_match = re.match(r'^(#{2,4})\s+(?:\d+\.?\d*\.?\d*\.?\s*)?(.+)$', line)
        if heading_match:
            level = len(heading_match.group(1))
            heading_text = heading_match.group(2).strip()
            # Generate section label from heading
            sec_label = re.sub(r'[^a-z0-9]+', '-', heading_text.lower()).strip('-')

            # Unnumbered sections (journal backmatter)
            UNNUMBERED = {'Highlights', 'CRediT authorship contribution statement',
                          'Declaration of competing interest', 'Data availability',
                          'Acknowledgements', 'References'}
            # Skip Highlights entirely (separate submission file for Elsevier)
            if heading_text == 'Highlights':
                i += 1
                # Skip everything until next ## heading
                while i < len(body_lines) and not re.match(r'^#{2}\s', body_lines[i]):
                    i += 1
                continue
            elif heading_text in UNNUMBERED:
                if level == 2:
                    output.append('\\section*{' + heading_text + '}')
                elif level == 3:
                    output.append('\\subsection*{' + heading_text + '}')
            elif level == 2:
                output.append('\\section{' + heading_text + '}')
                output.append('\\label{sec:' + sec_label + '}')
            elif level == 3:
                output.append('\\subsection{' + heading_text + '}')
                output.append('\\label{sec:' + sec_label + '}')
            elif level == 4:
                output.append('\\subsubsection{' + heading_text + '}')
            output.append('')
            i += 1
            continue

        # Table with caption: [Table:label] Caption
        table_cap_match = re.match(r'\[Table:([^\]]+)\]\s*(.*)', line)
        if table_cap_match:
            tab_label = table_cap_match.group(1)
            tab_caption = table_cap_match.group(2).strip()
            # Skip blank line, then parse table
            j = i + 1
            while j < len(body_lines) and body_lines[j].strip() == '':
                j += 1
            if j < len(body_lines) and '|' in body_lines[j]:
                tabular, consumed = convert_markdown_table(body_lines[j:])
                if tabular:
                    output.append('\\begin{table}[htbp]')
                    output.append('\\centering')
                    output.append('\\caption{' + tab_caption + '}')
                    output.append('\\label{tab:' + tab_label + '}')
                    output.append(tabular)
                    output.append('\\end{table}')
                    output.append('')
                    i = j + consumed
                    continue
            # If no table found, output as text
            i += 1
            continue

        # Standalone markdown table (no caption prefix)
        if '|' in line and i + 1 < len(body_lines) and re.match(r'\s*\|[-:\s|]+\|\s*$', body_lines[i + 1]):
            tabular, consumed = convert_markdown_table(body_lines[i:])
            if tabular:
                output.append(tabular)
                output.append('')
                i += consumed
                continue

        # Figure caption inline: [Fig:label] text
        fig_cap_match = re.match(r'^\[Fig:([^\]]+)\]\s*(.*)', line)
        if fig_cap_match:
            label = fig_cap_match.group(1)
            caption = fig_cap_match.group(2).strip()
            output.append(convert_figure_caption_block(label, caption))
            output.append('')
            i += 1
            continue

        # Regular paragraph line
        output.append(line)
        i += 1

    # Join and apply inline conversions
    body_text = '\n'.join(output)

    # Convert equations first (before other transformations mess with $$)
    body_text = convert_equations(body_text)

    # Convert references
    body_text = convert_citations(body_text)
    body_text = convert_figure_refs(body_text)
    body_text = convert_table_refs(body_text)
    body_text = convert_equation_refs(body_text)
    body_text = convert_section_refs(body_text)

    # Convert inline formatting
    body_text = convert_inline_formatting(body_text)

    # Escape % and convert Unicode
    body_text = escape_percent(body_text)
    body_text = convert_unicode(body_text)

    # Escape _ outside math/commands
    body_text = escape_underscores_outside_math(body_text)

    # Escape & in text lines (not in table rows)
    result_lines = []
    for line in body_text.split('\n'):
        result_lines.append(escape_ampersand_in_text(line))
    body_text = '\n'.join(result_lines)

    # Insert figure environments from figure_captions near first reference
    # Strategy: find first Fig.~\ref{fig:LABEL} occurrence, insert figure
    # block after the end of that paragraph (next blank line)
    if figure_captions:
        lines_out = body_text.split('\n')
        inserted = set()
        final_lines = []
        pending_figs = []  # figures to insert at next paragraph break

        for idx, line in enumerate(lines_out):
            final_lines.append(line)

            # Check if this line references any un-inserted figure
            for label, caption in figure_captions.items():
                if label in inserted:
                    continue
                ref_pattern = r'\\ref\{fig:' + re.escape(label) + r'\}'
                if re.search(ref_pattern, line):
                    pending_figs.append(label)
                    inserted.add(label)

            # At paragraph break (blank line), insert any pending figures
            if line.strip() == '' and pending_figs:
                for fig_label in pending_figs:
                    fig_caption = figure_captions[fig_label]
                    # Apply same escaping to caption
                    fig_caption = escape_percent(fig_caption)
                    fig_caption = convert_unicode(fig_caption)
                    fig_caption = convert_inline_formatting(fig_caption)
                    fig_caption = escape_underscores_outside_math(fig_caption)
                    fig_caption = convert_figure_refs_caption(fig_caption)
                    fig_caption = convert_table_refs(fig_caption)
                    filename = label_from_fig_ref(fig_label)
                    final_lines.append('\\begin{figure*}[htbp]')
                    final_lines.append('\\centering')
                    final_lines.append('\\includegraphics[width=\\textwidth]{' + filename + '}')
                    final_lines.append('\\caption{' + fig_caption + '}')
                    final_lines.append('\\label{fig:' + fig_label + '}')
                    final_lines.append('\\end{figure*}')
                    final_lines.append('')
                pending_figs = []

        # Any remaining un-inserted figures go at the end
        for label, caption in figure_captions.items():
            if label not in inserted:
                fig_caption = escape_percent(caption)
                fig_caption = convert_unicode(fig_caption)
                fig_caption = convert_inline_formatting(fig_caption)
                fig_caption = escape_underscores_outside_math(fig_caption)
                filename = label_from_fig_ref(label)
                final_lines.append('\\begin{figure*}[htbp]')
                final_lines.append('\\centering')
                final_lines.append('\\includegraphics[width=\\textwidth]{' + filename + '}')
                final_lines.append('\\caption{' + fig_caption + '}')
                final_lines.append('\\label{fig:' + label + '}')
                final_lines.append('\\end{figure*}')
                final_lines.append('')

        body_text = '\n'.join(final_lines)

    # Also convert abstract
    abstract = convert_citations(abstract)
    abstract = convert_inline_formatting(abstract)
    abstract = escape_percent(abstract)
    abstract = convert_unicode(abstract)
    abstract = escape_underscores_outside_math(abstract)

    # Assemble full document
    preamble = PREAMBLE % {'title': title, 'abstract': abstract, 'figdir': FIGURE_DIR}
    postamble = POSTAMBLE % {'bibfile': BIB_FILE}

    full_tex = preamble + body_text + '\n' + postamble

    # Write output
    out_path = Path(tex_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(full_tex, encoding='utf-8')
    print(f'Wrote {out_path} ({len(full_tex)} chars)')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f'Usage: {sys.argv[0]} input.md output.tex')
        sys.exit(1)
    convert_document(sys.argv[1], sys.argv[2])
