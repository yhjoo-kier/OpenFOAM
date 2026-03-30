#!/usr/bin/env python3
"""
Citation Verification Script
============================
.bbl 파일에서 참고문헌을 파싱하여 DOI resolve + CrossRef/Semantic Scholar API로
실제 존재 여부를 검증하고 근거 표를 생성한다.

Usage:
    python verify_citations.py [--bbl path/to/paper.bbl] [--output results.md]
"""

import argparse
import json
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Citation:
    """Parsed citation from .bbl file."""
    key: str
    raw_text: str
    doi: str = ""
    title: str = ""
    authors: str = ""
    year: str = ""
    journal: str = ""


@dataclass
class VerificationResult:
    """Result of verifying one citation."""
    citation: Citation
    doi_resolves: bool | None = None          # None = no DOI
    doi_status_code: int | None = None
    doi_redirect_url: str = ""
    crossref_found: bool | None = None        # None = no DOI to query
    crossref_title: str = ""
    crossref_authors: str = ""
    crossref_year: str = ""
    crossref_journal: str = ""
    crossref_type: str = ""
    title_match: bool | None = None
    year_match: bool | None = None
    semantic_scholar_found: bool | None = None
    semantic_scholar_title: str = ""
    semantic_scholar_year: str = ""
    semantic_scholar_venue: str = ""
    semantic_scholar_doi: str = ""
    verdict: str = ""                         # VERIFIED / LIKELY_REAL / UNVERIFIED / SUSPICIOUS
    evidence_summary: str = ""
    errors: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# .bbl parser
# ---------------------------------------------------------------------------

def parse_bbl(bbl_path: Path) -> list[Citation]:
    """Parse a .bbl file and extract citation entries."""
    text = bbl_path.read_text(encoding="utf-8")

    entries = []
    # Split by \bibitem
    parts = re.split(r"\\bibitem\{", text)
    for part in parts[1:]:  # skip preamble
        # Extract key
        key_match = re.match(r"([^}]+)\}", part)
        if not key_match:
            continue
        key = key_match.group(1)
        body = part[key_match.end():]

        # Clean LaTeX commands for text extraction
        clean = _clean_latex(body)

        # Extract DOI
        doi = ""
        doi_match = re.search(r"doi:(10\.[^\s}]+)", body)
        if doi_match:
            doi = doi_match.group(1).rstrip(".")

        # Extract year — look for (YYYY) pattern
        year = ""
        year_match = re.search(r"\((\d{4})\)", clean)
        if year_match:
            year = year_match.group(1)
        else:
            # fallback: any 4-digit year between 1900-2099
            year_match2 = re.search(r"\b(19\d{2}|20\d{2})\b", clean)
            if year_match2:
                year = year_match2.group(1)

        # Extract title — typically the text after authors, before journal
        # In .bbl, title is usually the first sentence after author names
        title = _extract_title(clean)

        # Extract first author(s)
        authors = _extract_authors(clean)

        entries.append(Citation(
            key=key,
            raw_text=clean.strip(),
            doi=doi,
            title=title,
            authors=authors,
            year=year,
        ))

    return entries


def _clean_latex(text: str) -> str:
    """Remove LaTeX commands, keeping text content."""
    t = text
    # Remove \href{...}{...} → keep second arg
    t = re.sub(r"\\href\s*\{[^}]*\}\s*\{([^}]*)\}", r"\1", t)
    # Remove \path{...} → keep content
    t = re.sub(r"\\path\{([^}]*)\}", r"\1", t)
    # Remove \newblock
    t = t.replace("\\newblock", "")
    # Remove common accents: {\"o} → o, etc.
    t = re.sub(r"\{\\\"(\w)\}", r"\1", t)
    t = re.sub(r"\{\\'(\w)\}", r"\1", t)
    t = re.sub(r"\{\\`(\w)\}", r"\1", t)
    t = re.sub(r"\{\\aa\}", "a", t)
    # Remove {{...}} → content
    t = re.sub(r"\{\{([^}]*)\}\}", r"\1", t)
    # Remove remaining {...}
    t = re.sub(r"\{([^}]*)\}", r"\1", t)
    # Remove ~
    t = t.replace("~", " ")
    # Collapse whitespace
    t = re.sub(r"\s+", " ", t)
    return t.strip()


def _extract_title(clean_text: str) -> str:
    """Heuristic title extraction from cleaned .bbl entry text."""
    # Title is typically after author list (ending with comma) and before journal/venue
    # Pattern: "Author1, Author2, Title here, Journal Name Vol (Year) Pages."
    # Try to get the longest capitalized phrase
    parts = clean_text.split(",")
    if len(parts) >= 3:
        # Title is usually the segment after author names
        # Authors typically have initials like "A." or "A.~B."
        for i, part in enumerate(parts):
            stripped = part.strip()
            # Skip parts that look like author names (short, have dots)
            if len(stripped) > 30 and not re.match(r"^[A-Z]\.\s", stripped):
                return stripped.split(".")[0].strip() if "." in stripped else stripped
    return ""


def _extract_authors(clean_text: str) -> str:
    """Extract author string (first ~80 chars)."""
    # Authors are at the beginning, before the title
    # Typically ends before a long phrase (the title)
    match = re.match(r"^(.{10,80}?),\s+[A-Z]", clean_text)
    if match:
        return match.group(1).strip()
    # Fallback: first 60 chars
    return clean_text[:60].strip()


# ---------------------------------------------------------------------------
# DOI resolution
# ---------------------------------------------------------------------------

def check_doi_resolves(doi: str) -> tuple[bool, int | None, str]:
    """Check if a DOI resolves via HTTP HEAD to doi.org."""
    if not doi:
        return False, None, ""
    url = f"https://doi.org/{doi}"
    req = Request(url, method="HEAD")
    req.add_header("User-Agent", "CitationVerifier/1.0 (academic-verification)")
    try:
        resp = urlopen(req, timeout=15)
        return True, resp.status, resp.url
    except HTTPError as e:
        if e.code in (301, 302, 303, 307, 308):
            return True, e.code, e.headers.get("Location", "")
        return False, e.code, ""
    except (URLError, TimeoutError):
        return False, None, ""


# ---------------------------------------------------------------------------
# CrossRef API
# ---------------------------------------------------------------------------

def query_crossref(doi: str) -> dict | None:
    """Fetch metadata from CrossRef API for a given DOI."""
    if not doi:
        return None
    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    req = Request(url)
    req.add_header("User-Agent", "CitationVerifier/1.0 (mailto:verify@example.com)")
    try:
        resp = urlopen(req, timeout=15)
        data = json.loads(resp.read().decode("utf-8"))
        return data.get("message", {})
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return None


def extract_crossref_metadata(cr: dict) -> dict:
    """Extract relevant fields from CrossRef response."""
    title = cr.get("title", [""])[0] if cr.get("title") else ""
    authors_list = cr.get("author", [])
    authors = ", ".join(
        f"{a.get('family', '')}, {a.get('given', '')}" for a in authors_list[:3]
    )
    if len(authors_list) > 3:
        authors += " et al."

    # Year from published-print or published-online or issued
    year = ""
    for date_field in ("published-print", "published-online", "issued"):
        parts = cr.get(date_field, {}).get("date-parts", [[]])
        if parts and parts[0] and parts[0][0]:
            year = str(parts[0][0])
            break

    journal = ""
    for jfield in ("container-title", "short-container-title"):
        jnames = cr.get(jfield, [])
        if jnames:
            journal = jnames[0]
            break

    return {
        "title": title,
        "authors": authors,
        "year": year,
        "journal": journal,
        "type": cr.get("type", ""),
    }


# ---------------------------------------------------------------------------
# Semantic Scholar API (fallback for DOI-less entries)
# ---------------------------------------------------------------------------

def query_semantic_scholar_by_title(title: str) -> dict | None:
    """Search Semantic Scholar by title."""
    if not title or len(title) < 10:
        return None
    url = f"https://api.semanticscholar.org/graph/v1/paper/search?query={quote(title)}&limit=3&fields=title,year,venue,externalIds"
    req = Request(url)
    req.add_header("User-Agent", "CitationVerifier/1.0")
    try:
        resp = urlopen(req, timeout=15)
        data = json.loads(resp.read().decode("utf-8"))
        papers = data.get("data", [])
        if papers:
            return papers[0]
        return None
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return None


def query_semantic_scholar_by_doi(doi: str) -> dict | None:
    """Look up a paper on Semantic Scholar by DOI."""
    if not doi:
        return None
    url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{quote(doi, safe='')}?fields=title,year,venue,externalIds"
    req = Request(url)
    req.add_header("User-Agent", "CitationVerifier/1.0")
    try:
        resp = urlopen(req, timeout=15)
        data = json.loads(resp.read().decode("utf-8"))
        return data
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return None


# ---------------------------------------------------------------------------
# Fuzzy title matching
# ---------------------------------------------------------------------------

def normalize_title(title: str) -> str:
    """Normalize title for comparison."""
    t = title.lower()
    t = re.sub(r"[^a-z0-9\s]", "", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def titles_match(t1: str, t2: str) -> bool:
    """Check if two titles are similar enough."""
    n1, n2 = normalize_title(t1), normalize_title(t2)
    if not n1 or not n2:
        return False
    # Check if one contains the other (handles truncation)
    if n1 in n2 or n2 in n1:
        return True
    # Word overlap ratio
    words1, words2 = set(n1.split()), set(n2.split())
    if not words1 or not words2:
        return False
    overlap = len(words1 & words2) / max(len(words1), len(words2))
    return overlap > 0.6


# ---------------------------------------------------------------------------
# Verification logic
# ---------------------------------------------------------------------------

def verify_citation(cit: Citation) -> VerificationResult:
    """Verify a single citation using all available methods."""
    result = VerificationResult(citation=cit)
    evidence_parts = []

    # --- Step 1: DOI resolution ---
    if cit.doi:
        ok, status, redir = check_doi_resolves(cit.doi)
        result.doi_resolves = ok
        result.doi_status_code = status
        result.doi_redirect_url = redir
        if ok:
            evidence_parts.append(f"DOI resolves (HTTP {status})")
        else:
            evidence_parts.append(f"DOI FAILED (HTTP {status})")
        time.sleep(0.3)  # rate limit

        # --- Step 2: CrossRef metadata ---
        cr = query_crossref(cit.doi)
        if cr:
            meta = extract_crossref_metadata(cr)
            result.crossref_found = True
            result.crossref_title = meta["title"]
            result.crossref_authors = meta["authors"]
            result.crossref_year = meta["year"]
            result.crossref_journal = meta["journal"]
            result.crossref_type = meta["type"]

            # Title match
            result.title_match = titles_match(cit.title, meta["title"]) if cit.title else None
            if result.title_match:
                evidence_parts.append("CrossRef title matches")
            elif result.title_match is False:
                evidence_parts.append("CrossRef title MISMATCH")

            # Year match
            if cit.year and meta["year"]:
                result.year_match = cit.year == meta["year"]
                if result.year_match:
                    evidence_parts.append(f"Year confirmed ({meta['year']})")
                else:
                    evidence_parts.append(f"Year mismatch: bbl={cit.year} vs CR={meta['year']}")

            if meta["journal"]:
                evidence_parts.append(f"Journal: {meta['journal']}")
            if meta["type"]:
                evidence_parts.append(f"Type: {meta['type']}")
        else:
            result.crossref_found = False
            evidence_parts.append("CrossRef: no data")
        time.sleep(0.3)

    else:
        result.doi_resolves = None
        result.crossref_found = None
        evidence_parts.append("No DOI in .bbl")

    # --- Step 3: Semantic Scholar (for DOI-less or as supplementary) ---
    ss = None
    if cit.doi:
        ss = query_semantic_scholar_by_doi(cit.doi)
    if not ss and cit.title:
        ss = query_semantic_scholar_by_title(cit.title)
    # For DOI-less: also try searching by key-derived title
    if not ss and not cit.doi:
        # Try raw text snippet as query
        snippet = cit.raw_text[:120]
        ss = query_semantic_scholar_by_title(snippet)

    if ss:
        result.semantic_scholar_found = True
        result.semantic_scholar_title = ss.get("title", "")
        result.semantic_scholar_year = str(ss.get("year", ""))
        result.semantic_scholar_venue = ss.get("venue", "")
        ext = ss.get("externalIds", {})
        result.semantic_scholar_doi = ext.get("DOI", "")
        evidence_parts.append(f"Semantic Scholar: found ('{result.semantic_scholar_title[:50]}...')")
    else:
        result.semantic_scholar_found = False
        evidence_parts.append("Semantic Scholar: not found")
    time.sleep(0.5)  # S2 rate limit is stricter

    # --- Verdict ---
    result.verdict, result.evidence_summary = _determine_verdict(result, evidence_parts)
    return result


def _determine_verdict(r: VerificationResult, evidence_parts: list[str]) -> tuple[str, str]:
    """Determine overall verdict based on collected evidence."""
    summary = "; ".join(evidence_parts)

    score = 0
    if r.doi_resolves:
        score += 3
    if r.crossref_found:
        score += 2
    if r.title_match:
        score += 2
    if r.year_match:
        score += 1
    if r.semantic_scholar_found:
        score += 2

    # Penalties
    if r.doi_resolves is False and r.citation.doi:
        score -= 3
    if r.title_match is False:
        score -= 2
    if r.year_match is False:
        score -= 1

    if score >= 6:
        verdict = "VERIFIED"
    elif score >= 3:
        verdict = "LIKELY_REAL"
    elif score >= 1:
        verdict = "WEAK"
    else:
        verdict = "UNVERIFIED"

    return verdict, summary


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(results: list[VerificationResult], output_path: Path) -> str:
    """Generate a detailed Markdown report."""
    lines = []
    lines.append("# Citation Verification Report")
    lines.append("")
    lines.append(f"**Generated**: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**Total citations**: {len(results)}")

    # Summary counts
    counts = {}
    for r in results:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    lines.append(f"**Verdicts**: {', '.join(f'{v}: {c}' for v, c in sorted(counts.items()))}")
    lines.append("")

    # Legend
    lines.append("## Verdict Legend")
    lines.append("")
    lines.append("| Verdict | Meaning |")
    lines.append("|---------|---------|")
    lines.append("| VERIFIED | DOI resolves + CrossRef metadata confirms title/year/journal |")
    lines.append("| LIKELY_REAL | Multiple independent sources confirm existence |")
    lines.append("| WEAK | Only partial confirmation (e.g., Semantic Scholar only) |")
    lines.append("| UNVERIFIED | Could not confirm existence through any automated method |")
    lines.append("")

    # Detailed table
    lines.append("## Detailed Results")
    lines.append("")
    lines.append("| # | Citekey | Year | DOI Resolves | CrossRef | Semantic Scholar | Verdict | Evidence |")
    lines.append("|---|---------|------|:------------:|:--------:|:----------------:|---------|----------|")

    for i, r in enumerate(results, 1):
        doi_col = _icon(r.doi_resolves)
        cr_col = _icon(r.crossref_found)
        ss_col = _icon(r.semantic_scholar_found)

        # Truncate evidence for table
        ev = r.evidence_summary
        if len(ev) > 120:
            ev = ev[:117] + "..."

        lines.append(
            f"| {i} | `{r.citation.key}` | {r.citation.year} "
            f"| {doi_col} | {cr_col} | {ss_col} "
            f"| **{r.verdict}** | {ev} |"
        )

    lines.append("")

    # Per-citation detail cards for non-VERIFIED
    flagged = [r for r in results if r.verdict != "VERIFIED"]
    if flagged:
        lines.append("## Entries Requiring Attention")
        lines.append("")
        for r in flagged:
            lines.append(f"### `{r.citation.key}` — {r.verdict}")
            lines.append("")
            lines.append(f"- **DOI**: `{r.citation.doi or 'N/A'}`")
            lines.append(f"- **Year (bbl)**: {r.citation.year}")
            lines.append(f"- **Title (bbl)**: {r.citation.title or '(extraction failed)'}")
            if r.crossref_found:
                lines.append(f"- **CrossRef title**: {r.crossref_title}")
                lines.append(f"- **CrossRef year**: {r.crossref_year}")
                lines.append(f"- **CrossRef journal**: {r.crossref_journal}")
            if r.semantic_scholar_found:
                lines.append(f"- **S2 title**: {r.semantic_scholar_title}")
                lines.append(f"- **S2 year**: {r.semantic_scholar_year}")
                lines.append(f"- **S2 venue**: {r.semantic_scholar_venue}")
            if r.errors:
                lines.append(f"- **Errors**: {'; '.join(r.errors)}")
            lines.append(f"- **Evidence**: {r.evidence_summary}")
            lines.append("")

    # CrossRef metadata comparison table
    lines.append("## CrossRef Metadata Comparison")
    lines.append("")
    lines.append("| # | Citekey | CrossRef Title | CrossRef Year | CrossRef Journal | Type |")
    lines.append("|---|---------|----------------|:-------------:|------------------|------|")
    for i, r in enumerate(results, 1):
        if r.crossref_found:
            cr_title = r.crossref_title[:60] + "..." if len(r.crossref_title) > 60 else r.crossref_title
            lines.append(
                f"| {i} | `{r.citation.key}` | {cr_title} "
                f"| {r.crossref_year} | {r.crossref_journal} | {r.crossref_type} |"
            )
        else:
            lines.append(f"| {i} | `{r.citation.key}` | — | — | — | — |")

    lines.append("")

    report = "\n".join(lines)
    output_path.write_text(report, encoding="utf-8")
    return report


def _icon(val: bool | None) -> str:
    """Return a text icon for boolean/None."""
    if val is True:
        return "YES"
    elif val is False:
        return "NO"
    return "N/A"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Verify citations in a .bbl file")
    parser.add_argument(
        "--bbl",
        type=Path,
        default=Path(__file__).parent.parent / "latex" / "paper.bbl",
        help="Path to .bbl file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent.parent / "results" / "citation_verification_report.md",
        help="Output report path",
    )
    args = parser.parse_args()

    if not args.bbl.exists():
        print(f"ERROR: .bbl file not found: {args.bbl}", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing {args.bbl}...")
    citations = parse_bbl(args.bbl)
    print(f"Found {len(citations)} citations.")

    print("\nVerifying citations (this may take 1-2 minutes due to API rate limits)...\n")
    results = []
    for i, cit in enumerate(citations, 1):
        print(f"  [{i:2d}/{len(citations)}] {cit.key} ... ", end="", flush=True)
        r = verify_citation(cit)
        print(r.verdict)
        results.append(r)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = generate_report(results, args.output)

    # Print summary to stdout
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    counts = {}
    for r in results:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    for v, c in sorted(counts.items()):
        print(f"  {v}: {c}")
    print(f"\nFull report: {args.output}")

    flagged = [r for r in results if r.verdict not in ("VERIFIED", "LIKELY_REAL")]
    if flagged:
        print(f"\n⚠ {len(flagged)} citation(s) need manual review:")
        for r in flagged:
            print(f"  - {r.citation.key}: {r.verdict}")

    sys.exit(0)


if __name__ == "__main__":
    main()
