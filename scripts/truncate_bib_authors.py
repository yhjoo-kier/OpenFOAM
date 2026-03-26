#!/usr/bin/env python3
"""Truncate long author lists in .bib file for specific entries."""
import re
import sys

TRUNCATE_KEYS = [
    'openai.etal2023GPT4Technical',
    'geminiteam.etal2023GeminiFamily',
]
MAX_AUTHORS = 3


def truncate_authors(bib, citekey, max_authors=3):
    pattern = r'(@\w+\{' + re.escape(citekey) + r',.*?author\s*=\s*\{)(.*?)(\},)'
    match = re.search(pattern, bib, re.DOTALL)
    if match:
        authors = match.group(2).split(' and ')
        if len(authors) > max_authors:
            short = ' and '.join(authors[:max_authors]) + ' and others'
            bib = bib[:match.start(2)] + short + bib[match.end(2):]
            print(f'  Truncated {citekey}: {len(authors)} -> {max_authors} authors')
    return bib


def main():
    bib_path = sys.argv[1]
    bib = open(bib_path, 'r', encoding='utf-8').read()
    for key in TRUNCATE_KEYS:
        bib = truncate_authors(bib, key, MAX_AUTHORS)
    open(bib_path, 'w', encoding='utf-8').write(bib)


if __name__ == '__main__':
    main()
