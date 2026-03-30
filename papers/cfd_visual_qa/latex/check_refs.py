import re, sys
with open('main.tex', 'r', encoding='utf-8') as f:
    text = f.read()
cites = set()
for m in re.findall(r'\\cite\{([^}]+)\}', text):
    for k in m.split(','):
        cites.add(k.strip())
bibs = set(re.findall(r'\\bibitem\{([^}]+)\}', text))
missing = cites - bibs
unused = bibs - cites
print(f'Cite keys: {len(cites)}, Bibitem keys: {len(bibs)}')
if missing:
    print(f'MISSING from bib ({len(missing)}): {sorted(missing)}')
else:
    print('All cite keys have matching bibitems.')
if unused:
    print(f'UNUSED bibitems ({len(unused)}): {sorted(unused)}')
else:
    print('All bibitems are cited.')
