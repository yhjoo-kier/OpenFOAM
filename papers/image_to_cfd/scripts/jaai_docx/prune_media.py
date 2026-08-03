# -*- coding: utf-8 -*-
"""Remove template media not referenced by document.xml, and rename output.
Usage: python prune_media.py <output_name.docx>
Reads Image-to-CFD_JAAI.docx (finalize output), prunes, writes <output_name>.docx.
"""
import zipfile, re, os, sys

src = 'Image-to-CFD_JAAI.docx'
dst = sys.argv[1]

z = zipfile.ZipFile(src); names = z.namelist()
data = {n: z.read(n) for n in names}; z.close()

rels = data['word/_rels/document.xml.rels'].decode('utf-8')
doc = data['word/document.xml'].decode('utf-8')
used = set(re.findall(r'r:(?:embed|id)="(rId\d+)"', doc))
removed = []

def drop(m):
    rid, target = m.group(1), m.group(2)
    if target.startswith('media/') and rid not in used:
        removed.append('word/' + target)
        return ''
    return m.group(0)

rels = re.sub(r'<Relationship Id="(rId\d+)" Type="[^"]*/image" Target="([^"]+)"\s*/>', drop, rels)
data['word/_rels/document.xml.rels'] = rels.encode('utf-8')

out = zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED)
for n in names:
    if n in removed:
        continue
    out.writestr(n, data[n])
out.close()
print(f'{dst}: removed {len(removed)} orphan media, size {round(os.path.getsize(dst)/1e6,1)} MB')
