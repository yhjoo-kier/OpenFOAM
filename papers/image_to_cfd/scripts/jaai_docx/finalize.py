# -*- coding: utf-8 -*-
import zipfile, re
src='build1.docx'; dst='Image-to-CFD_JAAI.docx'
zin=zipfile.ZipFile(src); names=zin.namelist()
data={n:zin.read(n) for n in names}; zin.close()

PAGE_RUN=('<w:r><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr><w:fldChar w:fldCharType="begin"/></w:r>'
'<w:r><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="18"/></w:rPr><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
'<w:r><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="18"/></w:rPr><w:fldChar w:fldCharType="end"/></w:r>')

def footer_xml():
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'+PAGE_RUN+'</w:p></w:ftr>').encode()

def running_header_xml():
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:hdr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:i/><w:sz w:val="16"/></w:rPr>'
    '<w:t>J. Korean Appl. Artif. Intell.</w:t></w:r></w:p></w:hdr>').encode()

data['word/footer2.xml']=footer_xml()
data['word/footer3.xml']=footer_xml()
data['word/header1.xml']=running_header_xml()

st=data['word/styles.xml'].decode()
st=re.sub(r'<w:numPr>.*?</w:numPr>','',st,flags=re.S)
data['word/styles.xml']=st.encode()

s=data['word/settings.xml'].decode()
s=re.sub(r'<w:evenAndOddHeaders\s*/>','',s)
data['word/settings.xml']=s.encode()

zout=zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED)
for n in names: zout.writestr(n, data[n])
zout.close()
print('finalized', dst)
