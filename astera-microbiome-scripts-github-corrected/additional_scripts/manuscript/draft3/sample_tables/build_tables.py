import csv, json
from pathlib import Path
from collections import defaultdict, Counter
from zipfile import ZipFile, ZIP_DEFLATED
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
rows = list(csv.DictReader((ROOT/'manuscript/draft3/audit/sample_manifest_with_read_totals.csv').open()))
metadata = {}
for marker in ['16S', 'ITS', 'AMF']:
    records = list(csv.DictReader((ROOT/f'manuscript/audit/{marker}_metadata.csv').open()))
    ids = json.loads((ROOT/f'manuscript/audit/{marker}_matrix_metadata.json').read_text())['dimnames'][0]
    assert len(records) == len(ids)
    metadata[marker] = dict(zip(ids, records))
by_uid = defaultdict(list)
for row in rows:
    by_uid[row['soil_uid']].append(row)
records = []
for uid, group in sorted(by_uid.items(), key=lambda x: int(x[0].split('.')[1])):
    assert len({(r['site'], r['Depth_cm'], r['sample_date'], r['sample_id']) for r in group}) == 1
    r = group[0]
    carbon = {}
    for item in group:
        m = metadata[item['marker']][item['sample_id']]
        assert (m.get('UID_code') or m.get('Sample_UID_code')) == uid
        carbon[item['marker']] = m['C']
    conflicting = max(map(float,carbon.values())) - min(map(float,carbon.values())) > 1e-8
    if conflicting:
        display = f"16S/ITS: {float(carbon['16S']):.4f}\nAMF: {float(carbon['AMF']):.4f} ‡"
    else:
        display = f'{float(next(iter(carbon.values()))):.4f}'
    records.append(dict(soil_uid=uid,sample_name=r['sample_id'],site=r['site'],depth_cm=r['Depth_cm'],date=r['sample_date'],C_display=display,carbon_conflict=conflicting,**{f'C_{m}':carbon.get(m,'') for m in ['16S','ITS','AMF']}))
assert len(records)==45
assert sum(r['carbon_conflict'] for r in records)==20
counts=Counter((r['site'],r['depth_cm'],r['date']) for r in records)
assert sum(counts.values())==45
with (OUT/'sample_details.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
with (OUT/'site_summary.csv').open('w',newline='') as f:
    writer=csv.writer(f);writer.writerow(['Site','Depth (cm)','Date','# samples'])
    for (site,depth,date),n in sorted(counts.items()):writer.writerow([site,depth,date,n])
NS='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
def p(text='',bold=False,size=20,page_break=False):
    runs=''
    for i,line in enumerate(str(text).split('\n')):
        if i:runs+='<w:r><w:br/></w:r>'
        runs+='<w:r><w:rPr>'+('<w:b/>' if bold else '')+f'<w:sz w:val="{size}"/></w:rPr><w:t xml:space="preserve">'+escape(line)+'</w:t></w:r>'
    return '<w:p><w:pPr>'+('<w:pageBreakBefore/>' if page_break else '')+'<w:spacing w:after="80"/></w:pPr>'+runs+'</w:p>'
def table(headers,data,widths):
    borders=''.join(f'<w:{e} w:val="single" w:sz="4" w:color="A0A0A0"/>' for e in ['top','left','bottom','right','insideH','insideV'])
    s='<w:tbl><w:tblPr><w:tblW w:w="9360" w:type="dxa"/><w:tblBorders>'+borders+'</w:tblBorders><w:tblLayout w:type="fixed"/><w:tblCellMar><w:top w:w="70" w:type="dxa"/><w:left w:w="90" w:type="dxa"/><w:bottom w:w="70" w:type="dxa"/><w:right w:w="90" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblGrid>'+''.join(f'<w:gridCol w:w="{w}"/>' for w in widths)+'</w:tblGrid>'
    for i,row in enumerate([headers]+data):
        s+='<w:tr><w:trPr><w:cantSplit/>'+('<w:tblHeader/>' if i==0 else '')+'</w:trPr>'
        for value,width in zip(row,widths):
            s+=f'<w:tc><w:tcPr><w:tcW w:w="{width}" w:type="dxa"/>'+('<w:shd w:fill="E8EDF2"/>' if i==0 else '')+'</w:tcPr>'+p(value,bold=i==0,size=19)+'</w:tc>'
        s+='</w:tr>'
    return s+'</w:tbl>'
body=p('Sampling tables — Astera draft 3',True,28)
body+=p('Table 1. Retained soil samples by site code, depth, and collection date.',True,22)
data=[[site,depth.replace('-','–'),date,str(n)] for (site,depth,date),n in sorted(counts.items())]
data.append(['Total','','','45'])
body+=table(['Site','Depth (cm)','Date','# samples'],data,[2400,2100,2700,2160])
body+=p('Counts represent 45 distinct soil sample identifiers retained in at least one marker dataset; samples shared across markers are counted once. Dataset sizes are 44 (16S), 43 (ITS), and 44 (18S/AMF), with 42 samples shared across all three. These counts describe the analytical cohort, not the broader field collection.',size=19)
body+=p('Co_MT and Co_NT are the manuscript site codes. Co_MT includes both NT and till samples; it must not be treated as a synonym for tillage. Depth intervals and dates follow the workbook-linked draft 3 manifest. Dates are YYYY-MM-DD.',size=19)
body+=p('Table 2. Sample identities and recorded carbon values.',True,22,page_break=True)
body+=p('† The manuscript records carbon as C. Its fraction and units are unverified, so these values cannot yet be labeled TOC. ‡ Conflicting records are shown separately for 16S/ITS and AMF; neither is selected as the verified value.',size=19)
data=[[r['sample_name'],r['site'],r['depth_cm'].replace('-','–'),r['C_display'],r['date']] for r in records]
body+=table(['Sample Name*','Site (code)','Depth (cm)','C (units unverified)†','Date'],data,[2850,1200,1150,2410,1750])
body+=p('* Sample names reproduce the analytical identifiers. MF denotes the fresh-frozen subset. Marker coverage: cs25_MF_10_vtill occurs in 16S and AMF; cs13_MF_10_vtill in 16S and ITS; cs62_MF_30_LNT in AMF only. All remaining samples occur in all three markers.',size=19)
body+=p('Depth correction: cs85_MF_30_LNT has a workbook depth of 0–10 cm despite the name suffix. CS.23 has a workbook depth of 10–30 cm despite a conflicting 16S metadata depth. The table uses workbook-linked depth intervals without renaming samples.',size=19)
body+=p('Carbon values are rounded to four decimal places. Twenty sample identifiers have conflicting carbon records. The original control workbook contains a separate “Total Carbon (TC) (%)” field for intact cores; those values are not substituted for the manuscript C values or relabeled as TOC.',size=19)
body+=p('Data sources',True,22)
body+=p('Sample scope, dates, site codes, and depths: manuscript/draft3/audit/sample_manifest_with_read_totals.csv (linked to CFARM_Co_2022_sample_control_list_2023.05.12.xlsx). Carbon: manuscript/audit/16S_metadata.csv, ITS_metadata.csv, and AMF_metadata.csv, matched through each marker’s matrix_metadata.json sample identifiers. Interpretation: Astera_manuscript_draft3.docx and Draft3_revision_report.docx. Full-precision carbon values are preserved in sample_details.csv.',size=18)
doc='<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document xmlns:w="'+NS+'"><w:body>'+body+'<w:sectPr><w:pgSz w:w="12240" w:h="15840"/><w:pgMar w:top="1080" w:right="1440" w:bottom="1080" w:left="1440"/></w:sectPr></w:body></w:document>'
with ZipFile(ROOT/'Sample_Tables.docx') as z:files={n:z.read(n) for n in z.namelist()}
files['word/document.xml']=doc.encode()
out=OUT/'Astera_Draft3_Filled_Sample_Tables.docx'
with ZipFile(out,'w',ZIP_DEFLATED) as z:
    for name,content in files.items():
        ET.fromstring(content);z.writestr(name,content)
with ZipFile(out) as z:
    assert z.testzip() is None
    tree=ET.fromstring(z.read('word/document.xml'));tables=tree.findall('.//{'+NS+'}tbl')
    assert len(tables)==2 and len(tables[1].findall('{'+NS+'}tr'))==46
print(out)
print(f'Validated: 45 distinct samples; {len(counts)} summary rows plus total; 20 carbon conflicts; all sample joins checked.')
