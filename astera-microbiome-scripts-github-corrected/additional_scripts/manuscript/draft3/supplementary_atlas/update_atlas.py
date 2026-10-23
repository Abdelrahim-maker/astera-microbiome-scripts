from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from copy import deepcopy
from lxml import etree as E
import pandas as pd
O=Path('manuscript/draft3/supplementary_atlas');src=Path('manuscript/supplementary_material.docx');out=O/'Astera_Supplementary_Atlas_Updated.docx'
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main';ns={'w':W};q=lambda x:'{'+W+'}'+x
with ZipFile(src) as z:files={n:z.read(n) for n in z.namelist()}
root=E.fromstring(files['word/document.xml']);body=root.find('w:body',ns);original=list(body);template=deepcopy(original[6]);cap=deepcopy(original[5]);note=deepcopy(original[7]);page=deepcopy(original[8])
def paragraph(t,base=note):
 p=deepcopy(base);runs=p.findall('w:r',ns);r=deepcopy(runs[0]) if runs else E.Element(q('r'))
 for c in list(p):
  if c.tag!=q('pPr'):p.remove(c)
 for c in list(r):
  if c.tag!=q('rPr'):r.remove(c)
 E.SubElement(r,q('t')).text=t;p.append(r);return p

def table(headers,rows,widths):
 t=deepcopy(template)
 for e in list(t):
  if e.tag in [q('tr'),q('tblGrid')]:t.remove(e)
 grid=E.Element(q('tblGrid'))
 for width in widths:E.SubElement(grid,q('gridCol'),{q('w'):str(width)})
 t.insert(1,grid)
 for i,values in enumerate([headers]+rows):
  tr=E.Element(q('tr'));pr=E.SubElement(tr,q('trPr'));E.SubElement(pr,q('cantSplit'))
  if i==0:E.SubElement(pr,q('tblHeader'))
  for j,value in enumerate(values):
   cell=deepcopy(template.findall('w:tr',ns)[0 if i==0 else 1].findall('w:tc',ns)[0]);p=cell.find('w:p',ns);new=paragraph(str(value),p);cell.replace(p,new)
   cp=cell.find('w:tcPr',ns);E.SubElement(cp,q('tcW'),{q('w'):str(widths[j]),q('type'):'dxa'})
   tr.append(cell)
  t.append(tr)
 return t
samples=pd.read_csv('manuscript/draft3/sample_tables/sample_details.csv');samples=samples[samples.C_16S.notna()].copy();assert len(samples)==44
s=samples.groupby(['site','depth_cm','date']).size().reset_index(name='# samples').rename(columns={'site':'Site','depth_cm':'Depth (cm)','date':'Date'});s['Depth (cm)']=s['Depth (cm)'].str.replace('-', '–');s=s.sort_values(['Site','Date','Depth (cm)'],ascending=[False,True,True]);assert s['# samples'].sum()==44
rows=s.astype(str).values.tolist()+[['','', 'Total','44']]
rows2=[]
for _,r in samples.iterrows():
 name=r.sample_name+('*' if r.soil_uid=='CS.25' else '')
 c=r.C_16S if pd.notna(r.C_16S) else r.C_AMF
 value=f'{c:.4f}'
 if r.soil_uid=='CS.64':value='2.3'
 if r.soil_uid=='CS.67':value='2.4'
 rows2.append([name,r.site,r.depth_cm.replace('-','–'),value,r.date])
new=[paragraph('Sampling tables',cap),paragraph('Table S1a. Sampling by site code, depth, and collection date; each site was sampled from multiple locations.',cap),table(['Site','Depth (cm)','Date','# samples'],rows,[1800,1600,2800,1720]),paragraph('These sampling tables describe the 44-sample 16S cohort. Dataset sizes are 44 (16S), 43 (ITS), and 44 (18S/AMF), with 42 samples shared across all three. These counts describe the analytical cohort, not the broader field collection.'),page,paragraph('Table S1b. Sample identities and recorded carbon values.',cap),table(['Sample Name*','Site (code)','Depth (cm)','C%†','Date'],rows2,[3100,1100,1100,900,1720]),paragraph('* cs25_MF_10_vtill is the sole sample in the 16S cohort that lacks associated ITS amplicon data (sequencing failure, as reported by the author). The AMF-only sample cs62_MF_30_LNT is excluded from these tables; cs13_MF_10_vtill lacks AMF data.'),paragraph('† Carbon values follow the supplied table, using the 16S-associated records. C% is the author-supplied heading; measurement type and percent units have not been independently verified. Values 2.3 and 2.4 are retained as supplied. Sample identifiers are preserved verbatim; depth follows the workbook, including the 0–10 cm assignment for cs85_MF_30_LNT.'),paragraph('All selected crop labels are soybean. Site and depth codes are not expanded into unverified physical definitions. Management and site are partially confounded.')]
for e in original[5:8]:body.remove(e)
for i,e in enumerate(new):body.insert(5+i,e)
# Keep S2 onward numbering stable using S1a/S1b for the two replacement tables.
files['word/document.xml']=E.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
with ZipFile(out,'w',ZIP_DEFLATED) as z:
 for n,data in files.items():z.writestr(n,data)
with ZipFile(out) as z:
 assert z.testzip() is None
 for n in files:
  if n!='word/document.xml':assert z.read(n)==files[n]
newroot=E.fromstring(files['word/document.xml']);assert len(newroot.findall('.//w:body/w:tbl',ns))==13
assert 'Sample distribution across site and depth codes' not in ''.join(newroot.itertext())
# All original document elements from S2 onward remain byte-identical.
remaining=list(body)[5+len(new)+1:];assert len(remaining)==len(original[9:])
for a,b in zip(remaining,original[9:]):assert E.tostring(a)==E.tostring(b)
(O/'EDIT_NOTES.txt').write_text('Replaced original S1 with S1a and S1b, preserving S2–S12 numbering and all remaining document content and embedded media. Matched original table style, font, captions, repeated headers, and no-split rows. Restricted both tables to the44-sample 16S cohort; excluded AMF-only cs62. Co_MT surface Nov28 is6; Co_NT deep Nov21 is8. Clarified marker-specific missingness. Retained supplied carbon values and noted unverified units. Source document preserved.\n')
print(out.resolve())
