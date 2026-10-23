import rdata,warnings,pandas as pd,numpy as np,json
from pathlib import Path
warnings.filterwarnings('ignore')
allmeta={}; orders={}
for domain,file in [('16S','combined_16S_ps.rds'),('ITS','combined_its.rds'),('AMF','combined_amf.rds')]:
 x=rdata.parser.parse_file(file); a=x.object.attributes;slots={}
 while a and a.value:
  if a.tag and a.tag.value: slots[a.tag.value.value.decode()]=a.value[0]
  a=a.value[1]
 z=rdata.conversion.convert(slots['sam_data']);m=pd.DataFrame(z)
 print(domain,'columns',list(m.columns)); print(m.loc[m['treatment'].str.lower()=='mf'].to_string(index=False)[:1600])
 m.to_csv(f'manuscript/audit/{domain}_metadata.csv',index=False)
 # sample rownames are attributes of sam_data
 a=slots['sam_data'].attributes;attrs={}
 while a and a.value:
  if a.tag and a.tag.value:attrs[a.tag.value.value.decode()]=rdata.conversion.convert(a.value[0])
  a=a.value[1]
 m.index=attrs['row.names'];mf=m[m.treatment.str.lower()=='mf'];orders[domain]=list(mf.index);allmeta[domain]=mf
 print('MF categories', {c:mf[c].value_counts().to_dict() for c in m.columns if mf[c].nunique()<12})
 otu=rdata.conversion.convert(slots['otu_table']);print('otu',type(otu),otu.shape)
 np.savez_compressed(f'manuscript/audit/{domain}_counts.npz',values=np.asarray(otu))
 # inspect dimnames
 a=slots['otu_table'].attributes;att={}
 while a and a.value:
  if a.tag and a.tag.value:att[a.tag.value.value.decode()]=rdata.conversion.convert(a.value[0])
  a=a.value[1]
 print('matrix attrs',[(k,str(v)[:180]) for k,v in att.items()])
 Path(f'manuscript/audit/{domain}_matrix_metadata.json').write_text(json.dumps({k:v.tolist() if hasattr(v,'tolist') else v for k,v in att.items()},default=lambda v:v.tolist() if hasattr(v,'tolist') else str(v)))
common=set.intersection(*(set(v) for v in orders.values()));print('common',len(common),'order matches',all([i for i in orders[d] if i in common]==[i for i in orders['16S'] if i in common] for d in orders));print(allmeta['16S'].loc[list(common),'management'].value_counts())
Path('manuscript/audit/sample_order.json').write_text(json.dumps(orders))
