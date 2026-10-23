import pandas as pd,numpy as np,json
from pathlib import Path
from scipy.stats import spearmanr
root=Path('manuscript/audit');meta={};report=[];discord=[]
for d in ['16S','ITS','AMF']:
 m=pd.read_csv(root/f'{d}_metadata.csv');a=json.loads((root/f'{d}_matrix_metadata.json').read_text());m.index=a['dimnames'][0];m=m[(m.treatment=='MF') & m.management.isin(['NT','vtill'])];meta[d]=m
 report.append(d+' site x management\n'+pd.crosstab(m.site,m.management).to_string()+'\ndepth x management\n'+pd.crosstab(m.depth,m.management).to_string())
 report.append('C/depth Spearman '+str(spearmanr(m.C,m.depth)))
 counts=np.load(root/f'{d}_counts.npz')['values'];ids=[a['dimnames'][0].index(i) for i in m.index];v=counts[ids];rel=v/v.sum(1)[:,None];alpha=pd.DataFrame({'Observed':(v>0).sum(1),'Shannon':-(np.where(rel>0,rel*np.log(np.maximum(rel,1e-300)),0)).sum(1),'reads':v.sum(1),'management':m.management.values},index=m.index);alpha.to_csv(root/f'{d}_sample_alpha.csv');report.append(alpha.groupby('management').median().to_string())
 res=pd.read_csv(f'analysis_outputs/{d}/{d}_maaslin2_nt_vs_till_results.csv');report.append('MaAsLin q05='+str((res.qval<=.05).sum())+' q25='+str((res.qval<=.25).sum())+' NT='+str(((res.qval<=.25)&(res.coef<0)).sum()))
for d in ['ITS','AMF']:
 ids=meta['16S'].index.intersection(meta[d].index)
 for i in ids:
  for c in ['management','depth','site','C','N','Soil.pH']:
   a,b=meta['16S'].loc[i,c],meta[d].loc[i,c]
   if (isinstance(a,str) and a!=b) or (not isinstance(a,str) and not np.isclose(a,b,equal_nan=True)):
    discord.append({'sample':i,'comparison':'16S vs '+d,'field':c,'16S_value':a,'other_value':b})
pd.DataFrame(discord).to_csv(root/'metadata_discrepancies.csv',index=False)
(root/'metadata_checks.txt').write_text('\n\n'.join(report));print('\n\n'.join(report));print('discrepancies',pd.DataFrame(discord).groupby(['comparison','field']).size())
