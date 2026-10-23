"""Reconstruct six marker-specific graphs from source counts."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reproduced/network_edges';OUT.mkdir(parents=True,exist_ok=True)
rows=[]
for marker in ['16S','ITS','AMF']:
 x=pd.read_csv(ROOT/'data'/f'{marker}_counts.csv.gz',index_col=0)
 meta=pd.read_csv(ROOT/'data'/f'{marker}_samples.csv').set_index('sample_id')
 for group in ['NT','till']:
  m=x.loc[meta.index[meta.management_group.eq(group)]]
  keep=np.flatnonzero((m>0).mean()>=.1)
  if len(keep)<10:keep=np.argsort(-m.sum().to_numpy(),kind='stable')[:min(100,m.shape[1])]
  m=m.iloc[:,keep]
  if m.shape[1]>250:m=m.iloc[:,np.argsort(-m.sum().to_numpy(),kind='stable')[:250]]
  z=np.log(m.to_numpy()+1);z-=z.mean(1)[:,None]
  r,p=spearmanr(z,axis=0);i,j=np.triu_indices(m.shape[1],1)
  q=multipletests(p[i,j],method='fdr_bh')[1];ok=(abs(r[i,j])>=.6)&(q<.05)
  names=np.array([marker+'|'+c for c in m.columns])
  e=pd.DataFrame({'from':names[i[ok]],'to':names[j[ok]],'rho':r[i[ok],j[ok]],'q_BH':q[ok]})
  e.to_csv(OUT/f'{marker}_{group}_edges.csv',index=False)
  n=len(set(e['from'])|set(e['to']));E=len(e)
  rows.append(dict(label=marker+' '+group,nodes=n,edges=E,density=2*E/n/(n-1),mean_degree=2*E/n))
a=pd.DataFrame(rows)
a.to_csv(OUT/'network_summary.csv',index=False)
print('Completed six marker-specific networks; outputs generated locally.')
