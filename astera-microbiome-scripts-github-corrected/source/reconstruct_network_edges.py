"""Reconstruct original thresholded networks from shared counts; verify archived node degrees/hubs."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/network_edges';OUT.mkdir(exist_ok=True)
markers=['16S','ITS','AMF'];mat={};meta={};checks=[]
for marker in markers:
 mat[marker]=pd.read_csv(ROOT/f'data/{marker}_counts.csv.gz',index_col=0)
 meta[marker]=pd.read_csv(ROOT/f'data/{marker}_samples.csv').set_index('sample_id')
common=set.intersection(*(set(m.index) for m in mat.values()))
def select(m,maximum=250,fallback=True):
 prev=(m>0).mean();keep=np.flatnonzero(prev>=.1)
 if fallback and len(keep)<10:keep=np.argsort(-m.sum().to_numpy(),kind='stable')[:min(100,m.shape[1])]
 m=m.iloc[:,keep]
 if m.shape[1]>maximum:m=m.iloc[:,np.argsort(-m.sum().to_numpy(),kind='stable')[:maximum]]
 return m
def network(m,label,node_path,prefix=None):
 m=select(m);x=np.log(m.to_numpy()+1);x-=x.mean(1)[:,None];r,p=spearmanr(x,axis=0);i,j=np.triu_indices(m.shape[1],1);q=multipletests(p[i,j],method='fdr_bh')[1];ok=(abs(r[i,j])>=.6)&(q<.05);names=np.array([prefix+'|'+c if prefix else c for c in m.columns]);edges=pd.DataFrame({'from':names[i[ok]],'to':names[j[ok]],'weight':r[i[ok],j[ok]],'qval':q[ok]});edges['sign']=np.where(edges.weight>0,'positive','negative')
 nodes=pd.read_csv(ROOT/node_path).set_index('node');degree=pd.concat([edges['from'],edges['to']]).value_counts();assert set(degree.index)==set(nodes.index),label
 assert np.array_equal(degree.loc[nodes.index],nodes.degree),label
 hubs=degree>=np.quantile(degree,.9);assert np.array_equal(hubs.loc[nodes.index],nodes.hub),label
 edges.to_csv(OUT/f'{label}_edges.csv',index=False);checks.append({'network':label,'nodes':len(degree),'edges':len(edges),'degrees_match':True,'hubs_match':True})
for marker in markers:
 for group in ['NT','till']:
  ids=meta[marker].index[meta[marker].management_group.eq(group)]
  network(mat[marker].loc[ids],f'{marker}_{group}',f'results/{marker}/{marker}_network_{group}_nodes_taxonomy.csv',marker)
for group in ['NT','till']:
 ids=[s for s in mat['16S'].index if s in common and meta['16S'].loc[s,'management_group']==group];blocks=[]
 for marker in markers:
  block=select(mat[marker].loc[ids],120,False);block.columns=[marker+'|'+c for c in block];blocks.append(block)
 network(pd.concat(blocks,axis=1),f'cross_domain_{group}',f'results/cross_domain_network_{group}_nodes_taxonomy.csv')
pd.DataFrame(checks).to_csv(OUT/'validation.csv',index=False);print(pd.DataFrame(checks).to_string(index=False))
