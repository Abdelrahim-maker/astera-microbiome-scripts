from pathlib import Path
import pandas as pd,numpy as np,json,shutil,hashlib,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from openpyxl import Workbook
R=Path('analysis_outputs');O=Path('manuscript');D=['16S','ITS','AMF'];colors=['#168a76','#d76b13'];plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
def read(p):return pd.read_csv(R/p)
def save(fig,name):
 fig.savefig(O/'figures'/f'{name}.png',bbox_inches='tight');fig.savefig(O/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(2,3,figsize=(10,6));rng=np.random.default_rng(7)
for j,d in enumerate(D):
 a=pd.read_csv(O/'audit'/f'{d}_sample_alpha.csv');p=read(f'{d}/{d}_alpha_wilcox.csv')
 for i,k in enumerate(['Observed','Shannon']):
  v=[a.loc[a.management==t,k].values for t in ['NT','vtill']];b=ax[i,j].boxplot(v,patch_artist=True,showfliers=False)
  for box,c in zip(b['boxes'],colors):box.set(facecolor=c,alpha=.25)
  for x,(vs,c) in enumerate(zip(v,colors),1):ax[i,j].scatter(x+rng.uniform(-.13,.13,len(vs)),vs,s=13,color=c,alpha=.7)
  ax[i,j].set_xticks([1,2]);ax[i,j].set_xticklabels(['NT','till']);ax[i,j].set_title(f'{chr(65+i*3+j)}  {d}: p = {p.loc[p.metric==k,"p_value"].iloc[0]:.3f}');ax[i,j].set_ylabel(k)
fig.tight_layout();save(fig,'figure1_alpha')
fig,ax=plt.subplots(figsize=(7,4));x=np.arange(3);counts=[]
for d in D:
 a=read(f'{d}/{d}_maaslin2_nt_vs_till_results.csv');counts.append([(a.qval<=.05).sum(),((a.qval>.05)&(a.qval<=.25)).sum()])
v=np.array(counts);ax.bar(x,v[:,0],color='#334e68',label='q ≤ 0.05');ax.bar(x,v[:,1],bottom=v[:,0],color='#a6bbc8',label='0.05 < q ≤ 0.25');ax.set_xticks(x);ax.set_xticklabels(D);ax.set_ylabel('Treatment-associated features');ax.legend()
for i in x:ax.text(i,v[i].sum()+.7,str(v[i].sum()),ha='center')
ax.set_ylim(0,52);fig.tight_layout();save(fig,'figure4_differential')
net=read('network_stats_summary.csv');fig,axs=plt.subplots(1,3,figsize=(10,3.6))
for ax,k,label in zip(axs,['nodes','mean_degree','density'],['Connected nodes','Mean degree','Density']):
 for i,(t,c) in enumerate(zip(['NT','till'],colors)):
  vals=[net.loc[net.label==f'{d} {t}',k].iloc[0] for d in D+['cross']];ax.bar(np.arange(4)+(i-.5)*.35,vals,width=.35,color=c,label=t)
 ax.set_xticks(np.arange(4));ax.set_xticklabels(D+['Joint']);ax.set_ylabel(label)
axs[0].legend();fig.tight_layout();save(fig,'figure5_networks')
b=read('imbalance_sensitivity/ordination_balanced_subsamples_permanova_permdisp.csv');j=read('imbalance_sensitivity/joint_balanced_subsamples_permanova.csv');fig,axs=plt.subplots(1,2,figsize=(9,4))
for ax,key,jkey,label in zip(axs,['permanova_r2','permanova_p'],['r2','p_value'],['PERMANOVA R²','Permutation p value']):
 vals=[b.loc[b.domain==d,key].values for d in D]+[j[jkey].values];ax.boxplot(vals,showfliers=False)
 for i,vv in enumerate(vals,1):ax.scatter(i+rng.uniform(-.11,.11,len(vv)),vv,s=16,color='#334e68',alpha=.7)
 ax.set_xticks([1,2,3,4]);ax.set_xticklabels(D+['Joint']);ax.set_ylabel(label)
axs[1].axhline(.05,color='#b32632',linestyle='--');fig.tight_layout();save(fig,'figure6_sensitivity')
# Copy every original figure unchanged and give each a permanent supplementary ID.
figs=[]
for i,p in enumerate(sorted(R.rglob('*.png')),1):
 dest=f'S{i:02d}_{p.stem}.png';shutil.copy2(p,O/'figures'/dest);figs.append({'number':i,'source':str(p),'file':'figures/'+dest})
(O/'source'/'figure_index.json').write_text(json.dumps(figs,indent=2))
# Preserve every CSV/TSV and workbook; one indexed workbook for all CSV outputs.
wb=Workbook();idx=wb.active;idx.title='Index';idx.append(['Table ID','Source','Sheet','Rows','Columns'])
manifest=[]
for i,p in enumerate(sorted(R.rglob('*.csv')),1):
 data=pd.read_csv(p);sid=f'D{i:03d}';sheet=wb.create_sheet(sid);sheet.append(list(data.columns))
 for row in data.itertuples(index=False,name=None):sheet.append([None if pd.isna(z) else z for z in row])
 sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
 for col in sheet.columns:sheet.column_dimensions[col[0].column_letter].width=min(55,max(13,len(str(col[0].value))+2))
 idx.append([sid,str(p),sid,len(data),len(data.columns)]);manifest.append({'id':sid,'source':str(p),'rows':len(data),'columns':len(data.columns)})
idx.freeze_panes='A2';idx.column_dimensions['B'].width=100;wb.save(O/'supplementary_tables.xlsx')
for p in R.rglob('*'):
 if p.is_file() and p.suffix in ['.csv','.tsv','.xlsx']:
  dest=O/'supplementary_data'/p.relative_to(R);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
(O/'source'/'data_index.json').write_text(json.dumps(manifest,indent=2))
inv=[]
for p in sorted(Path('.').glob('*')):
 if p.is_file():inv.append({'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for p in sorted(R.rglob('*')):
 if p.is_file():inv.append({'file':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(inv).to_csv(O/'audit'/'input_inventory.csv',index=False)
# Merge authoritative DOI metadata saved by both retrieval passes.
refs=[]
for p in sorted((O/'audit'/'reference_metadata').glob('*.json')):
 data=json.loads(p.read_text());k=p.stem
 def one(x):return x[0] if isinstance(x,list) else x
 title=re.sub('<[^>]+>','',one(data['title']));title=re.sub(r'\s+',' ',title).strip()
 if k=='Cotrufo2013':title='The Microbial Efficiency-Matrix Stabilization (MEMS) framework integrates plant litter decomposition with soil organic matter stabilization: do labile plant inputs form stable soil organic matter?'
 year=data.get('published',data.get('issued'))['date-parts'][0][0]
 if k=='Brewer2017':year=2017
 refs.append({'id':k,'type':'article-journal','title':title,'author':[{'family':a.get('family',''),'given':a.get('given','')} for a in data.get('author',[])],'issued':{'date-parts':[[year]]},'container-title':one(data.get('container-title','')),'volume':data.get('volume',''),'issue':data.get('issue',''),'page':data.get('page',data.get('article-number','')),'DOI':data.get('DOI',data.get('doi',''))})
for k,title,authors,year,journal,volume,page,doi in [
 ('Rottjers2018','From hairballs to hypotheses–biological insights from microbial networks',[('Röttjers','Lisa'),('Faust','Karoline')],2018,'FEMS Microbiology Reviews','42','761–780','10.1093/femsre/fuy030'),
 ('Kurtz2015','Sparse and Compositionally Robust Inference of Microbial Ecological Networks',[('Kurtz','Zachary D.'),('Müller','Christian L.'),('Miraldi','Emily R.'),('Littman','Dan R.'),('Blaser','Martin J.'),('Bonneau','Richard A.')],2015,'PLOS Computational Biology','11','e1004226','10.1371/journal.pcbi.1004226')]:
 if k not in [r['id'] for r in refs]:refs.append({'id':k,'type':'article-journal','title':title,'author':[{'family':f,'given':g} for f,g in authors],'issued':{'date-parts':[[year]]},'container-title':journal,'volume':volume,'page':page,'DOI':doi})
(O/'source'/'references.json').write_text(json.dumps(refs,indent=2,ensure_ascii=False));print('Prepared',len(figs),'original figures,',len(manifest),'data tables,',len(refs),'references')
