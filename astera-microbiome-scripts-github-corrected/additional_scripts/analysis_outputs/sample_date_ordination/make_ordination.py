"""Run from workspace root: python analysis_outputs/sample_date_ordination/make_ordination.py"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import json, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
warnings.filterwarnings('ignore', module='rdata')
OUT=Path('analysis_outputs/sample_date_ordination')
BOOK=Path('CFARM_Co_2022_sample_control_list_2023.05.12.xlsx')
b=pd.read_excel(BOOK,sheet_name='sample control list')
b=b[b.UID_code.astype(str).str.fullmatch(r'CS\.\d+')].copy()
assert not b.UID_code.duplicated().any()
b['sample_date']=pd.to_datetime(b.date_sampled,format='%Y.%m.%d',errors='raise').dt.strftime('%Y-%m-%d')
allmeta={}; counts={}; audit=[]
for domain,file in [('16S','combined_16S_ps.rds'),('ITS','combined_its.rds'),('AMF','combined_amf.rds')]:
 info=json.loads(Path(f'manuscript/audit/{domain}_matrix_metadata.json').read_text())
 v=np.load(f'manuscript/audit/{domain}_counts.npz')['values'].astype(float)
 assert not info['taxa_are_rows'][0]
 matrix_ids=info['dimnames'][0]
 m=pd.read_csv(f'manuscript/audit/{domain}_metadata.csv')
 assert v.shape==tuple(info['dim']) and len(m)==len(matrix_ids)
 m['sample_id']=matrix_ids
 uid='UID_code' if 'UID_code' in m else 'Sample_UID_code'
 assert all(str(u).lower().replace('.', '')==re.match(r'cs[0-9]+', sid.lower()).group() for u,sid in zip(m[uid],matrix_ids))
 assert all('_'+str(f).lower()+'_' in sid.lower() for f,sid in zip(m.treatment,matrix_ids))
 m=m.rename(columns={uid:'soil_uid'}).merge(b[['UID_code','sample_date','Site_Name','land_use','Depth_cm']],left_on='soil_uid',right_on='UID_code',how='left',validate='many_to_one',sort=False)
 assert list(m.sample_id)==matrix_ids
 assert np.isfinite(v).all() and (v>=0).all()
 keep=m.sample_date.notna().to_numpy() & (v.sum(1)>0)
 audit.append(dict(marker=domain,input_samples=len(m),missing_date=int(m.sample_date.isna().sum()),zero_count_samples=int((v.sum(1)==0).sum()),included=int(keep.sum())))
 m.loc[~keep].to_csv(OUT/f'{domain}_excluded_samples.csv',index=False)
 m=m.loc[keep].reset_index(drop=True);v=v[keep];v=v/v.sum(1)[:,None]
 allmeta[domain]=m;counts[domain]=v
DATES=sorted(set(pd.concat(allmeta.values()).sample_date)); palette=plt.get_cmap('turbo')(np.linspace(.05,.95,len(DATES)));colors=dict(zip(DATES,palette));shapes={'MF':'o','MS':'^','SD':'s'}
scores=[];summaries=[]
def panel(ax,domain,fraction=None):
 m=allmeta[domain];k=np.ones(len(m),bool) if fraction is None else m.treatment.eq(fraction).to_numpy();m=m.loc[k].reset_index(drop=True);v=counts[domain][k]
 d=squareform(pdist(v,metric='braycurtis'));n=len(d)
 A=-.5*d*d;G=A-A.mean(0)[None,:]-A.mean(1)[:,None]+A.mean();ev,U=np.linalg.eigh(G)
 # Lingoes correction gives Euclidean coordinates for non-Euclidean Bray-Curtis distances.
 c=max(0.,-ev[0]);dc2=d*d+2*c;np.fill_diagonal(dc2,0)
 A=-.5*dc2;G=A-A.mean(0)[None,:]-A.mean(1)[:,None]+A.mean();e,U=np.linalg.eigh(G);order=np.argsort(e)[::-1];e=e[order];U=U[:,order];xy=U[:,:2]*np.sqrt(np.maximum(e[:2],0));pct=e[:2]/e[e>1e-10].sum()*100
 for f in sorted(m.treatment.unique()):
  ix=m.treatment.eq(f).to_numpy();ax.scatter(xy[ix,0],xy[ix,1],c=[colors[t] for t in m.loc[ix,'sample_date']],marker=shapes.get(f,'o'),s=34,alpha=.85,edgecolors='white',linewidths=.35)
 ax.set(title=f'{domain}'+(f' · {fraction}' if fraction else '')+f' (n = {n})',xlabel=f'PCoA 1 ({pct[0]:.1f}%)',ylabel=f'PCoA 2 ({pct[1]:.1f}%)');ax.axhline(0,color='#dddddd',lw=.6,zorder=0);ax.axvline(0,color='#dddddd',lw=.6,zorder=0);ax.spines[['top','right']].set_visible(False);ax.set_aspect('equal',adjustable='datalim')
 scope=fraction or 'all fractions';s=m[['sample_id','soil_uid','sample_date','treatment']].copy();s['marker']=domain;s['scope']=scope;s['PCoA1']=xy[:,0];s['PCoA2']=xy[:,1];scores.append(s)
 summaries.append(dict(marker=domain,scope=scope,n=n,dates=m.sample_date.nunique(),axis1_percent=pct[0],axis2_percent=pct[1],lingoes_constant=c,raw_negative_inertia_fraction=float(-ev[ev<0].sum()/abs(ev).sum())))
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'figure.facecolor':'white','pdf.fonttype':42})
handles=[Line2D([],[],marker='o',ls='',color=colors[d],label=d,markersize=7) for d in DATES]
for split in [False,True]:
 fig,axs=plt.subplots(3 if split else 1,3,figsize=(15,12 if split else 5.5),squeeze=False)
 for j,domain in enumerate(allmeta):
  for i,f in enumerate(['MF','MS','SD'] if split else [None]):panel(axs[i,j],domain,f)
 fig.suptitle('Beta diversity by sample date'+(' — separate sample fractions' if split else ''),fontsize=18,y=.985)
 fig.text(.5,.94 if not split else .955,'Bray–Curtis dissimilarity on relative abundance · PCoA with Lingoes correction',ha='center',fontsize=11,color='#555555')
 hh=handles if split else handles+[Line2D([],[],marker=s,ls='',color='#555555',label=f,markersize=7) for f,s in shapes.items()]
 fig.legend(handles=hh,loc='lower center',ncol=5 if split else 7,frameon=False,bbox_to_anchor=(.5,.015),title='Sample date'+(' / sample fraction' if not split else ''))
 fig.tight_layout(rect=(.015,.13 if not split else .1,.995,.915 if not split else .94),h_pad=2)
 name='beta_diversity_by_date_by_fraction' if split else 'beta_diversity_by_sample_date'
 fig.savefig(OUT/f'{name}.png',dpi=220);fig.savefig(OUT/f'{name}.pdf');plt.close(fig)
pd.concat(scores).to_csv(OUT/'ordination_coordinates.csv',index=False)
pd.DataFrame(summaries).to_csv(OUT/'ordination_summary.csv',index=False)
pd.DataFrame(audit).to_csv(OUT/'sample_matching_audit.csv',index=False)
m=pd.concat([m.assign(marker=d) for d,m in allmeta.items()])
m.groupby(['marker','sample_date','treatment']).size().rename('n').reset_index().to_csv(OUT/'sample_counts_by_date.csv',index=False)
pd.crosstab(m.sample_date,m.site).to_csv(OUT/'date_by_site_counts.csv')
(OUT/'README.md').write_text('''# Beta diversity by sample date\n\nDates come from the sample control list sheet of CFARM_Co_2022_sample_control_list_2023.05.12.xlsx, joined by physical soil UID. Count matrices and sample metadata use the existing manuscript/audit exports of the three combined RDS objects, generated by manuscript/source/audit_rds.py. Matrix dimensions, row sample UIDs, and fraction labels are checked against the exported sample order before joining. The RDS extraction is not rerun; these plots use the retained audit snapshot. All dated, nonzero libraries are included, without rarefaction or additional feature filtering.\n\nEach library is normalized to relative abundance, then Bray–Curtis dissimilarities are calculated. Principal coordinates analysis uses the Lingoes correction (adding twice the absolute minimum negative eigenvalue to off-diagonal squared distances). Axis percentages refer to the corrected positive eigenvalues. Each marker, and each marker/fraction subset in the supplementary figure, has its own ordination. Coordinates across panels are not directly comparable.\n\nThe main figure colors points by collection date and uses shapes for the original MF, MS, and SD sample-fraction labels. The second figure independently ordinates each fraction. Matching exclusions, coordinates, sample counts, and correction constants are saved as CSV files.\n\nThese are descriptive plots: date is associated with sampling site and land use, and multiple fractions can come from the same physical soil sample. Apparent date clustering cannot by itself establish a temporal effect. No significance test is reported.\n\nReproduce from workspace root: `python analysis_outputs/sample_date_ordination/make_ordination.py`.\n''')
print(pd.DataFrame(audit).to_string(index=False));print(pd.DataFrame(summaries).to_string(index=False));print('Dates:',DATES)
