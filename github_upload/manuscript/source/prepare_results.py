"""Curate retained results and figures; does not refit archived MaAsLin2/R networks."""
from pathlib import Path
import shutil,json
import numpy as np,pandas as pd
from scipy.spatial.distance import pdist,squareform
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageOps,ImageDraw
B=Path(__file__).resolve().parents[2];O=Path(__file__).resolve().parents[1];P=B/'inputs';T=O/'tables';F=O/'figures'
D=['16S','ITS','AMF'];labels={'16S':'16S','ITS':'ITS','AMF':'18S (AMF)'};colors=['#168a76','#d76b13']
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(fig,n):
 for ext in ['png','pdf']:fig.savefig(F/f'{n}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
meta=pd.read_csv(B/'data/manuscript_metadata.csv');public=meta.copy();public.to_csv(T/'sample_metadata_long.csv',index=False)
u=public.drop_duplicates('sample_id').copy().sort_values('soil_uid',key=lambda x:x.str.extract(r'(\d+)',expand=False).astype(int))
for d in D:u[d]=u.sample_id.isin(public.loc[public.marker.eq(d),'sample_id']).map({True:'Yes',False:'No'})
u=u[['sample_id','profile_id','management_group','site','Depth_cm','sample_date']+D];u.to_csv(T/'S1_sample_membership.csv',index=False)
assert len(u)==45 and u[D].eq('Yes').all(axis=1).sum()==42
summ=[]
for d in D:
 x=pd.read_csv(B/'data'/(d+'_counts.csv.gz'),index_col=0);s=public[public.marker.eq(d)];assert set(x.index)==set(s.sample_id)
 summ.append(dict(Marker=labels[d],NT=int(s.management_group.eq('NT').sum()),Till=int(s.management_group.eq('till').sum()),ASVs=int((x.sum()>0).sum()),Median_reads=float(x.sum(axis=1).median())))
pd.DataFrame(summ).to_csv(T/'Table1_coverage.csv',index=False)
a=pd.read_csv(P/'final_reference/rarefaction/alpha_tests_rarefied.csv');a.to_csv(T/'alpha_all_scenarios.csv',index=False)
a[a.scenario.eq('primary_all_samples')].to_csv(T/'Table2_alpha.csv',index=False)
a[a.scenario.eq('sensitivity_95pct')].to_csv(T/'S2_alpha_sensitivity.csv',index=False)
for f in ['depth_selection.csv','unique_ASV_retention_summary.csv','sample_alpha_rarefied.csv','per_draw_test_summary.csv']:shutil.copy2(P/'final_reference/rarefaction'/f,T/f)
beta=pd.read_csv(P/'final_reference/restricted_permanova/model_statistics.csv');beta=beta[(beta.model=='management')&(beta.method=='PERMANOVA_uncorrected_Bray')&(beta.restriction=='site_and_depth_coverage')].copy()
assert np.allclose(beta.q_BH,multipletests(beta.p,method='fdr_bh')[1]);beta.to_csv(T/'Table3_beta_restricted.csv',index=False)
orig=[]
for d in D:
 z=pd.read_csv(P/f'results/{d}/{d}_nmds_permanova.csv').query("term=='Model'").iloc[0];disp=pd.read_csv(P/'results/imbalance_sensitivity/ordination_observed_permanova_permdisp.csv').set_index('domain').loc[d]
 orig.append(dict(marker=d,R2=z.R2,pseudo_F=z.F,p=z['Pr(>F)'],dispersion_F=disp.permdisp_f,dispersion_p=disp.permdisp_p))
pd.DataFrame(orig).to_csv(T/'S3_unrestricted_beta_dispersion.csv',index=False)
allres=[];discover=[]
for d in D:
 z=pd.read_csv(P/f'results/{d}/{d}_maaslin2_nt_vs_till_results.csv');z.insert(0,'marker',d);allres.append(z)
 discover.append(dict(marker=d,tested=len(z),q005=int((z.qval<=.05).sum()),q025=int((z.qval<=.25).sum()),higher_NT_q005=int(((z.qval<=.05)&(z.coef<0)).sum()),higher_till_q005=int(((z.qval<=.05)&(z.coef>0)).sum())))
res=pd.concat(allres,ignore_index=True);res.to_csv(T/'MaAsLin2_all_results.csv',index=False);res[res.qval<=.05].to_csv(T/'S4_MaAsLin2_q005.csv',index=False);pd.DataFrame(discover).to_csv(T/'MaAsLin2_discovery_counts.csv',index=False)
assert [x['q005'] for x in discover]==[5,14,1]
net=pd.read_csv(P/'results/network_stats_summary.csv');net=net[['label','nodes','edges','density','mean_degree']];net.to_csv(T/'S5_network_metrics.csv',index=False)
nb=pd.read_csv(P/'results/imbalance_sensitivity/network_balanced_summary.csv');nb=nb[['domain','label','median_nodes','median_edges','median_mean_degree','median_density']];nb.to_csv(T/'S6_network_balanced.csv',index=False)
bs=pd.read_csv(P/'results/imbalance_sensitivity/ordination_balanced_summary.csv');bs.to_csv(T/'S7_beta_balanced.csv',index=False)
for name in ['ordination_balanced_subsamples_permanova_permdisp.csv','network_balanced_subsamples_stats.csv']:
 z=pd.read_csv(P/'results/imbalance_sensitivity'/name)
 if name.startswith('network'):z=z.drop(columns=[c for c in z if 'wilcox' in c or 'ks_' in c or c in ['transitivity','modularity']],errors='ignore')
 z.to_csv(T/name,index=False)
# Copy high-resolution primary alpha and rarefaction figures, which already match the verified tables.
for old,new in [('figure1_alpha','Figure1_alpha'),('figure9_rarefaction','FigureS1_rarefaction')]:
 for ext in ['png','pdf']:shutil.copy2(P/'final_reference/manuscript_figures'/f'{old}.{ext}',F/f'{new}.{ext}')
# Replot only the retained three marker PCoAs, with restricted-management tests.
co=pd.read_csv(P/'results/ordination/ordination_coordinates.csv');ss=pd.read_csv(P/'results/ordination/ordination_summary.csv').set_index('marker');beta=beta.set_index('marker')
fig,axs=plt.subplots(1,3,figsize=(10,3.5))
for ax,d,l in zip(axs,D,'ABC'):
 z=co[co.marker.eq(d)];b=beta.loc[d];s=ss.loc[d]
 for mg,c,m in zip(['NT','till'],colors,['o','^']):
  w=z[z.management_group.eq(mg)];ax.scatter(w.PCoA1,w.PCoA2,c=c,marker=m,s=32,alpha=.85,label='NT' if mg=='NT' else 'Till')
 ax.set(title=f'{l}  {labels[d]}\nR² = {b.R2:.3f}; q = {b.q_BH:.4f}',xlabel=f'PCoA 1 ({s.axis1_percent:.1f}%)',ylabel=f'PCoA 2 ({s.axis2_percent:.1f}%)');ax.axhline(0,color='.85',lw=.5);ax.axvline(0,color='.85',lw=.5)
axs[0].legend(frameon=False,fontsize=9);fig.tight_layout();save(fig,'Figure2_composition')
co[['marker','sample_id','management_group','PCoA1','PCoA2']].to_csv(T/'PCoA_coordinates.csv',index=False)
# Clear summary of DA discovery counts; complete 20-feature coefficient table supplied separately.
fig,ax=plt.subplots(figsize=(6.6,3.5));x=np.arange(3)
for off,cut,c in [(-.17,.05,'#334e68'),(.17,.25,'#a6bbc8')]:
 v=[int(((res.marker==d)&(res.qval<=cut)).sum()) for d in D];bars=ax.bar(x+off,v,.29,color=c,label=f'q ≤ {cut:.2f}');ax.bar_label(bars,padding=3)
ax.set(xticks=x,xticklabels=[labels[d] for d in D],ylabel='Management-associated ASVs',ylim=(0,52));ax.legend(frameon=False);fig.tight_layout();save(fig,'Figure3_differential_abundance')
fig,axs=plt.subplots(1,3,figsize=(10,3.5))
for ax,col,label in zip(axs,['nodes','mean_degree','density'],['Connected ASVs','Mean degree','Density']):
 for i,(mg,c) in enumerate(zip(['NT','till'],colors)):
  vals=[net.loc[net.label.eq(f'{d} {mg}'),col].iloc[0] for d in D];ax.bar(np.arange(3)+(i-.5)*.34,vals,.3,color=c,label='NT' if mg=='NT' else 'Till')
 ax.set(xticks=range(3),xticklabels=[labels[d] for d in D],ylabel=label);ax.tick_params(axis='x',labelrotation=15)
axs[0].legend(frameon=False);fig.tight_layout();save(fig,'Figure4_network_metrics')
# Balance sensitivity retained for interpretation, without cross-marker or date-model panels.
z=pd.read_csv(T/'ordination_balanced_subsamples_permanova_permdisp.csv');fig,axs=plt.subplots(1,2,figsize=(8,3.8));rng=np.random.default_rng(7)
for ax,col,label in zip(axs,['permanova_r2','permanova_p'],['PERMANOVA R²','Unrestricted permutation p']):
 vals=[z.loc[z.domain.eq(d),col].values for d in D];ax.boxplot(vals,showfliers=False,medianprops={'color':'black'})
 for i,v in enumerate(vals,1):ax.scatter(i+rng.uniform(-.09,.09,len(v)),v,s=17,color='#334e68',alpha=.7)
 ax.set(xticks=range(1,4),xticklabels=[labels[d] for d in D],ylabel=label)
axs[1].axhline(.05,color='#ad2939',ls='--');fig.tight_layout();save(fig,'FigureS2_balance')
# Supporting spreadsheet keeps all 1,129 tested associations, including nonsignificant results.
with pd.ExcelWriter(O/'Supplementary_Data.xlsx',engine='openpyxl') as w:
 for f in T.glob('*.csv'):
  z=pd.read_csv(f);name=f.stem[:31];z.to_excel(w,sheet_name=name,index=False)
 for ws in w.book.worksheets:
  ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
  for col in ws.columns:
   width=min(55,max(12,max(len(str(c.value or '')) for c in list(col)[:60])+2));ws.column_dimensions[col[0].column_letter].width=width
# Contact sheet for visual inspection.
ims=[]
for p in sorted(F.glob('*.png')):
 im=Image.open(p).convert('RGB');im.thumbnail((700,430));tile=Image.new('RGB',(720,470),'white');tile.paste(im,((720-im.width)//2,30));ImageDraw.Draw(tile).text((10,5),p.stem,fill='black');ims.append(tile)
sheet=Image.new('RGB',(1440,470*((len(ims)+1)//2)),'#dddddd')
for i,im in enumerate(ims):sheet.paste(im,((i%2)*720,(i//2)*470))
sheet.save(O/'audit/figure_contact_sheet.jpg')
print('Prepared 4 main figures, 2 supplementary figures, curated statistics and complete DA workbook.')
