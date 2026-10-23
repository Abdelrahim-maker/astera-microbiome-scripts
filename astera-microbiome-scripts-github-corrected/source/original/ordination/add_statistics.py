"""Run from workspace root. Regenerates base ordination, then adds date-group statistics."""
import make_ordination as base
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.stats import f_oneway
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests
P=base.OUT; B=9999; SEED=20260922
rows=[]
def permanova(d2,g):
 n=len(g);levels=np.unique(g);total=d2.sum()/(2*n)
 within=sum(d2[np.ix_(g==v,g==v)].sum()/(2*np.sum(g==v)) for v in levels)
 between=total-within
 return (between/(len(levels)-1))/(within/(n-len(levels))),between/total

def dispersion(X,g):
 groups=[X[g==v] for v in np.unique(g)]
 z=[np.linalg.norm(x-x.mean(0),axis=1) for x in groups]
 means=np.array([a.mean() for a in z]);ns=np.array([len(a) for a in z]);grand=np.sum(means*ns)/ns.sum()
 between=np.sum(ns*(means-grand)**2);within=sum(np.sum((a-a.mean())**2) for a in z)
 return (between/(len(z)-1))/(within/(ns.sum()-len(z)))

for k,domain in enumerate(base.allmeta):
 m=base.allmeta[domain].copy();v=base.counts[domain];d=squareform(pdist(v,metric='braycurtis'));d2=d*d;g=pd.Categorical(m.sample_date,categories=base.DATES).codes
 A=-.5*d2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();ev=np.linalg.eigvalsh(G);c=max(0,-ev[0]);dc2=d2+2*c;np.fill_diagonal(dc2,0);A=-.5*dc2;G=A-A.mean(0)-A.mean(1)[:,None]+A.mean();ev,U=np.linalg.eigh(G);X=U[:,ev>1e-10]*np.sqrt(ev[ev>1e-10])
 assert np.allclose(squareform(pdist(X))**2,dc2,atol=1e-9)
 F,R2=permanova(d2,g);DF=dispersion(X,g)
 # Independent verification of pseudo-F by projection matrices.
 Z=np.eye(len(np.unique(g)))[g];H=Z@np.linalg.inv(Z.T@Z)@Z.T;J=np.ones_like(H)/len(g)
 A=-.5*d2;rawG=A-A.mean(0)-A.mean(1)[:,None]+A.mean()
 alt=np.trace((H-J)@rawG)/(Z.shape[1]-1)/(np.trace((np.eye(len(g))-H)@rawG)/(len(g)-Z.shape[1]))
 assert np.isclose(F,alt)
 assert np.isclose(DF,f_oneway(*[np.linalg.norm(X[g==v]-X[g==v].mean(0),axis=1) for v in np.unique(g)]).statistic)
 profiles=sorted(m.Site_Name.unique());pindex={p:i for i,p in enumerate(profiles)};ri=np.array([pindex[p] for p in m.Site_Name]);pg=[];sig=[]
 for p in profiles:
  z=m[m.Site_Name==p];assert z.sample_date.nunique()==1
  assert not z.Depth_cm.duplicated().any()
  pg.append(g[z.index[0]]);sig.append('|'.join(sorted(z.Depth_cm)))
 pg=np.array(pg);sig=np.array(sig);strata=[np.where(sig==s)[0] for s in np.unique(sig)]
 rng=np.random.default_rng(SEED+k);exceed=np.zeros(2,int)
 for _ in range(B):
  shuffled=pg.copy()
  for ids in strata:shuffled[ids]=rng.permutation(pg[ids])
  perm=shuffled[ri]
  assert np.array_equal(np.bincount(perm),np.bincount(g))
  exceed += [permanova(d2,perm)[0]>=F-1e-12,dispersion(X,perm)>=DF-1e-12]
 pv=(exceed+1)/(B+1)
 rows.append(dict(marker=domain,n=len(m),profiles=len(profiles),date_groups=len(np.unique(g)),df_date=len(np.unique(g))-1,df_residual=len(m)-len(np.unique(g)),pseudo_F=F,R2=R2,permanova_p=pv[0],dispersion_F=DF,dispersion_p=pv[1],permutations=B,seed=SEED+k))
 print(domain,rows[-1],flush=True)
stats=pd.DataFrame(rows);stats['permanova_q_BH']=multipletests(stats.permanova_p,method='fdr_bh')[1];stats['dispersion_q_BH']=multipletests(stats.dispersion_p,method='fdr_bh')[1];stats.to_csv(P/'date_group_statistics.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(16,7.6));fig.subplots_adjust(left=.065,right=.98,bottom=.40,top=.82,wspace=.30)
for ax,domain in zip(axs,base.allmeta):
 base.panel(ax,domain);r=stats.set_index('marker').loc[domain]
 verdict='Date groups differ (exploratory)' if r.permanova_q_BH<.05 else 'No significant date-group difference'
 txt=f'PERMANOVA: F = {r.pseudo_F:.2f}, R² = {r.R2:.3f}\np = {r.permanova_p:.4f}; BH q = {r.permanova_q_BH:.4f}\n{verdict}\nDispersion: p = {r.dispersion_p:.4f}; BH q = {r.dispersion_q_BH:.4f}'
 ax.text(.5,-.24,txt,transform=ax.transAxes,ha='center',va='top',fontsize=10.5,linespacing=1.6,bbox=dict(boxstyle='round,pad=.6',facecolor='#f2f5f7',edgecolor='#d5dce0'))
fig.suptitle('Manuscript samples: beta diversity by sample date',fontsize=19,y=.975)
fig.text(.5,.92,'Bray–Curtis PCoA · colors: sampling date · shapes: management',ha='center',fontsize=12)
handles=base.handles+[Line2D([],[],marker=s,ls='',color='#555555',label=f,markersize=7) for f,s in base.shapes.items()]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.09),ncol=6,frameon=False)
fig.text(.5,.055,'9,999 profile-level permutations, restricted by depth coverage; BH correction across the three markers.',ha='center',fontsize=10)
fig.text(.5,.025,'Date overlaps with site and management: these tests do not isolate a temporal effect.',ha='center',fontsize=10,color='#7b3d13')
for name in ['beta_diversity_by_sample_date','beta_diversity_by_sample_date_with_statistics']:
 for ext in ['png','pdf']:fig.savefig(P/f'{name}.{ext}',dpi=240)
plt.close(fig)
(P/'STATISTICS.md').write_text('''# Date-group tests for the manuscript cohort

The figure and date_group_statistics.csv report omnibus tests across four collection dates. PERMANOVA uses the full uncorrected Bray-Curtis distance matrix on relative abundances, not just the two plotted axes. R² is date-group between-sum-of-squares divided by total sum-of-squares. Dispersion uses ANOVA of distances to date-group centroids in all dimensions of Lingoes-corrected PCoA, with centroids recalculated under each label permutation (centroid PERMDISP). This differs from vegan's residual-permutation default.

Both tests use 9,999 random label permutations, with (exceedances + 1)/(9999 + 1). Soil profiles (workbook Site_Name) move together. Profile date labels are shuffled only among profiles with the same observed depth-interval set; this preserves date-group sample sizes and depth coverage. The 24 profiles per marker are treated as exchangeable within these depth strata. Broader field clustering is not controlled, so results remain exploratory. No date effect independent of site/management is claimed. All till samples belong to November 28; other dates are NT. Shuffling dates within profiles would be invalid for this question because each profile has only one collection date.

Benjamini-Hochberg correction is applied across the three markers, separately for PERMANOVA and dispersion. A significant omnibus result does not identify which pairs of dates differ. Lack of significance does not establish equivalence. A significant dispersion test complicates attribution of PERMANOVA significance to centroid differences alone.

Validation: sample IDs match manuscript cohorts; profile dates are unique; permutation group sizes are invariant; PERMANOVA pseudo-F matches an independent projection-matrix calculation; dispersion F matches scipy ANOVA; full corrected coordinates reconstruct corrected distances. Seeds and degrees of freedom are in the results CSV.

Method references: https://scikit.bio/docs/latest/generated/skbio.stats.distance.permanova.html ; https://scikit.bio/docs/dev/generated/skbio.stats.distance.permdisp.html ; https://search.r-project.org/CRAN/refmans/vegan/html/permutations.html

Reproduce: python analysis_outputs/manuscript_sample_date_ordination/add_statistics.py
''')
with (P/'README.md').open('a') as f:f.write('\n\nStatistics are now included in the main figure. See STATISTICS.md and date_group_statistics.csv; rerun add_statistics.py to regenerate annotated figures. This supersedes the earlier statement that no tests were reported.\n')
print(stats.to_string(index=False))
