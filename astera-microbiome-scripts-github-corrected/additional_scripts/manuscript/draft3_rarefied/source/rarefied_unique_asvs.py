"""Count each ASV once across samples within a rarefaction draw; reproduce original draws."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'manuscript/draft3_rarefied/rarefaction'
settings=pd.read_csv(OUT/'depth_selection.csv')
saved=pd.read_csv(OUT/'rarefaction_replicates.csv')
summary=[];draws=[]
for marker in ['16S','ITS','AMF']:
 info=json.loads((ROOT/f'manuscript/audit/{marker}_matrix_metadata.json').read_text())
 raw=np.load(ROOT/f'manuscript/audit/{marker}_counts.npz')['values']
 m=pd.read_csv(ROOT/f'manuscript/audit/{marker}_sample_alpha.csv',index_col=0).sort_index()
 lookup={s:i for i,s in enumerate(info['dimnames'][0])}
 v=raw[[lookup[s] for s in m.index]].astype(np.int64)
 original=int((v.sum(0)>0).sum())
 for r in settings[settings.marker.eq(marker)].itertuples():
  rng=np.random.default_rng(r.seed);present=np.zeros((r.repeats,v.shape[1]),dtype=bool)
  keep=v.sum(1)>=r.depth;retained_original=int((v[keep].sum(0)>0).sum())
  for i,sid in enumerate(m.index):
   if not keep[i]:continue
   ix=np.flatnonzero(v[i]);z=rng.multivariate_hypergeometric(v[i,ix],r.depth,size=r.repeats)
   expected=saved[saved.marker.eq(marker)&saved.scenario.eq(r.scenario)&saved.sample_id.eq(sid)].sort_values('replicate')
   assert np.array_equal((z>0).sum(1),expected.Observed.to_numpy())
   present[:,ix]|=z>0
  counts=present.sum(1)
  assert np.all(counts<=retained_original)
  for rep,count in enumerate(counts,1):draws.append(dict(marker=marker,scenario=r.scenario,replicate=rep,unique_ASVs=int(count),percent_original=100*count/original))
  summary.append(dict(marker=marker,scenario=r.scenario,depth=r.depth,samples=int(keep.sum()),original_unique_ASVs=original,retained_samples_original_unique_ASVs=retained_original,mean_unique_ASVs=counts.mean(),sd_unique_ASVs=counts.std(ddof=1),min_unique_ASVs=counts.min(),max_unique_ASVs=counts.max(),percent_retained=100*counts.mean()/original,percent_reduction=100*(1-counts.mean()/original),union_across_all_200_draws=int(present.any(0).sum())))
pd.DataFrame(summary).to_csv(OUT/'unique_ASV_retention_summary.csv',index=False)
pd.DataFrame(draws).to_csv(OUT/'unique_ASV_retention_per_draw.csv',index=False)
print(pd.DataFrame(summary).to_string(index=False))
