from common import ROOT
import pandas as pd
p=ROOT/'reproduced/full_dbrda/subset_models.csv'
d=pd.read_csv(p);rows=[]
# Original masks: management=1, date=2, site=4, depth=8.
for marker,g in d.groupby('marker',sort=False):
    a=g.set_index('mask').adjusted_R2;full=a[13]
    u={name:full-a[13^bit] for name,bit in [('Management',1),('Site',4),('Depth',8)]}
    rows.append(dict(marker=marker,management_unique=100*u['Management'],site_unique=100*u['Site'],depth_unique=100*u['Depth'],shared=100*(full-sum(u.values())),unexplained=100*(1-full)))
o=ROOT/'reproduced/variance_without_date';o.mkdir(parents=True,exist_ok=True)
pd.DataFrame(rows).to_csv(o/'variance_partition_summary.csv',index=False)
