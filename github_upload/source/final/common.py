from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
def load_counts(marker):
    a=pd.read_csv(ROOT/'data'/f'{marker}_counts.csv.gz',index_col=0)
    return {'dimnames':[list(a.index),list(a.columns)],'taxa_are_rows':[False],'dim':list(a.shape)},a.to_numpy()
def load_metadata():
    return pd.read_csv(ROOT/'data/sample_metadata.csv').rename(columns={'profile_id':'Site_Name'})
def load_alpha(marker):
    m=pd.read_csv(ROOT/'data'/f'{marker}_samples.csv').set_index('sample_id')
    m['reads']=m.retained_count_total
    m['management']=m.management_group.replace({'till':'vtill'})
    return m
