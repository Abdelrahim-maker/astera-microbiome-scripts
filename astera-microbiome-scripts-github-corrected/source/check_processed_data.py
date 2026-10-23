from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
for marker,expected in [('16S',44),('ITS',43),('AMF',44)]:
    counts=pd.read_csv(root/'data'/f'{marker}_counts.csv.gz',index_col=0)
    samples=pd.read_csv(root/'data'/f'{marker}_samples.csv').set_index('sample_id')
    assert len(counts)==expected and counts.index.equals(samples.index)
    assert (counts.sum(axis=1)==samples.retained_count_total).all()
    print(marker, len(counts), 'samples;', counts.shape[1], 'ASVs; read totals verified')
