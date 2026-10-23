from pathlib import Path
from PIL import Image
import re,json
O=Path('manuscript')
for f in ['manuscript.md','supplementary_material.md']:
 p=O/f;s=p.read_text()
 def dim(m):
  im=Image.open(O/m[1]);w=min(float(m[2]),float(m[3])*im.width/im.height)
  return f'![]({m[1]}){{width={w:.3f}in}}'
 s=re.sub(r'!\[\]\(([^)]+)\)\{width=([\d.]+)in height=([\d.]+)in\}',dim,s)
 if f=='manuscript.md':
  replacements={
  'Individual source plots and rarefaction curves are reproduced in the supplementary figure atlas.':'Individual source plots and rarefaction curves are reproduced in Figures S10–S12, S22–S24, S34–S36, S51, and S53.',
  'Full source ordinations are retained in the supplementary atlas.':'Full source ordinations are retained in Figures S9, S21, S33, and S52; the integrated image is also archived as Figure S50.',
  'A permissive discovery threshold and sensitivity to subsampling':'The original coefficient plots are Figures S1, S13, and S25. A permissive discovery threshold and sensitivity to subsampling',
  'These are contrasts between thresholded statistical association graphs':'Original graphs and treatment overlays are Figures S2, S4–S6, S8, S14, S16–S18, S20, S26, S28–S30, S32, and S46–S49. These are contrasts between thresholded statistical association graphs',
  'These are overlapping sample sets, not three independent confirmations':'The original carbon boxplots are Figures S37–S39. These are overlapping sample sets, not three independent confirmations',
  'These counts represent model terms, so':'Carbon term counts and community summaries are Figures S40–S41. These counts represent model terms, so',
  'These are supplied taxonomic assignments':'The original overlays and overlap summaries are Figures S3, S7, S15, S19, S27, S31, S42, and S43. These are supplied taxonomic assignments',
  'All post-selection feature–carbon correlations, bootstrap metric records, and original visualizations remain available in the supplement.':'All post-selection feature–carbon correlations and bootstrap metric records remain available in the supplement, with original visualizations in Figures S44–S45.'}
  for a,b in replacements.items():s=s.replace(a,b)
 p.write_text(s)
p=O/'source/references.json';rs=json.loads(p.read_text())
for r in rs:
 if r['id'] in {'Anderson2006','Bowles2017','Warton2012'}:r['issued']['date-parts'][0][0]={'Anderson2006':2006,'Bowles2017':2017,'Warton2012':2012}[r['id']]
p.write_text(json.dumps(rs,indent=2,ensure_ascii=False))
# Preserve source provenance in the CSV index without overlong unbreakable caption paths.
p=O/'supplementary_material.md';s=p.read_text();s=re.sub(r'Original source: `[^`]+`\. ', 'Original file identified by this figure number in figure_index.csv. ',s);s=s.replace('Fills/shapes encode treatment-only enrichment.','Fills encode treatment-only coefficient direction regardless of significance; squares indicate q ≤ 0.25.');p.write_text(s)
p=O/'manuscript.md';s=p.read_text().replace('Fill colors and square shapes indicate treatment-only associations as defined in the source legend;','Fill colors indicate treatment-only coefficient direction regardless of significance; square shapes indicate q ≤ 0.25;');p.write_text(s)
