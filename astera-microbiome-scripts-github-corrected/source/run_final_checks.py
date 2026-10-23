"""Run portable manuscript workflows; statistical choices remain as recorded."""
from pathlib import Path
import subprocess,sys,os,time,json
root=Path(__file__).resolve().parents[1];logs=root/'validation';logs.mkdir(exist_ok=True)
steps=['check_processed_data.py','final/rarefied_alpha.py','final/rarefied_unique_asvs.py','final/full_dbrda.py','final/variance_without_date.py','final/restricted_permanova.py','final/adjusted_carbon.py','final/rebuild_figures.py','final/compare_outputs.py']
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR=str(logs/'matplotlib'))
records=[]
for s in steps:
 print('RUN',s,flush=True);start=time.time()
 with (logs/(Path(s).stem+'.log')).open('w') as f:r=subprocess.run([sys.executable,str(root/'source'/s)],cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
 records.append(dict(script=s,returncode=r.returncode,seconds=round(time.time()-start,2)))
 (logs/'execution.json').write_text(json.dumps(records,indent=2)+'\n')
 print('END',s,'exit',r.returncode,'seconds',round(time.time()-start,1),flush=True)
 if r.returncode:print((logs/(Path(s).stem+'.log')).read_text()[-3000:],flush=True);sys.exit(r.returncode)
