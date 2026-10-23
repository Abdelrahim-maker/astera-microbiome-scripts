"""Run analyses from source data; no archived manuscript results or figures required."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parents[1]
for script in ['source/final/rarefied_alpha.py','source/final/rarefied_unique_asvs.py','source/final/restricted_management.py','source/reconstruct_marker_networks.py']:
 subprocess.run([sys.executable,str(root/script)],cwd=root,check=True)
