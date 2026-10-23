# Environment provenance
Python pins record the installed environment used successfully for the newly reconstructed networks and deposit checks. Install with `python -m pip install -r environment/requirements-python-validated.txt`. Python version is in python_environment.json. A fresh installation was not tested.

The original R analysis sessionInfo/renv lock was not found. install_R_dependencies.R lists dependencies extracted from the original scripts; it is unpinned and has not been executed here. It does not claim to reproduce historical package versions. Original MaAsLin2 parameter logs are under audit/model_logs. The raw-read processing environment is not reconstructed by this deposit.
