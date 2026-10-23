# Reproduction and provenance
Every result/figure is listed in WORKFLOW_INDEX.csv. Workflow assignments identify the analysis family; postprocessing steps whose original script was not recovered are explicitly marked. This is not a claim of an unbroken provenance chain for every historical result.

From the extracted ZIP root:
1. Install the Python dependencies in environment/requirements-python-validated.txt.
2. Run `python source/check_processed_data.py` to verify sample counts and retained totals.
3. Run `python source/reconstruct_network_edges.py` to regenerate all eight edge tables and validate degrees/hubs.

Original R scripts and ordination scripts are retained under source/original. They still require original phyloseq inputs or workbook/exports beyond this deposit; see README.md. Original NMDS plots and saved manuscript PCoA coordinates are included; the PCoA coordinates are not the coordinates of the NMDS figures. Historical NMDS coordinates were not exported. The original R package-version lock and some network annotation postprocessing scripts were not found; see environment/README.md and the index status column.
