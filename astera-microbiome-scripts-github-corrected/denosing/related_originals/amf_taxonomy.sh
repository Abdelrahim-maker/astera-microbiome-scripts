#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=PAS1212
#SBATCH --mem-per-cpu=2G
#SBATCH --job-name=tax_assignment


###train the maarjam classifier


#qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads unite_99_seqs.qza --i-reference-taxonomy unite_99_taxonomy.qza --o-classifier unite99_classifier.qza


##classify sequences

#qiime feature-classifier classify-sklearn --i-classifier  unite99_classifier.qza --i-reads rep-seqs-demux.qza --o-classification amf_cs_ASV_unit_taxonomy.qza


###16s taxonomy assignment using sliva 99%

qiime feature-classifier classify-sklearn --i-classifier  silva-138-99-nb-classifier.qza --i-reads merged_16S_rep_seqs.qza --o-classification 16s_merged_ASV_cs_99_taxonomy.qza


