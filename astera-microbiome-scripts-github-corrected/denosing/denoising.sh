#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=PAS1212
#SBATCH --job-name=denoising
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=abdelrahim.11@osu.edu
#SBATCH --partition=serial


##denoising 16s rrna sequences
#qiime dada2 denoise-paired --i-demultiplexed-seqs demuxed-sequences1.qza --p-trunc-len-f 151 --p-trunc-len-r 151 --o-table table_1.qza --o-representative-sequences rep_seqs_1.qza --o-denoising-stats denoising_demux_stats_1.qza

qiime dada2 denoise-paired --i-demultiplexed-seqs 16S_2_cs_pairedend_sequences.qza --p-trunc-len-f 142 --p-trunc-len-r 151 --o-table table_2.qza --o-representative-sequences rep_seqs_2.qza --o-denoising-stats denoising_demux_stats_2.qza
