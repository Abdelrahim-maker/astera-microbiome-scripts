#!/bin/bash
#SBATCH --time=28:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=PAS1212
#SBATCH --job-name=import
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=abdelrahim.11@osu.edu

###this code was used to import already demultiplexed sequences (I did not demux them, they were done by Argonne) > to import it we need a manifest file (it has a certain format based on each QIIME2 version) and manifest file is very case sensitive (https://docs.qiime2.org/2022.2/tutorials/importing/)
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /fs/ess/PAS1212/Afaf_rich/16srRNA_2/manifest_16s.csv --output-path 16S_2_paired_end_cs/16S_2_cs_pairedend_sequences.qza --input-format PairedEndFastqManifestPhred33V2



###for visualization of demultiplexed data
#qiime demux summarize --i-data AMF_cs_pairedend_sequences.qza --o-visualization amf_cs_demux.qzv

##Apply an initial quality filtering process based on quality scores. This method is an implementation of the quality filtering approach described by Bokulich et al. (Bokulich et al., 2013).
#qiime quality-filter q-score --i-demux  AMF_cs_pairedend_sequences.qza --o-filtered-sequences demux-filtered.qza --o-filter-stats demux-filter-stats.qza


#Apply the Deblur workflow using the denoise-16S action
#qiime deblur denoise-16S --i-demultiplexed-seqs demux-filtered.qza --p-trim-length 150 --p-sample-stats --p-jobs-to-start 4 --o-stats deblur-stats.qza --o-representative-sequences rep-seqs-deblur.qza --o-table table-deblur.qza (did not use this)

#qiime dada2 denoise-paired --i-demultiplexed-seqs demux-filtered.qza --p-trunc-len-f 214 --p-trunc-len-r 195 --o-table table.qza --o-representative-sequences rep-seqs-demux.qza --o-denoising-stats denoising-demux-stats.qza

#qiime dada2 denoise-paired --i-demultiplexed-seqs AMF_cs_pairedend_sequences.qza  --p-trunc-len-f 214 --p-trunc-len-r 195 --o-table table.qza --o-representative-sequences rep-seqs-demux.qza --o-denoising-stats denoising-demux-stats.qza

##useful link https://docs.qiime2.org/2022.2/tutorials/atacama-soils/

#useful link https://curr-protoc-bioinformatics.qiime2.org/

#useful link https://chmi-sops.github.io/mydoc_qiime2.html



##importing maarjam database (sequences and taxonomy file 2019) from https://maarjam.ut.ee/?action=bDownload 
#qiime tools import --type 'FeatureData[Sequence]' --input-path maarjam_database_SSU_TYPE.qiime.fasta  --output-path maarjam_sequences.qza 
#qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path maarjam_database_SSU_TYPE.qiime.txt --output-path maarjam_taxonomy.qza 

#import for 16s rrna

#qiime tools import --type 'EMPPairedEndSequences' --input-path emp-paired-end-sequences/ --output-path 16s_cs_multiplexed_1.qza (this worked)

#--input-format PairedEndFastqManifestPhred33V2

#qiime cutadapt demux-paired --i-seqs 16s_cs_multiplexed.qza --m-barcodes-file 231220_Rich_16S_MF.txt --m-barcodes-column BarcodeSequence --o-per-sample-sequences demuxed-16s-run1-seqs.qza --o-error-correction-details demux-run1-details.qza

#qiime cutadapt demux-paired --i-seqs 16s_cs_multiplexed.qza --m-forward-barcodes-file 231220_Rich_16S_MF.txt --m-forward-barcodes-column BarcodeSequence --output-dir 16S_demux_final 

#qiime cutadapt demux-single --i-seqs 16s_cs_multiplexed.qza --m-barcodes-file 231220_Rich_16S_MF.txt --m-barcodes-column BarcodeSequence --output-dir 16S_demux_single


#qiime cutadapt demux-paired --i-seqs 16s_cs_multiplexed.qza --m-forward-barcodes-file 231220_Rich_16S_MF.txt --m-forward-barcodes-column BarcodeSequence --output-dir demuxed-reads






#qiime demux emp-paired --i-seqs 16s_cs_multiplexed_1.qza  --m-barcodes-file 231220_Rich_16S_MF.txt --m-barcodes-column BarcodeSequence --o-per-sample-sequences demux_16s_1 --o-error-correction-details error-correction-details.qza


#qiime demux emp-paired --i-seqs 16s_cs_multiplexed_1.qza --m-barcodes-file 231220_Rich_16S_MF.txt  --m-barcodes-column BarcodeSequence --p-no-golay-error-correction --o-per-sample-sequences demuxed-sequences1.qza --o-error-correction-details error-correction-details.qza


#qiime dada2 denoise-paired --i-demultiplexed-seqs demuxed-sequences1.qza --p-trunc-len-f 151 --p-trunc-len-r 151 --o-table table_16s1.qza --o-representative-sequences rep-seqs-16s-demux1.qza --o-denoising-stats denoising-demux-16s-stats1.qza




##denoising 16s rrna sequences

qiime dada2 denoise-paired --i-demultiplexed-seqs demuxed-sequences1.qza --p-trunc-len-f 151 --p-trunc-len-r 151 --o-table table_1.qza --o-representative-sequences rep_seqs_1.qza --o-denoising-stats denoising_demux_stats_1.qza
