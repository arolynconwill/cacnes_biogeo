cd /scratch/mit_lieberman/projects/alex_aro_cacnes/QIIME_03_26_21_rerun_acnes/
source activate /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/qiime2-2020.1

echo "INITIALIZED"
# IF DEMULTIPLEXED
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /scratch/mit_lieberman/projects/alex_aro_cacnes/QIIME_03_26_21_rerun_acnes/path_to_16S_Alex_manifest.tsv \
  --output-path data.demux_se.qza \
  --input-format SingleEndFastqManifestPhred33V2

echo "IMPORTED"
# IF MULTIPLEXED
# import data
qiime tools import \
   --type EMPSingleEndSequences \
   --input-path /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/raw_data_links/160420Alm_D16-3717_1_sequence_link.fastq \
   --output-path data.qza 

# Demultiplex samples (using only single read from paired end)
qiime demux emp-single \
  --i-seqs data.qza \
  --m-barcodes-file metadata_qiime_v1_Alex.tsv \
  --m-barcodes-column BarcodeSeqS5N7 \
  --o-per-sample-sequences data.demux_se.qza

# check what we got
qiime demux summarize \
  --i-data data.demux_se.qza \
  --o-visualization data.demux_se.qzv

echo "CUTADAPT"
# remove primer (use only SE data here, thus need to remove only one)
qiime cutadapt trim-single \
  --p-cores 8 \
  --i-demultiplexed-sequences data.demux_se.qza \
  --p-front AGAGTTTGATCMTGGCTCAG \
  --p-discard-untrimmed \
  --o-trimmed-sequences data.demux_se.trim.qza

# check what we got
qiime demux summarize \
  --i-data data.demux_se.trim.qza \
  --o-visualization data.demux_se.trim.qzv


##################
###### Quality trim and denoise using deblur
##################

echo "TRIM AND DENOISE"

# remove low seq quality reads
qiime quality-filter q-score \
  --i-demux data.demux_se.trim.qza \
  --p-min-quality 4 \
  --o-filtered-sequences data.demux_se.trim.sq4.qza \
  --o-filter-stats filter_stats_sq4.qza

 # check what we got
qiime demux summarize \
  --i-data data.demux_se.trim.sq4.qza \
  --o-visualization data.demux_se.trim.sq4.qzv
 
 
 # # NOTE 8 cores
#qiime deblur denoise-16S \
#  --i-demultiplexed-seqs data.demux_se.trim.sq4.qza \
#  --p-trim-length 150 \
#  --p-jobs-to-start 8 \
#  --o-representative-sequences data.demux_se.trim.sq4.deblur.qza \
#  --o-table deblur-table.qza \
#  --p-sample-stats \
#  --o-stats deblur-stats.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs data.demux_se.trim.sq4.qza \
  --p-trunc-len 180 \
  --p-n-threads 8 \
  --o-representative-sequences data.demux_se.trim.sq4.deblur.qza \
  --o-table deblur-table.qza \
  --o-denoising-stats deblur-stats.qza

# visualize results 
qiime metadata tabulate \
  --m-input-file filter_stats_sq4.qza \
  --o-visualization filter_stats_sq4.qzv

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table deblur-table.qza \
  --o-visualization deblur-table.qzv \
  --m-sample-metadata-file metadata_summ_feature_table.tsv 

qiime feature-table tabulate-seqs \
  --i-data data.demux_se.trim.sq4.deblur.qza \
  --o-visualization data.demux_se.trim.sq4.deblur.qzv
 