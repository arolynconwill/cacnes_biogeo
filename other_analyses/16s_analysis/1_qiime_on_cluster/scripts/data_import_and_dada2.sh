# IF DEMULTIPLEXED
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path input_files/path_to_16S_Alex_manifest.tsv \
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
  --m-barcodes-file input_files/metadata_qiime_v1_Alex.tsv \
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

# export data files for SRA
qiime tools export \
  --input-path data.demux_se.trim.qza \
  --output-path exported-data-demux_se-trim


##################
###### Quality trim and denoise using dada2
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

qiime dada2 denoise-single \
  --i-demultiplexed-seqs data.demux_se.trim.sq4.qza \
  --p-trunc-len 180 \
  --p-n-threads 8 \
  --o-representative-sequences data.demux_se.trim.sq4.dada2.qza \
  --o-table dada2-table.qza \
  --o-denoising-stats dada2-stats.qza

# visualize results 
qiime metadata tabulate \
  --m-input-file filter_stats_sq4.qza \
  --o-visualization filter_stats_sq4.qzv

qiime deblur visualize-stats \
  --i-dada2-stats dada2-stats.qza \
  --o-visualization dada2-stats.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table dada2-table.qza \
  --o-visualization dada2-table.qzv \
  --m-sample-metadata-file metadata_summ_feature_table.tsv 

qiime feature-table tabulate-seqs \
  --i-data data.demux_se.trim.sq4.dada2.qza \
  --o-visualization data.demux_se.trim.sq4.dada2.qzv
 