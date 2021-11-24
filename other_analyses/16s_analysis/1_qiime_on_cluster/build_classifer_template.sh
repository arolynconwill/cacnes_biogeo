cd /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/
source activate /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/qiime2-2020.1

echo "BEGIN CLASSIFICATION FILE"

#SILVA_otus.qza
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/classifer_refseqs/classifer_input_seqs.fna \
  --output-path ref-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/classifer_refseqs/classifer_input_seqs.tax \
  --output-path ref-taxonomy.qza
 
 echo "IMPORTED DATA"
 
# If your 16S database is untrimmed, feature-classifier can be used to trim your reads at a chosen primer to an appropiate length

#qiime feature-classifier extract-reads \
#  --i-sequences SILVA_otus.qza \
#  --p-f-primer AGAGTTTGATCMTGGCTCAG \
#  --p-r-primer ATTACCGCGGCTGCTGG \
#  --p-trunc-len 150 \
#  --p-min-length 100 \
#  --p-max-length 400 \
#  --o-reads ref-seqs.qza

# echo "EXTRACTED READS"
	
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

echo "CREATED CLASSIFER - BEGINNING CLASSIFICATION"

#qiime feature-classifier classify-sklearn \
#  --i-classifier classifier.qza \
#  --i-reads /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2/data.demux_se.trim.sq4.deblur.qza \
#  --o-classification taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza

#qiime metadata tabulate \
#  --m-input-file taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza \
#  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qzv

#echo "END CLASSIFCATION"
#echo "BEGIN EXPORT"

# QIIME2 can automatically condense all similarly called reads into a tsv file (ex. condensing all C. acnes classifcations into one grouping)
# If you wish for this uncommend the barplot lines below

# # barplots
#qiime taxa barplot \
#  --i-table /scratch/path/QIIME_analysis_dir/deblur-table.qza \
#  --i-taxonomy taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qza \
#  --m-metadata-file /scratch/path/QIIME_analysis_dir/metadata_sample_desc.tsv \
#  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots.qzv

# get data
#qiime tools export \
#  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots.qzv \
#  --output-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.bar-plots_export

# In order to see all your raw OTUs, their classification, classification-confidence-level, and sample counts, the below code can be used. 
# feature-table.tsv will hold all your sample counts with QIIME2 generated OTU-labels (this will look like a long string of random numbers and letters)
# The directory taxonomic_OTU_classification will hold metadata.tsv, which will have QIIME2's OTU label, classification, and classification confidence value
# The directory SILVA_OTU_sequences will contain sequences.fasta that holds your OTU fasta sequences and their QIIME2 label
# All these files can be combined manually to make one descriptive CSV. I've found that QIIME2 has outputted all OTUs in the same order every time, which makes
# it easy to copy and paste everything into one big spreadsheet. I would check this, however, by copy and pasting the OTUs from different files
# into a string array in MATLAB and comparing them. 

# Note: QIIME2 includes a whole range of packages for diversity analysis, phylogenetics, etc. There's a lot more you can do with QIIME2 beyond this script! 

#qiime tools export \
#  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2/data.demux_se.trim.sq4.deblur.qzv \
#  --output-path SILVA_OTU_sequences
	
#qiime tools export \
#  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.deblur.qzv \
#  --output-path taxonomic_OTU_classification
	
#qiime tools export \
#  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2/deblur-table.qza \
#  --output-path deblur_OTU_export
  
#biom convert -i deblur_OTU_export/feature-table.biom -o deblur_OTU_export/feature-table.tsv --to-tsv

#echo "TAXA ASSIGNMENT DATA EXPORTED!!!"

