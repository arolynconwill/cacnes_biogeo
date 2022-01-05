echo "BEGIN CLASSIFICATION FILE"
echo "CREATED CLASSIFER - BEGINNING CLASSIFICATION"

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads data.demux_se.trim.sq4.dada2.qza \
  --o-classification taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza \
  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv

echo "END CLASSIFCATION"
echo "BEGIN EXPORT"

qiime tools export \
  --input-path data.demux_se.trim.sq4.dada2.qzv \
  --output-path SILVA_OTU_sequences
	
qiime tools export \
  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv \
  --output-path taxonomic_OTU_classification
	
qiime tools export \
  --input-path dada2-table.qza \
  --output-path dada2_OTU_export
  
biom convert -i dada2_OTU_export/feature-table.biom -o dada2_OTU_export/feature-table.tsv --to-tsv

echo "TAXA ASSIGNMENT DATA EXPORTED!!!"

