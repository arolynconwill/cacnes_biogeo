echo "BEGIN CLASSIFICATION FILE"

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path input_files/classifer_input_seqs.fna \
  --output-path ref-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path input_files/classifer_input_seqs.tax \
  --output-path ref-taxonomy.qza
 
echo "IMPORTED DATA"

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza