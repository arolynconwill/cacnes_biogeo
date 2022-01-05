====================================
== 16S processing pipeline steps: ==
====================================

1_DADA2_process.slurm:

Here, a barcode metadata sheet with the format name-tab-barcode is uploaded as metadata_qiime_v1_Alex.tsv for adaptor trimming purposes. DADA2_process.slurm then imports a multiplexed data file, demultiplexes it based on metadata_qiime_v1_Alex.tsv, trims the AGAGTTTGATCMTGGCTCAG primer sequence from the first paired end read using cutadapt, and then uses DADA2 to resolve single nucleotide differences between ASVs. This results in a variety of qza and qzv files, returning all filtered and processed sequence data as  data.demux_se.trim.sq4.dada2.qza. 


2_create_classifer.slurm: 

In order to built the 16S classifier used here, classifer_input_seqs.fna and classifer_input_seqs.tax, files containing SATIVA-cleaned sequences and taxonomy classifications from the SILVA 16S rRNA database are processed via QIIME2's fit-classifier-naive-bayes command to produce classifier.qza. The code to produce classifer_input_seqs.fna and classifer_input_seqs.tax is explained below.


3_classify_seqs.slurm:

In classify_seqs.slurm, classifier.qza and data.demux_se.trim.sq4.dada2.qza are processing using QIIME2's classify-sklearn command, resulting in the output taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza file. This is transformed into a QIIME2 visualizable table taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv. Additionally, raw ASV sequences are outputted as SILVA_OTU_sequences/sequences.fasta. Raw ASV counts are additionally outputted to dada2_OTU_export/feature-table.tsv and classifications to taxonomic_OTU_classification/metadata.tsv. 



=============================================
== Building a curated SILVA rRNA database: ==
=============================================

In order to clean dna-sequences.fasta, which holds the entirely of SILVA's v132 rRNA database, taxa to be classified down to the species level are first split into files separated by genus or family by pull_out_species_for_SATIVA.m. These included the genuses Cutibacterium, Acidipropionibacterium, Pseudopropionibacterium, Propionibacterium (which are grouped together into one file), as well as Corynebacteriaceae and the chronically-mislabelled Neisseriaceae families (Neisseriaceae families are again grouped together). Only genuses or families to be cleaned by SATIVA are split into individual files named SILVA_v132_{taxon here}_unmodified.fasta; all other sequences to remain in the classifier remain in SILVA_v132_not_specified_genus_unmodified. Sequences with taxonomies that contain multiple terms that would separate them under multiple taxa files (ex. "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Corynebacterium diphtheriae") are held in SILVA_v132_confusing_genus_unmodified.fasta and discarded from analysis. 

From here, split_taxa_for_SATIVA.m removes sequences missing a species classification or assigned to non-species taxa (ex. Corynebacterium sp.) as well as any sequence with a "Eukaryota" kingdom. This file also manually cleans some taxonomic designations, for example, replacing the "Corynebacterium 1" genus with "Corynebacterium", "Propionibacterium" with "Cutibacterium", and "Neis8seriaceae" with "Neisseriaceae". This file outputs "SILVA_v132_{cleaned taxon here}_CLEAN" .fasta and .tax files. 

Here, SATIVA was run locally using the command: 
../sativa.py -s SILVA_v132_{cleaned taxon here}_CLEAN.fasta -t SILVA_v132_{cleaned taxon here}_CLEAN.tax -x BAC -p 12345

The SATIVA output files SILVA_v132_{cleaned taxon here}_CLEAN.mis were manually renamed SILVA_v132_{cleaned taxon here}_CLEAN.mis.tsv. Using cleaning_post_sativa.m, these tsv files were then parsed such that taxonomically mislabeled sequences with greater than 90% confidence
could be relabeled and sequences with below 90% confidence could be removed. This outputs the files final_SATIVA_cleaned_{cleaned taxon here}.fasta and final_SATIVA_cleaned_{cleaned taxon here}.tax, as well as final_SATIVA_other_species.fasta and final_SATIVA_other_species.tax. Additionally, cleaning_post_sativa.m also intakes the file dna-sequences_noMislabel_StaphSpecLvl.taxonomy, which contains SATIVA-cleaned Staphylococcus species from Khadka et al., 2021. and splits it into .fasta and .tax files. 

From here, the commands in linux_commands_for_comb_files.txt are used to combine all final_SATIVA_cleaned_{taxon} files together into classifer_input_seqs.tax and classifer_input_seqs.fna, which are directly inputted into QIIME2's Naive Bayes classifier function. 

All the necessary files to build the database can be found in the subdirectory entitled files_for_building_classifier_database.

