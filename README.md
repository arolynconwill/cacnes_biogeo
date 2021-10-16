_C. acnes_ biogeography
=======================

This repository contains all code necessary to reproduce the analyses described in:

> [__Arolyn Conwill, Anne C. Kuan, Ravalika Damerla, Alexandra J. Poret, Jacob S. Baker, A. Delphine Tripp, Eric J. Alm, Tami D. Lieberman.__ Anatomy promotes neutral coexistence of strains in the human skin microbiome. _bioRxiv_ 2021.](https://www.biorxiv.org/content/10.1101/2021.05.12.443817v1)

This project is part of my postdoctoral research in the [Lieberman Lab](http://lieberman.science).



# Overview

The analyses are grouped into three parts:

* [_C. acnes_ genomic analysis](#c-acnes-genomic-analysis)
    * [Snakemake pipelines](#snakemake-pipelines)
    	* Mapping short reads to a reference genome
    	* Building and analyzing pan-genomes for each lineage
    	* Analyzing plasmid architecture
    * [Matlab analysis scripts](#matlab-analysis-scripts)
    	* SNV-based analyses
    	* Gene content analyses
* [Mathematical models](#mathematical-models)
    * [Pore colonization](#pore-colonization)
    * [Competition on skin](#competition-on-skin)
* [Other analyses](#other-analyses)
    * [_C. acnes_ growth curve analysis](#c-acnes-growth-curve-analysis)
    * [_C. granulosum_ analysis](#c-granulosum-analysis)
    * [16S data analysis](#16S-data-analysis)
    * [Reanalysis of metagenomic data from the literature](#reanalysis-of-metagenomic-data-from-the-literature)

See below for information on [Data availability](#data-availability) and for [Acknowledgments](#acknowledgments).



# _C. acnes_ genomic analysis
Directory: `/cacnes_genomic_analysis`

This directory includes all analyses derived from _Cutibacterium acnes_ sequencing data.

The first phase of analyses (Snakemake pipilines) uses well-known computational tools (e.g. cutadapt, sickle, bowtie2, samtools, SPAdes) to process raw sequencing reads into tables of candidate SNVs along with related metrics like read coverage at each position on the genome. The second phase of analyses (Matlab analysis scripts) performs finer filtering of pre-processed data in order to produce SNV tables, evolutionary analysis, etc. This second phase ultimately produces all figures presented in the manuscript.


### Snakemake pipelines
Directory: `/cacnes_genomic_analysis/snakemakes`

Computational pipelines are implemented in [snakemake](https://snakemake.readthedocs.io/en/stable/) and are designed for use on a [slurm](https://slurm.schedmd.com/documentation.html) cluster.

Mapping short reads to a reference genome and identifying candidate SNVs:
* `/cacnes_genomic_analysis/snakemakes/1_process_raw_data`: removes adapters and trims reads based on quality
* `/cacnes_genomic_analysis/snakemakes/2_bracken`: estimates _C. acnes_ abundance for reads from each colony
* `/cacnes_genomic_analysis/snakemakes/3_refgenome_mapping`: maps reads from each colony to a _C. acnes_ reference genome
* `/cacnes_genomic_analysis/snakemakes/4_refgenome_case`: generates data tables over preliminary variant positions

Assembling and analyzing pan-genomes for each lineage:
* `/cacnes_genomic_analysis/snakemakes/5_lineage_assembly_mapping`: assembly of pan-genomes for each lineage; mapping of reads onto assembled pan-genomes
* `/cacnes_genomic_analysis/snakemakes/6_lineage_assembly_case`: generates table of read coverage over assembled pan-genomes for each member colony
* `/cacnes_genomic_analysis/snakemakes/7_lineage_assembly_cdhit`: identifies gene clusters present across lineage pan-genomes

Analyzing plasmid architecture:
* `/cacnes_genomic_analysis/snakemakes/8_plasmid_hybridassembly`: hybrid assemblies (short read and long read data) of a subset of colonies containing plasmids
* `/cacnes_genomic_analysis/snakemakes/9_plasmid_alignments`: alignments of short reads onto plasmid scaffolds produced by the hybrid asembly


### Matlab analysis scripts
Directory: `/cacnes_genomic_analysis/matlab`

Note: Matlab analyses depend on data from their respective snakemake pipeline(s) as indicated below.

SNV-based analyses: `/cacnes_genomic_analysis/matlab/1_snv_analysis/genomic_analysis_main.m`
* Clustering of colonies into lineages
* Within-lineage SNV identification
* Spatial biogeography analysis (across skin regions and within pores)
* Evolutionary analysis (search for parallel evolution; inference of ancestral alleles; treemaking)
* and more!
Snakemake dependencies: `/cacnes_genomic_analysis/snakemakes/4_refgenome_case`

Gene content analyses: 
* Within-lineage gains/losses: `/cacnes_genomic_analysis/matlab/2_within-lineage-gain-loss`
    * Snakemake dependencies: `/cacnes_genomic_analysis/snakemakes/6_lineage_assembly_case`
* Between-lineage gene cluster differences: `/cacnes_genomic_analysis/matlab/3_btwn-lineage-gene-content`
    * Snakemake dependencies: `/cacnes_genomic_analysis/snakemakes/7_lineage_assembly_cdhit`
* Plasmid presence and evolutionary analysis: `/cacnes_genomic_analysis/matlab/4_plasmid-presence`
    * Snakemake dependencies: `/cacnes_genomic_analysis/snakemakes/9_plasmid_alignments`



# Mathematical Models
Directory: `/mathematical_models`


### Pore colonization
Directory: `/mathematical_models/pore_colonization_simulation`

Simulations of pore colonization can be found here.


### Competition on skin
Directory: `/mathematical_models/skin_competition_model`

Mathematical models of competition between _C. acnes_ strains on the skin can be found here.



# Other analyses
Directory: `/other_analyses`


### _C. acnes_ growth curve analysis
Directory: `/other_analyses/cacnes_growth_curves`

Analysis of growth rates of _C. acnes_ isolates measured _in vitro_. 


### _C. granulosum_ analysis
Directory: `/other_analyses/c_granulosum`

_This section is under construction._

Evolutionary analysis of _C. granulosum_ colonies:
* `/other_analyses/c_granulosum/refgenome`: alignments of short reads to a _C. granulosum_ reference genome, SNV calling, and evolutionary analysis
* `/other_analyses/c_granulosum/assemblies`: assembly of sample-specific _C. granulosum_ genomes, SNV calling, and evolutionary analysis


### 16S data analysis
Directory: `/other_analyses/16s_analysis`

_This section is under construction._

Analysis of 16S data:
* `/other_analyses/16s_analysis/snakemake`: classifies 16S reads with QIIME2/DADA2
* `/other_analyses/16s_analysis/matlab`: analyzes 16S data and makes plots


### Reanalysis of metagenomic data from the literature
Directory: `/other_analyses/reanalysis_oh_2014`

Re-analysis of spatial biogeography from [Oh et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185404/), with reads characterized in terms of strain-types proposed by [Scholz et al., 2014](https://pubmed.ncbi.nlm.nih.gov/25111794/).



# Data availability

Data relating to these analyses are available from the following sources:

* Raw sequencing reads: SRA Accession [PRJNA771717](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA771717) 
* Assembled pan-genomes for each lineage: _not yet available_ <!-- [TODO](TODO). -->
* Matlab data files: _not yet available_ <!-- [TODO](TODO). -->



# Acknowledgments

* Ravalika Damerla wrote the plasmid hybrid assembly code.
* Alex Poret wrote the _C. granulosum_ analysis code and the 16S data analysis pipeline.
* A. Delphine Tripp wrote the growth curve analysis code.


