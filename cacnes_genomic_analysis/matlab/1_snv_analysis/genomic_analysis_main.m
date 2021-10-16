%%%%%%%%%%%%%
%% SUMMARY %%
%%%%%%%%%%%%%

% This script performs the C. acnes SNV-based analysis


%% Setup

% Add scripts path
path(path,'scripts')



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% CLUSTERING COLONIES INTO LINEAGES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Loads data from alignments to reference genome (snakemake step: 4_refgenome_case)
% 2. Loads sample metadata
% 3. Applies quality filters to samples
% 4. Performs preliminary SNV calling and generates a distance matrix
% 5. Clusters colonies into lineages
% 6. Generates summary plots
% 7. Creates input files for lineage pan-genome assemblies

mkdir( '1_clustering' )
identify_clusters_aro


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IDENTIFYING DE NOVO SNVs IN EACH LINEAGE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Identifies SNVs differentiating colonies within the same lineage
% 2. Infers the ancestral genotype for each lineage using an outgroup
% 2. Builds a tree for each lineage

mkdir( '2_snvs' )
identify_de_novo_muts_aro


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIOGEOGRAPHY ANALYSIS: STRAIN TYPES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Makes strain-type and lineage bar charts (Figs 2 and 3)
mkdir( '3_straintype_biogeo' )
make_lineage_barcharts_1
make_lineage_barcharts_2

% 2. Make visualization of which lineages are in which pores from pore
% strip data (Fig 3 and Supp)
mkdir( '3_straintype_porestripvis' )
make_pore_coordinate_plots

% 3. Compare lineages found on faces vs backs on individual subjects (Supp)
mkdir( '3_straintype_faceback' )
compare_faceback


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BIOGEOGRAPHY ANALYSIS: PORES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir( '4_pore_biogeo' )

% Characterize low lineage-level diversity inside pores (Fig 3)
investigate_lineages_in_pores

% Investigate intra-pore and inter-pore dMRCAs (Fig 5)
analyze_pore_mrcas


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARALLEL EVOLUTION ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir( '5_parallel_evo' )

% 1. Compute observed mutational spectrum (needed for parallel evo
% analysis)
mutspec_general % observed mutation spectrum
mutspec_hypers % also look at mutation spectrum of hypermutator clade (Supp)

% 2. Look for signatures of parallel evolution (Supp)
parallel_evo_analysis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MISCELLANEOUS ANALYSES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Molecular clock analysis
mkdir( '6_misc_molec_clock' )
molecular_clock_main

% Spatial confinement analysis
mkdir( '6_misc_spatial_confinement' )
spatial_confinement_coordinates
spatial_confinement_faceback
spatial_confinement_facezone


%%%%%%%%%%%%%%%%
%% TREEMAKING %%
%%%%%%%%%%%%%%%%

% Phylogeny for each lineage
mkdir('7_trees_lineages')
treemaking_lineages

% Phylogeny for each strain type
mkdir('7_trees_straintypes')
treemaking_straintypes

% Phylogeny for lineage ancestors
mkdir('7_trees_lineageancestors')
treemaking_lineage_ancestors
% Note: This script generates a directory that can be run on the cluster.


%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPLEMENTAL TABLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Makes supplemental tables for SNV-based analyses
mkdir( '8_tables' )

% Sampling summary
make_sampling_table
make_sampling_table_fig

% Lineage info
make_lineage_table

% Colony metadata and mutation table
make_lineage_mut_tables

% Recombination regions
analyze_recombo_positions

