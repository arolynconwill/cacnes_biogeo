function [ contigs_to_keep, assembled_genome_length, eff_genome_length, locustags_to_keep ] = ...
    contig_filter_lineage( this_clade_num, this_clade_name, dir_assemblies, ...
    min_contig_copynum_to_keep, min_aa_length_to_keep )

%% SUMMARY

% This script filters contigs for looking at between lineage gene content
% differences.


%%

fprintf(1,['Filtering contigs for ' this_clade_name '!' '\n'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SETUP %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Transfer data for this cluster

% Paths to downloaded files
path_genome = [ dir_assemblies '/' this_clade_name '/genome.fasta' ];
dir_genome = [ dir_assemblies '/' this_clade_name ];
path_annotations = [ dir_assemblies '/' this_clade_name '/prokka_out.gbk' ];
path_plasmidblast =  [ dir_assemblies '/' this_clade_name '/' 'clade_' num2str(this_clade_num) '_plasmids.csv' ];
path_covmat =  [ dir_assemblies '/' this_clade_name '/coverage_matrix.mat' ];

% Warn if files missing
if ~exist( path_genome, 'file' ) || ~exist( path_annotations, 'file' ) || ~exist( path_plasmidblast, 'file' ) || ~exist( path_covmat, 'file' )
    fprintf(1, 'Warning! File(s) missing.\n' )
end


%% Load data for this cluster
fprintf(1,['Loading data...' '\n'])

% Get genome info
% Genome stats
[ChrStarts, GenomeLength, ~, ~] = genomestats(dir_genome); assembled_genome_length = GenomeLength;
contig_num_by_index = p2chrpos(1:1:GenomeLength,ChrStarts); contig_num_by_index = contig_num_by_index(:,1);
contig_start_indices = ChrStarts+1; % since ChrStarts starts with 0
contig_end_indices = [ ChrStarts(2:end), GenomeLength ];
% Fasta file
genome_info = fastaread( path_genome );
genome_contig_names = { genome_info.Header };
genome_contig_seqs = { genome_info.Sequence };
genome_contig_lengths = cellfun(@(x) numel(x), genome_contig_seqs );
clear genome_info
% Annotations
if ~exist([dir_genome '/cds_sorted.mat'], 'file')
    read_gb_prokka_aro( dir_genome )
end % makes cds_sorted.mat if it doesn't already exist
load([dir_genome '/cds_sorted.mat']) 

% Load coverage matrix
load( path_covmat ) % all_coverage_per_bp, SampleNames
clear all_maf_per_bp % don't need
clear cov_modes % don't need
num_samples = numel(SampleNames);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% GENERATE COVERAGE MATRICES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute coverage statistics and normalized coverage matrix
fprintf(1,['Generating coverage matrices...' '\n'])

% Option to ignore ends of genome
start_pos=1; 
end_pos=size(all_coverage_per_bp,2)-start_pos+1;
% Coverage matrix
test_mat_cov = all_coverage_per_bp(:,start_pos:end_pos);

% Raw overage stats by sample (w/ plasmid masking)
avg_cov_by_sample = mean(test_mat_cov,2);
stdev_cov_by_sample = std(double(test_mat_cov),0,2);

% Also compute copy number matrix (w/ plasmid masking)
median_cov_by_sample = median(test_mat_cov,2);
all_coverage_per_bp_copynum = double(all_coverage_per_bp)./double(median_cov_by_sample);
test_mat_copynum = all_coverage_per_bp_copynum(:,start_pos:end_pos);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FIND CONTIGS WITH LOW COVERAGE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mean copy number of each contig
fprintf(1,['Filtering contigs...' '\n'])

% Number of contigs in this assembly
num_contigs = numel( contig_start_indices );

% Mean copy number of each contig
contig_mean_copynums = zeros( num_contigs,1 );
for c=1:num_contigs
    contig_mean_copynums(c) = mean(mean( all_coverage_per_bp_copynum(:,contig_start_indices(c):contig_end_indices(c)) ));
end

% Contigs to keep
contigs_to_keep = find( contig_mean_copynums >= min_contig_copynum_to_keep );

% "Genome length" of contigs to keep
eff_genome_length = sum( genome_contig_lengths( contigs_to_keep ) );
eff_genome_length_rel_to_CacnesP1 = eff_genome_length/2519002;

% Plot of contig copy numbers
fig=figure(1);
clf(1)
hold on
box on
histogram( contig_mean_copynums, 0:0.1:3.5 )
line( [min_contig_copynum_to_keep min_contig_copynum_to_keep], ylim, 'Color', 'k', 'LineWidth', 2 )
xlabel('contig copy number')
ylabel('num of contigs')
title(['lineage: ' this_clade_name])
text( 1.5, 0.8*max(ylim), [ 'Effective genome length: ' num2str(eff_genome_length) ] )
text( 1.5, 0.75*max(ylim), [ '...relative to Cacnes_P1: ' num2str(eff_genome_length_rel_to_CacnesP1) ], 'Interpreter', 'none' )
hold off

% Save plot
dir_plots = 'figures';
if ~exist( [pwd '/' dir_plots], 'dir' )
    mkdir( [pwd '/' dir_plots] )
end
print(fig,[ dir_plots '/Contig_Copynums_' this_clade_name '.png' ],'-r400','-dpng')

% Plot of contig copy numbers
fig=figure(2);
clf(2)
hold on
box on
scatter( genome_contig_lengths, contig_mean_copynums )
xlabel('contig length')
ylabel('contig copy number')
title(['lineage: ' this_clade_name])
axes('Position',[.5 .5 .4 .4])
box on
scatter( genome_contig_lengths, contig_mean_copynums )
xlim([500 5000])
hold off

% Save plot
dir_plots = 'output_figures';
if ~exist( [pwd '/' dir_plots], 'dir' )
    mkdir( [pwd '/' dir_plots] )
end
print(fig,[ dir_plots '/Contig_CopynumsLengths_' this_clade_name '.png' ],'-r400','-dpng')



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% GET LIST OF GENES NOT REMOVED %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Concatenate locus tags across contigs to keep
locustags_to_keep = {};
for c=1:numel(contigs_to_keep)
    next_contig_to_keep = contigs_to_keep(c);
    next_CDS = CDS{next_contig_to_keep};
    if ~isempty(next_CDS) % some contigs might not have any genes
        next_locustags_on_contig = {next_CDS.locustag};
        next_locustags_on_contig_aas = {next_CDS.translation};
        next_locustags_on_contig_length = cellfun(@(x) numel(x), next_locustags_on_contig_aas );
        next_locustags_to_keep = next_locustags_on_contig( next_locustags_on_contig_length >= min_aa_length_to_keep ); % keep if above 50 aa
        locustags_to_keep = horzcat( locustags_to_keep, next_locustags_to_keep );
    end
end

% Save data
dir_locustags = 'output_locustags';
if ~exist( [pwd '/' dir_locustags], 'dir' )
    mkdir( [pwd '/' dir_locustags] )
end
save( [ dir_locustags '/LineageGenomeInfo_' this_clade_name '.mat' ], 'contigs_to_keep', 'eff_genome_length', 'locustags_to_keep' );


end

