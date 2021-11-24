function [ all_regions_table, num_regions_final ] = ...
    gainloss_lineage( this_clade_num, this_clade_name, ...
    dir_data_assemblies, dir_tree_info, dir_save_regions )

%% Summary

% This function detects gain/loss regions in lineage

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SETUP %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Where to find data on lineage pan-genome assemblies

% Paths to downloaded files
path_genome = [ dir_data_assemblies '/' this_clade_name '/genome.fasta' ];
dir_genome = [ dir_data_assemblies '/' this_clade_name ];
path_annotations = [ dir_data_assemblies '/' this_clade_name '/prokka_out.gbk' ];
path_plasmidblast =  [ dir_data_assemblies '/' this_clade_name '/' 'clade_' num2str(this_clade_num) '_plasmids.csv' ];
path_covmat =  [ dir_data_assemblies '/' this_clade_name '/coverage_matrix.mat' ];


%% Load data for this cluster
fprintf(1,['Loading data for ' this_clade_name '...' '\n'])

% Get genome info
% Genome stats
[ChrStarts, GenomeLength, ~, ~] = genomestats(dir_genome);
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
% Note: this ONLY includes samples that passed strict purtiy filtering
% % i.e. not all colonies in the cluster



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PRE-PROCESSING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Mask plasmid positions in coverage matrix
fprintf(1,['Finding plasmid regions to mask...' '\n'])

plasmid_pos_on_genome = find_plasmid_pos_on_genome( path_plasmidblast, path_genome );
%p2chrpos(plasmid_pos_on_genome,ChrStarts) %chrpos2index
%plasmid_pos_on_genome = []; % for testing: turn off plasmid filtering

% Set coverage to zero in these regions
all_coverage_per_bp(:,plasmid_pos_on_genome) = 0;


%% Reorder according to SNP tree

% Where to find file
path_treeorder = [ dir_tree_info '/' 'clade_' num2str(this_clade_num) '_treeorder.txt' ];

% If exists, load file and reorder
if exist( path_treeorder, 'file' )
    % Get sample names in treeorder from text file
    fid = fopen( path_treeorder );
    temp=textscan(fid,'%s'); temp=temp{1,1}; temp=cellfun(@(x) [x(1) ':' x(3:end)], temp, 'UniformOutput', false );
    fclose(fid);
    temp = flipud(temp); % reverse order so first sample at top of plot not at 1
    % Reorder SampleNames and coverage matrix
    [ ~, reorder ] = sort( cellfun(@(x) find(ismember(temp,x)), SampleNames ) );
    if numel(reorder)<numel(SampleNames)
        fprintf(1,'Missing sample names!\n')
    end
    SampleNames = SampleNames(reorder)';
    all_coverage_per_bp = all_coverage_per_bp( reorder,: );
else
    fprintf(1,'Warning! No treeorder file for this clade!\n' )
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FIND CANDIDATE GAINS/LOSSES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Compute coverage statistics and normalized coverage matrix
fprintf(1,['Generating coverage matrices...' '\n'])

% Option to ignore ends of genome
start_pos=1; 
end_pos=size(all_coverage_per_bp,2)-start_pos+1;
% Coverage matrix
test_mat_cov = all_coverage_per_bp(:,start_pos:end_pos);

% Raw overage stats by sample (w/ plasmid masking)
test_mat_cov_noplasmid = test_mat_cov(:,~ismember(1:1:GenomeLength,plasmid_pos_on_genome));
avg_cov_by_sample = mean(test_mat_cov_noplasmid,2);
stdev_cov_by_sample = std(double(test_mat_cov_noplasmid),0,2);

% Double-normalized coverage matrix
% First normalize by sample
difference_mean_in_sdunits = (double(test_mat_cov) - repmat(avg_cov_by_sample,1,size(test_mat_cov,2))) ./ repmat(stdev_cov_by_sample,1,size(test_mat_cov,2)) ;
% Then normalize by position normalize by sample
avg_cov_bypos = mean(difference_mean_in_sdunits);
test_mat_doublenorm = (double(difference_mean_in_sdunits) - repmat(avg_cov_bypos,size(test_mat_cov,1),1)) ;
% Make a version with genome indexing
buffer_mat = zeros( num_samples, start_pos-1 );
test_mat_doublenorm_buffered = horzcat( horzcat( buffer_mat, test_mat_doublenorm ), buffer_mat );

% Also compute copy number matrix (w/ plasmid masking)
median_cov_by_sample = median(test_mat_cov_noplasmid,2);
all_coverage_per_bp_copynum = double(all_coverage_per_bp)./double(median_cov_by_sample);
test_mat_copynum = all_coverage_per_bp_copynum(:,start_pos:end_pos);


%% Find candidate gain/loss events
fprintf(1,['Finding candidate gain/loss events...' '\n'])

% Parameters: regions absent
min_size_loss = 500; % bp
loose_threshold_for_loss = 0.25; % copy number
max_cutoff_for_loss = 0.15; % copy number
max_avg_cov_for_loss = 2.5; % reads
min_avg_copynum_for_loss_high_control = 0.85; % copy number

% Parameters: regions present
min_size_gain = min_size_loss; % bp
loose_threshold_for_gain = 0.33; % copy number % 0.5
min_cutoff_for_gain = min_avg_copynum_for_loss_high_control;  % copy number
min_avg_cov_for_gain = 20; % reads
max_avg_copynum_for_gain_low_control = max_cutoff_for_loss; % copy number
max_avg_cov_for_gain_low_control = max_avg_cov_for_loss; % copy number

% Parameters: testing contigs
max_contig_len_to_test = GenomeLength; % test all contigs (even though obviously long ones won't show up as gain/losses)

% Sample filter
min_cov_to_eval_sample = 20; % reads

% Find positions/samples with high/low coverage
has_high_coverage = test_mat_copynum > loose_threshold_for_gain; 
has_low_coverage = test_mat_copynum < loose_threshold_for_loss;
has_high_coverage(:,1:2)=0;  has_high_coverage(:,end-1:end)=0; %force the first and last positions to be 0 so that there are the same number of starts and ends
has_low_coverage(:,1:2)=0;  has_low_coverage(:,end-1:end)=0; %force the first and last positions to be 0 so that there are the same number of starts and ends

% Initialize
copy_number_variants_strains=[];
copy_number_variants_ends=[];
copy_number_variants_starts=[];
isdel=[];

% Find candidate gains and losses
for i=1:num_samples
    fprintf(1,['Sample ' num2str(i) '/' num2str(num_samples) ': ' SampleNames{i} '\n'])
    
    if avg_cov_by_sample(i)>min_cov_to_eval_sample % only looking at strains with high enough coverage

        % Candidate gain/loss regions in this sample: 
        del_starts_0=find(diff(has_low_coverage(i,:))>0)+1; % diff is one shorter 
        del_ends_0=find(diff(has_low_coverage(i,:))<0);
        % Impose contig boundaries
        [ del_starts, del_ends ] = break_at_contig_boundaries( del_starts_0, del_ends_0, ...
            contig_num_by_index, contig_start_indices, contig_end_indices, start_pos );
        % Also test all contigs with length under min_contig_len_to_test
        del_starts = [ del_starts, contig_start_indices( genome_contig_lengths <= max_contig_len_to_test ) ];
        del_ends = [ del_ends, contig_end_indices( genome_contig_lengths <= max_contig_len_to_test ) ];
        % Filter candidate losss
        for j=1:numel(del_starts)
            if (del_ends(j) - del_starts(j) + 1) >= min_size_loss ... %length
                    && mean(test_mat_copynum(i,del_starts(j):del_ends(j))) <= max_cutoff_for_loss ... %normalized 
                    && mean(test_mat_cov(i,del_starts(j):del_ends(j))) <= max_avg_cov_for_loss ...%raw data
                    && max( mean( test_mat_copynum(:,del_starts(j):del_ends(j)),2 ) ) >= min_avg_copynum_for_loss_high_control % require another sample to have the region
                
                % Record candidate loss
                copy_number_variants_strains(end+1)=i;
                copy_number_variants_starts(end+1)=del_starts(j);
                copy_number_variants_ends(end+1)=del_ends(j);
                isdel(end+1)=1;
                
            end
        end
        
    end
    
end

fprintf(1,['Number of candidate gain/loss regions: ' num2str(numel(copy_number_variants_strains)) '.\n'] )


%% Generate clickable plot

% Convert back to genome position
if numel(copy_number_variants_strains)>0
    div_clickable_scatter_coverage( ...
        copy_number_variants_strains, copy_number_variants_starts+start_pos-1, copy_number_variants_ends+start_pos-1, isdel, ...
        SampleNames, all_coverage_per_bp, all_coverage_per_bp_copynum, test_mat_doublenorm_buffered, ...
        this_clade_name, GenomeLength, ChrStarts, genome_contig_seqs, CDS, ...
        dir_save_regions )
end


%% Collapse overlapping regions across samples
fprintf(1,['Collapsing gains/loss regions...' '\n'] )

% Boolean matrix of all regions
all_regions_bool = zeros( num_samples, GenomeLength, 'logical' ); % initialize
start_pos_genome = copy_number_variants_starts+start_pos-1;
end_pos_genome = copy_number_variants_ends+start_pos-1; 
for r0=1:numel(start_pos_genome)
    all_regions_bool( copy_number_variants_strains(r0), start_pos_genome(r0):end_pos_genome(r0) ) = 1;
end
% Detect all collapsed regions
all_regions_bool_combined = sum(all_regions_bool);
all_regions_bool_combined( all_regions_bool_combined>1 ) = 1;
all_regions_starts_0 = find( diff(all_regions_bool_combined) > 0 ) + 1;
all_regions_ends_0 = find( diff(all_regions_bool_combined) < 0 );
if all_regions_bool_combined(end)==1 % case where region includes last position on genome
    all_regions_ends_0 = [ all_regions_ends_0, GenomeLength ];
end
% Split across contig boundaries if regions on adjacent contigs happened to merge
[ all_regions_starts, all_regions_ends ] = break_at_contig_boundaries( all_regions_starts_0, all_regions_ends_0, ...
    contig_num_by_index, contig_start_indices, contig_end_indices, 1 );

% Print number of total regions
num_regions_final = numel(all_regions_starts);
fprintf(1,['Number of collapsed gain/loss regions: ' num2str(num_regions_final) '.\n'] )


%% Save info for all regions found
fprintf(1,['Recording info on gains/loss regions...' '\n'] )


% Make a table make a info on all regions
% Make a fasta with all region sequences
% Make a summary plot for each region

% Initialize
all_regions_table = struct;

% Loop through regions
for r=1:numel(all_regions_starts)
    
    % Basic info on region
    region_name = [ 'C-' num2str(this_clade_num) '_reg-' num2str(r) ];
    start_index = all_regions_starts(r);
    end_index = all_regions_ends(r);
    contig_index = contig_num_by_index(start_index);
    start_contig_pos = p2chrpos(start_index,ChrStarts);
    end_contig_pos = p2chrpos(end_index,ChrStarts);
    
    % Make a table
    % Names
    all_regions_table(r).Name = region_name;
    all_regions_table(r).Cluster = this_clade_name;
    % Compute how many colonies have this region present
    bool_samples_pos = ( mean( all_coverage_per_bp_copynum(:,start_index:end_index),2 ) >= min_cutoff_for_gain );
    bool_samples_neg = ( mean( all_coverage_per_bp_copynum(:,start_index:end_index),2 ) <= max_cutoff_for_loss );
    bool_samples_ambig = ~bool_samples_pos & ~bool_samples_neg; 
    num_samples_pos = sum( bool_samples_pos );
    num_samples_neg = sum( bool_samples_neg );
    all_regions_table(r).NumColonies_Pos = num_samples_pos;
    all_regions_table(r).NumColonies_Neg = num_samples_neg;
    all_regions_table(r).NumColonies_Ambig = num_samples - num_samples_pos - num_samples_neg;
    all_regions_table(r).NamesColonies_Pos = SampleNames( bool_samples_pos );
    all_regions_table(r).NamesColonies_Neg = SampleNames( bool_samples_neg );
    all_regions_table(r).NamesColonies_Ambig = SampleNames( bool_samples_ambig );
    % Length and position
    all_regions_table(r).Region_Length = end_index-start_index+1;
    all_regions_table(r).Genome_IndexStart = start_index;
    all_regions_table(r).Genome_IndexEnd = end_index;
    all_regions_table(r).Contig_Num = contig_index;
    all_regions_table(r).Contig_Length = genome_contig_lengths( contig_index );
    all_regions_table(r).Contig_PosStart = start_contig_pos(2);
    all_regions_table(r).Contig_PosEnd = end_contig_pos(2);
    % Gene content
    contig_seq = genome_contig_seqs{ contig_index };
    region_seq = contig_seq( start_contig_pos(2):end_contig_pos(2) );
    all_regions_table(r).Sequence = region_seq;
    [ genes_annotations, genes_translations ] = get_gene_annnotations( all_regions_starts(r), all_regions_ends(r), ChrStarts, CDS );
    all_regions_table(r).Gene_Annotations = genes_annotations;
    all_regions_table(r).Gene_AAs = genes_translations;
    % Alignment info
    all_regions_table(r).Align_PacnesC1 = 'TBD'; % fill in later
    all_regions_table(r).Align_WrongTaxa = 'TBD'; % fill in later

    % Make a fasta file
    fasta_filename = [ dir_save_regions '/' region_name '.fasta' ];
    if exist( fasta_filename, 'file' ) % remove if already exists
        delete( fasta_filename )
    end
    fastawrite( fasta_filename, region_name, region_seq );

    % Make a summary plot
    fig_num = 10;
    figure(fig_num); clf(fig_num); hold on;
    st=suptitle([ 'region: ' region_name ] );
    st.Interpreter = 'none';
    st.FontSize=20;
    subplot_num = 3;
    color_pos = rgb('Blue');
    color_neg = rgb('Red');
    color_ambig = 0.15*[ 1 1 1 ] ;
    % 1. Absolute coverage plot
    subplot_index = 1;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        all_coverage_per_bp_copynum, 'copy number', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % 2. Copy number coverage plot
    subplot_index = 2;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        all_coverage_per_bp, 'raw coverage', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % 3. Normalized coverage plot
    subplot_index = 3;
    make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
        bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
        test_mat_doublenorm_buffered, 'normalized coverage', ...
        start_index, end_index, contig_start_indices, GenomeLength, ...
        color_pos, color_neg, color_ambig )
    % Save summary plot
    print([dir_save_regions '/' region_name '.png'],'-dpng')

end

% Save table
region_set_name = [ 'C-' num2str(this_clade_num) '_regions' ];
save( [ dir_save_regions '/' region_set_name '.mat' ], 'all_regions_table' )



%%

fprintf(1,['Done!' '\n\n\n'])


end

