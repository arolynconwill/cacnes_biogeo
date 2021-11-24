%% SUMMARY

% This script determines which samples have plasmids, evaluates the plasmid
% type (based on alilgnments to different plasmid scaffolds from hybrid
% assemblies), and infers evolutionary relationships among plasmid
% genotypes.



%%

%%%%%%%%%
% SETUP %
%%%%%%%%%

% set up directories, import basic info about samples, etc.


%% SET UP DIRECTORIES

% Local directory names
dir_hybrid_scaffolds = 'data_hybrid_scaffolds'; % hybrid assembled scaffolds
dir_case_step = 'data_case_step'; % short reads mapped onto hybrid assembled plasmids
dir_samples = 'data_samples'; % coverage info from alignments to the reference genome
dir_lineages = 'data_lineage_trees'; % info about lineage phylogenies

% Add scripts to path
path(path,[ pwd '/' 'miniscripts' ]) % functions for this analysis
path(path,'../lab_scripts') % functions from shared lab folder

% Where to find info on each cluster
dir_clusters = '/Users/arolyn/Dropbox (MIT)/Postdoc/Pacnes_Biogeo/ANALYSIS_GITHUB/cacnes_genomic_analysis/matlab/1_snv_analysis/2_snvs';



%% OTHER BASICS

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';


%% LINEAR PLASMID TYPES

% Query list of names of hybrid asssembled scaffolds
hybrid_scaffold_dir = dir( dir_hybrid_scaffolds );
hybrid_scaffold_list = { hybrid_scaffold_dir(:).name };
hybrid_scaffold_list = hybrid_scaffold_list( cellfun(@(x) x(1)=='A', hybrid_scaffold_list) );

% Number of hybrid assembled scaffolds
num_hybrid_scaffolds = numel( hybrid_scaffold_list );


%% NUMBER OF COLONIES

% Query list of sample names 
load( [ dir_case_step '/' 'case_' hybrid_scaffold_list{1} '/' 'candidate_mutation_table.mat' ], 'SampleNames' )
num_samples_case = numel( SampleNames ); % includes all colonies that passed filtering


%% LOAD COVERAGE DATA FROM REFERENCE GENOME ALIGNMENTS

% Load simple sample names
load([ dir_samples '/' 'sample_names.mat'],'SampleNamesSimple_all' );
% Load cluster info
load([ dir_samples '/' 'cluster_names.mat'],'cluster_names_new' );
% Load data (chromosome)
load( [ dir_samples '/' 'chromosomal_coverage.mat' ] )
coverage_chromosome_all = coverage_chromosome_median; clear coverage_chromosome_median; % rename
% Load data (regB)
load( [ dir_samples '/' 'regB_coverage.mat' ] )
coverage_regB_all = coverage_regB_median_norm; clear coverage_regB_median_norm; % rename

% Save coverage data and long names in same order as case step data
% Initialize 
SampleNames_keep = cell( num_samples_case,1 ); % short sample names
SampleNamesSimple_keep = cell( num_samples_case,1 ); % simple sample names
SampleNamesLong_keep = cell( num_samples_case,1 ); % long sample names
chromosomal_coverage = zeros(num_samples_case,1) ; % median coverage of whole chromosome
regB_coverage = zeros(num_samples_case,1); % median coverage of "region B" chromosomal region
sample_indices_scaffcase = zeros(num_samples_case,1); % index of sample in candidate mutation table from plasmid data
sample_indices_snpcase = zeros(num_samples_case,1); % index of sample in candidate mutation table from reference genome data
% Fill in info for each sample
for n=1:num_samples_case
    next_name = SampleNames{n};
    next_name_index = find( ismember( SampleNames_all, next_name ) );
    if ~isempty( next_name_index ) % option to add additional criteria here (i.e. only sample A)
        sample_indices_scaffcase(n) = n;
        sample_indices_snpcase(n) = next_name_index;
        SampleNames_keep{n} = next_name;
        SampleNamesSimple_keep{n} = SampleNamesSimple_all{ next_name_index };
        SampleNamesLong_keep{n} = SampleNamesLong_all{ next_name_index };
        chromosomal_coverage(n) = coverage_chromosome_all( next_name_index );
        regB_coverage(n) = coverage_regB_all( next_name_index );
    else % print names that aren't found
        fprintf(1,[ ' Alert! Missing sample: ' next_name '\n' ] )
    end
end

% Number of samples to analyze
num_samples = numel( SampleNames_keep );

% Subjects
subject_membership = cellfun(@(x) x(1), SampleNamesLong_keep );

% Grab lineage info for each sample from long sample names
cluster_membership = cellfun(@(x) str2num(x(end-1:end)), SampleNamesLong_keep );
SampleNamesSimpleLong_keep = cell( size( SampleNamesSimple_keep ) );
for i=1:num_samples
    temp_cluster_index = round(cluster_membership(i));
    if temp_cluster_index>0
        SampleNamesSimpleLong_keep{i} = [ strrep( cluster_names_new{round(cluster_membership(i))}, '_', '-' ) '_' SampleNamesSimple_keep{i} ];
    else
        SampleNamesSimpleLong_keep{i} = [ 'Unclustered' '_' SampleNamesSimple_keep{i} ];
    end
end
cluster_list_all = unique( cluster_membership );

% Reorder samples by lineage and by lineage tree
reorder_by_cluster = []; % initialize
for c=1:numel(cluster_list_all)
    next_lineage = cluster_list_all(c);
    if exist( [ dir_lineages '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ], 'file' )
        % reorder large Subject A clades
        next_tree_order = importdata([ dir_lineages '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ]);
    else
        % don't bother with other clades
        next_tree_order = SampleNamesSimple_keep( cluster_membership == next_lineage );
    end
    % Find index of each name in SampleNames_keep
    for i=1:numel(next_tree_order)
        next_name = next_tree_order{i};
        next_index = find( ismember( SampleNamesSimple_keep, next_name ) );
        if ~isempty( next_index )
            reorder_by_cluster(end+1) = next_index;
        end
    end
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE SHORT READ ALIGNMENTS TO PLASMID SCAFFOLDS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% EVALUATE ALIGNMENTS OF EACH SAMPLE TO EACH SCAFFOLD
% Metrics: percent of scaffold covered, mean copy number of scaffold, and percent identity of short reads relative to hybrid assembled scaffold

% Parameters
coverage_cutoff = 0.33;
bin_size = 100; % bp

% Coverage of each scaffold for each sample
mat_mean_cov = zeros( num_hybrid_scaffolds, num_samples ); % initialize
mat_perc_cov = zeros( num_hybrid_scaffolds, num_samples ); % initialize
mat_bp_cov = zeros( num_hybrid_scaffolds, num_samples ); % initialize
% Percent identity of each scaffold for each sample
mat_perc_id = zeros( num_hybrid_scaffolds, num_samples ); % initialize
% Heatmap data
coverage_scaffold_norm_binned_all = {};

% Scaffold length
scaffold_length = []; % initialize

% Parameters
min_qual_for_call = 30; % qual = quality, specifically FQ 
min_maf_for_call = .66; % maf = major allele frequency 
min_cov_each_strand_for_call = 2; % cov = coverage 

% Directory for saving figures
dir_preliminary = 'plasmid_typing_preliminary';
if ~exist( dir_preliminary, 'dir')
    mkdir( dir_preliminary )
end

for s = 1:num_hybrid_scaffolds % next scaffold

    % Load case step data and get coverage, calls, positions
    filename = [ dir_case_step '/' 'case_' hybrid_scaffold_list{s} '/' 'candidate_mutation_table.mat' ];
    [ Calls, coverage, p ] = get_scaff_cov_and_calls( filename, min_qual_for_call, min_maf_for_call, min_cov_each_strand_for_call );
    
    %%

    % Grab calls for samples to keep
    calls_scaffold = Calls(:,sample_indices_scaffcase);
    % Compute coverage for samples to keep
    coverage_scaffold = coverage(:,sample_indices_scaffcase);
    coverage_scaffold_norm = coverage_scaffold./chromosomal_coverage';

    %%

    % Load reference genome (i.e. scaffold)
    % Define directory with reference genome
    dir_refgenomes = [ dir_hybrid_scaffolds '/' hybrid_scaffold_list{s} ];
    % Load reference genome
    [ChrStarts, GenomeLength, ~, ScafNames] = genomestats(dir_refgenomes);
    refnt_all = extract_outgroup_mutation_positions(dir_refgenomes, p2chrpos(p,ChrStarts));
    [~,refnti_all]=ismember(refnt_all,NTs); 
    % Record
    scaffold_length(s) = GenomeLength;

    %%

    % Bin normalized coverage
    num_bins = ceil(GenomeLength/bin_size);
    coverage_scaffold_norm_binned = zeros( num_samples, num_bins );
    for b=1:num_bins
        coverage_scaffold_norm_binned(:,b) = mean( coverage_scaffold_norm( bin_size*(b-1)+1:min(bin_size*b,GenomeLength), : ) );
    end
    scaffold_mean_cov = mean(coverage_scaffold_norm);

    % Percent coverage
    scaffold_bp_cov = sum( coverage_scaffold_norm > coverage_cutoff );
    scaffold_perc_cov = sum( coverage_scaffold_norm > coverage_cutoff )/GenomeLength;
    [ ~,reorder_by_cov ] = sort( scaffold_perc_cov, 'descend' );

    % Percent identity
    scaffold_perc_id = zeros( size( scaffold_perc_cov ) );
    for i=1:num_samples
        next_calls = calls_scaffold( :,i );
        next_cov = coverage_scaffold_norm(:,i);
        pos_to_compare = next_cov>coverage_cutoff & next_calls>0; % only positions with coverage and with calls
        if sum(pos_to_compare)/GenomeLength > 0.5
            num_pos_to_compare = sum( pos_to_compare );
            next_perc_id = sum( next_calls(pos_to_compare) == refnti_all(pos_to_compare ) )/num_pos_to_compare;
            scaffold_perc_id(i) = next_perc_id;
        else
            scaffold_perc_id(i) = -0.1;
        end
    end

    %% Save data in matrices
    
    mat_mean_cov(s,:) = scaffold_mean_cov;
    mat_perc_cov(s,:) = scaffold_perc_cov;
    mat_bp_cov(s,:) = scaffold_bp_cov;
    mat_perc_id(s,:) = scaffold_perc_id;
    
    coverage_scaffold_norm_binned_all{s} = coverage_scaffold_norm_binned;
    
    %% Summary plot for this scaffold
    
    for v=1:2
        
        % Coverage heatmap and percent coverage bar chart
        figure(1);
        clf(1)
        hold on
        n_cols = 6;
        % Reorder samples
        if v==1
            reorder = reorder_by_cov;
        elseif v==2
            reorder = reorder_by_cluster;
        end
        %
        % Heatmap: binned coverage over scaffold
        subplot( 2,n_cols, [ 1 2 3 7 8 9 ] )
        cov_norm_max_for_colormap = 3;
        plot_title = [ 'plasmid scaffold: ' hybrid_scaffold_list{s} ];
        plot_heatmap( coverage_scaffold_norm_binned, cov_norm_max_for_colormap, reorder, plot_title, num_bins, bin_size )
        %
        % Bar chart: mean scaffold coverage
        subplot( 2,n_cols, [ 4 10 ] )
        bar_data = scaffold_mean_cov;
        x_limits = [0 2*cov_norm_max_for_colormap];
        x_label = 'mean coverage (norm)';
        plot_title = 'coverage (mean,norm)';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        % Bar chart: percent of scaffold covered
        subplot( 2,n_cols, [ 5 11 ] )
        bar_data = scaffold_perc_cov;
        x_limits = [0 1];
        x_label = 'percent scaffold covered';
        plot_title = '% covered';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        % Bar chart: percent of scaffold identity
        subplot( 2,n_cols, [ 6 12 ] )
        bar_data = scaffold_perc_id;
        x_limits = [0.975 1];
        x_label = 'average nucleotide identity';
        plot_title = '% identity';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        hold off

        % Save plot
        if v==1
            print([dir_preliminary '/' 'preliminary_summary_' hybrid_scaffold_list{s} '_bycov.png'],'-dpng')
        elseif v==2
            print([dir_preliminary '/' 'preliminary_summary_' hybrid_scaffold_list{s} '_bylineage.png'],'-dpng')
        end
        
    end
    
end


%% PLASMID TYPING

% Approach
% % 1. Known plasmid scaffolds: "Greedy" approach evaluates match to each
% scaffold (yes/no) and chooses the longest match.
% % 2. Unknown plasmid scaffolds: If there are no good "matches" but still 
% significant scaffold coverage, then call it an unknown plasmid type.
% % 3. No plasmid: Assume everything else doesn't have a plasmid.

% Thresholds
cutoff_perc_cov = 0.9;
cutoff_perc_id = 0; 
cutoff_per_cov_transposon = 0.75;

% Initialize: 0 = no match; 1 = match
scaffold_match = zeros( num_samples, num_hybrid_scaffolds+1 );
scaffold_match_greedy = zeros( num_samples, num_hybrid_scaffolds+1 );

% Sort scaffolds by length
[ scaffold_length_sorted, scaffold_length_sorted_indices ] = sort(scaffold_length, 'descend');
% Add placeholder for unknown plasmid type
scaffold_length_sorted_indices = [ scaffold_length_sorted_indices, num_hybrid_scaffolds+1 ];
scaffold_length_sorted = [ scaffold_length_sorted, 0 ];

% Loop through scaffolds: find any match
for s=1:num_hybrid_scaffolds+1
    next_scaff = scaffold_length_sorted_indices(s);
    if s<=num_hybrid_scaffolds
        is_match = mat_perc_cov(next_scaff,:) >= cutoff_perc_cov & ...
            mat_perc_id(next_scaff,:) >= cutoff_perc_id;
    elseif s==num_hybrid_scaffolds+1
        is_match = max(mat_perc_cov(1:4,:)) > 0.75;
    end
    scaffold_match(:,next_scaff) = is_match;
end

% Loop through scaffolds: find longest match
for s=1:num_hybrid_scaffolds+1
    next_scaff = scaffold_length_sorted_indices(s);
    is_match = scaffold_match(:,next_scaff);
    already_matched = sum( scaffold_match_greedy, 2);
    scaffold_match_greedy(:,next_scaff) = is_match & ~already_matched;
end
scaffold_match_greedy = logical( scaffold_match_greedy ); % convert to booleans

% Boolean for whether or not the sample has any plasmid (anything except
% scaffold 5)
has_plasmid = (sum(scaffold_match_greedy(:,[1 2 3 4 6]),2)>0);
has_plasmid_othertype = scaffold_match_greedy(:,6);
% Boolean for whether or not the sample has the transposon (scaffold 5)
has_transposon = scaffold_match(:,5); % can be on a plasmid or not


%% Save plasmid presence and types

% Directory for saving figures and data
dir_typing = 'plasmid_typing_final';
if ~exist( dir_typing, 'dir')
    mkdir( dir_typing )
end

% Save CSV
fid = fopen( [dir_typing '/' 'data_plasmid_presence.csv'], 'w' );
fprintf( fid, [ 'lineage_name, sample_name, has_plasmid, has_plasmid_othertype, has_transposon, ' '\n' ] );
for i=1:num_samples
    i_reorder = reorder_by_cluster(i);
    if cluster_membership(i_reorder)>0
        fprintf( fid, [ cluster_names_new{cluster_membership(i_reorder)} ',' SampleNamesSimple_keep{i_reorder} ',' num2str(has_plasmid(i_reorder)) ',' num2str(has_plasmid_othertype(i_reorder)) ',' num2str(has_transposon(i_reorder)) '\n' ] );
    else
        fprintf( fid, [ 'Unclustered' ',' SampleNamesSimple_keep{i_reorder} ',' num2str(has_plasmid(i_reorder)) ',' num2str(has_plasmid_othertype(i_reorder)) ',' num2str(has_transposon(i_reorder)) '\n' ] );
    end
end
fclose( fid );


%% POST-TYPING EVALUATION

% loop through each scaffold plus NO MATCHES set

% Choose which set of matches to examine
for match_index = 1:num_hybrid_scaffolds+2 % 1-7
    
    if match_index <= num_hybrid_scaffolds
        this_set = scaffold_match_greedy(:,match_index);
        match_name = hybrid_scaffold_list{match_index};
        match_header = [ hybrid_scaffold_list{match_index} ' (scaffold P' num2str(match_index) ')'];
    elseif match_index == num_hybrid_scaffolds+1 % samples with no matches
        this_set = scaffold_match_greedy(:,match_index);
        match_name = 'new_plasmid';
    elseif match_index == num_hybrid_scaffolds+2 % samples with no matches
        this_set = sum(scaffold_match_greedy, 2)<1;
        match_name = 'no_plasmid_detected';   
    end
    num_samples_set = sum( this_set );
    if num_samples_set == 0
        fprintf(1,[ 'No matches to: ' match_name '\n' ])
        continue
    end
    this_set_names = SampleNames_keep( this_set );
    this_set_names_long = SampleNamesLong_keep( this_set );
    cluster_membership_set = cellfun(@(x) str2num(x(end-1:end)), this_set_names_long );

    reorder_set_by_cluster = []; % initialize
    for c=1:numel(cluster_list_all)
        next_lineage = cluster_list_all(c);
        if exist( [ dir_samples '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ], 'file' )
            next_tree_order = importdata([ dir_samples '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ]);
        else
            next_tree_order = this_set_names_long( cluster_membership_set == next_lineage );
        end
        % Find index of each name in SampleNames_keep
        for i=1:numel(next_tree_order)
            next_name = next_tree_order{i};
            next_index = find( ismember( this_set_names_long, next_name ) );
            if ~isempty( next_index )
                reorder_set_by_cluster(end+1) = next_index;
            end
        end
    end
    
    % Start figure
    figure(10)
    clf(10)
    hold on
    
    n_cols = 9;

    st = suptitle( [ 'SHOWING SAMPLES MATCHED TO: ' match_header ' (n= ' num2str(sum(this_set)) ')'] );
    st.Interpreter = 'none';

    for s = 1:num_hybrid_scaffolds % next scaffold

        query_header = [ hybrid_scaffold_list{s} ' (scaffold P' num2str(s) ')'];

        % Grab data for this scaffold
        coverage_scaffold_norm_binned = coverage_scaffold_norm_binned_all{s};
        scaffold_mean_cov = mat_mean_cov(s,:);
        scaffold_perc_cov = mat_perc_cov(s,:);
        scaffold_perc_id = mat_perc_id(s,:); scaffold_perc_id(this_set);
        num_bins = size( coverage_scaffold_norm_binned ); num_bins = num_bins(2);

        % Coverage heatmap and percent coverage bar chart
        % Reorder samples
        reorder = reorder_set_by_cluster;
        %
        % Heatmap: binned coverage over scaffold
        subplot( num_hybrid_scaffolds, n_cols, [ 1 2 3 ] + n_cols*(s-1) )
        cov_norm_max_for_colormap = 3;
        if match_index == s
            plot_title = [ '* BEST MATCH * ' 'plasmid scaffold: ' query_header ' * BEST MATCH *' ];
        else
            plot_title = [ 'plasmid scaffold: ' query_header ];
        end
        plot_heatmap_mini( coverage_scaffold_norm_binned(this_set,:), cov_norm_max_for_colormap, reorder, plot_title, num_bins, bin_size, this_set_names_long )
        %
        % Bar chart: mean scaffold coverage
        subplot( num_hybrid_scaffolds, n_cols, [ 4 ] + n_cols*(s-1) )
        bar_data = scaffold_mean_cov(this_set);
        x_limits = [0 2];
        x_label = '';
        plot_title = 'coverage (norm)';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        % Bar chart: percent of scaffold covered
        subplot( num_hybrid_scaffolds, n_cols, [ 5 ] + n_cols*(s-1) )
        bar_data = scaffold_perc_cov(this_set);
        x_limits = [0.5 1];
        x_label = '';
        plot_title = '% covered';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        % Bar chart: percent of scaffold identity
        subplot( num_hybrid_scaffolds, n_cols, [ 6 ] + n_cols*(s-1) )
        bar_data = scaffold_perc_id(this_set);
        x_limits = [0.99 1];
        x_label = '';
        plot_title = '% identity';
        plot_barh( bar_data, reorder, x_limits, x_label, plot_title )
        %
        % Scatter plot: % cov vs % identity
        subplot( num_hybrid_scaffolds, n_cols, [ 7 ] + n_cols*(s-1) )
        hold on
        box on
        scatter( scaffold_perc_cov, scaffold_perc_id, 20, 'k' )
        scatter( scaffold_perc_cov(this_set), scaffold_perc_id(this_set), 20, 'r' )
        xlabel('% covered')
        xlim([0.5 1])
        ylabel('% identity')
        ylim([0.975 1])
        if s<=match_index
            line( [0.5 1], [ cutoff_perc_id cutoff_perc_id ], 'Color', rgb('DeepPink'))
            line([ cutoff_perc_cov cutoff_perc_cov ], [0.975 1], 'Color', rgb('DeepPink'))        
        end
        hold off
        %
        % Line plot: normalized coverage by sample
        subplot( num_hybrid_scaffolds, n_cols, [ 8 9 ] + n_cols*(s-1) )
        plot_title = '';
        plot_linescov_mini( coverage_scaffold_norm_binned(this_set,:), cov_norm_max_for_colormap, coverage_cutoff, plot_title, num_samples_set, num_bins, bin_size );

    end

    hold off

    % Save plot
    print([dir_typing '/' 'posttyping_summary_' match_name '_matches_long.png'],'-dpng')
    
    %pause

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY FIGURES AND DATA ABOUT WHICH SAMPLES HAVE PLASMIDS AND WHICH TYPE OF PLASMID %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RECORD FINAL ASSIGNMENTS SAMPLE NAMES

% Add annotations with plasmid type to sample names
SampleNames_keep_annotated = cell(size(SampleNames_keep)); % initialize
SampleNamesLong_keep_annotated = cell(size(SampleNames_keep)); % initialize
SampleNamesSimple_keep_annotated = cell(size(SampleNames_keep)); % initialize
SampleNamesSimple_keep_annotated_short = cell(size(SampleNames_keep)); % initialize
SampleNamesSimpleLong_keep_annotated = cell(size(SampleNames_keep)); % initialize
plasmid_type = zeros(size(SampleNames_keep));
for n=1:num_samples
    next_type = find( scaffold_match_greedy(n,:) );
    if isempty( next_type )
        next_annotation = '_plasmid-none';
        next_annotation_simple = '_plasmid-0';
        next_annotation_simple_mini = '';
        plasmid_type(n) = 0;
    else
        next_annotation = [ '_plasmid-' num2str(next_type) ];
        next_annotation_simple = '_plasmid-1';
        next_annotation_simple_mini = '_*';
        plasmid_type(n) = next_type;
    end
    SampleNames_keep_annotated{n} = [ SampleNames_keep{n} next_annotation ];
    SampleNamesLong_keep_annotated{n} = [ SampleNamesLong_keep{n} next_annotation ];
    SampleNamesSimple_keep_annotated{n} = [ SampleNamesSimple_keep{n} next_annotation_simple ];
    SampleNamesSimple_keep_annotated_short{n} = [ SampleNamesSimple_keep{n} next_annotation_simple_mini ];
    SampleNamesSimpleLong_keep_annotated{n} = [SampleNamesSimpleLong_keep{n} next_annotation_simple ];
end
plasmid_type_list = unique(plasmid_type);

% Save matlab data
save([ dir_typing '/' 'data_plasmid_presence' ],'has_plasmid','has_plasmid_othertype','has_transposon',...
    'SampleNamesSimple_keep','SampleNamesSimpleLong_keep','SampleNames_keep',...
    'SampleNamesSimple_keep_annotated','SampleNamesSimple_keep_annotated_short','SampleNamesSimpleLong_keep_annotated','SampleNames_keep_annotated')


%% HEATMAPS SHOWING PLASMID ABUNDANCE IN EACH LINEAGE

% Load cluster info
load([ dir_samples '/' 'cluster_names.mat'],'cluster_names','cluster_names_new','clusters_all_slst' );
% Get super SLST for clusters
clusters_all_slst_super = cellfun(@(x) x(end-1), clusters_all_slst );
cluster_names_slst = arrayfun(@(x) [cluster_names{x} '_SLST-' clusters_all_slst{i}], 1:1:numel(cluster_names), 'UniformOutput', false);

% Add unclustered
cluster_names = [ {'Unclustered'}, cluster_names ];
cluster_names_slst = [ {'Unclustered'}, cluster_names_slst ]; 
clusters_all_slst = [ {'Unclustered'}, clusters_all_slst];
clusters_all_slst_super = [ 'U', clusters_all_slst_super ]; % U for unclustered

% Filter (toggle Subj _ vs all)
%subject_tag = 'A'; cluster_list_for_table = unique( cluster_membership( subject_membership==subject_tag ) ); 
% 'ABCDEFHIJKLMNOQR'
subject_tag = 'all'; cluster_list_for_table = cluster_list_all; % all subjects

% Toggle to reorder by SLST (approximate tree order)
sort_slst = true;
cluster_list_for_table_slst = clusters_all_slst( cluster_list_for_table+1 );
cluster_list_for_table_names = cluster_names_slst( cluster_list_for_table+1 );
if sort_slst
    [~,reorder] = sort( cluster_list_for_table_slst );
    cluster_list_for_table_names(reorder)'
else
    if ismember(0,cluster_list_for_table)
        reorder = 2:1:numel(cluster_list_for_table);
        reorder = [reorder, 1];
    else
        reorder = 1:1:numel(cluster_list_for_table);
    end
end

% Heatmap of plasmid types in each lineage
% Prep data
mat_lineage_plasmid_types = zeros( numel(cluster_list_for_table), numel(plasmid_type_list) );
% Prep and make heatmap
for c=1:numel(cluster_list_for_table)
    for p=1:numel(plasmid_type_list)
        if subject_tag=='A'
            mat_lineage_plasmid_types(c,p) = sum( ...
                cluster_membership==cluster_list_for_table(c) & ...
                plasmid_type==plasmid_type_list(p) & ...
                subject_membership == 'A' ...
                );
        elseif numel(subject_tag)==1
            mat_lineage_plasmid_types(c,p) = sum( ...
                cluster_membership==cluster_list_for_table(c) & ...
                plasmid_type==plasmid_type_list(p) & ...
                subject_membership == subject_tag ...
                );
        elseif isequal( subject_tag, 'all' )
            mat_lineage_plasmid_types(c,p) = sum( ...
                cluster_membership==cluster_list_for_table(c) & ...
                plasmid_type==plasmid_type_list(p) ...
                );
        else
            fprintf(1,['Subject tag not recognized...\n'])
        end
    end
end

% Actual heatmap
if numel(subject_tag)==1
    figure(15)
    clf(15)
else
    figure(16)
    clf(16)
end

% Colormap % colorbrewer blue
my_colors = [...
    256,256,256;...
    247,251,255;...
    222,235,247;...
    198,219,239;...
    158,202,225;...
    107,174,214;...
    66,146,198;...
    33,113,181;...
%    8,81,156;...
%    8,48,107 ...
    ]/256;
hold on
% Image
box on
mat_lineage_plasmid_types_norm = mat_lineage_plasmid_types./sum(mat_lineage_plasmid_types,2);
imagesc( flipud(mat_lineage_plasmid_types_norm(reorder,:)), [0 1] ); % normalized for each lineage
colormap( my_colors ); %colormap( summer(100) )
colorbar( 'XTick', 0:.25:1 );
% Axes and labels
set(gca, 'FontSize', 18 )
xlabel( 'plasmid type' )
xlim( [ 0.5 numel(plasmid_type_list)+0.5 ] )
xticks( 1:1:numel(plasmid_type_list) )
labels_xticks = arrayfun(@(x) num2str(x), plasmid_type_list, 'UniformOutput', false );
labels_xticks{1} = 'none'; labels_xticks{end} = 'other';
xticklabels( labels_xticks )
ylabel( 'lineage' )
ylim( [ 0.5 numel(cluster_list_for_table)+0.5 ] )
yticks( 1:1:numel(cluster_list_for_table) )
labels_yticks = cluster_list_for_table_names(reorder);
yticklabels( fliplr(labels_yticks) )
set(gca,'TickLabelInterpreter','none')
title(['subject ' subject_tag])
% Text with numbers of colonies
for c=1:numel(cluster_list_for_table)
    c_reorder = reorder(c);
    for p=1:numel(plasmid_type_list)
        text( p, numel(cluster_list_for_table)-c+1, num2str(mat_lineage_plasmid_types(c_reorder,p)), 'FontSize', 14, 'HorizontalAlignment', 'center' )
    end
end
hold off

% Save plot
if sort_slst
    print([dir_typing '/' 'overall_summary_typetable_subj-' subject_tag '_sort-slst.png'],'-dpng')
else
    print([dir_typing '/' 'overall_summary_typetable_subj-' subject_tag '.png'],'-dpng')
end

%% HEATMAPS SHOWING PLASMID ABUNDANCE IN EACH LINEAGE

% Filter (toggle Subj _ vs all)
%subject_tag = 'A'; cluster_list_for_table = unique( cluster_membership( subject_membership=='A' ) ); % subj A only
%subject_tag = 'O'; cluster_list_for_table = unique( cluster_membership( subject_membership==subject_tag ) ); % 'ABCDEFHIJKLMNOQR'
subject_tag = 'all'; cluster_list_for_table = cluster_list_all; % all subjects

% Toggle to reorder by SLST (approximate tree order)
sort_slst = true;
cluster_list_for_table_slst = clusters_all_slst( cluster_list_for_table+1 );
cluster_list_for_table_names = cluster_names_slst( cluster_list_for_table+1 );
if sort_slst
    [~,reorder] = sort( cluster_list_for_table_slst );
    cluster_list_for_table_names(reorder)';
else
    reorder = 1:1:numel(cluster_list_for_table);
end

% Heatmap of plasmid types in each lineage
% Prep data
mat_lineage_plasmid_frac = zeros( numel(cluster_list_for_table), 3 );
% Prep and make heatmap
for c=1:numel(cluster_list_for_table)
    if numel(subject_tag)==1
        num_plasmid_pos = sum( ...
            cluster_membership==cluster_list_for_table(c) & ...
            plasmid_type>0 & ...
            subject_membership == subject_tag ...
            );
        num_plasmid_either = sum( ...
            cluster_membership==cluster_list_for_table(c) & ...
            subject_membership == subject_tag ...
            );
    elseif isequal( subject_tag, 'all' )
        num_plasmid_pos = sum( ...
            cluster_membership==cluster_list_for_table(c) & ...
            plasmid_type>0 ...
            );
        num_plasmid_either = sum( ...
            cluster_membership==cluster_list_for_table(c) ...
            );
    else
        fprintf(1,['Subject tag not recognized...\n'])
    end
    mat_lineage_plasmid_frac(c,1) = num_plasmid_pos;
    mat_lineage_plasmid_frac(c,2) = num_plasmid_either;
    mat_lineage_plasmid_frac(c,3) = num_plasmid_pos/num_plasmid_either;

end

% Actual heatmap
if numel(subject_tag)==1
    figure(15)
    clf(15)
else
    figure(16)
    clf(16)
end

% Colormap % colorbrewer blue
my_colors = [...
    252,251,253; ...
    239,237,245; ...
    218,218,235; ...
    188,189,220; ...
    158,154,200; ...
    128,125,186; ...
    106,81,163; ...
    84,39,143; ...
    63,0,125; ...
    ]/256;
hold on
% Image
box on
imagesc( flipud(mat_lineage_plasmid_frac(reorder,3)), [0 1] ); % normalized for each lineage
colormap( my_colors ); %colormap( summer(100) )
colorbar( 'XTick', 0:.25:1 );
% Axes and labels
set(gca, 'FontSize', 18 )
xlabel( 'plasmid abundance' )
xlim( [ 0.5 1+0.5 ] )
xticks( [1] )
xticklabels( {''} )
ylabel( 'lineage' )
ylim( [ 0.5 numel(cluster_list_for_table)+0.5 ] )
yticks( 1:1:numel(cluster_list_for_table) )
labels_yticks = cluster_list_for_table_names(reorder);
yticklabels( fliplr(labels_yticks) )
set(gca,'TickLabelInterpreter','none')
title(['subject ' subject_tag])
% % Text with numbers of colonies
% for c=1:numel(cluster_list_for_table)
%     c_reorder = reorder(c);
%     for p=1:numel(plasmid_type_list)
%         text( p, numel(cluster_list_for_table)-c+1, num2str(mat_lineage_plasmid_types(c_reorder,p)), 'FontSize', 14, 'HorizontalAlignment', 'center' )
%     end
% end
hold off

% Save plot
if sort_slst
    print([dir_typing '/' 'overall_summary_fractable_subj-' subject_tag '_sort-slst.png'],'-dpng')
else
    print([dir_typing '/' 'overall_summary_fractable_subj-' subject_tag '.png'],'-dpng')
end

%%

if isequal( subject_tag, 'all' )
    figure(17)
    clf(17)
    hold on
    % Image
    box on
    histogram( mat_lineage_plasmid_frac(:,3), 0:0.1:1, 'FaceColor', [188,189,220]/256 ); % normalized for each lineage
    % Axes and labels
    set(gca, 'FontSize', 18 )
    xlabel( 'fraction of colonies in lineage with plasmid' )
    ylabel( 'number of lineages' )
    title(['subject all' ])
    print([dir_typing '/' 'overall_summary_fractable_subj-' subject_tag '_histogram.png'],'-dpng')
end



%%

%%%%%%%%%%%%%%%%%%%%
% MAKE PHYLOGENIES %
%%%%%%%%%%%%%%%%%%%%

% annotates lineage phylogenies with plasmid presence/absence
% make phylogeny over core plasmid architecture


%% ANNOTATE LINEAGE TREES

% Choose cluster
lineages_for_treemaking = { 'Cluster-1-A-216', 'Cluster-2-A-125', 'Cluster-3-A-63', 'Cluster-5-F-45', ...
    'Cluster-6-B-33', 'Cluster-8-F-22', 'Cluster-10-Q-20', 'Cluster-12-B-14', 'Cluster-13-B-12', ...
    'Cluster-14-R-11', 'Cluster-16-Q-10', 'Cluster-20-B-7', 'Cluster-21-F-7', 'Cluster-29-B-5', ...
    'Cluster-32-Q-5', 'Cluster-35-B-4', 'Cluster-37-F-4', 'Cluster-40-R-4', 'Cluster-48-B-3' };

% Directory for trees
dir_trees = 'trees_lineages_annotated';
if ~exist(dir_trees,'dir') % ARO: Groups subfolder
    mkdir(dir_trees)
end
copyfile( 'dnapars.app', [ dir_trees '/dnapars.app' ] )
copyfile( 'dnapars', [ dir_trees '/dnapars' ] )

% Make annotated tree for this cluster
for c=1:numel(lineages_for_treemaking)
    this_cluster_name = lineages_for_treemaking{c};
    make_tree_lineage_annotated( this_cluster_name, SampleNamesSimple_keep, SampleNamesSimple_keep_annotated_short, NTs, dir_trees, dir_clusters )
end

% Manually color-code and format trees


%% TREES FOR EACH PLASMID SCAFFOLD

% Calls filters
min_qual_for_call = 30; % qual = quality, specifically FQ 
min_maf_for_call = .66; % maf = major allele frequency 
min_cov_each_strand_for_call = 3; % cov = coverage 
% Position filters
max_fraction_ambiguous_samples = .5;

% Directory for tree
dir_trees = 'trees_plasmids_annotated';
if ~exist(dir_trees,'dir') % ARO: Groups subfolder
    mkdir(dir_trees)
end

% Toggle
make_tree = 1;

for match_index = 1:3%1:num_hybrid_scaffolds % next scaffold

    % Plasmid scaffold to examine
    this_set = scaffold_match_greedy(:,match_index);
    match_name = hybrid_scaffold_list{match_index};
    SampleNames_this_set = SampleNames_keep( this_set );
    SampleNamesSimple_this_set = SampleNamesSimple_keep( this_set );
    SamplesNamesSimple_this_set_long = SampleNamesSimpleLong_keep( this_set );
    SamplesNamesSimple_this_set_annotated = SampleNamesSimpleLong_keep_annotated( this_set );
    Nsample = numel( SampleNamesSimple_this_set );
    cluster_membership_set = cluster_membership( this_set );

    % Load case step data and get coverage, calls, positions
    filename = [ dir_case_step '/' 'case_' hybrid_scaffold_list{match_index} '/' 'candidate_mutation_table.mat' ];
    [ Calls_this_set, ~, p, counts, Quals ] = get_scaff_cov_and_calls_this_set( filename, SampleNames_this_set, ...
        min_qual_for_call, min_maf_for_call, min_cov_each_strand_for_call );
    
    % Lazy SNP calling
    % Ignore positions that don't vary
    nonvariablep=(sum(Calls_this_set==1,2)==sum(Calls_this_set~=0,2)) ...
        | (sum(Calls_this_set==2,2)==sum(Calls_this_set~=0,2)) ...
        | (sum(Calls_this_set==3,2)==sum(Calls_this_set~=0,2)) ...
        | (sum(Calls_this_set==4,2)==sum(Calls_this_set~=0,2)) ; % all samples with calls have the same base
    % Ignore positions where too many samples have N's
    ambiguousp = sum(Calls_this_set==0,2)/Nsample > max_fraction_ambiguous_samples;
    % SNP positions
    goodpos = ~ ( nonvariablep | ambiguousp ); 
    
    % Check if SNPs exist
    if sum(goodpos) == 0
        fprintf(1, ['No SNPs for Plasmid ' num2str( match_index) '!\n' ] )
        continue
    end
    
    sum(goodpos) 
    
    % Calls for analysis
    Calls_for_analysis = Calls_this_set( goodpos,: );
    
    % Locations of SNPs
    figure(30); 
    clf(30)
    hold on
    box on
    histogram(p(goodpos),50, 'FaceColor', [.75 .75 .75]); xlabel('position (bp)'); ylabel('# SNPs'); title(['plasmid ' num2str(match_index)])
    text(max(p)/2, 0.975*max(ylim), ['(n=' num2str(sum(goodpos)) ')'], 'HorizontalAlignment', 'center')
    hold off
    print([dir_trees '/' 'snppos_plasmid-' num2str(match_index) '.png'],'-dpng')
    
    % Reference genome (ie plasmid scaffold)
    dir_refgenome = [ dir_hybrid_scaffolds '/' hybrid_scaffold_list{match_index} ];
    [ChrStarts, ~, ~, ~] = genomestats(dir_refgenome);
    refnt_all = extract_outgroup_mutation_positions(dir_refgenome, p2chrpos(p,ChrStarts));
    [~,refnti_all]=ismember(refnt_all,NTs); 
    % Set "ancestor" as reference (not really the ancestor
    anc_nti = refnti_all; 
    fixedmutation = ( (Calls_for_analysis~=repmat(anc_nti(goodpos),1,Nsample)) & Calls_for_analysis>0 );
    [MutQual, ~] = ana_mutation_quality( Calls_for_analysis,Quals(goodpos,:) );
    % Reorder
    reorder_set_by_cluster = []; % initialize
    for c=1:numel(cluster_list_all)
        next_lineage = cluster_list_all(c);
        if exist( [ dir_lineages '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ], 'file' )
            next_tree_order = importdata([ dir_lineages '/' 'sample_order_lineage_' num2str(next_lineage) '.txt' ]);
        else
            next_tree_order = SamplesNamesSimple_this_set_long( cluster_membership_set == next_lineage );
        end
        % Find index of each name in SampleNames_keep
        for i=1:numel(next_tree_order)
            next_name = next_tree_order{i};
            next_index = find( ismember( SamplesNamesSimple_this_set_long, next_name ) );
            if ~isempty( next_index )
                reorder_set_by_cluster(end+1) = next_index;
            end
        end
    end
    
    if ~make_tree
        continue
    end
    
    cd(dir_trees);

    % Option to change set of samples included in tree or to change order
    samplestoplot = 1:numel(SampleNamesSimple_this_set); % currently everything

    % Collect calls for tree making and stores them as characters
    calls_for_tree=zeros(size(Calls_for_analysis));
    calls_for_tree(Calls_for_analysis>0)=NTs(Calls_for_analysis(Calls_for_analysis>0));
    calls_for_tree(Calls_for_analysis==0)='N';
    calls_for_tree=calls_for_tree(:,samplestoplot); % only grabs samples to plot, as defined above

    % Make the tree
    [treefilename, UsedTreeNames] = generate_parsimony_tree_old_aro(calls_for_tree, SamplesNamesSimple_this_set_long(samplestoplot), [ 'Plasmid- ' num2str(match_index) '_' match_name]);

    cd('..')
    
end


%% TREES OF PLASMID COMMON REGION SNPs
% Make trees

% Compare samples with any plasmid type
this_set_name = 'all'; 
this_set = has_plasmid; % all samples that have a plasmid
match_index = 3; % use alignments to plasmid 3
max_fraction_ambiguous_samples = .1; % want common regions to both sets of samples
% note: using plasmid 3 alignments because has a region common the the
% greatest number of plasmid positive samples and also these samples
% represent the greatest diversity across lineages

% Directory for tree
dir_trees = 'trees_plasmids_common';
if ~exist(dir_trees,'dir')
    mkdir(dir_trees)
end

% Position filter
pos_filter_min_copynum = 0.75;
pos_filter_frac_colonies = 0.85; %0.67; 
% want positions with where at least XX% of colonies have a copy number of
% at least 0.XX

% Sample filter
sample_filter_min_copynum = 0.75;
sample_filter_min_frac_pos = 0.75;
% want samples that have at least 75% of positions with a copy number of at
% least 0.75

% Basecalling filters
% Calls filters
min_qual_for_call = 30; % qual = quality, specifically FQ 
min_maf_for_call = .66; % maf = major allele frequency 
min_cov_each_strand_for_call = 3; % cov = coverage 
% Positions filter
max_fraction_ambiguous_samples = 1; % not using this filter here

% Plasmid scaffold to examine
match_name = hybrid_scaffold_list{match_index};
SampleNamesSimple_this_set = SampleNames_keep( this_set );
SamplesNamesSimple_this_set_long = SampleNamesSimpleLong_keep( this_set );
SamplesNamesSimple_this_set_annotated = SampleNamesSimpleLong_keep_annotated( this_set );
Nsample = numel( SampleNamesSimple_this_set );

% Load case step data and get coverage, calls, positions
filename = [ dir_case_step '/' 'case_' hybrid_scaffold_list{match_index} '/' 'candidate_mutation_table.mat' ];
[ Calls_this_set, coverage, p, counts, Quals ] = get_scaff_cov_and_calls_this_set( filename, SampleNamesSimple_this_set, ...
    min_qual_for_call, min_maf_for_call, min_cov_each_strand_for_call );
coverage_copynum = coverage./chromosomal_coverage(this_set)';

% Filter positions
pos_frac_colonies_covered = sum( coverage_copynum > pos_filter_min_copynum,2 )/Nsample;
pos_to_keep = ( pos_frac_colonies_covered > pos_filter_frac_colonies );
%
sum(pos_to_keep) % 46457 (prev 33423)
figure(10); histogram( pos_frac_colonies_covered )

% Filter samples
samp_frac_pos_covered = sum( coverage_copynum(pos_to_keep,:) > sample_filter_min_copynum )/sum(pos_to_keep);
samp_to_keep = ( samp_frac_pos_covered > sample_filter_min_frac_pos );
%
sum(samp_to_keep) % 215 of 216 (prev 215 of 291 )
figure(11); histogram( samp_frac_pos_covered, 0:0.05:1 )

% Find variable positions
Calls_this_set_filtered = Calls_this_set( pos_to_keep, samp_to_keep );
SampleNames_filtered = SamplesNamesSimple_this_set_long( samp_to_keep );
Nsample_filtered = sum( samp_to_keep );
nonvariablep=(sum(Calls_this_set_filtered==1,2)==sum(Calls_this_set_filtered~=0,2)) ...% all samples with calls have the same base
    | (sum(Calls_this_set_filtered==2,2)==sum(Calls_this_set_filtered~=0,2)) ...
    | (sum(Calls_this_set_filtered==3,2)==sum(Calls_this_set_filtered~=0,2)) ...
    | (sum(Calls_this_set_filtered==4,2)==sum(Calls_this_set_filtered~=0,2)) ; % 769 positions are variable
% Ignore positions where too many samples have N's
ambiguousp = sum(Calls_this_set_filtered==0,2)/Nsample_filtered > max_fraction_ambiguous_samples; % 0

% SNP positions
goodpos = ~ ( nonvariablep | ambiguousp ); 
sum(goodpos) % 868 (prev 769)
% Calls for analysis
Calls_for_analysis = Calls_this_set_filtered( goodpos,: );

% Locations of SNPs
figure(30); 
clf(30)
hold on
box on
histogram(p(goodpos),50, 'FaceColor', [.75 .75 .75]); xlabel('position (bp)'); ylabel('# SNPs'); title(['plasmid ' num2str(match_index)])
text(max(p)/2, 0.975*max(ylim), ['(n=' num2str(sum(goodpos)) ')'], 'HorizontalAlignment', 'center')
hold off
print([dir_trees '/' 'snppos_plasmid3-common-' this_set_name '.png'],'-dpng')

%%

% Treemaking
cd(dir_trees);

% Option to change set of samples included in tree or to change order
samplestoplot = 1:Nsample_filtered; % currently everything

% Collect calls for tree making and stores them as characters
calls_for_tree=zeros(size(Calls_for_analysis));
calls_for_tree(Calls_for_analysis>0)=NTs(Calls_for_analysis(Calls_for_analysis>0));
calls_for_tree(Calls_for_analysis==0)='N';
calls_for_tree=calls_for_tree(:,samplestoplot); % only grabs samples to plot, as defined above

% Make the tree
[treefilename, UsedTreeNames] = generate_parsimony_tree_old_aro(calls_for_tree, SampleNames_filtered(samplestoplot), [ 'CommonRegionTree_Plasmid-' num2str(match_index) '_' match_name '_' this_set_name]);

% Distance matrix for tree scale bar

% Generate distance matrix
dist_mat = zeros( Nsample_filtered, Nsample_filtered );
for i=1:Nsample_filtered
    for j=1:Nsample_filtered
        dist_mat(i,j) = sum( Calls_for_analysis(:,i)~=Calls_for_analysis(:,j) & Calls_for_analysis(:,i)>0 & Calls_for_analysis(:,j)>0 );
    end
end

% Write csv file
fid = fopen( 'dist_mat.csv' , 'w' );
% first line
fprintf(fid, ',');
for i=1:Nsample_filtered
    fprintf(fid, [ SampleNames_filtered{i} ',' ] );
end
fprintf(fid,'\n');
% next lines
for i=1:Nsample_filtered
    fprintf(fid, [ SampleNames_filtered{i} ',' ] );
    for j=1:Nsample_filtered
        fprintf(fid, [ num2str(dist_mat(i,j)) ',' ] );
    end
    fprintf(fid,'\n');
end
fclose(fid);

cd('..')


