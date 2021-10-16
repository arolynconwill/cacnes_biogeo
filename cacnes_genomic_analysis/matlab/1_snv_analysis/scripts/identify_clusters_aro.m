%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SET PARAMETERS, LOAD DATA, PROCESS METADATA %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters: IMPORTANT TO ADJUST !!!!
% Much of this will need to vary according to your reference genome,
% coverage, and particular samples
fprintf(1,'Setting parameters...\n')

% Parameters for filtering Calls; base set to N if does not pass filter
max_frac_reads_supp_indels = 0.5; % maximum fraction of reads supporting indels 
min_qual_for_call = 30; % qual = quality, specifically FQ  
min_maf_for_call = .9; % maf = major allele frequency
min_cov_each_strand_for_call = 3; % cov = coverage

% Parameters for basic filtering of samples
% Filter 0: Bracken percent reads assigned to C. acnes
min_frac_cacnes_bracken = 0.9; 
% Filter 1: Coverage
min_median_coverage_to_include_sample = 10;
% Filter 2: Purity
max_low_maf_positions_to_include_sample = 0.01; % purity threshold 
min_maf_for_purity = .65; % used in purity filter
min_cov_for_purity = 4; % used in purity filter

% Parameters for within-cluster contamination detection
min_size_to_examine_cluster = 5; % minimum size to examine cluster
min_pos_to_examine_cluster = 10; % minimum number of positions required for filtering
filtering_mean_maf_per_sample = 0.95; % minimum allowed mean MAF per sample over candidate SNPs with good coverage
filtering_min_cov_to_consider_pos = 8; % minimum coverage for position's MAF to be included in filtering

% Parameters for suspicious subject filtering
min_num_colonies_from_minority_subject = 3;

% Clustering parameters
clustering_dist_cutoff = 35; % dbscan epsilon
clustering_minpts = 3; % dbscan min pts
clustering_merging_threshold = 80; % mean distance between clusters for merging

% Load plasmid positions to mask
load( 'data/plasmid_pos_to_mask.mat' ) % plasmid_pos_on_genome

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs='ATCG';


%% Set up directories and environment
fprintf(1,'Setting up directories...\n')

% Directory for reference genome folder:
dir_ref_genome = [ pwd '/reference_genomes/Pacnes_C1'];
% Directory for my scripts:
dir_scripts_myscripts = [ pwd '/scripts/myscripts_clustering'];
path(dir_scripts_myscripts,path);
dir_scripts_myscripts = [ pwd '/../lab_scripts'];
path(dir_scripts_myscripts,path);

% Clustering directory
dir_clustering = '1_clustering';
if ~exist(dir_clustering,'dir')
    mkdir(dir_clustering)
end


%% Load data from previous steps
fprintf(1,'Loading data...\n')

% From case step (ran on cluster):
load('data/candidate_mutation_table') % includes SampleNames, counts, Quals, indel_counter, and p

% Rename variables
counts_unfiltered = counts; clear counts;
Quals_unfiltered = -1*Quals; clear Quals; % use -Quals because this way higher numbers are more confident
p_all = p; 
SampleNames_unfiltered = SampleNames;
%
%
%
%
% FIX NONSTANDARD SAMPLE NAME FORMAT %
SampleNames_unfiltered = cellfun(@(x) strrep(x,'LcFo','lcFo'), SampleNames, 'UniformOutput', false )
%
%
%


% Remove controls
is_control = cellfun(@(x) contains(x,'Control'), SampleNames_unfiltered );
SampleNames_unfiltered = SampleNames_unfiltered(~is_control);
counts_unfiltered = counts_unfiltered(:,:,~is_control);
Quals_unfiltered = Quals_unfiltered(:,~is_control);
indel_counter_unfiltered = indel_counter(:,:,~is_control);


%% Parse sample names
fprintf(1,'Parsing sample names...\n')

% Keys
type_names={'Extract', 'Strip', 'Scrape'};
type_names_short={'Ext','Stp','Scr'};
zone_list={{'Forehead','Fo','Tz'},{'Chin','Ch','Cn','Ce'},{'RightCheek','Rc'},{'LeftCheek','Lc'},{'Nose','No','no'},{'NearMouth','Mo','nc'},{'Back','Ba'},{'Neck','Nk'},{'Shoulder','Sh'}};
zones_list_short = cellfun(@(x) x{2}, zone_list, 'UniformOutput', false );
zonecodes='FCRLNMBKHU';
typecodes='XYZU';
zonearea_list={{'TopRight','tr','Tr'},{'TopCenter','tc','Tc'},{'TopLeft','tl','Tl'},...
    {'CenterRight','cr','Cr','Mr'},{'CenterCenter','cc','Cc','C'},{'CenterLeft','cl','Cl'},... % R and L can also mean center right and center left but need to be dealt with in special ways because they are not unique tags
    {'BottomRight','br','Br','Lr'},{'BottomCenter','bc','Bc'},{'BottomLeft','bl','Bl','Ll'}};
mult_list = {'Ws','W-','Bs','Mx'}; % characters in names that indicate extract may include >1 pore

% Extracts info from sample name
% Function to parse sample names

% Subject IDs:
subjects_unfiltered = cellfun(@(x) x(1)-64, SampleNames_unfiltered );
% Sebaceous skin region:
zones_unfiltered_str = cellfun(@(x) strsplit(x,'-'), SampleNames_unfiltered, 'UniformOutput', false );
zones_unfiltered_str = cellfun(@(x) x{3}, zones_unfiltered_str, 'UniformOutput', false );
zones_unfiltered = cell2mat( cellfun(@(x) find( cellfun(@(y) contains(x,y), zones_list_short ) ), zones_unfiltered_str, 'UniformOutput', false ) );
% Sample type:
types_unfiltered_str = cellfun(@(x) x(7:9), SampleNames_unfiltered, 'UniformOutput', false );
types_unfiltered = cell2mat( cellfun(@(x) find(ismember(type_names_short,x)), types_unfiltered_str, 'UniformOutput', false ) );
% Sampling time:
times_unfiltered = cellfun(@(x) str2double(x(end-1:end)), SampleNames_unfiltered );
% Sample number:
specimen_number_unfiltered = cellfun(@(x) str2double(x(3:5)), SampleNames_unfiltered);
% Pore multiplicity:
load('data/spec_mult.mat')
multiples_unfiltered = zeros( 1,numel(SampleNames_unfiltered) );
for i=1:numel(multiples_unfiltered)
    temp_spec = specimen_number_unfiltered(i);
    if ~isempty(find(temp_spec==spec_mult_specnums))
        multiples_unfiltered(i) = spec_mult_porenums(spec_mult_specnums==temp_spec);
    elseif types_unfiltered(i) == 3
        multiples_unfiltered(i) = 5; % scrapes
    elseif sum( cellfun(@(x) contains( SampleNames_unfiltered{i}, x), mult_list ) )
        multiples_unfiltered(i) = 5; % name indicates mixes of pores
    else
        fprintf( 1, [ 'Warning! Pore sample multiplicity not available for ' SampleNames_unfiltered{i} '.\n' ] )
        multiples_unfiltered(i) = 5; % assume not a single pore
    end
end


%% Load bracken data
fprintf(1,'Getting bracken info...\n')

% From kraken/bracken:
load('data/data_bracken_cacnes')
% SampleNames_bracken, cacnes_fracs, cacnes_reads
%
%
%
%
% FIX NONSTANDARD SAMPLE NAME FORMAT %
SampleNames_bracken = cellfun(@(x) strrep(x,'LcFo','lcFo'), SampleNames_bracken, 'UniformOutput', false )
%
%
%

cacnes_frac_unfiltered = zeros( numel(SampleNames_unfiltered),1 ); % initialize
cacnes_reads_unfiltered = zeros( numel(SampleNames_unfiltered),1 ); % initialize
for i=1:numel(SampleNames_unfiltered)
    temp_loc = find( ismember( SampleNames_bracken,SampleNames_unfiltered{i} ) );
    if isempty(temp_loc)
        fprintf(['Warning! No bracken data for sample ' SampleNames_unfiltered{i} '.\n'])
        % Note: known that six samples did not have sufficient reads for
        % kraken/bracken pipeline.
    else
        cacnes_frac_unfiltered( i ) = cacnes_fracs( temp_loc );
        cacnes_reads_unfiltered( i ) = cacnes_reads( temp_loc );
    end
end
clear temp_loc



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PRELIMINARY FILTERING STEPS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Keep track of filtering results
fprintf(1,'Initializing filter tracker...\n')

% Initialize a cell array to keep track of which samples are filtered at which stage
sample_filter_types = {}; % keep track of types of filters
sample_filter_failed = {}; % keep track of samples that didn't pass each filter


%% Calculate coverage, call major alleles, filter calls...
fprintf(1,'Calculating coverage and major alleles...\n')

% Calcualte coverage
coverage_unfiltered=squeeze(sum(counts_unfiltered(1:8,:,:)));
coverage_forward_strand_unfiltered=squeeze(sum(counts_unfiltered(1:4,:,:)));
coverage_reverse_strand_unfiltered=squeeze(sum(counts_unfiltered(5:8,:,:)));
[maf_unfiltered, maNT_unfiltered, ~, ~] = div_major_allele_freq(counts_unfiltered);

% Call major alleles
Calls_unfiltered = maNT_unfiltered;
% Set Calls to zero (N) if they are low quality, have low coverage, or have evidence of indels
Calls_unfiltered( Quals_unfiltered < min_qual_for_call | maf_unfiltered < min_maf_for_call | coverage_forward_strand_unfiltered < min_cov_each_strand_for_call | coverage_reverse_strand_unfiltered < min_cov_each_strand_for_call ) = 0; 
Calls_unfiltered( squeeze(indel_counter_unfiltered(1,:,:)) >= max_frac_reads_supp_indels*coverage_unfiltered ) = 0; 
% Mask plasmid regions
Calls_unfiltered( ismember(p_all, plasmid_pos_on_genome),: ) = 0; 


%% Preliminary filtering
fprintf(1,'Preliminary filtering of samples...\n')

% Filter 0: Fraction reads assigned to C. acnes by bracken
% uses: cacnes_frac_unfiltered
fprintf(1,'...bracken filter...\n')
goodsamples_bracken = ( cacnes_frac_unfiltered' >= min_frac_cacnes_bracken ); 

% Filter 1: Coverage of C. acnes reference genome
% uses: min_median_coverage_to_include_sample
fprintf(1,'...coverage filter...\n')
coverage_unfiltered_median = median(coverage_unfiltered(~ismember(p_all,plasmid_pos_on_genome),:)); % coverage over non-plasmid positions only
goodsamples_coveragefilter = coverage_unfiltered_median >= min_median_coverage_to_include_sample;

% Filter 2: Purity filter across all samples
% uses: max_low_maf_positions_to_include_sample, min_maf_for_purity, min_cov_for_purity
fprintf(1,'...purity (global maf) filter...\n')
maf_unfiltered_fraclowmaf = sum(maf_unfiltered<min_maf_for_purity & coverage_unfiltered>= min_cov_for_purity)./sum(coverage_unfiltered>= min_cov_for_purity);
goodsamples_purityfilter = ( maf_unfiltered_fraclowmaf < max_low_maf_positions_to_include_sample); %purity filter
% Looking only at positions with at least the minimum coverage
    % Find what percent of them have a major allele frequency that is too low (<0.65)
    % And make sure that percentage is under the threshold (< 0.01)

% Logical matrix indicating whether or not each sample is good
fprintf(1,'...manual removal of samples...\n')
% Remove all Subject X data (contaminated)
contaminated_samples = SampleNames_unfiltered( cellfun(@(x) x(1)=='X', SampleNames_unfiltered ) );
goodsamples_manual = ~ismember(SampleNames_unfiltered,contaminated_samples);

% Combine all filters together
goodsamples_all = goodsamples_bracken & ... % enough reads are C. acnes
    goodsamples_coveragefilter & ... % enough coverage of ref genome
    goodsamples_purityfilter & ... % purity filter
    goodsamples_manual; % ignore samples manually removed


%% Update filter tracker and make plots summarizing initial filtering

% Update filter tracker (sequential)
sample_filter_types{end+1} = 'Bracken';
sample_filter_failed{end+1} = find( ~goodsamples_bracken );
sample_filter_types{end+1} = 'Coverage';
sample_filter_failed{end+1} = setdiff( find( ~goodsamples_coveragefilter ), find( ~goodsamples_bracken ));
sample_filter_types{end+1} = 'Purity';
sample_filter_failed{end+1} = setdiff( find( ~goodsamples_purityfilter ), find( ~goodsamples_bracken | ~goodsamples_coveragefilter ));
sample_filter_types{end+1} = 'ContamFlag';
sample_filter_failed{end+1} = setdiff( find( ismember(SampleNames_unfiltered,contaminated_samples) ), find( ~goodsamples_bracken | ~goodsamples_coveragefilter | ~goodsamples_purityfilter ));

% Directory for recording initial filtering
dir_initialfilt = [dir_clustering '/' 'InitialFilters'];
if ~exist(dir_initialfilt,'dir')
    mkdir(dir_initialfilt)
end

% Make plots
plot_initial_filtering( dir_initialfilt, ...
    cacnes_frac_unfiltered, min_frac_cacnes_bracken, goodsamples_bracken, ...
    coverage_unfiltered_median, min_median_coverage_to_include_sample, goodsamples_coveragefilter, ...
    maf_unfiltered_fraclowmaf, max_low_maf_positions_to_include_sample, min_maf_for_purity, min_cov_for_purity, goodsamples_purityfilter, ...
    goodsamples_manual, goodsamples_all )


%% Record data only for samples that passed filters so far
fprintf(1,'Downsizing data for only samples that passed preliminary filtering...\n')

% Record old indexing 
goodsamples_oldindexing = find(goodsamples_all);

% Record info for samples to consider (using _all suffix on variable names)
% Sample info:
SampleNames_all=SampleNames_unfiltered(goodsamples_all);
subjects_all=subjects_unfiltered(goodsamples_all); clear subjects_unfiltered;
zones_all=zones_unfiltered(goodsamples_all); clear zones_unfiltered;
types_all=types_unfiltered(goodsamples_all); clear types_unfiltered;
times_all=times_unfiltered(goodsamples_all); clear times_unfiltered; 
specimen_number_all=specimen_number_unfiltered(goodsamples_all); clear specimen_number_unfiltered;
multiples_all=multiples_unfiltered(goodsamples_all); clear multiples_unfiltered;
% Calls and major alleles:
Calls_all=Calls_unfiltered(:,goodsamples_all); clear Calls_unfiltered;
maNT_all=maNT_unfiltered(:,goodsamples_all); clear maNT_unfiltered;
maf_all=maf_unfiltered(:,goodsamples_all); clear maf_unfiltered;
% Coverage:
coverage_all=coverage_unfiltered(:,goodsamples_all); clear coverage_unfiltered;
coverage_forward_strand_all=coverage_forward_strand_unfiltered(:,goodsamples_all); clear coverage_forward_strand_unfiltered;
coverage_reverse_strand_all=coverage_reverse_strand_unfiltered(:,goodsamples_all); clear coverage_reverse_strand_unfiltered;
% Counts and quality:
counts_all=counts_unfiltered(:,:,goodsamples_all); clear counts_unfiltered;
Quals_all=Quals_unfiltered(:,goodsamples_all); clear Quals_unfiltered;
indel_counter_all=indel_counter_unfiltered(:,:,goodsamples_all); clear indel_counter; 
% Bracken stuff (to use for more filtering for samples to use for clade assemblies)
cacnes_frac_all = cacnes_frac_unfiltered( goodsamples_all ); clear cacnes_frac_unfiltered;
cacnes_reads_all = cacnes_reads_unfiltered( goodsamples_all ); clear cacnes_reads_unfiltered;


%% Preliminary SLST-typing of each colony

[ slst_all, SampleNamesLong_all ] = type_slst_preliminary( ...
    dir_ref_genome, p_all, maNT_all, NTs, SampleNames_all );


%% Make a distance matrix using all positions
% Note: This step takes a while! Run on cluster for large sampleset.

% Number of samples included in distance matrix
numsamples=numel(SampleNames_all);

% Only generate the distance matrix if the file doesn't already exist
if ~(exist([ dir_clustering '/isolate_distance_matrix.mat' ],'file'))
    fprintf(1,'Generating distance matrix...\n')
    % Initialize distance matrix variable:
    distance_matrix=zeros(numsamples);
    % Specify calls to use for distance calculation:
    Calls_dist=Calls_all;
    % Calculate distances: 
    for i=1:numsamples
        fprintf(1,['Sample progress: ' num2str(i) '/' num2str(numsamples) '. \n'])
        % Number of positions where calls differ, only over positions that have calls in both samples
        distance_matrix(i,:) = sum( Calls_dist~=repmat(Calls_dist(:,i),1,numsamples) & Calls_dist>0 & repmat(Calls_dist(:,i),1,numsamples)>0 ); % note: not normalized
    end
    % Save matrix:
    save( [ dir_clustering '/isolate_distance_matrix' ],'distance_matrix')
else
    fprintf(1,'Loading distance matrix...\n')
    load( [ dir_clustering '/isolate_distance_matrix' ])
end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% WITHIN-CLUSTER CONTAMINATION FILTERING %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Setup for tracking which samples get filtered

% Keep track of original subject assignments
subjects_all_original = subjects_all;

% Fake subject number to assign to samples we removed during filtering
subject_fake = 100;


%% Within-cluster filtering: catch colonies contaminated from the same subject
fprintf(1,'Within-cluster filtering by subject (step 1)...\n')

% See above for clustering and filtering parameters

% Clustering options
allow_clustering_btwn_subjects = 0; % does not allow clusters to contain samples from more than one subject
allow_cluster_merging = 0; % does not allow close clusters to merge

% Filtering options
figures_boolean = true; % whether or not to save figures 
pause_boolean = false; % whether or not to pause after each cluster
dir_figs = [dir_clustering '/Step1']; % directory suffix for this step

% Initialize subjects tracker (set filtered samples to fake_subject)
subjects_filtering_current = subjects_all_original;

% Do within cluster filtering iteratively until no more samples can be filtered
filter_keep_going = true;
filter_tally = 1;
while filter_keep_going
    
    % Track which round of filtering we're on
    filter_tally_letter = char( filter_tally+64 );
    fprintf(1,['\n' 'Within-subject clustering: Step 1, Round ' filter_tally_letter '\n'] )
    
    % Make preliminary clusters
    clusters_all_current = do_clustering ( ...
        allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
        subjects_filtering_current, subject_fake, distance_matrix, SampleNames_all );
    % Sort clusters by size and find subject
    [ clusters_all_current, clusters_all_current_subject ] = ...
        sort_cluster_order( clusters_all_current, SampleNames_all );
    % Check number of samples clustered
    fprintf(1,['Number of clusters: ' num2str(numel(clusters_all_current)) '\n'] )
    fprintf(1,['Number of samples clustered: ' ...
        num2str(sum(cellfun(@(x) numel(x), clusters_all_current))) '/' ...
        num2str(numel(SampleNames_all)) '\n'] )
    
    % Check which clusters have changed
    if filter_tally>1
        % Flag only clusters that were not there previously for filtering
        clusters_all_for_filtering = ~cellfun(@(x) sum( cellfun(@(y) isequal(x,y), clusters_all_previous) ), clusters_all_current);
    else
        % Flag all clusters for filtering
        clusters_all_for_filtering = ones(numel(clusters_all_current),1);
    end
    fprintf(1,['Number of clusters for filtering: ' ...
        num2str(sum(clusters_all_for_filtering)) '/' num2str(numel(clusters_all_for_filtering)) '. \n'])
    clusters_all_for_filtering_indices = find(clusters_all_for_filtering);
    for q=1:numel(clusters_all_for_filtering_indices)
        fprintf(1,[ num2str(clusters_all_for_filtering_indices(q)) ' ' ])
    end
    fprintf(1,'\n')

    % Do within-cluster filtering
    subjects_filtering_next = do_within_cluster_filtering( ...
        subjects_filtering_current, clusters_all_current_subject, ...
        clusters_all_current, clusters_all_for_filtering, ...
        SampleNames_all, Calls_all, coverage_all, maf_all, Quals_all, ...
        distance_matrix, p_all, ...
        figures_boolean, pause_boolean, subject_fake, ...
        min_size_to_examine_cluster, min_pos_to_examine_cluster, ...
        filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
        [dir_figs filter_tally_letter]);
    
    % Decide if more filtering is necessary
    if ~isequal( subjects_filtering_current, subjects_filtering_next )
        % Keep going if anything changed
        filter_tally = filter_tally+1;
        subjects_filtering_current = subjects_filtering_next;
        clusters_all_previous = clusters_all_current;
    else
        % Finish if there was no change
        subjects_all_filtered_1 = subjects_filtering_next; 
        filter_keep_going = false;
    end
    
end

%%

% Update filter tracker
sample_filter_types{end+1} = 'WithinCluster_1';
sample_filter_failed{end+1} = goodsamples_oldindexing( subjects_all_filtered_1 == subject_fake );



%% Within-cluster filtering: catch colonies contaminated from another subject
fprintf(1,'Within-cluster filtering among subjects (step 2)...\n')

% See above for clustering and filtering parameters

% Clustering options
allow_clustering_btwn_subjects = 1; % allows clusters to contain samples from more than one subject
allow_cluster_merging = 1; % allows cluster merging

% Filtering options
figures_boolean = true; % whether or not to save figures 
pause_boolean = false; % whether or not to pause after each cluster
dir_figs = [dir_clustering '/Step2']; % directory suffix for this step

% Initialize subjects tracker (set filtered samples to fake_subject)
subjects_filtering_current = subjects_all_filtered_1;

% Do within cluster filtering iteratively until no more samples can be filtered
filter_keep_going = true;
filter_tally = 1;
while filter_keep_going
    
    % Track which round of filtering we're on
    filter_tally_letter = char( filter_tally+64 );
    fprintf(1,['Within-subject clustering: Step 2, Round ' filter_tally_letter '\n'] )
    
    % Make preliminary clusters
    clusters_all_current = do_clustering ( ...
        allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
        subjects_filtering_current, subject_fake, distance_matrix, SampleNames_all );
    % Sort clusters by size and find subject
    [ clusters_all_current, clusters_all_current_subject ] = ...
        sort_cluster_order( clusters_all_current, SampleNames_all );
    % Check number of samples clustered
    fprintf(1,['Number of clusters: ' num2str(numel(clusters_all_current)) '\n'] )
    fprintf(1,['Number of samples clustered: ' ...
        num2str(sum(cellfun(@(x) numel(x), clusters_all_current))) '/' ...
        num2str(numel(SampleNames_all)) '\n'] )
    
    % Check which clusters have changed
    if filter_tally>1
        % Flag only clusters that were not there previously for filtering
        clusters_all_for_filtering = ~cellfun(@(x) sum( cellfun(@(y) isequal(x,y), clusters_all_previous) ), clusters_all_current);
    else
        % Flag all clusters for filtering
        clusters_all_for_filtering = ones(numel(clusters_all_current),1);
    end
    fprintf(1,['Number of clusters for filtering: ' ...
        num2str(sum(clusters_all_for_filtering)) '/' num2str(numel(clusters_all_for_filtering)) '. \n'])
    clusters_all_for_filtering_indices = find(clusters_all_for_filtering);
    for q=1:numel(clusters_all_for_filtering_indices)
        fprintf(1,[ num2str(clusters_all_for_filtering_indices(q)) ' ' ])
    end
    fprintf(1,'\n')
    
    % Do within-cluster filtering
    subjects_filtering_next = do_within_cluster_filtering( ...
        subjects_filtering_current, clusters_all_current_subject, ...
        clusters_all_current, clusters_all_for_filtering, ...
        SampleNames_all, Calls_all, coverage_all, maf_all, Quals_all, ...
        distance_matrix, p_all, ...
        figures_boolean, pause_boolean, subject_fake, ...
        min_size_to_examine_cluster, min_pos_to_examine_cluster, ...
        filtering_mean_maf_per_sample, filtering_min_cov_to_consider_pos, ...
        [dir_figs filter_tally_letter]);   
    % Decide if more filtering is necessary
    if ~isequal( subjects_filtering_current, subjects_filtering_next )
        % Keep going if there was a change
        filter_tally = filter_tally+1;
        subjects_filtering_current = subjects_filtering_next;
        clusters_all_previous = clusters_all_current;
    else
        % Finish if there was no change
        subjects_all_filtered_2 = subjects_filtering_next; 
        filter_keep_going = false;
    end
    
end


%%

% Update filter tracker
sample_filter_types{end+1} = 'WithinCluster_2';
sample_filter_failed{end+1} = goodsamples_oldindexing( ~(subjects_all_filtered_1 == subject_fake) & ...
    (subjects_all_filtered_2 == subject_fake) );


%% Re-cluster now that within-cluster filtering is done
fprintf(1,'Re-clustering samples after contamination filtering...\n')

% See above for clustering and filtering parameters

% Clustering options
allow_clustering_btwn_subjects = 1; % allows clusters to contain samples from more than one subject
allow_cluster_merging = 1; % allows cluster merging

% Cluster
clusters_all = do_clustering( ...
    allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
    subjects_all_filtered_2, subject_fake, distance_matrix, SampleNames_all );

% Sort clusters by size and find subject
[ clusters_all, clusters_all_subjects ] = ...
    sort_cluster_order( clusters_all, SampleNames_all );


%% Suspicious subject filtering
fprintf(1,'Checking for suspicious subjects...\n')

% This section checks which subjects are represented in each cluster. If a
% minority subject has fewer than min_num_colonies_from_minority_subject 
% colonies in the cluster, then remove them.

[ subjects_all_filtered_all, removed_list ] = do_suspicious_subject_filtering( ...
    min_num_colonies_from_minority_subject, clusters_all, ...
    subjects_all_filtered_2, subject_fake, ...
    SampleNames_all,  ...
    dir_clustering );

%%

% Update filter tracker
sample_filter_types{end+1} = 'SuspiciousSubject';
sample_filter_failed{end+1} = goodsamples_oldindexing( sort(removed_list) );


%% Make final clusters using dbscan
fprintf(1,'Final clustering...\n')

% See above for clustering and filtering parameters

% Clustering options
allow_clustering_btwn_subjects = 1; % allows clusters to contain samples from more than one subject
allow_cluster_merging = 1; % allows cluster merging

% Final clustering
clusters_all = do_clustering( ...
    allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
    subjects_all_filtered_all, subject_fake, distance_matrix, SampleNames_all );

% Sort clusters by size and find subject
[ clusters_all, clusters_all_subjects ] = ...
    sort_cluster_order( clusters_all, SampleNames_all );

% Check number of samples clustered
fprintf(1,['Number of clusters: ' num2str(numel(clusters_all)) '\n'] )
fprintf(1,['Number of samples clustered: ' ...
    num2str(sum(cellfun(@(x) numel(x), clusters_all))) '/' ...
    num2str(numel(SampleNames_all)) '\n'] )

% Print cluster info
for c=1:numel(clusters_all)
    fprintf(1,['Cluster-' num2str(c) '-' clusters_all_subjects{c} '-' num2str(numel(clusters_all{c})) '\n' ])
end


%% Determine which samples are unclustered
fprintf(1,'Finding unclustered samples...\n')

% Update filter tracker with clustered samples
sample_filter_types{end+1} = 'Clustered';
sample_filter_failed{end+1} = goodsamples_oldindexing( cell2mat(clusters_all) );

% Update filter tracker with unclustered samples (these samples passed
% filtering but did not join into a cluster)
sample_filter_types{end+1} = 'Non-Clustered';
sample_filter_failed{end+1} = goodsamples_oldindexing( ...
    setdiff( find(subjects_all_filtered_all~=subject_fake), cell2mat(clusters_all) ) ...
    );

% Record nonclustered samples
unclustered_all = setdiff( find(subjects_all_filtered_all~=subject_fake), cell2mat(clusters_all) );


%% Record filtering outcomes for all samples in a csv file

% Write file
fid=fopen([ dir_clustering '/Log_SampleFilterReport.txt' ],'w');
fprintf(fid,['SampleName_new, FilterOutcome \n']);
for i=1:numel(sample_filter_types)
    if i~=8 % Anything that wasn't clustered
        temp_indices = sample_filter_failed{i};
        for j=1:numel(temp_indices)
            fprintf(fid,[ ...
                SampleNames_unfiltered{temp_indices(j)} ', ' ...
                sample_filter_types{i} ' \n']);
        end
    else % Anything that was clustered
        for c=1:numel(clusters_all) 
            temp_indices = goodsamples_oldindexing(clusters_all{c});
            temp_clusterstring = num2str(c);
            if numel(temp_clusterstring)==1
                temp_clusterstring = ['0' temp_clusterstring ];
            end
            for j=1:numel(temp_indices)
                fprintf(fid,[ ...
                    SampleNames_unfiltered{temp_indices(j)} ', ' ...
                    'Clustered-C' temp_clusterstring ' \n']);
            end
        end
    end
end
clear temp_indices
clear temp_clusterstring
fclose(fid);


%% Save filter tracker information and make a figure summarizing all filtering steps

% sample_filter_types
%     {'Bracken'          }
%     {'Coverage'         }
%     {'Purity'           }
%     {'ContamFlag'       }
%     {'WithinCluster_1'  }
%     {'WithinCluster_2'  }
%     {'SuspiciousSubject'}
%     {'Clustered'        }
%     {'Non-Clustered'    }

% Save data
save([dir_clustering '/data_clusterfiltering.mat'],'-v7.3',...
    'sample_filter_types','sample_filter_failed', ...
     'SampleNames_unfiltered');

dir_finalfilt = [dir_clustering '/' 'FinalFilters'];
if ~exist(dir_finalfilt,'dir')
    mkdir(dir_finalfilt)
end

% Make figure for Part I of filtering
filt_num_original = numel( SampleNames_unfiltered );
filt_num_pass_bracken = filt_num_original - numel( sample_filter_failed{1} );
filt_num_pass_bracken_coverage = filt_num_original - numel( sample_filter_failed{1} ) - numel( sample_filter_failed{2} );
filt_num_pass_bracken_coverage_purity = filt_num_original - numel( sample_filter_failed{1} ) - numel( sample_filter_failed{2} ) - numel( sample_filter_failed{3} );
% Figure
figure(20)
clf(20)
box on
bar([filt_num_original,filt_num_pass_bracken,filt_num_pass_bracken_coverage,filt_num_pass_bracken_coverage_purity])
xticklabels({'None','+Bracken','+Coverage','+Purity'})
xlabel('Filter')
ylabel('Number of samples')
title('Initial filtering')
set(gca,'FontSize',14)
% Save
print([dir_finalfilt '/' 'Summary_Filter_Bar_Purity+Coverage+Bracken'],'-dpng')

% Make figure for Part II of filtering
filt_num_pass_clusterfilt1 = filt_num_pass_bracken_coverage_purity - numel( sample_filter_failed{5} );
filt_num_pass_clusterfilt2 = filt_num_pass_clusterfilt1 - numel( sample_filter_failed{6} );
filt_num_pass_suspicioussubject = filt_num_pass_clusterfilt2 - numel( sample_filter_failed{7} );
filt_num_pass_clustered = numel( sample_filter_failed{8} );
% Figure
figure(21)
clf(21)
box on
bar([filt_num_original,filt_num_pass_bracken,filt_num_pass_bracken_coverage,filt_num_pass_bracken_coverage_purity,filt_num_pass_clusterfilt1,filt_num_pass_clusterfilt2,filt_num_pass_suspicioussubject,filt_num_pass_clustered])
xticklabels({'None','+Bracken','+Coverage','+Purity','+ClusterFilter1','+ClusterFilter2','+SuspiciousSubject','+Clustered'})
xtickangle(45)
xlabel('Filter')
ylabel('Number of samples')
title('All filtering')
set(gca,'FontSize',14)
% Save
print([dir_finalfilt '/' 'Summary_Filter_Bar_PartI+PartII'],'-dpng')


%% Find outgroups
fprintf(1,'Finding outgroups...\n')

% This section sets the outgroup of a cluster to be the closest other
% cluster (measured by mean pairwise distance) from a different subject.

% Get majority subject for each cluster (#s)
clusters_all_subjects_mode = cellfun(@(x) mode(subjects_all(x)), clusters_all);

% Find outgroups
for i=1:numel(clusters_all) % for each cluster
    % Calculate mean pairwise distance to other clusters
    distance_each_cluster=zeros(1,numel(clusters_all));
    for j=1:numel(clusters_all)
        btwn_cluster_diff=distance_matrix(clusters_all{i},clusters_all{j});
        distance_each_cluster(j)=mean(btwn_cluster_diff(:));
    end
    distance_each_cluster(i)=100000; % to avoid self-finding
    distance_each_cluster(clusters_all_subjects_mode(i)==clusters_all_subjects_mode)=100000; % to avoid finding cluster from the same subject
    % Find closest cluster
    closest_cluster=find(distance_each_cluster==min(distance_each_cluster));
    fprintf(1,['Cluster subject: ' char(clusters_all_subjects_mode(i)+64) '\n'])
    fprintf(1,['Outgroup subject: ' char(clusters_all_subjects_mode(closest_cluster)+64) '\n'])
    % Record closest cluster samples as outgroup (max of 10 samples)
    if numel(clusters_all{closest_cluster}) > 10
        % Change_Scratch: Grab random indices rather than first 10 indices
        outgroup_rand_indices = randperm(numel(clusters_all{closest_cluster}),10);
        outgroup_full = clusters_all{closest_cluster};
        outgroup{i}=outgroup_full(1:10);
    else
        outgroup{i}=clusters_all{closest_cluster};
    end
end


%% Add clade number to SampleNamesLong!

% Figure out which cluster each sample is in
SampleClusters_final = -ones(numel(SampleNamesLong_all),1);
% New tree names that include cluster assignment
SampleNamesLong_all_withClade = {};
% Loop through all samples and check if they are in a cluster
for nextNameIndex=1:numel(SampleNamesLong_all)
    nextName = SampleNamesLong_all{nextNameIndex}; % name of the next sample
    nextNameClusterFind = cell2mat(cellfun(@(x) sum(ismember(x,nextNameIndex)), clusters_all, 'UniformOutput', false)); % boolean membership of each cluster
    if sum(nextNameClusterFind)>0 % if sample is in a cluster
        nextClusterNum = find(nextNameClusterFind); % cluster number
        SampleClusters_final(nextNameIndex) = nextClusterNum; % save cluster number
        % Update name to indicate cluster number
        if nextClusterNum < 10
            SampleNamesLong_all_withClade{end+1} = [ nextName '_C-0' num2str(nextClusterNum) ];
        else
            SampleNamesLong_all_withClade{end+1} = [ nextName '_C-' num2str(nextClusterNum) ];
        end
    elseif ismember(nextNameIndex,unclustered_all)
        % Update name to indicate this sample is not in a cluster
        SampleNamesLong_all_withClade{end+1} = [ nextName '_C-00' ];
        SampleClusters_final(nextNameIndex) = 0;
    else
        % Update name to indicate this sample was filtered
        SampleNamesLong_all_withClade{end+1} = [ nextName '_filtered' ];
    end
end

% Check naming
for i=1:numel(SampleNamesLong_all_withClade); fprintf(1,[SampleNamesLong_all_withClade{i} '\n']); end



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% POST-CLUSTERING DIAGNOSTIC PLOTS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Make diagnostic plots for each cluster
% Includes all samples in cluster plus outgroup samples :)

% Combine clusters with outgroups for making plots
clusters_all_with_outgroups = {};
for c=1:length(clusters_all)
    clusters_all_with_outgroups{end+1} = union( clusters_all{c}, outgroup{c} );
end

% Filtering output options
figures_boolean = true; % whether or not to save figures 
pause_boolean = false; % whether or not to pause after each cluster
dir_figs = [dir_clustering '/StepFinal']; % directory suffix for this step

% Examine any cluster that has >= 3 samples and >= 3 candidate SNPs
min_size_to_examine_cluster_diagnostic = 3; % minimum size to examine cluster 
min_pos_to_examine_cluster_diagnostic = 3; % minimum number of positions allowed for filtering samples

% Examine final clusters to check filtering
do_within_cluster_examination_woutgroup( ...
    subjects_all, clusters_all_subjects, clusters_all, outgroup, ...
    SampleNames_all, Calls_all, coverage_all, maf_all, distance_matrix, p_all, ...
    figures_boolean, pause_boolean, subject_fake, ...
    min_size_to_examine_cluster_diagnostic, min_pos_to_examine_cluster_diagnostic, ...
    dir_figs);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SAVE CLUSTER STEP DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Save data only from samples that are in clades
fprintf(1,['Saving ' num2str(sum(cellfun(@(x) numel(x), clusters_all))) ' samples in ' num2str(numel(clusters_all)) ' clusters...\n'])
fprintf(1,['...and saving ' num2str(sum(SampleClusters_final==0)) ' unclustered samples...\n'])

% Only save samples that are in a cluster
samples_to_save = (SampleClusters_final>-1); % isequal((SampleClusters_final>-1),(subjects_all_filtered_all~=subject_fake))

% Get variables with only samples to save
fprintf(1,'\n...new vars...\n')
SampleNames_final = SampleNames_all(samples_to_save); % grab subset
SampleNamesLong_final = SampleNamesLong_all_withClade(samples_to_save);
clusters_final = {}; % re-index clusters_all
for c=1:max(SampleClusters_final) % loop through clusters
   clusters_final{end+1} = find( SampleClusters_final(samples_to_save) == c ); 
end
outgroup_final = {}; % re-index outgroup
for c=1:max(SampleClusters_final) % loop through clusters
    indicesOld = outgroup{c}; % get old indices
    outgroupNames = SampleNames_all(indicesOld); % get sample names in outgroup
    [~,indicesNew] = ismember( outgroupNames, SampleNames_final ); % get new indices
    outgroup_final{end+1} = indicesNew;
end
unclustered_final = find( SampleClusters_final(SampleClusters_final>-1) == 0 );
clusters_final_subjects = clusters_all_subjects; % okay as is
subjects_final = subjects_all(samples_to_save); % grab subset % isequal(subjects_all(samples_to_save),subjects_all_filtered_all(samples_to_save))
zones_final = zones_all(samples_to_save); % grab subset
types_final = types_all(samples_to_save); % grab subset
times_final = times_all(samples_to_save); % grab subset
specimen_number_final = specimen_number_all(samples_to_save); % grab subset
multiples_final = multiples_all(samples_to_save); % grab subset
slst_final = slst_all(samples_to_save); % grab subset; 
Calls_final = Calls_all(:,samples_to_save); % grab subset
counts_final = counts_all(:,:,samples_to_save); % grab subset
Quals_final = Quals_all(:,samples_to_save); % grab subset
coverage_final = coverage_all(:,samples_to_save); % grab subset
maNT_final = maNT_all(:,samples_to_save); % grab subset
maf_final = maf_all(:,samples_to_save); % grab subset
%p_all; % okay as is
indel_counter_final = indel_counter_all(:,:,samples_to_save);
cacnes_frac_final = cacnes_frac_all( samples_to_save );
cacnes_reads_final = cacnes_reads_all( samples_to_save );


%% Make simple sample names

SampleNamesSimple_final = cell( size( SampleNames_final ) );

% Reindex subjects
subjects_list = unique(subjects_final);
subjects_list_num_cols = arrayfun(@(x) sum(x==subjects_final), subjects_list );
[~, sort_subjects] = sort( subjects_list_num_cols, 'descend' );
subjects_list_sorted = subjects_list( sort_subjects );

% Loop through specimen numbers
list_specimens = unique( specimen_number_final );
for s=1:numel(list_specimens)
    
    % Specimen info
    next_spec = list_specimens(s);
    next_spec_indices_bool = ( specimen_number_final == next_spec );
    next_spec_indices = find( next_spec_indices_bool );
    next_spec_num_cols = numel( next_spec_indices );
    next_spec_type = mode( types_final( next_spec_indices ) );
    next_spec_subj = mode( subjects_final( next_spec_indices ) );
    next_spec_subj_char = char( next_spec_subj+64 );
    next_spec_subj_num = find( next_spec_subj_char == char(64+subjects_list_sorted) );

    % Update sample name based on sample type and pore multiplicity
    if next_spec_type == 1 || next_spec_type == 2 % if pore specimen (extract or strip specimen)
        % Get pore multiplicity
        nextspec_mult = spec_mult_porenums( spec_mult_specnums==next_spec );
        if nextspec_mult == 1 % if single pores
            for j=1:next_spec_num_cols
                SampleNamesSimple_final{ next_spec_indices(j) } = [ 'subj-' num2str(next_spec_subj_num) '_pore-' num2str(next_spec) '_col-' num2str(j) ];
            end
        else % not single pore
            for j=1:next_spec_num_cols
                SampleNamesSimple_final{ next_spec_indices(j) } = [ 'subj-' num2str(next_spec_subj_num) '_multipore-' num2str(next_spec) '_col-' num2str(j) ];
            end
        end
    else % not a pore (scrape specimen)
        for j=1:next_spec_num_cols
            SampleNamesSimple_final{ next_spec_indices(j) } = [ 'subj-' num2str(next_spec_subj_num) '_scrape-' num2str(next_spec) '_col-' num2str(j) ];
        end
    end

end

% Make directory for storing csvs
mkdir('keys')

% Subjects key
fid = fopen( [ dir_clustering '/' 'key_subjects.csv' ], 'w');
fprintf(fid,'subject_number,subject_letter,\n');
for i=1:numel(subjects_list_sorted)
    fprintf(fid, [num2str(i) ',' char(64+subjects_list_sorted(i)) ',\n' ]);
end
fclose(fid);

% Sample names key
fid = fopen( [ dir_clustering '/' 'key_samples.csv' ], 'w');
fprintf(fid,'sample_name,sample_name_simple,\n');
for i=1:numel(SampleNames_final)
    fprintf(fid, [ SampleNames_final{i} ',' SampleNamesSimple_final{i} ',\n' ]);
end
fclose(fid);


%% Make cluster names


%% Save data
% Takes forever...

fprintf(1,'\n...saving...\n')
save('data/cluster_step_variables.mat','-v7.3',...
    'clusters_final','outgroup_final','unclustered_final','clusters_final_subjects',... 
    'SampleNames_final','SampleNamesLong_final','SampleNamesSimple_final', ...
    'subjects_final','zones_final','types_final','times_final',...
    'specimen_number_final','multiples_final','slst_final',...
    'Calls_final','counts_final','Quals_final','coverage_final',...
    'maNT_final','maf_final','p_all','indel_counter_final',... 
    'cacnes_frac_final','cacnes_reads_final',...
    'subjects_list','subjects_list_sorted'...
    ); 

%load('data/cluster_step_variables.mat')



%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% MAKE SAMPLES.CSV WITH NEW NAMES %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Writes generic samples.csv for all samples that passed filtering

% Set up directory
dir_samplescsv = 'AssemblyCSVs';
if ~exist( [dir_clustering '/' dir_samplescsv], 'dir' )
    mkdir( [dir_clustering '/' dir_samplescsv] )
end

% Get original samples.csv
csv_original = table2cell(readtable(['data/samples.csv']));
original_paths = csv_original(:,1);
original_samples = csv_original(:,2);
original_refgenomes = csv_original(:,3);
original_providernames = csv_original(:,4);
original_subjects = csv_original(:,5);

% Make table to keep track of info for new samples.csv
new_table = {};
new_table{1,1} = 'Path';
new_table{1,2} = 'Sample';
new_table{1,3} = 'ReferenceGenome';
new_table{1,4} = 'ProviderName';
new_table{1,5} = 'Subject';
new_table_row = 2;

% Loop through samples
for s=1:numel(SampleNames_final)
    
    % Next sample
    next_name = SampleNames_final{s};

    % Find sample in original csv
    next_index_csv = find( cellfun(@(x) isequal(x,next_name),original_samples) );
    
    % Add to table
    new_table{new_table_row,1} = original_paths{next_index_csv}; % 'Path';
    new_table{new_table_row,2} = next_name; % 'Sample';
    new_table{new_table_row,3} = original_refgenomes{next_index_csv}; % 'ReferenceGenome';
    new_table{new_table_row,4} = original_providernames{next_index_csv}; % 'ProviderName';
    new_table{new_table_row,5} = next_name(1); % 'Subject'; actually clade
    new_table_row = new_table_row+1;

end

% Write file
fid = fopen([ dir_clustering '/' dir_samplescsv '/' 'samples_allpassedfiltering.csv'],'w');
for i=1:length(new_table)
    fprintf(fid,[new_table{i,1} ',' ...
        new_table{i,2} ',' ...
        new_table{i,3} ',' ...
        new_table{i,4} ',' ...
        new_table{i,5} ',' '\n' ]);
end
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PREPARE FILES FOR CLUSTER GENOME ASSEMBLIES %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Writes files for genome assembly (99% bracken C. acnes filter)

% Note: If no samples pass the purity filter, grab any colonies that are at
% least 0.95*min_frac_cacnes_for_assembly purity.

% Filters
min_frac_cacnes_for_assembly = 0.99;
min_frac_cacnes_for_assembly_alternative = 0.95; % if no colonines with 99% purity
min_clade_size_for_assembly = 1;

% Set up directory
dir_cladeassembly = 'AssemblyCSVs';
if ~exist( [dir_clustering '/' dir_cladeassembly], 'dir' )
    mkdir( [dir_clustering '/' dir_cladeassembly] )
end
% Directory for clade text files
dir_txtfiles = '2-clades';
if ~exist( [dir_clustering '/' dir_cladeassembly '/' dir_txtfiles], 'dir' )
    mkdir( [dir_clustering '/' dir_cladeassembly '/' dir_txtfiles] );
end

% Get original samples.csv
csv_original = table2cell(readtable(['data/samples.csv']));
original_paths = csv_original(:,1);
original_samples = csv_original(:,2);
original_refgenomes = csv_original(:,3);
original_providernames = csv_original(:,4);
original_subjects = csv_original(:,5);

% Make table to keep track of info for new samples.csv
new_table = {};
new_table{1,1} = 'Path';
new_table{1,2} = 'Sample';
new_table{1,3} = 'ReferenceGenome';
new_table{1,4} = 'ProviderName';
new_table{1,5} = 'Subject';
new_table_row = 2;

% Loop through clades
for c=1:numel(clusters_final)

    % Get clade
    next_clade = clusters_final{c};
    next_clade_names = SampleNames_final(next_clade);
    next_clade_bracken_fracs = cacnes_frac_final(next_clade);
    next_clade_bracken_reads = cacnes_reads_final(next_clade);

    % Apply bracken filter
    next_clade_keep = ( next_clade_bracken_fracs >= min_frac_cacnes_for_assembly );
    if sum(next_clade_keep)==0
        fprintf(1,[ 'Warning! Needed to lower purity filter for Cluster ' num2str(c) '.' '\n' ])
        next_clade_keep = ( next_clade_bracken_fracs >= min_frac_cacnes_for_assembly_alternative );
    end

    % Save clade info
    if sum(next_clade_keep) >= min_clade_size_for_assembly
        
        % Samples that passed strict bracken filter
        next_clade_keep_indices = find( next_clade_keep );
        % Write text file that indicates which samples to use in assembly
        fid = fopen([ dir_clustering '/' dir_cladeassembly '/' dir_txtfiles '/' 'clade' num2str(c) '_samples.txt'],'w');
        for i=1:numel(next_clade_keep_indices) 
            fprintf(fid,[ next_clade_names{next_clade_keep_indices(i)} '\n' ]); % fprintf(fid,[ next_cluster_old_names{i} '\n' ]);
        end
        fclose(fid);

        % All samples in clade even if they didn't pass strict bracken filter
        next_clade_all_indices = 1:1:numel(next_clade);
        % Add to table that will become samples.csv 
        for s=1:numel(next_clade_all_indices)
            % Get info
            next_index = next_clade_all_indices(s);
            next_index_csv = find( cellfun(@(x) isequal(x,next_clade_names{next_index}),original_samples) );
            % Add to table
            new_table{new_table_row,1} = original_paths{next_index_csv}; % 'Path';
            new_table{new_table_row,2} = next_clade_names{next_index}; % 'Sample';
            new_table{new_table_row,3} = original_refgenomes{next_index_csv}; % 'ReferenceGenome';
            new_table{new_table_row,4} = original_providernames{next_index_csv}; % 'ProviderName';
            new_table{new_table_row,5} = num2str(c); % 'Subject'; actually clade
            new_table_row = new_table_row+1;
        end
    end

end

% Write file
fid = fopen([ dir_clustering '/' dir_cladeassembly '/' 'samples.csv'],'w');
for i=1:length(new_table)
    fprintf(fid,[new_table{i,1} ',' ...
        new_table{i,2} ',' ...
        new_table{i,3} ',' ...
        new_table{i,4} ',' ...
        new_table{i,5} ',' '\n' ]);
end
fclose(fid);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CLUSTER DISTANCE HISTOGRAMS COLLECTOR'S CURVE FIGURES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up directory
dir_collcurve = [ dir_clustering '/' 'SummaryHistCollec' ];
if ~exist( dir_collcurve, 'dir' )
    mkdir( dir_collcurve )
end

%% Make histograms summarizing distance between clusters 
% 
% % Min distance histograms
% make_distance_histograms_by_subj( clusters_all, unclustered_all, distance_matrix, clusters_all_subjects, subjects_all, subject_fake, dir_clustering )


%% Collector's curve to examine % w/o clade vs number of samples
fprintf(1,'Making collector''s curves...\n')

% Make distance matrix for these samples only

% Only generate the distance matrix if the file doesn't already exist
if ~(exist( [ dir_collcurve '/mini_distance_matrix.mat' ],'file'))
    fprintf(1,'Generating distance matrix...\n')
    % Initialize distance matrix variable:
    numsamples = numel(SampleNames_final);
    distance_matrix_mini=zeros(numsamples);
    % Specify calls to use for distance calculation:
    Calls_dist=Calls_final;
    % Calculate distances: 
    for i=1:numsamples
        fprintf(1,['Sample progress: ' num2str(i) '/' num2str(numsamples) '. \n'])
        % Number of positions where calls differ, only over positions that have calls in both samples
        distance_matrix_mini(i,:) = sum( Calls_dist~=repmat(Calls_dist(:,i),1,numsamples) & Calls_dist>0 & repmat(Calls_dist(:,i),1,numsamples)>0 ); % note: not normalized
    end
    % Save matrix:
    save( [ dir_collcurve '/mini_distance_matrix.mat' ],'distance_matrix_mini')
else
    fprintf(1,'Loading distance matrix...\n')
    load( [ dir_collcurve '/mini_distance_matrix.mat' ] )
end


%% Make histogram of pairwise distance between samples

make_distance_histograms_pairwise( distance_matrix_mini, subjects_final, clusters_final, unclustered_final, SampleNamesLong_final, dir_collcurve )


%% Make histograms summarizing distance between clusters 

% Min distance histograms
subject_fake = 100;
make_distance_histograms_by_subj_end( clusters_final, unclustered_final, distance_matrix_mini, clusters_final_subjects, subjects_final, subject_fake, dir_clustering, dir_collcurve )


%% Clustering collector's curves

% Clustering parameters (same as final clusters)
allow_clustering_btwn_subjects = 1; % allows clusters to contain samples from more than one subject
allow_cluster_merging = 1; % allows cluster merging

%%

% Collector's curves
make_clustering_collectors_curves_downsample_colonies( ...
    allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
    subjects_final, distance_matrix_mini, clusters_final, unclustered_final, dir_collcurve )

%%

% Collector's curves by specimen
make_clustering_collectors_curves_downsample_specimens( ...
    allow_clustering_btwn_subjects, clustering_dist_cutoff, clustering_minpts, allow_cluster_merging, clustering_merging_threshold, ...
    subjects_final, specimen_number_final, distance_matrix_mini, clusters_final, unclustered_final, dir_collcurve )

