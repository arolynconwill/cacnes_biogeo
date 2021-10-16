%%%%%%%%%%%%%
%% SUMMARY %%
%%%%%%%%%%%%%

%% Summary

% This script identifies SNPs within each cluster and saves information
% about them.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 0: SETUP (PARAMETERS, DIRECTORIES, ETC) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set parameters
fprintf(1,'Setting parameters...\n')

% Used to filter maNT for identifying SNPs
% For filtering individual calls
max_frac_reads_supp_indels = 0.25; % maximum fraction of reads supporting indels 
max_frac_reads_supp_dels = 0.25; % maximum fraction of reads supporting dels 
min_qual_for_call = 30; % qual = quality, specifically FQ % 60 in identify_clusters
min_maf_for_call = .75; % maf = major allele frequency % 0.9 in identify_clusters % originally 0.66
min_cov_each_strand_for_call = 3; % cov = coverage % 3 in identify_clusters
% For filtering positions
min_median_coverage_position = 12; %15; % across samples per position % originally 10 % in future would be good to do this relative to other positions...
max_fraction_ambiguous_samples = .34; % 0.33; % across samples per position % originally .3
% For filtering recombination regions
recombination_block_size=500;
correlation_cutoff=0.75;
% For filtering plasmid positions
load('data/plasmid_pos_to_mask.mat','plasmid_pos_on_genome'); % positions on reference genome with homology to plasmid gene content

% Used to filter maNT for analysis
min_maf_for_analysis = 0.67; % .65; 
min_cov_each_strand_for_analysis = 1; % added
min_coverage_for_analysis = 3; % for each call

% Promoter mutations: how far upstream of the nearest gene to annotate something a promoter mutation
promotersize = 150; % bp

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';


%% Set up directories and environment
fprintf(1,'Setting up directories...\n')

% Directory for my scripts:
dir_scripts_myscripts = [ pwd '/scripts/myscripts_denovomuts'];
path(dir_scripts_myscripts,path);
dir_scripts_myscripts = [ '../lab_scripts'];
path(dir_scripts_myscripts,path);
dir_scripts_aro_sh = '/Users/arolyn/Dropbox\ \(MIT\)/Postdoc/Pacnes_Biogeo/ANALYSIS_GITHUB/cacnes_genomic_analysis/matlab/1_snv_analysis/scripts/myscripts_denovomuts';
% Reference genome folder on Dropbox:
dir_ref_genome = [ pwd '/reference_genomes/Pacnes_C1'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1: LOAD DATA AND METADATA %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load data from clustering step
fprintf(1,'Loading data (from clustering step)...\n')

% Load cluster data (rename when necessary)
clusters_all_subjects = load('data/cluster_step_variables','clusters_final_subjects'); clusters_all_subjects = clusters_all_subjects.clusters_final_subjects;
SampleNames_all = load('data/cluster_step_variables','SampleNames_final'); SampleNames_all = SampleNames_all.SampleNames_final;
SampleNamesSimple_all = load('data/cluster_step_variables','SampleNamesSimple_final'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_final;
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
clusters_all = load('data/cluster_step_variables','clusters_final'); clusters_all = clusters_all.clusters_final;
unclustered_all = load('data/cluster_step_variables','unclustered_final'); unclustered_all = unclustered_all.unclustered_final;
outgroup = load('data/cluster_step_variables','outgroup_final'); outgroup = outgroup.outgroup_final; outgroup = cellfun(@(x) transpose(x), outgroup, 'UniformOutput', false); % need to transpose to match direction of clusters_all
load('data/cluster_step_variables','subjects_list');
load('data/cluster_step_variables','subjects_list_sorted');

% Load data about candidate SNP positions
p_all = load('data/cluster_step_variables','p_all'); p_all = p_all.p_all;
counts_all = load('data/cluster_step_variables','counts_final'); counts_all = counts_all.counts_final;
Quals_all = load('data/cluster_step_variables','Quals_final'); Quals_all = Quals_all.Quals_final;
coverage_all = load('data/cluster_step_variables','coverage_final'); coverage_all = coverage_all.coverage_final;
maNT_all = load('data/cluster_step_variables','maNT_final'); maNT_all = maNT_all.maNT_final;
maf_all = load('data/cluster_step_variables','maf_final'); maf_all = maf_all.maf_final;
indel_counter_all = load('data/cluster_step_variables','indel_counter_final'); indel_counter_all = indel_counter_all.indel_counter_final;

% Load metadata
subjects_all = load('data/cluster_step_variables','subjects_final'); subjects_all = subjects_all.subjects_final;
zones_all = load('data/cluster_step_variables','zones_final'); zones_all = zones_all.zones_final;
types_all = load('data/cluster_step_variables','types_final'); types_all = types_all.types_final;
times_all = load('data/cluster_step_variables','times_final'); times_all = times_all.times_final;
specimen_number_all = load('data/cluster_step_variables','specimen_number_final'); specimen_number_all = specimen_number_all.specimen_number_final;
slst_all = load('data/cluster_step_variables','slst_final'); slst_all = slst_all.slst_final; slst_all_original = slst_all; % save because will be updated with assembley info later
cacnes_frac_all = load('data/cluster_step_variables','cacnes_frac_final'); cacnes_frac_all = cacnes_frac_all.cacnes_frac_final;


%% Load additional metadata
fprintf(1,'Loading other metadata...\n')

% Load strip specimen coordinates
load('data/spec_coor.mat') % spec_coor_specnums, spec_coor_stripnum, spec_coor_x, spec_coor_y

% Load mutliplicity of pore specimens (strips and extracts)
load('data/spec_mult.mat') % spec_mult_specnums, spec_mult_porenums

% Load SLST from assemblies
load('data/assemblies_straintypes.mat') % slst_all, clusters_all_slst
% adjusts any missing SLSTs that were not identified from reference genome alignment

% Update sample names with metadata
SampleNamesLong_all = update_names_with_metadata( SampleNames_all, slst_all, clusters_all, unclustered_all );



%% Read in genome information
fprintf(1,'Reading in reference genome...\n')

[ChrStarts, GenomeLength, ~, ScafNames] = genomestats(dir_ref_genome);
refnt_all = extract_outgroup_mutation_positions(dir_ref_genome, p2chrpos(p_all,ChrStarts));
[~,refnti_all]=ismember(refnt_all,NTs); 


%% Save very basic coverage data

mkdir('data_other')

% Save average coverage over chromosome
coverage_chromosome_mean = mean(coverage_all(~ismember(p_all,plasmid_pos_on_genome),:));
coverage_chromosome_median = median(coverage_all(~ismember(p_all,plasmid_pos_on_genome),:));
save('data_other/chromosomal_coverage','SampleNames_all','SampleNamesLong_all','SampleNamesSimple_all','coverage_chromosome_mean','coverage_chromosome_median');
% mean(coverage_chromosome_mean) % 76
% median(coverage_chromosome_mean) % 66
% mean(coverage_chromosome_median) % 76
% median(coverage_chromosome_median) % 66


% "region B": 1398569 to 1410901
coverage_regB = coverage_all( p_all >=1398569 & p_all <=1410901, : );
coverage_regB_median = median( coverage_regB );
coverage_regB_median_norm = coverage_regB_median./coverage_chromosome_median;
save('data_other/regB_coverage','SampleNames_all','SampleNamesLong_all','coverage_regB_median','coverage_regB_median_norm');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: IDENTIFY DE NOVO SNPS IN EACH CLUSTER %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Make clusters directory

% SNP-calling info from each cluster will be saved in clusters/name_of_cluster
if ~exist('2_snvs','dir') % Clusters subfolder inside main directory
    mkdir('2_snvs')
end
% Copy dnapars in Clusters folder if not already there
if ~exist('2_snvs/dnapars','file')
    copyfile( 'tools/dnapars', '2_snvs/dnapars' )
    copyfile( 'tools/dnapars.app', '2_snvs/dnapars.app' )
end
% Move to clusters subdirectory
cd('2_snvs') 


%% Save cluster names
fprintf(1,'Saving cluster names...\n')

% Old cluster names
cluster_names = {};
for k=1:numel(clusters_all) % Loop through a bunch of clusters
    this_cluster_name = [ 'Cluster-' num2str(k) '-' clusters_all_subjects{k} '-' num2str(numel(clusters_all{k}))];
    cluster_names{end+1} = this_cluster_name;
end
% New cluster names and cluster subjects by number
cluster_names_new = cell( size( cluster_names ) ); % initialize 
cluster_names_intermediate = cell( size( cluster_names ) ); % initialize 
clusters_all_subjects_nums = cell(size( clusters_all_subjects ) ); % initialize 
subjects_list_letters = unique( clusters_all_subjects );
for s=1:numel(subjects_list_letters)
    this_subject_clusters = find( cell2mat(clusters_all_subjects) == subjects_list_letters{s} );
    for i=1:numel(this_subject_clusters)
        subject_index = num2str(find( subjects_list_letters{s} == char(64+subjects_list_sorted) ));
        lineage_index = char( 96+i );
        cluster_names_new{ this_subject_clusters(i) } = ['Lineage_' subject_index lineage_index];
        clusters_all_subjects_nums{this_subject_clusters(i)} = subject_index;
        % inermediate version
        cluster_names_intermediate{ this_subject_clusters(i) } = ['Lineage_' subjects_list_letters{s} '-' num2str(i)];
    end
end

% Clusters key
fid = fopen( '../keys/key_clusters.csv', 'w');
fprintf(fid,'clusters_og,clusters_v1,clusters_v2,cluster_slst,\n');
for i=1:numel(cluster_names)
    fprintf(fid, [ cluster_names{i} ',' cluster_names_intermediate{i} ',' cluster_names_new{i} ',' clusters_all_slst{i} ',\n' ]);
end
fclose(fid);


% Save
save('cluster_names','cluster_names','cluster_names_new','cluster_names_intermediate','clusters_all','clusters_all_subjects','clusters_all_subjects_nums','clusters_all_slst'); % save cluster names
save('sample_names','SampleNames_all','SampleNamesLong_all','SampleNamesSimple_all','slst_all'); % save sample names


%% Display options for SNP-calling

% Clickable table??
show_clickable_tables = 0; % set to 1 to make interactive tables for investigating raw data at each variant position
% Plotting??
save_plots = 1; % save plots generated (can be slow for large clusters)
% Save data??
save_data = 1; % data file for each clade which can be reloaded later (can be slow for large clusters)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ...ANALYZING ONE CLADE AT A TIME... %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Detecting mutations in each cluster
fprintf(1,'Detecting mutations within each cluster...\n')

% Keep track of number of SNPs detected in each cluster
num_snps_per_lineage = zeros( numel(clusters_all),1 );

%%

% Evaluate cluster-by-cluster
for k=1:numel(clusters_all)

    % Variables to re-set for each cluster
    fig_num = 1; % Keep track of the number of figures

    %% Basic cluster information
    this_cluster_name = [ 'Cluster-' num2str(k) '-' clusters_all_subjects{k} '-' num2str(numel(clusters_all{k}))];
    fprintf(1,['\nExamining ' this_cluster_name '...\n'])

    %% Make a directory for this cluster
    if ~exist(this_cluster_name,'dir')
        mkdir(this_cluster_name)
    end
    % Add dnapars
    if ~exist( [ this_cluster_name '/dnapars' ],'file')
        copyfile( 'dnapars', [this_cluster_name '/dnapars'] )
        copyfile( 'dnapars.app', [this_cluster_name '/dnapars.app'] )
    end
    cd(this_cluster_name) % go into subdirectory for cluster

    % Directory for diagnostic plots
    dir_diagnostic = 'Diagnostic';
    if save_plots
        if ~exist(dir_diagnostic,'dir') % ARO: Groups subfolder
            mkdir(dir_diagnostic)
        end
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 0: GET INFORMATION FOR THIS CLUSTER %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Get info on samples in this cluster

    % Identify samples in this cluster
    clustersamples = [clusters_all{k}; outgroup{k}]; % in cluster and in outrgoup
    fprintf(1,['Number of samples (incl. outgroup): ' num2str(sum(clustersamples>0)) '\n']);
    SampleNames=SampleNames_all(clustersamples);
    SampleNamesLong=SampleNamesLong_all(clustersamples);
    SampleNamesSimple=SampleNamesSimple_all(clustersamples);
    outgroup_isolates=ismember(SampleNames,SampleNames_all(outgroup{k}));
    Nsample=numel(SampleNames);

    % Add outgroup tag to a new version of SampleNames
    SampleNames_tagoutgroup = SampleNames;
    SampleNamesLong_tagoutgroup = SampleNamesLong;
    for s=length(clusters_all{k})+1:length(clustersamples)
        SampleNames_tagoutgroup{s} = [ 'OUTGROUP_', SampleNames{s} ];
        SampleNamesLong_tagoutgroup{s} = [ 'OUTGROUP_', SampleNamesLong{s} ];
    end

    % Get sample metadata for this cluster
    subjects = subjects_all(clustersamples);
    zones = zones_all(clustersamples);
    types = types_all(clustersamples);
    times = times_all(clustersamples);
    specimen_numbers = specimen_number_all(clustersamples);
    specimen_time_months = times_all(clustersamples);
    slst = slst_all(clustersamples);

    % Get quals and counts for filtering
    Quals = Quals_all(:,clustersamples); % quals at each position
    indel_counter_both = squeeze(indel_counter_all(1,:,clustersamples)); % num reads supporting indels
    indel_counter_delonly = squeeze(indel_counter_all(2,:,clustersamples)); % num reads supporting deletions only
    maf = maf_all(:,clustersamples); % maf at each position
    counts = counts_all(:,:,clustersamples); % counts at each position
    cov_fwd_strand = squeeze(sum(counts(1:4,:,:)));
    cov_rev_strand = squeeze(sum(counts(5:8,:,:)));

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 1: IDENTIFY SNPs IN THIS CLUSTER %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% Determine calls for identifying SNPs
    
    % Initialize from maNT
    Calls = maNT_all(:,clustersamples); % all positions in p, but only for this cluster (plus outgroup)
    % Calls filtering (part 1/2): for high quality calls
    Calls( Quals < min_qual_for_call ...
        | indel_counter_both >= max_frac_reads_supp_indels*(cov_fwd_strand+cov_rev_strand) ...
        | indel_counter_delonly >= max_frac_reads_supp_dels*(cov_fwd_strand+cov_rev_strand) ...
        | maf < min_maf_for_call ...
        | cov_fwd_strand < min_cov_each_strand_for_call ...
        | cov_rev_strand < min_cov_each_strand_for_call ...
        ) = 0; % set to N

    % Save Calls across p_all for "modifying" the reference genome to be cluster-specific
    Calls4annotations = Calls;
    Calls4annotations( ...
        sum(Calls<1,2) >= (Nsample*max_fraction_ambiguous_samples) ... % too many N's
        | median(cov_fwd_strand+cov_rev_strand,2) < min_median_coverage_position ... % not enough coverage
        ) = 0; % remove positions non reference genome that aren't covered well or that are ambiguous
    
    % Find nonvariable positions
    nonvariablep = ( ...
        sum(Calls==1,2)==Nsample ...
        | sum(Calls==2,2)==Nsample ...
        | sum(Calls==3,2)==Nsample ...
        | sum(Calls==4,2)==Nsample ...
        | sum(Calls==0,2)==Nsample );
    fprintf(1,['Number of nonvariable positions: ' num2str(sum(nonvariablep)) '\n']);

    % Downsize data for nonvariable positions only
    p = p_all(~nonvariablep); % list of nonvariable positions relative to reference genome
    Calls = Calls(~nonvariablep,:); % nucleotide identity at each position
    counts = counts(:,~nonvariablep,:); % counts at each position
    Quals = Quals(~nonvariablep,:); % quals at each position
    indel_counter_both = indel_counter_both(~nonvariablep,:);
    indel_counter_delonly = indel_counter_delonly(~nonvariablep,:);
    maf = maf(~nonvariablep,:); % maf at each position
    cov_fwd_strand = cov_fwd_strand(~nonvariablep,:); % fwd cov at each position
    cov_rev_strand = cov_rev_strand(~nonvariablep,:); % rev cov at each position
    coverage = coverage_all(~nonvariablep,clustersamples); % coverage at each position
    maNT = maNT_all(~nonvariablep,clustersamples); % maNT at each position

    % Calls filtering (part 2/2): only keep high quality positions 
    Nsample_ingroup = sum(~outgroup_isolates); 
    Calls_ingroup = Calls(:,~outgroup_isolates);
    coverage_ingroup = coverage(:,~outgroup_isolates);
    low_quality_positions = ( ...
        sum( Calls_ingroup<1, 2 ) >= Nsample_ingroup*max_fraction_ambiguous_samples ... % too many N's
        | median( coverage_ingroup,2 ) < min_median_coverage_position ... % not enough coverage
        );
    Calls( low_quality_positions,: ) = 0; 
    
    
    %% Infer ancestor of cluster
    fprintf(1,'Determining outgroup nucleotides...\n');

    % Define ancestor based on outgroup
    outgroup_calls = Calls(:,outgroup_isolates);
    outgroup_calls( outgroup_calls==0 ) = nan;  % since mode takes the lowest value when there is a tie, setting 0's to nan allows a call when the reads are 50% N and 50% some base
    anc_nti = mode(outgroup_calls,2); % set ancestor nucleotides to mode of outgroup sample calls at that position

    % Replace with reference if no value exists from outgroup
    refnt = refnt_all(~nonvariablep);
    refnti = refnti_all(~nonvariablep); 
    anc_nti(isnan(anc_nti)) = refnti(isnan(anc_nti)); 

    % Use mode of ingroup calls for locations without outgroup or reference (would only happen if ref genome has an N at that position)
    locs_without_ancestor = find(sum(Calls(:,~outgroup_isolates)==repmat(anc_nti,1,sum(~outgroup_isolates)),2)==0);
    ingroup_calls = Calls(locs_without_ancestor,~outgroup_isolates);
    ingroup_calls(ingroup_calls==0) = nan; 
    anc_nti(locs_without_ancestor) = mode(ingroup_calls,2);

    
    %% Identify positions with fixed mutations
    fprintf(1,'Filtering positions and identifying fixed mutations...\n');

    % Determine which samples and positions have fixed mutations
    % Metric for the quality of evidence of a SNP
    [MutQual, MutQualIsolates] = ana_mutation_quality( Calls(:,~outgroup_isolates),Quals(:,~outgroup_isolates) );  % assumes quals already inverted
    % Boolean for where calls differ from the inferred ancestor
    fixedmutation_preliminary = ( (Calls~=repmat(anc_nti,1,Nsample)) & Calls>0 ); 
    fixedmutation_preliminary(MutQual < 1,:) = 0; % don't include positions where there is only one type of non-N call among cluster members
    fixedmutation_preliminary(:,outgroup_isolates)=0; % don't consider mutations in outgroup samples

    % Preliminary SNP positions
    goodpos_preliminary = find( sum( fixedmutation_preliminary,2 ) > 0 );
    p_goodpos_preliminary = p(goodpos_preliminary);
    
    
    %% Filter preliminary SNP positions where there is suspected recombination
    
    % Identify preliminary SNP positions on the reference genome where there is suspected recombination
    maNT_temp = maNT(goodpos_preliminary,:); % maNT only at preliminary goodpos
    maf_temp = maf(goodpos_preliminary,:); % maf only at preliminary goodpos
    anc_nti_temp = anc_nti(goodpos_preliminary); % anc nti only at preliminary SNPs
    mutantAF_temp=zeros(size(maNT_temp));
    mutantAF_temp(maNT_temp~=anc_nti_temp)=maf_temp(maNT_temp~=anc_nti_temp);
    mutantAF_temp = mutantAF_temp( :,~outgroup_isolates ); % remove outgroup isolates bcs we only care about within-lineage correlation
    % Find pairs of nearby preliminary SNP positions where mutant patterns are highly correlated
    p_involved_in_non_snp_event = identify_non_snp_events( recombination_block_size, correlation_cutoff, goodpos_preliminary, p, mutantAF_temp, GenomeLength, save_plots & show_clickable_tables, dir_diagnostic, this_cluster_name );
    if k==13
        % Manually remove a few more positions
        p_involved_in_non_snp_event = [ p_involved_in_non_snp_event; 2093929; 2093930 ]; % ancestral call is probably wrong for one of these positions (since outgroup doesn't have good calls here), which makes mutant allele frequency correlation negative
    end
    % Clear temp variables
    clear maNT_temp; clear maf_temp; clear anc_nti_temp; clear mutantAF_temp; 
    
    % Identify preliminary SNP positions on the reference genome where there is homology with the plasmid
    p_on_plasmid = plasmid_pos_on_genome; % imported above
    
    % Identify preliminary SNP positions on the reference genome in region recombination region L11
    p_on_regionL11 = p(goodpos_preliminary( find((p(goodpos_preliminary)>=1793873 & p(goodpos_preliminary)<=1793881) ) )); % Arolyn, 2019.07.01 (much shorter; basically only relevant to one clade)
    
    % Remove suspected recombination positions from fixedmutation and update goodpos
    p_from_recombination = unique( union( p_involved_in_non_snp_event, union( p_on_plasmid, p_on_regionL11 ) ) );
    % Update fixedmutation
    fixedmutation = fixedmutation_preliminary;
    fixedmutation(ismember(p,p_from_recombination),:) = 0; % remove mutations at recombination positions    
    % Indices relative to p
    goodpos = find( sum( fixedmutation,2 ) > 0 ); 
    p_goodpos = p(goodpos);
    
    % Save info
    if save_data
        save(['data_' this_cluster_name '_recombination'], ...
            'p_involved_in_non_snp_event', 'p_on_plasmid', 'p_on_regionL11', ...
            'p_from_recombination', ...
            'p_goodpos_preliminary', 'p_goodpos', ...
            'SampleNames', 'SampleNamesSimple', 'outgroup_isolates', ...
            'fixedmutation_preliminary', 'p', ...
            'Calls', 'anc_nti' )
    end

%     cd('..')
%     continue

    % Report number of mutations
    fprintf(1,['Number of positions with fixed mutations: ' num2str(length(goodpos)) '\n']);

    % Ancestral NTs only for SNP positions
    anc_nti_goodpos = anc_nti(goodpos); 

    
    %% Infer ancestor of this cluster
    
    % Initialize with mode of filtered calls over p_all from cluster members only
    Calls_lineage_ancestor = mode( Calls4annotations(:,~outgroup_isolates), 2 ); 
    % Set SNP positions to inferred ancestral allele
    Calls_lineage_ancestor( ismember(p_all,p_goodpos) ) = anc_nti_goodpos; 
    % not "Calls_lineage_ancestor( ismember(p_all,p) ) = anc_nti" because don't want to put calls at positions that don't exist for the ingroup
    % Set recombination positions to N's
    Calls_lineage_ancestor( ismember(p_all,p_from_recombination) ) = 0;

    
    %% No more analysis if there are no SNPs

    if isempty(goodpos)
        fprintf(1,'WARNING! There are no good positions with fixed mutations!\n');
        fid = fopen( 'goodpos_0.txt', 'wt' ); % way to leave record of warning
        fprintf(fid,'WARNING! goodpos is empty!!\n');
        fclose(fid);
        if save_data
            fprintf(1,'Saving data for this cluster...\n');
            annotation_full = struct;
            save(['data_' this_cluster_name '.mat'],'-v7.3',...
                'SampleNames','SampleNamesLong','SampleNamesSimple','outgroup_isolates','Nsample',...
                'subjects','zones','types','times','specimen_numbers','specimen_time_months','slst',... % removed 'grids'
                'Calls','p','counts','Quals','coverage','maf','maNT',...
                'indel_counter_delonly','indel_counter_both',...
                'anc_nti','goodpos',...
                'fixedmutation',...
                'Calls_lineage_ancestor',... 
                'anc_nti_goodpos','annotation_full' ... % these last two are empty if no SNPs
                );
        end        
        cd('..')
        continue % loop to next cluster
    end
    

    %% Get information about mutations
    
    % Pull information from genome at SNP positions (this version adapts the reference genome locally to have lineage-specific codons
    NTs4annotations = getNTs4annotations(Calls4annotations,p_all,Nsample,goodpos,outgroup_isolates,refnti_all);
    p4annotations = p_all;
    annotations = annotate_mutations_gb_lineage(p2chrpos(p(goodpos),ChrStarts),dir_ref_genome,p4annotations, NTs4annotations) ;
    annotation_full = append_annotations(annotations, anc_nti(goodpos), Calls(goodpos,:), counts(:,goodpos,:), fixedmutation(goodpos,:), promotersize) ; %% adds information about particular mutations observed, based on
    
    
    %% Make a clickable table to evaluate SNPs

    % Record number of SNPs found in this lineage
    num_snps_per_lineage(k) = numel(goodpos);

    % Metric for quality of position
    MutQual_fracambig = sum((Calls_ingroup<1),2)/Nsample_ingroup;
    % Metric for quality of position based on mac
    MutQual_maf = sum( maf(:,~outgroup_isolates)<0.85, 2)/Nsample_ingroup;
    % Metric for quality of position based on coverage vs expected coverage
    MutQual_cov = median( coverage_ingroup,2 )/median( median( coverage_ingroup,2 ) );
    
    % Option to reorder samples
    order = 1:numel(SampleNames);
    % Make clickable SNP table if requseted by the user
    if show_clickable_tables && numel(goodpos)>1
        QualSort=0; % set 1 to show mutations with lowest FQ scores up top, 0 to show in order on the genome
        clickable_snp_table_aro(annotation_full, Calls(goodpos,order), counts(:,goodpos,order), SampleNames(order), ScafNames, MutQual(goodpos), MutQual_fracambig(goodpos), QualSort);
        % If no annotations, use the next line instead:
        % clickable_snp_table_no_annotations(Calls(goodpos,order), counts(:,goodpos,order), SampleNames(order), ScafNames, MutQual(goodpos), QualSort);
        fig_num = fig_num + 1;
    end

%     cd('..')
%     continue
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 2: ANALYZE SNPs IN THIS CLUSTER %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% Get calls for analysis
    % Now that we have high-quality positions, call filters are a little looser in order to avoid unnecessary N's

    % Some filtering on maNT to get calls for analysis of goodpos
    Calls_for_analysis = maNT(goodpos,:); % Note that this includes the outgroup
    Calls_for_analysis( maf(goodpos,:) < min_maf_for_analysis ) = 0; % maf filter
    Calls_for_analysis( cov_fwd_strand(goodpos,:) < min_cov_each_strand_for_analysis | cov_rev_strand(goodpos,:) < min_cov_each_strand_for_analysis ) = 0; % stranded coverage filter
    Calls_for_analysis( coverage(goodpos,:) < min_coverage_for_analysis ) = 0; % coverage filter
    Calls_for_analysis( ( indel_counter_delonly(goodpos,:) >= (max_frac_reads_supp_dels*coverage(goodpos,:)) ) ) = 0; % del filter (but keep base from insertions)

    % Generate distance matrix
    Calls_for_analysis_ingroup = Calls_for_analysis(:,~outgroup_isolates);
    SampleNamesSimple_ingroup = SampleNamesSimple(~outgroup_isolates);
    dist_mat = zeros( Nsample_ingroup, Nsample_ingroup );
    for i=1:Nsample_ingroup
        for j=1:Nsample_ingroup
            dist_mat(i,j) = sum( Calls_for_analysis_ingroup(:,i)~=Calls_for_analysis_ingroup(:,j) & Calls_for_analysis_ingroup(:,i)>0 & Calls_for_analysis_ingroup(:,j)>0 );
        end
    end
    
    % Write csv file
    fid = fopen( ['dist_mat.csv' ], 'w' );
    % first line
    fprintf(fid, ',');
    for i=1:Nsample_ingroup
        fprintf(fid, [ SampleNamesSimple_ingroup{i} ',' ] );
    end
    fprintf(fid,'\n');
    % next lines
    for i=1:Nsample_ingroup
        fprintf(fid, [ SampleNamesSimple_ingroup{i} ',' ] );
        for j=1:Nsample_ingroup
            fprintf(fid, [ num2str(dist_mat(i,j)) ',' ] );
        end
        fprintf(fid,'\n');
    end
    fclose(fid);    
    
%     cd('..')
%     continue
    
    %% Make diagnostic plots summarizing mutations for this cluster
    fprintf(1,'Making figures summarizing mutations for this cluster...\n');

    % Make a calls plot
    plot_sampleset_calls_v2_diagnostic( Calls_for_analysis, 1:1:numel(goodpos), p(goodpos), ...
        SampleNames_tagoutgroup, 1:1:numel(SampleNames), this_cluster_name, dir_diagnostic, fig_num, save_plots )
    fig_num = fig_num+1; % update
    
    % Make a plot of mutation positions along genomes
    plot_sampleset_muts_on_genome( GenomeLength, p(goodpos), this_cluster_name, dir_diagnostic, fig_num, save_plots )
    fig_num = fig_num+1; % update
    
%     cd('..')
%     continue


    %% Make a tree of all samples and make additional annotated trees showing identity of each SNP 
    fprintf(1,'Making a tree...\n');

    % Directory for tree
    if ~exist('Tree','dir') % ARO: Groups subfolder
        mkdir('Tree')
    end

    cd('Tree');

    % Option to change set of samples included in tree or to change order
    samplestoplot = 1:numel(SampleNames); % currently everything

    % Collect calls for tree making and stores them as charactersw
    calls_for_tree=zeros(size(Calls_for_analysis));
    calls_for_tree(Calls_for_analysis>0)=NTs(Calls_for_analysis(Calls_for_analysis>0));
    calls_for_tree(Calls_for_analysis==0)='N';
    calls_for_tree=calls_for_tree(:,samplestoplot); % only grabs samples to plot, as defined above

    % Make the tree
    [treefilename_shortnames, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree, SampleNames_tagoutgroup(samplestoplot), this_cluster_name);
    [treefilename_longnames, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree, SampleNamesLong_tagoutgroup(samplestoplot), [this_cluster_name '-LongNames']);

    % Run treecounting
    run_treecounting_for_cluster( SampleNames_tagoutgroup, outgroup_isolates, ...
        calls_for_tree, anc_nti_goodpos, NTs, this_cluster_name, ...
        p, ChrStarts, goodpos, dir_scripts_aro_sh );

    cd('..')
    
    
    %% Save data for this cluster for later use

    if save_data
        fprintf(1,'Saving data for this cluster...\n');
        save(['data_' this_cluster_name '.mat'],'-v7.3',...
            'SampleNames','SampleNamesLong','SampleNamesSimple','outgroup_isolates','Nsample',...
            'subjects','zones','types','times','specimen_numbers','specimen_time_months','slst',... % removed 'grids'
            'Calls','p','counts','Quals','coverage','maf','maNT',...
            'indel_counter_delonly','indel_counter_both',...
            'anc_nti','goodpos',...
            'fixedmutation',...
            'Calls_lineage_ancestor',...
            'anc_nti_goodpos','annotation_full',...
            'Calls_for_analysis'...
            ); 
    end
    
    %% Return back to the parent Clusters directory
    cd('..')
    
end % end loop through clusters

save('data_num_muts.mat','num_snps_per_lineage')
fprintf(1,['Total SNPs detected: ' num2str(sum(num_snps_per_lineage)) '\n'])

cd('..')

