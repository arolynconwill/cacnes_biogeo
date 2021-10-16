%% SUMMARY

% This script makes strain (superSLST) trees with easily interpretable sample names.


%% Set up directories and environment

% Directory for trees
dir_treemaking = '7_trees_straintypes';
if ~exist( dir_treemaking, 'dir' )
    mkdir( dir_treemaking )
end

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Directory for my scripts:
dir_scripts_myscripts = [ pwd '/scripts/myscripts_denovomuts'];
path(dir_scripts_myscripts,path);

% Add dnapars (parsimony tool)
if ~exist('dnapars.app','file')
    copyfile( 'tools/dnapars.app', [ dir_treemaking '/dnapars.app' ] )
    copyfile( 'tools/dnapars', [ dir_treemaking '/dnapars' ] )
end


%% Set parameters

% Implementation notes
% % same as in identify_de_novo_muts for finding within-lineage SNVs
% % filtering regions homologous to plasmids
% % not doing additional recombination detection

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
% For filtering plasmid positions
load('data/plasmid_pos_to_mask.mat','plasmid_pos_on_genome'); % positions on reference genome with homology to plasmid gene content

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';


%% Load data

% Load data for all samples
p_all = load('data/cluster_step_variables','p_all'); p_all = p_all.p_all;
counts_all = load('data/cluster_step_variables','counts_final'); counts_all = counts_all.counts_final;
Quals_all = load('data/cluster_step_variables','Quals_final'); Quals_all = Quals_all.Quals_final;
coverage_all = load('data/cluster_step_variables','coverage_final'); coverage_all = coverage_all.coverage_final;
maNT_all = load('data/cluster_step_variables','maNT_final'); maNT_all = maNT_all.maNT_final;
maf_all = load('data/cluster_step_variables','maf_final'); maf_all = maf_all.maf_final;
indel_counter_all = load('data/cluster_step_variables','indel_counter_final'); indel_counter_all = indel_counter_all.indel_counter_final;

% Load lineage names
dir_clusters = '2_snvs';
load( [ dir_clusters '/'  'cluster_names.mat' ]) % lineage names
lineage_names = cluster_names; clear cluster_names; 
cluster_sizes = cellfun(@(x) numel(x), clusters_all);
clusters_all_superslst = cellfun(@(x) x(1), clusters_all_slst);

% Load sample names
SampleNamesSimple_all = load('data/cluster_step_variables','SampleNamesSimple_final'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_final;
SampleNamesLong_all = load('data/cluster_step_variables','SampleNamesLong_final'); SampleNamesLong_all = SampleNamesLong_all.SampleNamesLong_final;
% make simple names with lineage tag
SampleNamesSimpleLineage_all = cell( size( SampleNamesSimple_all ) );
for k=1:numel(clusters_all)
    temp_indices = clusters_all{k};
    for i=1:numel(temp_indices) 
        temp_index = temp_indices(i);
        SampleNamesSimpleLineage_all{temp_index} = [ cluster_names_new{k} '_' SampleNamesSimple_all{temp_index} ];
    end
end


%% Define straintypes

strain_list = 'ACDFHK'; %unique( clusters_all_superslst );
strain_list_outgroup = [ strain_list(1:end-1) strain_list(end-1) ];


%% Make a tree for each superSLST

make_tree = true; % takes as long time so can switch on/off

for s=1:numel(strain_list)
    
    %% Get samples in this strain type
    my_strain = strain_list(s);
    my_clusters = find( my_strain == clusters_all_superslst );
    my_colonies = [];
    for i=1:numel(my_clusters)
        my_colonies = [ my_colonies; clusters_all{my_clusters(i)} ];
    end
    SampleNamesSimple=SampleNamesSimpleLineage_all(my_colonies);
    Nsample = numel(my_colonies);
    
    %% Get outgroup calls
    out_strain = strain_list_outgroup(s);
    out_clusters = find( my_strain == clusters_all_superslst );
    Nout = numel(out_clusters);
    Calls_outgroup = zeros( Nout,numel(p_all) );
    for i=1:Nout
        load( [ dir_clusters '/' lineage_names{i} '/data_' lineage_names{i} '.mat' ], 'Calls_lineage_ancestor' );
        Calls_outgroup(i,:) = Calls_lineage_ancestor;
    end
    Calls_outgroup = mode( Calls_outgroup );
    
    %% Find SNVs within this strain type
    
    % Get quals and counts for filtering
    Quals = Quals_all(:,my_colonies); % quals at each position
    indel_counter_both = squeeze(indel_counter_all(1,:,my_colonies)); % num reads supporting indels
    indel_counter_delonly = squeeze(indel_counter_all(2,:,my_colonies)); % num reads supporting deletions only
    maf = maf_all(:,my_colonies); % maf at each position
    counts = counts_all(:,:,my_colonies); % counts at each position
    cov_fwd_strand = squeeze(sum(counts(1:4,:,:)));
    cov_rev_strand = squeeze(sum(counts(5:8,:,:)));

    % Initialize from maNT
    Calls = maNT_all(:,my_colonies); % all positions in p, but only for this cluster (plus outgroup)
    % Calls filtering (part 1/2): for high quality calls
    Calls( Quals < min_qual_for_call ...
        | indel_counter_both >= max_frac_reads_supp_indels*(cov_fwd_strand+cov_rev_strand) ...
        | indel_counter_delonly >= max_frac_reads_supp_dels*(cov_fwd_strand+cov_rev_strand) ...
        | maf < min_maf_for_call ...
        | cov_fwd_strand < min_cov_each_strand_for_call ...
        | cov_rev_strand < min_cov_each_strand_for_call ...
        ) = 0; % set to N

    % Find nonvariable positions
    nonvariablep = ( ...
        sum(Calls==1,2)==sum(Calls~=0,2) ... % all non-ambiguous calls are A
        | sum(Calls==2,2)==sum(Calls~=0,2) ... % ... %
        | sum(Calls==3,2)==sum(Calls~=0,2) ... % ... C
        | sum(Calls==4,2)==sum(Calls~=0,2) ... % ... G
        | sum(Calls==0,2)==Nsample ); % all calls are ambiguous
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
    coverage = coverage_all(~nonvariablep,my_colonies); % coverage at each position
    maNT = maNT_all(~nonvariablep,my_colonies); % maNT at each position
    % And downsize outgroup calls
    Calls_outgroup_nonvar = Calls_outgroup( ~nonvariablep );

    % Calls filtering (part 2/3): only keep high quality positions 
    low_quality_positions = ( ...
        sum( Calls<1, 2 ) >= Nsample*max_fraction_ambiguous_samples ... % too many N's
        | median( coverage,2 ) < min_median_coverage_position ... % not enough coverage
        );
    Calls( low_quality_positions,: ) = 0; 
    fprintf(1,['Number of low quality positions: ' num2str(sum(low_quality_positions)) '\n']);

    % Calls filtering (part 3/3)
    Calls( ismember( p,plasmid_pos_on_genome' ),: ) = 0; % do not consider plasmid positions

    % SNV positions =  remaining variable positions
    goodpos = find( ~( ...
        sum(Calls==1,2)==sum(Calls~=0,2) ... % all non-ambiguous calls are A
        | sum(Calls==2,2)==sum(Calls~=0,2) ... % ... %
        | sum(Calls==3,2)==sum(Calls~=0,2) ... % ... C
        | sum(Calls==4,2)==sum(Calls~=0,2) ... % ... G
        | sum(Calls==0,2)==Nsample ) ); % all calls are ambiguous
    p_goodpos = p(goodpos);
    fprintf(1,['Number of nonvariable positions: ' num2str(numel(goodpos)) '\n']);
    
    % Downsize calls again
    Calls_goodpos = Calls( goodpos,: );
    Calls_outgroup_goodpos = Calls_outgroup_nonvar( goodpos );
    
    %% Make a tree
    
    % Directory setup
    dir_save = [ dir_treemaking '/' 'SLST-' my_strain ]; 
    if ~exist( dir_save, 'dir' )
        mkdir( dir_save )
    end

    % Initialize
    SampleNames_tree = horzcat( {'outgroup'}, SampleNamesSimple );
    Calls_goodpos_with_outgroup = [ Calls_outgroup_goodpos', Calls_goodpos ];
    calls_for_tree=zeros(size(Calls_goodpos_with_outgroup));
    calls_for_tree(Calls_goodpos_with_outgroup>0)=NTs(Calls_goodpos_with_outgroup(Calls_goodpos_with_outgroup>0));
    calls_for_tree(Calls_goodpos_with_outgroup==0)='N';

    % Make a tree
    if make_tree
        cd( dir_save );
        [treefilename, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree, SampleNames_tree, ['SLST-' my_strain]);
        cd( '../..' )
    end
    
    %% Distance matrix (for scale bars)
    
    % Generate distance matrix
    dist_mat = zeros( Nsample, Nsample );
    for i=1:Nsample
        for j=1:Nsample
            dist_mat(i,j) = sum( Calls_goodpos(:,i)~=Calls_goodpos(:,j) & Calls_goodpos(:,i)>0 & Calls_goodpos(:,j)>0 );
        end
    end
    
    % Write csv file
    fid = fopen( [dir_save '/dist_mat.csv' ], 'w' );
    % first line
    fprintf(fid, ',');
    for i=1:Nsample
        fprintf(fid, [ SampleNamesSimple{i} ',' ] );
    end
    fprintf(fid,'\n');
    % next lines
    for i=1:Nsample
        fprintf(fid, [ SampleNamesSimple{i} ',' ] );
        for j=1:Nsample
            fprintf(fid, [ num2str(dist_mat(i,j)) ',' ] );
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %% Other
    
    % Extra: plot p_goodpos
    plot(p_goodpos); xlabel('SNV index'); ylabel('pos on genome');
    print( [ dir_save '/' 'plot_goodpos.png'], '-dpng')

end



