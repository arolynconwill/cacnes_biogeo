%% SUMMARY

% This script makes lineage trees with easily interpretable sample names.


%% Set up directories and environment

% Directory for trees
dir_treemaking = '7_trees_lineages';
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

% Nucleotides: 1=A, 2=T, 3=C, 4=G
NTs = 'ATCG';


%% Load lineage info

% Load lineage names
dir_clusters = '2_snvs';
load( [ dir_clusters '/'  'cluster_names.mat' ],'cluster_names','clusters_all') % lineage names
lineage_names = cluster_names; clear cluster_names; 
cluster_sizes = cellfun(@(x) numel(x), clusters_all);

% List of path to files containing annotation_full for each lineage
paths_to_files = cellfun(@(lineage_name) [ dir_clusters '/' lineage_name '/' 'data_' lineage_name '.mat' ], lineage_names, 'UniformOutput', false );


%% Loop through lineages

for c=1:numel(lineage_names)
        
    % Load lineage data
    this_lineage_name = lineage_names{c};
    fprintf(1,['Treemaking for ' this_lineage_name '...\n'])
    load(['2_snvs/' this_lineage_name '/data_' this_lineage_name '.mat']);

    % Skip lineages if no SNVs
    if isempty(goodpos) 
        fprintf(1,['No SNPs in ' this_lineage_name '...\n'])
        continue
    end

    % Directory setup
    dir_save = [ dir_treemaking '/' this_lineage_name '_trees' ]; 
    if ~exist( dir_save, 'dir' )
        mkdir( dir_save )
    end

    %% Process sample names and outgroup

    % Implementation notes
    % % For outgroup: only keeps one; renames it "inferred lineage ancestor"
    % % For others: rename in a standardized way
    % % % Single pore specimens: "pore-[spec#]_colony-[col#]"
    % % % Multi pore specimens: "multipore-..."
    % % % Scrape specimens:  "mixedsurface-[spec#]_colony-[col#]"

    % Initialize
    SampleNames_tree = SampleNamesSimple;
    calls_for_tree=zeros(size(Calls_for_analysis));
    calls_for_tree(Calls_for_analysis>0)=NTs(Calls_for_analysis(Calls_for_analysis>0));
    calls_for_tree(Calls_for_analysis==0)='N';

    % Replace outgroup colonies with the inferred lineage ancestor
    calls_for_tree_ancestor = zeros( size(anc_nti_goodpos) );
    calls_for_tree_ancestor( anc_nti_goodpos>0 ) = NTs( anc_nti_goodpos );
    % Get ingroup calls
    calls_for_tree_ingroup = calls_for_tree(:,~outgroup_isolates);
    % Combine together with ancestral calls first
    calls_for_tree_new = [ calls_for_tree_ancestor, calls_for_tree_ingroup ];
    % Reorder 
    SampleNames_new = { 'inferred_lineage_ancestor', SampleNames_tree{~outgroup_isolates} };

    
    %% Make a tree

    cd( dir_save );

    % Make the tree
    [treefilename, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree_new, SampleNames_new, [this_lineage_name '_clean-names']);

    cd( '../..' )
    
end


%% Special tree for Lineage A-1

for c=1

    % Load lineage data
    this_lineage_name = lineage_names{c};
    load(['../Clusters/' this_lineage_name '/data_' this_lineage_name '.mat']);

    % Skip lineages if no SNVs exist or if no pore samples
    %     if isempty(goodpos) || sum( types(~outgroup_isolates) == 3 )==sum(~outgroup_isolates)
    %         continue
    %     else
    %         fprintf(1,['Making a tree for lineage ' this_lineage_name '...\n'])
    %     end

    % Directory setup
    dir_save = [ this_lineage_name '_trees' ]; 
    if ~exist( dir_save, 'dir' )
        mkdir( dir_save )
    end

    %% Make super simple "simplenames"

    SampleNamesSimple_short = cell( size( SampleNamesSimple ) );
    for s=1:numel(SampleNamesSimple)
        next_name = SampleNamesSimple{s};
        next_indices = strfind(next_name,'_');
        new_name = next_name(next_indices(1)+1:end);
        SampleNamesSimple_short{s} = new_name;
    end

    %% Process sample names and outgroup

    % Initialize
    SampleNames_tree = SampleNamesSimple_short;
    calls_for_tree=zeros(size(Calls_for_analysis));
    calls_for_tree(Calls_for_analysis>0)=NTs(Calls_for_analysis(Calls_for_analysis>0));
    calls_for_tree(Calls_for_analysis==0)='N';

    % Replace outgroup colonies with the inferred lineage ancestor
    calls_for_tree_ancestor = zeros( size(anc_nti_goodpos) );
    calls_for_tree_ancestor( anc_nti_goodpos>0 ) = NTs( anc_nti_goodpos );
    % Get ingroup calls
    calls_for_tree_ingroup = calls_for_tree(:,~outgroup_isolates);
    % Combine together with ancestral calls first
    calls_for_tree_new = [ calls_for_tree_ancestor, calls_for_tree_ingroup ];
    % Reorder 
    SampleNames_new = { 'inferred_lineage_ancestor', SampleNames_tree{~outgroup_isolates} };


    %% Make a tree

    cd( dir_save );

    % Make the tree
    [treefilename, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree_new, SampleNames_new, [this_lineage_name '_clean-names']);

    cd( '..' )

end
