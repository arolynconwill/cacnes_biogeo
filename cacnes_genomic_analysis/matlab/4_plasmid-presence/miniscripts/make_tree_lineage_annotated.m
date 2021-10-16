function make_tree_lineage_annotated( this_cluster_name, ...
    SampleNames_regular, SampleNamesSimple_plasmid, ...
    NTs, dir_trees, dir_clusters )

%% Load data

fprintf(1,[this_cluster_name '\n'])
load([ dir_clusters '/' this_cluster_name '/data_' this_cluster_name '.mat']);


%% Get calls for tree

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
    

%% Start updating SampleNames for annotated tree

% Updated samples names
SampleNames_annotated = SampleNames_new;

% Loop through and rename
for n=1:numel(SampleNames_annotated)
    next_name = SampleNames_annotated{n};
    next_index = find( ismember( SampleNames_regular, next_name ) );
    if ~isempty( next_index )
        SampleNames_annotated{n} = SampleNamesSimple_plasmid{ next_index };
    else
        SampleNames_annotated{n} = [ SampleNames_annotated{n} ];
    end
end


%% Make the tree

cd(dir_trees);

if ~exist(this_cluster_name,'dir') % ARO: Groups subfolder
    mkdir(this_cluster_name)
end
cd(this_cluster_name);

% Make the tree
[treefilename, UsedTreeNames] = generate_parsimony_tree_old_aro(calls_for_tree_new, SampleNames_annotated, [this_cluster_name '_with-plasmid-info']);

cd('..')
cd('..')

end


