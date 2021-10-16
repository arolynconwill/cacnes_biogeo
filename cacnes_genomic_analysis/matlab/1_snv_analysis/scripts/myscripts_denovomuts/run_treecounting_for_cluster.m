function run_treecounting_for_cluster( SampleNames, outgroup_isolates, ...
    calls_for_tree, anc_nti_goodpos, NTs, this_cluster_name, ...
    p, ChrStarts, goodpos, dir_myscripts_sh )

%% Notes

% Python script requires that LCA samples come first in csv.


%% Reorder samples and combine outgroup into one sample

% Get ancestral calls
calls_for_tree_ancestor = zeros( size(anc_nti_goodpos) );
calls_for_tree_ancestor( anc_nti_goodpos>0 ) = NTs( anc_nti_goodpos );

% Get ingroup calls
calls_for_tree_ingroup = calls_for_tree(:,~outgroup_isolates);

% Combine together with ancestral calls first
calls_for_tree_new = [ calls_for_tree_ancestor, calls_for_tree_ingroup ];

% Reorder 
SampleNames_new = { 'inferred_ancestor', SampleNames{~outgroup_isolates} };


%% Directory for these trees

dir_for_trees = [ 'treecounting' ];
if ~exist(dir_for_trees,'dir') 
    mkdir(dir_for_trees)
end

% Add dnapars
copyfile( '../dnapars', 'dnapars' )
copyfile( '../dnapars.app', 'dnapars.app' )



%% Modify sample names

% python script requires all labels to start with the same character
SampleNames_tree = SampleNames_new;
for s=1:length(SampleNames_new)
    SampleNames_tree{s} = [ 'Z_' SampleNames_new{s} ]; % need to all start with the same letter...
end

% Remove colon from names
SampleNames_tree_carrot = SampleNames_tree; % put ^ instead of : so that these match actual leaf names 
for s=2:length(SampleNames_new)
    % find colon
    oldname = SampleNames_tree{s};
    i_colon = find( oldname == ':' );
    SampleNames_tree_carrot{s} = [ oldname(1:i_colon-1) '^' oldname(i_colon+1:end) ]; % need to all start with the same letter...
end


%% Make the regular tree using all SNPs

% Option to change set of samples included in tree or to change order
samplestoplot = 1:numel(SampleNames_new); % currently everything

% Make the tree
cd(dir_for_trees);
[treefilename, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree_new(:,samplestoplot), SampleNames_tree(samplestoplot), this_cluster_name);
cd('..')


%% Make a tree for each SNP location

cd(dir_for_trees);

dir_for_snp_trees = 'tree_counting_mini';
if exist(dir_for_snp_trees,'dir')
    rmdir(dir_for_snp_trees,'s') 
end
mkdir(dir_for_snp_trees)
cd(dir_for_snp_trees)

% Generate csv necessary for treecounting
fid=fopen('for_tree_labeling.csv','w');
fprintf(fid,'chr,pos');
contig_positions=p2chrpos(p, ChrStarts);
for i=1:numel(SampleNames_tree_carrot)
    fprintf(fid,[',' SampleNames_tree_carrot{i}]);
end
for i=1:numel(goodpos)
    fprintf(fid,['\n' num2str(contig_positions(goodpos(i),1)) ',' num2str(contig_positions(goodpos(i),2))]);
    for j=1:size(calls_for_tree_new,2)
        fprintf(fid,[',' calls_for_tree_new(i,j)]);
    end
end
fprintf(fid,'\n');
fclose(fid);

% Run treecounting
eval(['! python2.7 ' dir_myscripts_sh '/countMutations_aro.py ../' treefilename ' for_tree_labeling.csv'])

cd('..')
cd('..')


end