%% Set up directories

% Directory for my scripts:
workingdir = char(pwd); % working directory
MYSCRIPTSDIRECTORY = [workingdir '/myscripts' ];
path(MYSCRIPTSDIRECTORY,path);


%% Save to make on cluster

load('clade_tree_cluster_data.mat')


%%

% ARO: Make Trees subfolder to hold all files about to be generated
if ~exist('Trees','dir') % ARO: Groups subfolder
    mkdir(['Trees'])
end
% ARO: Put dnapars in Trees directory if not already there
if ~exist('Trees/dnapars','file')
    copyfile dnapars Trees
end

% Put tree in own directory
cd('Trees') % Put all of this in Trees subfolder

% Make tree
fprintf('\nStarting dnapars...\n') 
[treefilename, UsedTreeNames] = generate_parsimony_tree_old(calls_for_tree',names_tree);

cd('..')
