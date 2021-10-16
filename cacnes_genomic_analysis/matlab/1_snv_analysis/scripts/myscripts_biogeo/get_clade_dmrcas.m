function [ dmrca_clade_max, dmrca_clade_median, dmrca_intrapore, dmrca_interpore_pairs_list, dmrca_interpore_min_list, num_pores, pore_specnum ] = ...
    get_clade_dmrcas( this_cluster_name, ...
    mult_spec_nums, mult_spec_mults )

%% General info

type_extract = 1;
type_strip = 2;
type_scrape = 3;



%%
%% Case 1: no SNVs in this lineage
%%

if exist(['2_snvs/' this_cluster_name '/goodpos_0.txt'], 'file' )
    fprintf(1,['Warning! No SNVs for this cluster!\n'])
    dmrca_clade_max = 0;
    dmrca_clade_median = 0;
    load(['2_snvs/' this_cluster_name '/data_' this_cluster_name '.mat']);

    % List of unique specimens and types
    [ specimens_list, specimens_list_indices ] = unique( specimen_numbers(~outgroup_isolates) );
    hypers_boolean = specimens_list == 31 | specimens_list == 41 ; % for removing hypermutator specimens
    specimens_list = specimens_list( ~hypers_boolean );
    specimens_list_indices = specimens_list_indices( ~hypers_boolean );
    types_list = types( specimens_list_indices );
    % Find out which ones are pores
    types_list_pores_bool = ( types_list == type_strip ) | ( types_list == type_extract );
    types_list_pores_indices = find( types_list_pores_bool );
    specimens_list_pores = specimens_list( types_list_pores_indices );
    fprintf(1,['Number of pore specimens: ' num2str(numel(specimens_list_pores)) '\n'])
    % Find out which pore specimens are actually one pore
    specimens_list_mults = arrayfun(@(x) mult_spec_mults(x==mult_spec_nums) , specimens_list_pores);
    specimens_list_pores_keep = specimens_list_pores( specimens_list_mults==1 ); % only keep singles
    fprintf(1,['Number of singlet pore specimens only: ' num2str(numel(specimens_list_pores_keep)) '\n'])
    % Find out which pores have more than one sample
    pore_specimen_sample_num = arrayfun(@(x) sum(x==specimen_numbers), specimens_list_pores_keep);
    specimens_list_pores_keep = specimens_list_pores_keep( pore_specimen_sample_num>1 );
    % Number of pores for analysis
    num_pores = numel( specimens_list_pores_keep );
    fprintf(1,['Number of singlet pore specimens with more than one sample only: ' num2str(num_pores) '\n'])
 
    if num_pores>1 % needs pair for interpore
        dmrca_intrapore = zeros( num_pores,1 );
        num_pairs = nchoosek( num_pores,2 );
        dmrca_interpore_pairs_list = zeros( num_pairs,1 );
        dmrca_interpore_min_list = zeros( num_pores,1 );
        pore_specnum = specimens_list_pores_keep;
    elseif num_pores==1
        dmrca_intrapore = zeros( num_pores,1 );
        dmrca_interpore_pairs_list = nan;
        dmrca_interpore_min_list = nan;
        pore_specnum = specimens_list_pores_keep;
    else
        dmrca_intrapore = nan;
        dmrca_interpore_pairs_list = nan;
        dmrca_interpore_min_list = nan;
        pore_specnum = nan;
    end
    
    return 
end



%%
%% Case 2: SNVs exist in this lineage
%%

% Load data
load(['2_snvs/' this_cluster_name '/data_' this_cluster_name '.mat']);

% Using filtered maNT from above, but witout the outgroup isolates
calls_clade=Calls_for_analysis(:,~outgroup_isolates);
diff_mrca=((calls_clade~=repmat(anc_nti_goodpos,1,sum(~outgroup_isolates))) & calls_clade>0);


%% Get lineage dMRCA

% Define most recent common ancestor
mrca=mode(Calls(:,~outgroup_isolates),2); % to include nonvariable positions
mrca(goodpos)=anc_nti_goodpos; % setting proper values for variable positions
mrca_goodpos = mrca(goodpos); % save MRCA at variable positions in lineage only

%Some useful printouts
fprintf(1,['Number of samples (w/o outgroup): ' num2str(numel(SampleNames)-sum(outgroup_isolates)) '\n']);
fprintf(1,['De novo muts: ' num2str(numel(goodpos)) '\n']); % Aro_Change: numel goodpos rather than numel annotation_full
fprintf(1,['Min dMRCA: ' num2str(min(sum(diff_mrca))) '\n']);
fprintf(1,['Median dMRCA: ' num2str(median(sum(diff_mrca))) '\n']);
fprintf(1,['Max dMRCA: ' num2str(max(sum(diff_mrca))) '\n']);

% Clade dMRCA = maximum distance from inferred ancestor
% But remove hypermutators for this part
hypers_bool = specimen_numbers(~outgroup_isolates) == 31 | specimen_numbers(~outgroup_isolates) == 41 ; % for removing hypermutator specimens
% Get max dMRCA of non-hypermutator sample
dmrca_clade_max = max(sum(diff_mrca(:,~hypers_bool)));
dmrca_clade_median = median(sum(diff_mrca(:,~hypers_bool)));


%% Find out how many valid pores there are (single follicle; at least 2 colonies)

% List of unique specimens and types
[ specimens_list, specimens_list_indices ] = unique( specimen_numbers(~outgroup_isolates) );
hypers_boolean = specimens_list == 31 | specimens_list == 41 ; % for removing hypermutator specimens
specimens_list = specimens_list( ~hypers_boolean );
specimens_list_indices = specimens_list_indices( ~hypers_boolean );
types_list = types( specimens_list_indices );

% Find out which ones are pores
types_list_pores_bool = ( types_list == type_strip ) | ( types_list == type_extract );
types_list_pores_indices = find( types_list_pores_bool );
specimens_list_pores = specimens_list( types_list_pores_indices );
fprintf(1,['Number of pore specimens: ' num2str(numel(specimens_list_pores)) '\n'])

% Find out which pore specimens are actually one pore
specimens_list_mults = arrayfun(@(x) mult_spec_mults(x==mult_spec_nums) , specimens_list_pores);
specimens_list_pores_keep = specimens_list_pores( specimens_list_mults==1 ); % only keep singles
fprintf(1,['Number of singlet pore specimens only: ' num2str(numel(specimens_list_pores_keep)) '\n'])

% Find out which pores have more than one sample
pore_specimen_sample_num = arrayfun(@(x) sum(x==specimen_numbers), specimens_list_pores_keep);
specimens_list_pores_keep = specimens_list_pores_keep( pore_specimen_sample_num>1 );

% Number of pores for analysis
num_pores = numel( specimens_list_pores_keep );
fprintf(1,['Number of singlet pore specimens with more than one sample only: ' num2str(num_pores) '\n'])

% Check to see if there are a sufficient number of pores to continue analysis
if num_pores == 0 % need a pair in order to get at least one interpore dMRCA
    fprintf(1,['No valid pores in clade.\n\n'])
    dmrca_intrapore = nan;
    dmrca_interpore_pairs_list = nan;
    dmrca_interpore_min_list = nan;
    pore_specnum = nan;
    return
elseif num_pores == 1 % can only to intrapore
    fprintf(1,['Only one valid pore in clade.\n\n'])
    this_pore = specimens_list_pores_keep(1);
    calls_pore = Calls_for_analysis( :,( specimen_numbers == this_pore ) );
    % Find number of positons that differ
    varpos_pore = calls_pore ~= repmat( mode(calls_pore,2), 1, sum( specimen_numbers == this_pore ) ) & calls_pore>0 ;
    goodpos_pore = ( sum(varpos_pore,2)>0 ); % which positions are variable
    % Infer ancestor of pore
    mrca_pore =  mode(calls_pore,2); % fine for positions that don't differ
    mrca_pore(goodpos_pore) = mrca_goodpos( goodpos_pore );
    % Find mean dMRCA of the pore
    diffs_pore = calls_pore ~= repmat( mrca_pore, 1, sum( specimen_numbers == this_pore ) ) & calls_pore>0 ;
    % Record data
    dmrca_intrapore = mean(sum(diffs_pore)); 
    dmrca_interpore_pairs_list = nan;
    dmrca_interpore_min_list = nan;
    pore_specnum = specimens_list_pores_keep; 
    return
end

% Get intrapore dMRCA and record pore MRCA
dmrca_intrapore = -ones( num_pores,1 ); % initialize as -1 (won't update for pores that don't have more than one sample)
pore_mrcas = zeros( num_pores, numel(goodpos) );
pore_specnum = [];
for p=1:num_pores % loop through pore specimens
    
    this_pore = specimens_list_pores_keep(p); %OLD: specimens_list( types_list_pores_indices(p) );
    pore_specnum(end+1) = this_pore;
    
    if sum( specimen_numbers == this_pore ) > 1 % look at specimens where there is more than one sample
        calls_pore = Calls_for_analysis( :,( specimen_numbers == this_pore ) );
        % Find number of positons that differ
        varpos_pore = calls_pore ~= repmat( mode(calls_pore,2), 1, sum( specimen_numbers == this_pore ) ) & calls_pore>0 ;
        goodpos_pore = ( sum(varpos_pore,2)>0 ); % which positions are variable
        % Infer ancestor of pore
        mrca_pore =  mode(calls_pore,2); % fine for positions that don't differ
        mrca_pore(goodpos_pore) = mrca_goodpos( goodpos_pore );
        % Find mean dMRCA of the pore
        diffs_pore = calls_pore ~= repmat( mrca_pore, 1, sum( specimen_numbers == this_pore ) ) & calls_pore>0 ;
        dmrca_pore = mean(sum(diffs_pore)); 
        dmrca_intrapore(p) = dmrca_pore;
    else
        mrca_pore = Calls_for_analysis( :,( specimen_numbers == this_pore ) );
    end
    pore_mrcas( p,: ) = mrca_pore;
end

dmrca_intrapore = dmrca_intrapore( dmrca_intrapore > -1 );

%% Get interpore pairwise dMRCA

% Inferred ancestor of whole clade
mrca_clade_goodpos = mrca(goodpos);

% Find pairwise mean interpore dMRCA
dmrca_interpore = -ones( num_pores );
for p1=1:num_pores
    for p2=1:num_pores
        p1_mrca = pore_mrcas(p1,:);
        p2_mrca = pore_mrcas(p2,:);
        % Infer ancestor between two pore
        p1p2_mrca = mode( [ p1_mrca; p2_mrca ] ); % for nonvariable
        p1p2_mrca( p1_mrca ~= p2_mrca ) = mrca_clade_goodpos( p1_mrca ~= p2_mrca ); % anc for diffs
        % Find dMRCA
        p1_dmrca = sum( ( p1_mrca ~= p1p2_mrca ) & ( p1_mrca > 0 ) );
        p2_dmrca = sum( ( p2_mrca ~= p1p2_mrca ) & ( p2_mrca > 0 ) );
        % Record interpore dMRCA
        dmrca_interpore( p1,p2 ) = (p2_dmrca+p1_dmrca)/2;
    end
end

% List of all pairwise mean interpore dMRCAs
dmrca_interpore_pairs_list = dmrca_interpore( 1:1:num_pores>[1:1:num_pores]' );

% List of min interpore dMRCAs (i.e. find closest pore)
dmrca_interpore_min_list = -ones( num_pores,1 );
for p=1:num_pores
    dmrca_interpore_min_list(p) = min( dmrca_interpore( p, setdiff(1:1:num_pores,p) ) );
end

fprintf(1,['\n\n'])

end

