function [ annotation_pores, annotation_other ] = ...
    filter_snps_by_dist( this_lineage_name, dir_clusters )

%% Summary

% This function considered one lineage and determines which set of
% mutations happened within-pore. It then splits annotation_full into
% annotation_pore and annotation_other. Assumes positions are only mutated
% once within the lineage (though there are a handful of exceptions to
% this; for exceptions, the pore will *not* be counted as monophyletic*). 

% * Example: lineage 1, pore 245, not monophyletic because of pos
% 196/942726, which is 0/4 in two pore colonies and 2/4 in other colonies,
% with 2 as the inferred ancestor


%% Version history

% 2020.01.25, Arolyn


%% Load data

load([ dir_clusters '/' this_lineage_name '/' 'data_' this_lineage_name '.mat' ] )


%% Find specimen numbers for monophyletic pores

% Get extract and strip specimen numbers
type_extract = 1;
type_strip = 2;
specimen_numbers_pores = unique( specimen_numbers( ( types == type_extract | types == type_strip ) & ~outgroup_isolates ) );
if isequal( this_lineage_name,'Cluster-1-A-212' )
    specimen_numbers_pores = setdiff( specimen_numbers_pores, [31 41] ); % remove hypermutators 
    fprintf(1,['Note: excluding hypermutator samples in lineage 1 pore specimen list...\n'])
end

% Get boolean showing which positions differ from the inferred lineage MRCA for each colony in that lineage
if isequal( this_lineage_name,'Cluster-1-A-212' )
    hypermutators = specimen_numbers == 31 | specimen_numbers == 41;
    samples_to_include = ~hypermutators' & ~outgroup_isolates;
    Calls_lineage = Calls_for_analysis(:,samples_to_include); % exclude hypermutators and outgroup isolates
    specimen_numbers_lineage = specimen_numbers( samples_to_include );
    num_samples_in_lineage = sum(samples_to_include);
    fprintf(1,['Note: excluding hypermutator samples in lineage 1 diff_mrca...\n'])
else
    Calls_lineage = Calls_for_analysis(:,~outgroup_isolates); % exclude outgroup isolates
    specimen_numbers_lineage = specimen_numbers( ~outgroup_isolates );
    num_samples_in_lineage = sum(~outgroup_isolates);
end

% Boolean showing which positions differ from the inferred lineage MRCA for each colony in that lineage
diff_mrca = ( ...
    ( Calls_lineage ~= repmat(anc_nti_goodpos,1,num_samples_in_lineage) ) ...
    & Calls_lineage > 0 ...
    );
% Determine which pore specimens are monophyletic
monophyletic_pores = [];
within_pore_muts_all = [];
for n=1:numel(specimen_numbers_pores)
    next_pore = specimen_numbers_pores(n);
    pore_indices = ( specimen_numbers_lineage == next_pore );
    [ is_mono, within_pore_muts ] = test_if_pore_is_monophyletic( diff_mrca, pore_indices, p, goodpos );
    if is_mono % if pore monophyletic
        monophyletic_pores(end+1) = next_pore;
        if numel(within_pore_muts)>0 % if within-pore SNPs exist
            within_pore_muts_all = [ within_pore_muts_all; within_pore_muts ];
        end
    end
end


%% Split annotation_full into annotation_pore and annotation_other

% Find indices in annotation_full for within-pore mutations
within_pore_muts_indices = ismember( p(goodpos), within_pore_muts_all );
% Split annnotation_full
annotation_pores = annotation_full( within_pore_muts_indices );
annotation_other = annotation_full( ~within_pore_muts_indices );


end