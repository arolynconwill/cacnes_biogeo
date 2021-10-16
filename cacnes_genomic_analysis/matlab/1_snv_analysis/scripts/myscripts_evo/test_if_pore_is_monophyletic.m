function [ is_mono, within_pore_muts ] = test_if_pore_is_monophyletic( diff_mrca, pore_indices, p, goodpos )

%% Summary

% This function determines if a pore specimen is monophyletic. It allows
% for other colonies to identical to the inferred pore ancestor. It returns
% a boolean (is_mono) indicating if the pore is monophyletic, as well as a
% list of positions on the genome on which there are within-pore mutations.
% This script assumes a given position on the genome is only mutated once
% within a lineage.

%% Version history

% 2020.01.26, Arolyn: switched to nested ifs
% 2020.01.25, Arolyn: original version


%% Determine if pore colonies are monophyletic

% Number of colonies in pore
num_colonies = sum(pore_indices);

% Cases to determine if pore is monophyletic
if num_colonies == 1 % Checks if only one pore colony

    is_mono = 1;
    within_pore_muts = [];

else 

    % Checks for within-pore mutations
    diff_mrca_pores_only = diff_mrca( :,pore_indices );
    
%     size(diff_mrca_pores_only)
%     %conditional breakpoint numel(diff_mrca_pores_only)<1
%     temp1=diff_mrca_pores_only(:,1) 
%     temp2=repmat( temp1, 1,num_colonies )
%     within_pore_muts_bool = sum( diff_mrca_pores_only ~= temp2,2 )
%     within_pore_muts_bool = ( within_pore_muts_bool<0 )    
    
    within_pore_muts_bool = sum( diff_mrca_pores_only ~= repmat( diff_mrca_pores_only(:,1), 1,num_colonies ),2 ) > 0;

    if sum( within_pore_muts_bool ) == 0 % no within pore mutations

        is_mono = 1;
        within_pore_muts = [];
        
    else
        
        % Checks that all within-pore SNPs are exclusive to that pore
        diff_mrca_others = diff_mrca( within_pore_muts_bool, ~pore_indices );
        within_pore_muts_shared = ( sum( diff_mrca_others,2 ) > 0 );

        if sum( within_pore_muts_shared ) > 0 % not monophyletic

            is_mono = 0;
            within_pore_muts = [];
            
        else % is monophyletic
            
            is_mono = 1;
            within_pore_muts = p( goodpos( within_pore_muts_bool ) );
            
        end
        
    end


end


end