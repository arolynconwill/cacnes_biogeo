function [ mutationmatrix, mut_observed, typecounts, prob_nonsyn ] = ...
    mutation_spectrum_module( goodpos, anc_nti_goodpos, Calls_for_analysis, annotation_full, p, ...
    clade_name, directory_save, save_data )

%% Summary

% This script computes the mutation spectrum for a given lineage



%% Execution

% Initialize vectors
mutationmatrix=zeros(4,4); % ATCG -> ATCG matrix; gathered into mut_observed (6 entry vector)
typecounts=zeros(4,1); %NSPI % used to calculate prob_nonsyn

% Count mutations
for i=1:numel(goodpos) % loop through positions at which there is a SNP
    anc = anc_nti_goodpos(i); % ancestor NT at this position
    new = unique(Calls_for_analysis(i,:)); %all NTs found at this position
    if anc > 0 & sum(anc == new) % if ancestor is known and ancestor NT is in annotation_full
        % Remove ancestor from new
        new = new( new ~= anc );
        % Remove N from new
        new = new( new ~= 0 );
        if numel(new)==0 % if nothing left
            fprintf(1,['Warning: No mutation found at p=' num2str(p(goodpos(i))) '\n'])
        elseif numel(new)==1 % if one mutation
            % Update mutation matrix
            mutationmatrix(anc,new)=mutationmatrix(anc,new)+1;
            % Count type (use annotation_full)
            j = find([annotation_full.pos]==p(goodpos(i))); % in case this position is at a different index in annotation_full
            if annotation_full(j).type=='N' % nonsyn mut
                typecounts(1)=typecounts(1)+1;
            elseif annotation_full(j).type=='S' % syn mut
                typecounts(2)=typecounts(2)+1;
            elseif annotation_full(j).type=='P' % promoter mut
                typecounts(3)=typecounts(3)+1;
            elseif annotation_full(j).type=='I' % intergenic mut
                typecounts(4)=typecounts(4)+1;
            else
                error('unrecognized type')
            end
        elseif numel(new) > 1 % if more than one mutation
            fprintf(1,['Warning: Multiple mutations found at p=' num2str(p(goodpos(i))) '\n'])
            % Update mutation matrix
            for m=1:numel(new) % once for each mutation
                mutationmatrix(anc,new(m))=mutationmatrix(anc,new(m))+1;
            end
        end
    else
        fprintf(1,['Warning: Ancestor not found at p=' num2str(p(goodpos(i))) '\n'])
    end
end

% Count six types of mutations
mut_observed = div_matrix2_6types(mutationmatrix); % tally of mutations, NOT normalised

% Calculate fraction of nonsynonymous mutations (only N and S)
prob_nonsyn = typecounts(1)/(typecounts(1)+typecounts(2));


%% Save mutation spectrum information
if save_data
    save( [directory_save '/' clade_name '_mutation-spectrum.mat'], ...
        'mutationmatrix', 'mut_observed', 'typecounts', 'prob_nonsyn' ...
        ) 
end


end