function probability = prob_nonsyn_codon_mutation(codon,mut)
% Calculates the probability that a given mutation is nonsynonymous on a 
% given codon 

aa0=nt2aa(codon,'GeneticCode',11); % amino acid before mutation

% Find the positions on the codon at which mutation could occur
possiblemuts=find(codon==mut(1)); 

if isempty(possiblemuts) % if the mutation cannot occur on this codon
    probability = nan;
else % if the mutation can occur on this codon
    tally = 0; % initialize variable to count mutations that are nonsynonymous
    for i=1:numel(possiblemuts) % loop through possible mutations
        newcodon=codon; % save original codon
        newcodon(possiblemuts(i))=mut(2); % perform mutation
        aaf = nt2aa(newcodon,'GeneticCode',11); % amino acid after mutation
        if aa0 ~= aaf % if amino acid changed
            tally = tally + 1; % then add one to the tally
        end
    end
    probability = tally / numel(possiblemuts); % fraction of mutations that were nonsynonymous
end

end