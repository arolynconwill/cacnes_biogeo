function table = codon_mutation_table( allmuts, allcodons )
% Generates table of probabilities that a given mutation on a codon is
% nonsynonymous

table = zeros( numel(allcodons), numel(allmuts) ); % Initialize table

for i=1:numel(allcodons) % loops through all possible codons
    for j=1:numel(allmuts) % loops through all possible mutations
        % Calculates the probability that a given mutation is nonsynonymous
        % on a given codon and then updates the table
        table(i,j) = prob_nonsyn_codon_mutation( allcodons{i}, allmuts{j} );
    end
end

end
