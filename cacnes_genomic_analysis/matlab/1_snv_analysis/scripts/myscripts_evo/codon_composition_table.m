function table = codon_composition_table( allbases, allcodons )
% Makes a table of codon composition by base

table = zeros( numel(allcodons), numel(allbases) ); % Initialize table

for i=1:numel(allcodons) % loop through all codons
    for j=1:numel(allbases) % loop through all bases
        % Count the number of times a base occurs in a given codon
        % Then update the count in the appropriate position in the table
        table(i,j) = length( strfind( allcodons{i}, allbases{j} ) );
    end
end

end
