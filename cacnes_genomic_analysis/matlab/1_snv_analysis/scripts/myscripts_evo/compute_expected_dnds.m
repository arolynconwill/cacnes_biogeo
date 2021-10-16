function probnonsyn = compute_expected_dnds( dir_ref_genome, mut_spec_prob, gene_nums_of_interest )

%% Summary

% This script examines the probability of a nonsynonymous mutation on a 
% reference genome given some mutational spectrum

% The third argument (gene_nums_of_interest) is optional. If it is
% provided, this function will return the neutral model dN/dS computed over
% the genes of interest only. Otherwise, this function will return the
% neutral model dN/dS compuated over all protein coding reigons.


%% Versions

% Arolyn, 2018.08.29: Original version
% Arolyn, 2020.07.06: Uses mutation spectrum from experiments
% % Incorporates observed mutation spectrum across all lineages
% Arolyn, 2021.01.20: Reference genome subset of genes
% % Option to only consider a subset of genes on the reference genome
% % Turned into a function


%% Define DNA bases and possible mutations

% All possible mutations to be considered:
allbases = { 'A', 'T', 'G', 'C' }; % four DNA bases
allmuts = { 'AT', 'AG', 'AC', 'TA', 'TG', 'TC', 'GA', 'GT', 'GC', 'CA', 'CT', 'CG' }; % 'AT' denotes A to T
mut_types_names={ 'AT/TA', 'AC/TG', 'AG/TC', 'GC/CG', 'GT/CA', 'GA/CT' };


%% Define codons

% All standard codons:
allcodons = fieldnames( geneticcode(11) );
allcodons = allcodons(2:end-1); % remove header and start codons

% Generate table of codon composition by base
codoncompositiontable = codon_composition_table( allbases, allcodons );
% Rows (64) = all possible codons
% Columns (4) = number of A's/T's/G's/C's in each codon


%% Generate table of probabilities of nonsynonymous mutations
% Rows (64) = all possible codons
% Columns (12) = all possible mutations
% Entries = probability of nonsynonymous mutation (given a mutation and a
% codon)
% Note: If the given mutation cannot be applied to a given codon, the entry
% is nan. 

% Generate table of probabilities that a given mutation on a codon is 
% nonsynonymous
codonnonsyntable = codon_mutation_table( allmuts, allcodons );


%% Calculate mutational spectrum

% Convert into format of allmuts
mutationalspectrum = zeros( numel(allmuts), 1 );
for i=1:numel(allmuts)
    next_mut = allmuts{i};
    next_index = find(cellfun(@(x) contains(x,next_mut), mut_types_names));
    mutationalspectrum(i) = mut_spec_prob(next_index)/2;
end


%% Calculate codon distribution in reference genome
% Import a reference genome and determine the codon distribution across all
% proteins

if nargin == 3 % if genes of interest are provided
    
    % Tally codons over subset of genome (gene specified in input)
    codondistribution = codons_in_genome( dir_ref_genome, allcodons, gene_nums_of_interest );

else % if not, consider all genes
    
    % Tally codons over the whole genome (proteins only)
    codondistribution = codons_in_genome( dir_ref_genome, allcodons );

end


%% Calculate probability of nonsynonymous mutation
% Takes into account: mutation spectrum, abundance of codons on genome,
% abundance of bases in codons, probability of a given mutation being
% nonsynonymous on a givne codon...

probnonsyn = 0; % probability of nonsynonymous mutations over all possible mutations

for thismutindex=1:length(allmuts) % loop through all possible mutations

    % The mutation being considered:
    thismutation=allmuts{thismutindex}; % ex. AT
    
    % Probability that this mutation occurs
    probthismutation=mutationalspectrum(thismutindex);
    
    % Find the codons that can undergo this mutation:
    thisbase=thismutation(1); % base that gets mutated; ex. A
    thisbaseindex=find(cellfun(@(x)~isempty(x), strfind(allbases,thisbase)),1); % ex. A is indexed in position 1
    thisbasecodonoccurrences=codoncompositiontable(:,thisbaseindex); % ex. how many A's in each codon
    % Indices of codons that have the relevant initial base:
    thismutcodons=(thisbasecodonoccurrences>0); % ex. AAT has A's but GGC does not
    
    % Probability that this mutation occurs on a given relevant codon
    % Take into account base composition of codons
    probmutoncodonbases=thisbasecodonoccurrences(thismutcodons);
    % Take into account codon abundance on reference genome
    probcodonongenome=codondistribution(thismutcodons);
    % Combine these two probabilities
    probmutoncodon=probmutoncodonbases.*probcodonongenome;
    % Renormalize (sum = 1 over all relevant codons)
    probmutoncodon=probmutoncodon/sum(probmutoncodon);
   
    % Probability that this mutation is nonsynonymous at each relevant
    % codon
    thismutnonsynoncodon=codonnonsyntable(:,thismutindex);
    probmutnonsynoncodon=thismutnonsynoncodon(thismutcodons);
    
    % Overall probability that this mutation is nonsynonymous over all
    % possible codons
    probmutnonsyn=probthismutation*sum(probmutoncodon.*probmutnonsynoncodon);

    % Add contribution of this mutation to the total probability of a nonsynonymous mutation
    probnonsyn = probnonsyn+probmutnonsyn;  
    
end

% Print the result to the console
%disp(['Probability of nonsynonymous mutation: ' num2str(probnonsyn)])

end