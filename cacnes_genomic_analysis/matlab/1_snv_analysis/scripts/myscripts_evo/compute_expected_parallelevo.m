function [ sim_num_muts_per_gene, sim_num_mutsper1kb_per_gene, gene_lengths, gene_cdsindex, prob_gene_mutated ] = ...
    compute_expected_parallelevo( dir_ref_genome, mut_spec_prob, num_muts_to_sim, num_trials, gene_nums_of_interest )

%% Summary

% This function simulates random mutations on a genome based the provided
% mutation spectrum and codon composition of each gene.


%% Define DNA bases, codons, and possible mutations

% Define bases
allbases = { 'A', 'T', 'G', 'C' }; % four DNA bases

% Define all standard codons:
allcodons = fieldnames( geneticcode(11) );
allcodons = allcodons(2:end-1); % remove header and start codons
% Generate table of codon composition by base
codoncompositiontable = codon_composition_table( allbases, allcodons );
% Rows (64) = all possible codons
% Columns (4) = number of A's/T's/G's/C's in each codon

% Define possible mutations
allmuts = { 'AT', 'AG', 'AC', 'TA', 'TG', 'TC', 'GA', 'GT', 'GC', 'CA', 'CT', 'CG' }; % 'AT' denotes A to T
mut_types_names={ 'AT/TA', 'AC/TG', 'AG/TC', 'GC/CG', 'GT/CA', 'GA/CT' };


%% Calculate mutational spectrum

% Break down mutation spectrum according to allmuts
mutationalspectrum_types = zeros( numel(allmuts), 1 );
for i=1:numel(allmuts)
    next_mut = allmuts{i};
    mutationalspectrum_types(i) = mut_spec_prob( cellfun(@(x) contains(x,next_mut), mut_types_names) )/2;
end

% Then collapse to find probability that each base is mutated (regardless of what that base turns into)
mutationalspectrum_bases = zeros( numel(allbases), 1 );
for i=1:numel(allbases)
    next_base = allbases{i};
    next_mut_types = cellfun(@(x) x(1)==next_base, allmuts );
    mutationalspectrum_bases(i) = sum( mutationalspectrum_types( next_mut_types ) );
end

% Finally determine the probability that a codon is mutated given its base composition
mutationalspectrum_codons = sum( mutationalspectrum_bases.*codoncompositiontable' ); % not normalized
mutationalspectrum_codons = mutationalspectrum_codons/sum(mutationalspectrum_codons); % normalized


%% Calculate codon distribution in reference genome
% Import a reference genome and determine the codon distribution across all
% proteins

if nargin == 5 % if genes of interest are provided
    % Tally codons over subset of genome (gene specified in input)
    fprintf(1,'Codon counts for subset of genome...\n')
    [ codon_counts_per_gene, gene_lengths, gene_cdsindex, ~, ~ ] = codons_per_gene_in_genome( dir_ref_genome, allcodons, gene_nums_of_interest );
else % if not, consider all genes
    % Tally codons over the whole genome (proteins only)
    [ codon_counts_per_gene, gene_lengths, gene_cdsindex, ~, ~ ] = codons_per_gene_in_genome( dir_ref_genome, allcodons );
end

% Number of genes
num_genes = length(codon_counts_per_gene);


%% Compute probability that each gene is mutated
% given mutation probability for each codon and given codon content of each gene 

prob_gene_mutated = sum(mutationalspectrum_codons'.*codon_counts_per_gene); % not normalized
prob_gene_mutated = prob_gene_mutated/sum(prob_gene_mutated);

prob_gene_mutated_cumsum = [0, cumsum(prob_gene_mutated)];


%% Simulate expected number of mutations per gene

% Generate random numbers to simulate random mutations
my_rand_nums = rand( num_trials, num_muts_to_sim );

% Table to store the number of mutations on each gene in each trial
sim_num_muts_per_gene = zeros( num_trials, num_genes ); % initialize

for i=1:num_trials
    for j=1:num_genes
        sim_num_muts_per_gene(i,j) = sum( ...
            prob_gene_mutated_cumsum(j) <= my_rand_nums(i,:) & ...
            my_rand_nums(i,:) < prob_gene_mutated_cumsum(j+1) ...
            );
    end
end

gene_lengths = gene_lengths';
sim_num_mutsper1kb_per_gene = 1000*sim_num_muts_per_gene./gene_lengths;


end