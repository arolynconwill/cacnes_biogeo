function [ codon_counts_per_gene, gene_lengths, gene_cdsindex, codon_counts_raw, codon_counts_norm ] = codons_per_gene_in_genome( genomedirectory, allcodons, subset_gene_nums )

%% Summary

% Takes a reference genome and finds the abundance of each codon per gene
% Optional third argument allows user to specify genes of interest; if none
% provided, function considers all genes

% Note that the CDS index is the same as the "gene number" in annotation_full
% % [cds_indices,locus_tag_nums] = genomic_position_all(cds_struct,genomelength, ChrStarts)
% % genenum=genomic_position_all(CDS, GLength, ChrStarts); % gene numbers for each position on genome


%% Get genes from reference genome

% Get CDS
if ~exist([genomedirectory '/cds_sorted.mat'], 'file')
    read_gb_modified(genomedirectory)
end % makes cds_sorted.mat if it doesn't already exist
load([genomedirectory '/cds_sorted.mat']) 

% Get structure representing genes of interest
all_genes=[CDS{:}]; % copies CDS
if nargin == 3 % if genes of interest are provided
    % Take only specified gene numbers
    my_genes = all_genes( subset_gene_nums );
    gene_cdsindex = subset_gene_nums;
else % otherwise use all coding regions
    % Remove tRNAs and rRNAs
    my_genes = all_genes([all_genes.tagnumber]>0.5); % removes tRNAs and rRNAs
    gene_cdsindex = find([all_genes.tagnumber]>0.5);
end


%% Tally codon occurrences over all proteins in genome

codon_counts_raw = zeros( numel(allcodons),1 ); % initialize % tally number of each codon across all genes
codon_counts_per_gene = zeros( numel(allcodons), length(my_genes) ); % initialize % tally number of each codon across each gene
gene_lengths = zeros( length(my_genes),1 ); % initialize % keep track of length of each gene

for i=1:length(my_genes) 
    % Get gene sequence and codon distribution
    nextSequence = upper(my_genes(i).Sequence); % get next sequence
    nextCodonCounts = cell2mat(struct2cell( codoncount(nextSequence) )); % tally codons in sequence
    % Get gene length
    nextGeneLength = my_genes(i).loc2 - my_genes(i).loc1 + 1;
    % Save
    codon_counts_raw = codon_counts_raw + nextCodonCounts; % update tally of codons
    codon_counts_per_gene(:,i) = nextCodonCounts;
    gene_lengths(i) = nextGeneLength;
end

% Renormalize to get probability of each codon over all genes
codon_counts_norm = codon_counts_raw / sum(codon_counts_raw);


end
