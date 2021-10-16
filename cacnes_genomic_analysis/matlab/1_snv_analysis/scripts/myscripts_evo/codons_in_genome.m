function codon_counts = codons_in_genome( genomedirectory, allcodons, subset_gene_nums )

%% Summary

% Takes a reference genome and finds the abundance of all codons 
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
else % otherwise use all coding regions
    % Remove tRNAs and rRNAs
    my_genes = all_genes([all_genes.tagnumber]>0.5); % removes tRNAs and rRNAs
end


%% Tally codon occurrences over all proteins in genome

codon_counts = zeros( numel(allcodons),1 ); % for storing tally of codons

for i=1:length(my_genes)
    nextSequence = upper(my_genes(i).Sequence); % get next sequence
    nextCodonCounts = cell2mat(struct2cell( codoncount(nextSequence) )); % tally codons in sequence
    codon_counts = codon_counts + nextCodonCounts; % update tally of codons
end

% Renormalize to probability of each codon over whole genome
codon_counts = codon_counts / sum(codon_counts);


end
