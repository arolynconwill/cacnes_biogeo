function mut_annotations = annotate_mutations_gb_lineage(Positions,REFGENOMEFOLDER,p_lineage,calls_lineage)

%% Version History

% Original: 
% % Written and modified by Idan Yelin, Roy Kishony, Tami Lieberman, and Hattie Chung
% % May need revisions for your reference genome, especially subfunction read_gb_modified

% 2020-05-18, Arolyn and Jacob
% Old version takes gene sequences containing SNPs directly from a reference genome
% New version updates gene sequences from the reference genome to be closer to your set of isolates of interest. 
% % Two new inputs (lineage specific positions that vary from the reference genome and lineage specific nucleotide identities at those positions)
% Warning: using this feature is NOT advised for species with a lot of genome variability; but it can be used to make a closely related reference genome slightly more lineage specific.
% Other minor changes: 
% % Prints warning if two SNPs occur on the same codon.
% % Previously "p" was used to represent two different variables; changed one of them to p_list

% 2021-01-28, Arolyn
% % Fixed usage of nt2aa such that alternative start codons are allowed for
% the first codon in a gene!

% % Companion function: getNTs4annotations


%% Intialize

nts='atcg';
rc='tagc';

% Reference genome info
[ChrStarts, GLength, ~, ~]= genomestats(REFGENOMEFOLDER); % gets basic info about reference genome
if ~exist([REFGENOMEFOLDER '/cds_sorted.mat'], 'file')
    read_gb_modified(REFGENOMEFOLDER) % processes genbank file if not already done
end
load([REFGENOMEFOLDER '/cds_sorted.mat']) % Loads information from genbank file for all genes


%% Process lineage-specific information

% p_lineage = lineage specific positions to update in the reference genome
% calls_lineage = lineage specific nucleotide identities at those positions
Positions_lineage = p2chrpos( p_lineage,ChrStarts ); % get scaffolds and scaffold positions
calls_lineage_nts = nts( calls_lineage ); % turn into nucleotides
calls_lineage_nts_rcs = rc( calls_lineage ); % turn into complements for grabbing nts for reverse genes


%%  Fill in annotation for mutated positions

% Find which genes the positions of interest are part of
p_list=chrpos2index(Positions,ChrStarts);
genenum=genomic_position_all(CDS, GLength, ChrStarts); % gene numbers for each position on genome
mut_genenum=genenum(p_list); % gene number for each position of interest

% Record gene number and position on gene for all mutations so that we can check if any appear on the same codon later
list_of_genenums = zeros( size(Positions,1),1 );
list_of_genepos = zeros( size(Positions,1),1 );

% Create data structure with details on these genes and these positions in these genes
mut_annotations=[] ;
for i=1:size(Positions,1) % loops through positions of interest
    
    nGene = mut_genenum(i) ; % which gene 
    nScf = Positions(i,1) ; % which scaffold
    
    list_of_genenums(i) = nGene;
    
    mut_annotations(i).gene_num = nGene ;
    mut_annotations(i).scaffold = nScf ;
    mut_annotations(i).pos = Positions(i,2) ;
    
    if nGene==round(nGene) % intragenic
        cdf = CDS{nScf}(nGene) ; % get info about this gene from CDS
        
        if size(cdf.product,1)>1
            cdf.product=reshape(cdf.product',size(cdf.product,1)*size(cdf.product,2),1)';
        end
        
        if isfield(cdf,'gene')
            mut_annotations(i).gene       = cdf.gene ;
        else
            mut_annotations(i).gene = '';
        end
        
        if isfield(cdf,'product')
            mut_annotations(i).protein    = cdf.product ;
        else
            mut_annotations(i).protein = '';
        end
        
        mut_annotations(i).protein_id = cdf.protein_id ;
        mut_annotations(i).strand     = cdf.strand ; %1 indicates reverse strand, 0 forward strand
        mut_annotations(i).loc1       = cdf.loc1 ;
        mut_annotations(i).loc2       = cdf.loc2 ;
        mut_annotations(i).note       = cdf.note ;
        mut_annotations(i).locustag   = cdf.locustag ;
        if isfield(cdf,'oldlocustag')
            mut_annotations(i).oldlocustag   = cdf.oldlocustag ;
        end

        % Get gene sequence from CDS and update with lineage-specific information
        gene_sequence_reference = cdf.Sequence; % from genbank file of reference genome
        gene_sequence_lineage = gene_sequence_reference; % initialize
        % Figure out which positions need updating
        pos_to_update_bool = Positions_lineage(:,1) == nScf & Positions_lineage(:,2) >= cdf.loc1 & Positions_lineage(:,2) <= cdf.loc2;
        n_pos_to_update = sum(pos_to_update_bool);
        if n_pos_to_update > 0 % if there is anything to update
            pos_to_update_p = Positions_lineage(pos_to_update_bool,2);
            pos_to_update_nt = calls_lineage_nts(pos_to_update_bool);
            pos_to_update_nt_rc = calls_lineage_nts_rcs(pos_to_update_bool);
            for n=1:n_pos_to_update
                if ~mut_annotations(i).strand % fwd gene
                    gene_sequence_lineage( pos_to_update_p(n) - cdf.loc1 + 1 ) = pos_to_update_nt(n);
                else % rev gene
                    gene_sequence_lineage( cdf.loc2 - pos_to_update_p(n) + 1 ) = pos_to_update_nt_rc(n);
                end
            end
        end
        % Record lineage-specific gene sequence in mut_annotations
        mut_annotations(i).Sequence = gene_sequence_lineage;
        
        mut_annotations(i).translation = nt2aa(gene_sequence_lineage, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS', true, 'GENETICCODE', 11);
        
        if mut_annotations(i).strand
            p = double(mut_annotations(i).loc2) - double(Positions(i,2)) + 1;
        else
            p = double(Positions(i,2)) - double(mut_annotations(i).loc1) + 1;
        end
        mut_annotations(i).nt_pos = p ;
        
        list_of_genepos(i) = p;
        
        aan = floor((p-1)/3) + 1 ;
        ncn = p-(aan-1)*3 ;
        codons=cell(4,1);
        AA='';
        mut_annotations(i).aa_pos=aan;
        
        if numel(mut_annotations(i).Sequence) >= aan*3 & mut_annotations(i).translation > 1;
            codon = mut_annotations(i).Sequence(aan*3-2:aan*3) ;
            for j=1:4 %for each nts
                if mut_annotations(i).strand
                    codon(ncn)=rc(j);
                else
                    codon(ncn)=nts(j);
                end
                codons{j}=codon;
                if aan==1 % first amino acid
                    AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS', true, 'GENETICCODE', 11) ;
                else
                    AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS', false, 'GENETICCODE', 11) ;
                end
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from being called as a start codon all the time
            end
        end
        
        mut_annotations(i).codons   = codons ;
        mut_annotations(i).AA   = AA ;
        mut_annotations(i).NonSyn = length(unique(AA))>1 ;
    else %intergenic
        if floor(nGene)>=1 % gene before position
            cdf = CDS{nScf}(floor(nGene)) ;
            if isfield(cdf,'gene')
                mut_annotations(i).gene1       = cdf.gene ;
            else
                mut_annotations(i).gene1 = '';
            end
            mut_annotations(i).locustag1 = cdf.locustag ;
            mut_annotations(i).distance1 = Positions(i,2) - cdf.loc2;
            mut_annotations(i).protein1    = cdf.product ;
            if cdf.strand==1
                mut_annotations(i).distance1 = mut_annotations(i).distance1 * -1;
            end
        end
        if ceil(nGene)<=length(CDS{Positions(i,1)}) % gene after position
            cdf = CDS{nScf}(ceil(nGene)) ;
            if isfield(cdf,'gene')
                mut_annotations(i).gene2       = cdf.gene ;
            else
                mut_annotations(i).gene2 = '';
            end
            mut_annotations(i).locustag2 = cdf.locustag;
            mut_annotations(i).distance2 = cdf.loc1 - Positions(i,2);
            mut_annotations(i).protein2    = cdf.product ;
            if cdf.strand==0
                mut_annotations(i).distance2 = mut_annotations(i).distance2 * -1;
            end
        end
        mut_annotations(i).NonSyn = nan;
        mut_annotations(i).AA   = nan ;
    end
    
end

% Check if any SNPs occur on the same codon of the same gene
list_of_geneaan = floor((list_of_genepos-1)/3) + 1 ; % get codon number for each mutation
list_of_genenums_unique = unique(list_of_genenums);
for i=1:numel(list_of_genenums_unique) % check each gene
    same_gene_bool = ( list_of_genenums_unique(i) == list_of_genenums );
    this_gene_aans = list_of_geneaan( same_gene_bool ); % list of which codon each mutation corresponds to
    if numel(unique(this_gene_aans)) < numel(this_gene_aans)
        this_gene_aans_unique = unique( this_gene_aans );
        this_gene_aans_unique_counts = arrayfun(@(x) sum(x==this_gene_aans), this_gene_aans_unique );
        this_gene_aans_multiples = this_gene_aans_unique( this_gene_aans_unique_counts>1 );
        for j=1:numel(this_gene_aans_multiples) 
            fprintf(1,['Warning! Multiple SNPs on gene ' num2str(list_of_genenums_unique(i)) ' codon ' num2str(this_gene_aans_multiples(j)) '. \n'])
        end
    end
end

return
