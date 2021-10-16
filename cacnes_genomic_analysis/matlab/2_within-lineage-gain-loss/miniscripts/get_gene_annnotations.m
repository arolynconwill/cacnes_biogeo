function [ genes_annotations, genes_translations ] = get_gene_annnotations( start_index, end_index, chr_starts, CDS )

% Don't need to check all positions for new genes
pos_chr = p2chrpos(start_index:1:end_index,chr_starts);
pos_skip_bp = 20; % check for new gene every xx bp in region
pos_chr_check = pos_chr(1:50:size(pos_chr,1),:);

% Initialize
genes_annotations = {};
genes_translations = {};
last_cds_contig_index = 0;
% Loop through positions to check
for p=1:size(pos_chr_check,1)
    contig_num = pos_chr_check(p,1);
    contig_pos = pos_chr_check(p,2);
    cds_contig = CDS{contig_num};
    if ~isempty( cds_contig ) % if contig has genes
        cds_contig_index = find( ( [cds_contig.loc1] <= contig_pos ) & ( contig_pos <= [cds_contig.loc2] ) );
        if ~isempty( cds_contig_index ) % if position is on a coding region
            if cds_contig_index~= last_cds_contig_index % if this coding region is a new gene that hasn't already been recorded
                last_cds_contig_index = cds_contig_index;
                % Record annotation
                product = cds_contig(cds_contig_index).product;
                if size(product,1)>1
                    product_str = []; % initialize
                    for i=1:size(product,1)
                        product_str = [ product_str, strtrim(product(i,:)) ];
                        if i~=size(product,1)
                            product_str = [ product_str, ' ' ];
                        end
                    end
                    genes_annotations{end+1} = strtrim( product_str );
                else
                    genes_annotations{end+1} = strtrim( cds_contig(cds_contig_index).product );
                end
                % Record translationn
                genes_translations{end+1} = cds_contig(cds_contig_index).translation2;
            end
        end
    end
end
% pad(product,size(product,1)+1,'right',' ')

% Remove duplicates
[ ~, unique_genes ] = unique( genes_translations );
genes_annotations = genes_annotations(unique_genes);
genes_translations = genes_translations(unique_genes);


end