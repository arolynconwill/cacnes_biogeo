function cdsindex_proteins = get_coding_genenums( dir_ref_genome, positions_to_mask )

% Get CDS
if ~exist([dir_ref_genome '/cds_sorted.mat'], 'file')
    read_gb_modified(dir_ref_genome)
end % makes cds_sorted.mat if it doesn't already exist
load([dir_ref_genome '/cds_sorted.mat']) 
all_genes=[CDS{:}]; % copies CDS

% Get CDS indices for proteins only
cdsindex_proteins = find([all_genes.tagnumber]>0.5); % excludes tRNA and rRNA

% Remove genes containing masked positions
if nargin>1
    gene_loc1 = [all_genes.loc1];
    gene_loc2 = [all_genes.loc2];
    cdsindex_mask = [];
    for p=1:numel(positions_to_mask)
        cds_index = find( positions_to_mask(p)>=gene_loc1 & positions_to_mask(p)<=gene_loc2 );
        if ~isempty( cds_index )
            cdsindex_mask = union( cdsindex_mask, cds_index );
        end
    end
    fprintf(1,['Masking ' num2str(numel(intersect(cdsindex_proteins,cdsindex_mask))) ' of ' num2str(numel(cdsindex_proteins)) ' proteins in genome.\n'  ])
    cdsindex_proteins = setdiff( cdsindex_proteins, cdsindex_mask );
end

end