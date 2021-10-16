function [ list_taxa, list_taxa_ids, list_fracs, list_reads ] = read_bracken_output( filename, min_frac )

% Read tab separated bracken file
my_data = tdfread(filename);

% Get components of data
list_taxa_all = ( my_data.name );
list_taxa_ids_all = ( my_data.taxonomy_id );
list_fracs_all = ( my_data.fraction_total_reads );
list_reads_all = ( my_data.new_est_reads );

% Filter only taxa that are >=1% of reads
taxa_to_keep = ( list_fracs_all >= min_frac );

% Return filtered lists
list_taxa = list_taxa_all( taxa_to_keep,: );
list_taxa_ids = list_taxa_ids_all( taxa_to_keep );
list_fracs = list_fracs_all( taxa_to_keep );
list_reads = list_reads_all( taxa_to_keep );

end