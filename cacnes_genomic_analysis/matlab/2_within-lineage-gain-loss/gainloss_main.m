%% MOBILE ELEMENT GAIN/LOSS ANALYSIS


%% Summary

% This script uses a coverage matrix over a lineage pan-genome assembly in
% order to to identify candidate gain/loss regions.


%% Directory setup

% Main directory:
dir_main = char(pwd);
path( dir_main,path );
% Directory for Lieberman Lab scripts:
dir_lab_scripts = '../lab_scripts';
path(dir_lab_scripts,path);
% Directory for my scripts:
dir_scripts_aro = [dir_main '/miniscripts' ];
path(dir_scripts_aro,path);

% Where to find data from assemblies
dir_data_assemblies = '../data/data_assemblies/';

% Where to find data on lineage trees
dir_lineage_tree_info = 'data_lineage_snvtrees'; % just for making nice figures


%% Load information on lineages

% Lineage names
load( '../1_snv_analysis/2_snvs/cluster_names.mat' )


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFY GAIN/LOSS REGIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output directory
dir_gainloss_output = 'output_gainloss_regions';
if ~exist( dir_gainloss_output, 'dir' )
    mkdir( dir_gainloss_output )
end


%% Evaluate cluster-by-cluster

for this_clade_num = 1:numel(cluster_names)
    
    this_clade_name = cluster_names{this_clade_num};
    all_regions_table = gainloss_lineage( ...
        this_clade_num, this_clade_name, ...
        dir_data_assemblies, dir_lineage_tree_info, dir_gainloss_output );

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BLAST GAIN/LOSS REGIONS FOUND %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import data from blast

% Part I: Align gain/loss regions to C. ances reference genome

% 1. Run blast on the cluster to align regions to the C acnes reference genome
% 2. Put output csv's in this directory for analysis

% Directory for Pacnes_C1 reference genome blast
dir_blast_refgenome = 'blast_refgenome';

% Part II: Blast gain/loss regions against NCBI database

% Directory for blast against NCBI
dir_blast_ncbi = 'blast_ncbi';

% NCBI blast data processing
filename_blast_ncbi = 'region_blast_ncbi_tophit.csv';
table_blast_ncbi = readtable( [ dir_blast_ncbi '/' filename_blast_ncbi ], 'Delimiter', ',' );
table_blast_ncbi_names = table_blast_ncbi.Region;
table_blast_ncbi_tophit = table_blast_ncbi.TopHit_Species;
table_blast_ncbi_tophit_accession = table_blast_ncbi.TopHit_Accession;
table_blast_ncbi_querycover = table_blast_ncbi.Query_Cover; table_blast_ncbi_querycover = cell2mat( cellfun(@(x) str2double(x(1:end-1)), table_blast_ncbi_querycover, 'UniformOutput', false ) );
table_blast_ncbi_percentid = table_blast_ncbi.Percent_Identity; table_blast_ncbi_percentid = cell2mat( cellfun(@(x) str2double(x(1:end-1)), table_blast_ncbi_percentid, 'UniformOutput', false ) );
table_blast_ncbi_evalue = table_blast_ncbi.Evalue;
% Params
min_query_cover = 20; % percent


%%

%%%%%%%%%%%%%%%%%%%%%%
% SUMMARIZE FINDINGS %
%%%%%%%%%%%%%%%%%%%%%%

%% Directory
dir_supp_table = 'output_supp_table';
if ~exist( dir_supp_table, 'dir' )
    mkdir( dir_supp_table )
end

%% Concatenate tables from all clusters and add blast info

% Initialize
full_regions_table = struct;
num_reg_total=0;

% Loop through clusters and regions
for c=1:numel(cluster_names)
    
    load( [ dir_gainloss_output '/' 'C-' num2str(c) '_regions.mat' ] )
    if isempty( fieldnames( all_regions_table ) )
        continue
    else
        num_reg = numel({all_regions_table.Name});
    end
    
    % Get full list of fieldnames from cluster 1
    if c==1
        table_fields = fieldnames( all_regions_table ); % all_regions_table
        table_fields = table_fields(1:end-2); % remove blast placeholders
    end
    
    if num_reg>0
        for r=1:num_reg
            
            % Update region number
            num_reg_total = num_reg_total+1;
            % Region name
            next_region_name = all_regions_table(r).Name;
            
            % Get fields from cluster table
            for f=1:numel(table_fields)
                full_regions_table(num_reg_total).(table_fields{f}) = all_regions_table(r).(table_fields{f});
            end
            
            % Add new cluster name
            full_regions_table(num_reg_total).ClusterNew = cluster_names_new{c};
            
            % Get info on NCBI blast
            blast_ncbi_index = find( ismember( table_blast_ncbi_names, next_region_name ) );
            if table_blast_ncbi_querycover(blast_ncbi_index) >= min_query_cover
                full_regions_table(num_reg_total).BlastNCBI_BestHit = table_blast_ncbi_tophit{blast_ncbi_index};
                full_regions_table(num_reg_total).BlastNCBI_BestHitAcession = table_blast_ncbi_tophit_accession{blast_ncbi_index};
            else
                full_regions_table(num_reg_total).BlastNCBI_BestHit = 'no good hits';
                full_regions_table(num_reg_total).BlastNCBI_BestHitAcession = 'n/a';
            end
            
            % Get info on refgenome blast
            blast_refgenome_filename = [ 'blast_refgenome_' next_region_name '.csv' ];
            region_fasta_filename = [ next_region_name '.fasta' ];
            [ refgenome_perc_region_covered, refgenome_positions ] = get_refgenome_blast_info( [ dir_blast_refgenome '/' blast_refgenome_filename ], [ dir_gainloss_output '/' region_fasta_filename ] );
            full_regions_table(num_reg_total).BlastPacnesC1_PercentRegFound = refgenome_perc_region_covered;
            full_regions_table(num_reg_total).BlastPacnesC1_GenomePositions = refgenome_positions;
            
        end
    end
end

% Save
save( [ dir_gainloss_output '/' 'full_regions_table.mat' ], 'full_regions_table' )


%% Load sample name conversions

dir_snp_analysis = '../1_snv_analysis';
SampleNames_all = load([ dir_snp_analysis '/2_snvs/sample_names.mat' ],'SampleNames_all'); SampleNames_all = SampleNames_all.SampleNames_all;
SampleNamesSimple_all = load([ dir_snp_analysis '/2_snvs/sample_names.mat' ],'SampleNamesSimple_all'); SampleNamesSimple_all = SampleNamesSimple_all.SampleNamesSimple_all;


%% Make supplemental table (csv file)

% Open file
fid = fopen( [dir_supp_table '/supp_table_mobile_elements.csv'], 'w' );
% Headers 
fprintf(fid, 'Region name,Lineage name,Region length (bp),Region sequence,Present in reference genome,Coordinates in reference genome,Top blast hit (species),Top blast hit (accession number),Number positive colonies,Number negative colonies,Number ambiguous colonies,Names of positive colonies (regular),Names of positive colonies (simple),Names of negative colonies (regular),Names of negative colonies (simple),Names of ambiguous colonies (regular),Names of ambiguous colonies (simple), \n');
% One row per mobile element
for i=1:length( full_regions_table )
    % Region name
    fprintf(fid, full_regions_table(i).Name );
    fprintf(fid, ',' );
    % Lineage name
    fprintf(fid, full_regions_table(i).ClusterNew );
    fprintf(fid, ',' );
    % Region length (bp)
    fprintf(fid, num2str(full_regions_table(i).Region_Length) );
    fprintf(fid, ',' );
    % Region sequence
    fprintf(fid, full_regions_table(i).Sequence );
    fprintf(fid, ',' );
    % Present in reference genome
    if full_regions_table(i).BlastPacnesC1_PercentRegFound > 90
        fprintf(fid, 'yes' );
    else
        fprintf(fid, 'no' );
    end
    fprintf(fid, ',' );
    % Coordinates in reference genome
    next_coor = full_regions_table(i).BlastPacnesC1_GenomePositions;
    for j=1:size(next_coor,1)
        fprintf(fid, [ num2str(next_coor(j,1)) '-' num2str(next_coor(j,2)) ' ' ] );
    end
    fprintf(fid, ',' );
    % Top blast hit
    fprintf(fid, full_regions_table(i).BlastNCBI_BestHit );
    fprintf(fid, ',' );
    fprintf(fid, full_regions_table(i).BlastNCBI_BestHitAcession );
    fprintf(fid, ',' );
    % Number and names of colonies;
    fprintf( fid, num2str( full_regions_table(i).NumColonies_Pos ) ) ;
    fprintf(fid, ',' );
    fprintf( fid, num2str( full_regions_table(i).NumColonies_Neg ) ) ;
    fprintf(fid, ',' );
    fprintf( fid, num2str( full_regions_table(i).NumColonies_Ambig ) );
    fprintf(fid, ',' );
    next_names = full_regions_table(i).NamesColonies_Pos;
    for j=1:numel(next_names)
        fprintf(fid, [ next_names{j} ' ' ] );
    end
    fprintf(fid, ',' );
    for j=1:numel(next_names)
        fprintf(fid, [ SampleNamesSimple_all{ismember( SampleNames_all, next_names{j})} ' ' ] );
    end
    fprintf(fid, ',' );
    next_names = full_regions_table(i).NamesColonies_Neg;
    for j=1:numel(next_names)
        fprintf(fid, [ next_names{j} ' ' ] );
    end
    fprintf(fid, ',' );
    for j=1:numel(next_names)
        fprintf(fid, [ SampleNamesSimple_all{ismember( SampleNames_all, next_names{j})} ' ' ] );
    end
    fprintf(fid, ',' );
    next_names = full_regions_table(i).NamesColonies_Ambig;
    for j=1:numel(next_names)
        fprintf(fid, [ next_names{j} ' ' ] );
    end
    fprintf(fid, ',' );
    for j=1:numel(next_names)
        fprintf(fid, [ SampleNamesSimple_all{ismember( SampleNames_all, next_names{j})} ' ' ] );
    end
    fprintf(fid, ',' );
    fprintf(fid, '\n' );
end
fclose( fid );




