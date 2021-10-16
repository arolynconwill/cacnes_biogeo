%% FILTER CONTIGS FROM ASSEMBLIES FOR INPUT INTO CDHIT


%% Summary

% This script uses a coverage matrix over a genome to identify suspicious
% contigs in lineage assembled genomes. Genes on suspicious content are not
% included in between-lineage gene cluster analysis.


%% Directory setup

% Main directory:
dir_main = char(pwd);
path( dir_main,path );
% Directory for Lieberman Lab scripts:
dir_lab_scripts = '../../lab_scripts';
path(dir_lab_scripts,path);
% Directory for my scripts:
dir_scripts_aro = [dir_main '/miniscripts' ];
path(dir_scripts_aro,path);

% Where to find data from assemblies
dir_data_assemblies = '../../data/data_assemblies/';


%% Load information on lineages

% Lineage names
load( '../../1_snv_analysis/2_snvs/cluster_names.mat' )


%% Evaluate lineage-by-lineage

% Filters
min_contig_copynum_to_keep = 0.5; % keeps contigs with average copy number >= 0.5
min_aa_length_to_keep = 50; % removes plasmid contigs

% Initialize
all_contigs_to_keep = {};
all_assembled_genome_lenghts = [];
all_eff_genome_lenghts = [];
all_locustags_to_keep = {};

% Loop through lineages
for this_clade_num = 1:numel(cluster_names)
    
    this_clade_name = cluster_names{this_clade_num};
    [ contigs_to_keep, assembled_genome_length, eff_genome_length, locustags_to_keep ] = ...
        contig_filter_lineage( this_clade_num, this_clade_name, dir_data_assemblies, ...
        min_contig_copynum_to_keep, min_aa_length_to_keep);

    all_contigs_to_keep{end+1} = contigs_to_keep;
    all_assembled_genome_lenghts(end+1) = assembled_genome_length;
    all_eff_genome_lenghts(end+1) = eff_genome_length;
    all_locustags_to_keep{end+1} = locustags_to_keep;
    
end

%%

% Save data
dir_locustags = 'locustags';
if ~exist( [pwd '/' dir_locustags], 'dir' )
    mkdir( [pwd '/' dir_locustags] )
end
save( [ dir_locustags '/LineageGenomeInfo_All' '.mat' ], 'all_contigs_to_keep', 'all_assembled_genome_lenghts', 'all_eff_genome_lenghts', 'all_locustags_to_keep' );


%%

% Plot of genome length (histograms)
fig=figure(2);
clf(2)
hold on
box on
histogram( all_eff_genome_lenghts, 2200000:50000:3000000 )
histogram( all_assembled_genome_lenghts, 2200000:50000:3000000 )
line( [2519002 2519002], ylim, 'Color', 'k', 'LineWidth', 2 )
xlabel('genome length')
ylabel('num of lineages')
title(['lineage: all'])
text( 2519002*1.01, 0.8*max(ylim), 'Pacnes_C1', 'Interpreter', 'none' )
legend({'filtered','assembled'})
hold off

% Save plot
dir_plots = 'figures';
if ~exist( [pwd '/' dir_plots], 'dir' )
    mkdir( [pwd '/' dir_plots] )
end
print(fig,[ dir_plots '/GenomeLengths_all_hist.png' ],'-r400','-dpng')

% Plot of genome length (scatter)
figure(3)
clf(3)
hold on
box on
scatter(all_assembled_genome_lenghts,all_eff_genome_lenghts,100)
ylabel('assembled genome length')
xlabel('filtered genome length')
xlim( [ 2400000 2800000 ] )
ylim( [ 2400000 2800000 ] )
hold off

% Save plot
dir_plots = 'figures';
if ~exist( [pwd '/' dir_plots], 'dir' )
    mkdir( [pwd '/' dir_plots] )
end
print(fig,[ dir_plots '/GenomeLengths_all_scatter.png' ],'-r400','-dpng')


%% Filter fasta files of genes

% Fasta file from prokka
dir_fastas_prokka = 'input_fastas_prokka';
fid0 = fopen( [ dir_fastas_prokka '/' 'lineages_all.fasta' ], 'w' );

% Where to put filtered fasta
dir_fastas_filtered = 'output_fastas_filtered';
if ~exist( dir_fastas_filtered, 'dir' )
    mkdir( dir_fastas_filtered )
end

for this_clade = 1:numel(cluster_names)

    % Import original fasta file with amino acid sequennces
    fasta_og = fastaread( [ dir_fastas_prokka '/clade_' num2str(this_clade) '_prokka_out.faa' ] );
    fasta_og_headers = { fasta_og.Header };
    fasta_og_seqs = { fasta_og.Sequence };

    % Filter sequences
    fasta_og_headers_locustags = cellfun(@(x) strtok(x), fasta_og_headers, 'UniformOutput', false );
    entries_to_keep = find(ismember( fasta_og_headers_locustags, all_locustags_to_keep{this_clade} ));

    % Write new file
    fid = fopen( [ dir_fastas_filtered '/' 'lineage_' num2str(this_clade) '.fasta' ], 'w' );
    for i=1:numel(entries_to_keep)
        next_index = entries_to_keep(i);
        if this_clade<10
            this_clade_str = [ '0' num2str(this_clade) ];
        else
            this_clade_str = [ num2str(this_clade) ];
        end
        fprintf( fid, [ '>' 'L_' this_clade_str '_' fasta_og_headers_locustags{next_index} '\n' ] );
        fprintf( fid, [ fasta_og_seqs{next_index} '\n' ] );
        fprintf( fid0, [ '>' 'L_' this_clade_str '_' fasta_og_headers_locustags{next_index} '\n' ] );
        fprintf( fid0, [ fasta_og_seqs{next_index} '\n' ] );
    end
    fclose(fid);

end

fclose(fid0);

% Concatenate filtered fasta files from each lineage together
! cat output_fastas_filtered/lineage_*.fasta > lineages_all.fasta
