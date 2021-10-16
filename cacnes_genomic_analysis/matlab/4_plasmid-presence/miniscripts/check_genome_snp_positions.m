function check_genome_snp_positions( this_cluster_name, filename_blast_alignments )

% Load SNP data for this cluster
fprintf(1,[this_cluster_name '\n'])
load(['/Users/Aro/Dropbox (MIT)/Postdoc/Pacnes/Analysis_SNP/analysis_evo_v9_2020-10/clusters/' this_cluster_name '/data_' this_cluster_name '.mat']);

% Load plasmid positions on reference genome
blast_table = readtable( filename_blast_alignments );
c1_pos_1 = [ blast_table.Var9 ];
c1_pos_2 = [ blast_table.Var10 ];
plasmid_pos_on_genome = [];
for r=1:numel(c1_pos_1)
    p1 = c1_pos_1(r);
    p2 = c1_pos_2(r); 
    if p1<p2
        ps = p1:1:p2;
    else
        ps = p2:1:p1;
    end
    plasmid_pos_on_genome = union( plasmid_pos_on_genome, ps );
end
fprintf(1,['Number of positions on reference genome that align to plasmid scaffolds: ' num2str(numel(plasmid_pos_on_genome)) '\n' ])

% Save plasmid positions on genome
i_slash = find(ismember( filename_blast_alignments, '/' ));
my_dir = filename_blast_alignments(1:1:i_slash(end));
save([my_dir 'plasmid_pos_on_C1_genome'],'plasmid_pos_on_genome')

% Compare to regions already masked
pos_already_masked = [ 1398569:1:1418224, 1793873:1:1793881 ];
pos_overlap = intersect( pos_already_masked, plasmid_pos_on_genome );
fprintf(1,['Number of these positions on reference genome already masked: ' num2str(numel(pos_overlap)) '\n' ])

% Count how many SNPs in this clade overlap with these positions
snps_on_plasmid = intersect( p(goodpos), plasmid_pos_on_genome);
fprintf(1,['Number of clade SNPs on potential plasmid regions: ' num2str(numel(snps_on_plasmid)) '\n' ])


end