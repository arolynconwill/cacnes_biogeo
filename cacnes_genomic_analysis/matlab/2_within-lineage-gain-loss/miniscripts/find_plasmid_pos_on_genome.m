function plasmid_pos_on_genome = find_plasmid_pos_on_genome( path_plasmidblast, path_genome )


%% Summary

% This script identifies positions on the assembled genome that likely
% represent plasmid regions


%% Filter parameters

% Which alignments are "good"
min_aln_len = 200; % minimum alignment length
max_eval = -3; % max evalue
min_pid = 0; % min percent identity % not currently using

% Threshold for removing an entire contig
max_frac_contig_hit = 0.5;


%% Load genome info

% Genome stats
dir_genome = path_genome(1:end-12);
[chr_starts, genome_length, ~, ~] = genomestats(dir_genome);
all_pos=p2chrpos(1:1:genome_length,chr_starts);
all_pos_chr = all_pos(:,1); % chr for each position on genome
all_pos_index = 1:1:genome_length;

% From fasta
genome_info = fastaread( path_genome );
genome_contig_names = { genome_info.Header };
genome_contig_seqs = { genome_info.Sequence };
genome_contig_lengths = cellfun(@(x) numel(x), genome_contig_seqs );
num_contigs = numel(genome_contig_names);


%% Load blast data

% Read csv file
blast_table = readtable( path_plasmidblast );
% Fields: query id   subject id  % identity  alignment length    mismatches  gap opens   q. start    q. end  s. start    s. end  evalue  bit score

% Get column of interest
gen_contig = [ blast_table.Var2 ];
gen_pos_1 = [ blast_table.Var9 ];
gen_pos_2 = [ blast_table.Var10 ];
aln_len = [ blast_table.Var4 ];
aln_eval = log10( [ blast_table.Var11 ] );
inf_alt = -200;
aln_eval( aln_eval == -Inf ) = inf_alt; % -Inf annoying for plots
aln_pid = [ blast_table.Var3 ];


%% Make some plots to determine filters

% Histograms
%histogram(aln_len,1:100:max(aln_len))
%histogram(aln_eval,-100:1:max(aln_eval)); xlabel( log10(eval) )

% figure(1)
% clf(1)
% hold on
% box on
% scatter( log10(aln_len), aln_eval, 100 )
% xlabel('alignment length (bp, log10)')
% ylabel('evalue (log10)')
% yticks(inf_alt:10:0)
% yticklabels( [ -Inf, inf_alt+10:10:0 ] )
% title('alignments of plasmids to reference genome')
% set(gca,'FontSize',12)
% hold off
% print('scatter_alnlen-alneval','-dpng')
% 
% 
% figure(2)
% clf(2)
% hold on
% box on
% scatter( log10(aln_len), aln_pid, 100 )
% xlabel('alignment length (bp, log10)')
% ylabel('percent identity')
% title('alignments of plasmids to reference genome')
% set(gca,'FontSize',12)
% hold off
% print('scatter_alnlen-alnpid','-dpng')


%% Grab positions for alignments that meet the filter criteria

% Initialize
plasmid_pos_on_genome = [];

% Go through rows
for r=1:numel(gen_pos_1)
    
    % Get alignment info
    next_len = aln_len(r);
    next_eval = aln_eval(r);
    next_pid = aln_pid(r);
    
    % Filter
    if next_len >= min_aln_len && next_eval <= max_eval && next_pid > min_pid

        % Which contig
        contig_index = find( ismember( genome_contig_names, gen_contig{r} ) );

        % Get positions on contig
        p1 = gen_pos_1(r);
        p2 = gen_pos_2(r); 
        if p1<p2
            ps = p1:1:p2;
        else
            ps = p2:1:p1;
        end
        
        % Get contig and length of all contigs before
        prev_contigs_lengths_summed = sum( genome_contig_lengths( 1:contig_index-1 ) );
        
        % Record positions
        plasmid_pos_on_genome = [ plasmid_pos_on_genome, prev_contigs_lengths_summed+ps ];
    
    end

end

% Remove duplicates
plasmid_pos_on_genome = unique( plasmid_pos_on_genome );


%% Additionally remove all contigs where >=50% are plasmid hits

frac_contig_hit_all = zeros( num_contigs,1 );
for c=1:num_contigs
    this_chr_len = genome_contig_lengths(c); % length of contig
    this_chr_pos = all_pos_index(all_pos_chr==c); % all position indices on contig
    frac_contig_hit = sum(ismember( plasmid_pos_on_genome, this_chr_pos ))/this_chr_len;
    frac_contig_hit_all(c) = frac_contig_hit;
    if frac_contig_hit >= max_frac_contig_hit
        fprintf(1,[ 'Contig ' num2str(c) ': ' num2str( 100*frac_contig_hit ) '%% hit...' '\n' ])
        plasmid_pos_on_genome = [ plasmid_pos_on_genome, this_chr_pos ]; % add all contig pos to list
    end
end
%figure(2); histogram(frac_contig_hit_all,0:0.05:1)

% Remove duplicates
plasmid_pos_on_genome = unique( plasmid_pos_on_genome );


%% Report

% Print
fprintf(1,['Number of positions on reference genome that align to plasmid scaffolds: ' num2str(numel(plasmid_pos_on_genome)) '\n' ])


end
