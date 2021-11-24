function make_summary_coverage_plot( subplot_num, subplot_index, num_samples, ...
    bool_samples_pos, bool_samples_neg, bool_samples_ambig, ...
    coverage_matrix, matrix_type, ...
    start_genome_index, end_genome_index, chr_starts, genome_length, ...
    color_pos, color_neg, color_ambig )

% Absolute coverage plot
subplot(subplot_num,1,subplot_index)
hold on
% Plot parameters
extraonsides=2500; % bp
plotstart=max(1,start_genome_index-extraonsides);
plotend=min(genome_length,end_genome_index+extraonsides);
% Plot appearance
box on
xlim([plotstart plotend])
xlabel('position on genome')
if isequal( matrix_type, 'raw coverage')
    ylim([ 0 1.1*max(max(coverage_matrix(:,plotstart:plotend))) ])
elseif isequal( matrix_type, 'copy number' )
    ylim([ 0 1.1*max(max(coverage_matrix(:,plotstart:plotend))) ])
elseif isequal( matrix_type, 'normalized coverage' )
    abs_max = max(max(abs(coverage_matrix(:,plotstart:plotend))));
    ylim([ -1.1*abs_max 1.1*abs_max ])
else
    fprintf(1,'Warning! Unrecognized ylabel; ylim may not be optimized.\n')
end
yl=ylabel(matrix_type);
yl.Position(1)=plotstart-(plotend-plotstart)/15; % ylabel always in same relative position to axes
set(gca,'FontSize',16)
% Plot lines in background for contig boundaries
chr_starts_to_plot = chr_starts( chr_starts >= plotstart & chr_starts <= plotend );
for i=1:numel( chr_starts_to_plot )
    line( [chr_starts_to_plot(i) chr_starts_to_plot(i)], ylim, 'Color', [ rgb('Purple') 0.5 ] )
end
% Plots raw coverage of each strain individually
line_alpha = min(max(5/num_samples,0.1),1);
line_width = 1;
for j=1:num_samples
    if bool_samples_pos(j)
        plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',line_width, 'Color', [ color_pos,line_alpha ] );
    elseif bool_samples_neg(j)
        plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',line_width, 'Color', [ color_neg,line_alpha ] );
    else
        plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',line_width, 'Color', [ color_ambig,line_alpha ] );
    end
end
% This plots a the average coverage of all samples
plot(plotstart:plotend,mean(coverage_matrix(:,plotstart:plotend)),'LineWidth',2.5, 'Color', rgb('Black'))
% This adds vertical lines on the region boundaries
temp_ylim=ylim;
plot([start_genome_index start_genome_index], [temp_ylim(1) temp_ylim(2)], '--', 'Color',rgb('Green'), 'LineWidth',1.5 );
plot([end_genome_index end_genome_index],[temp_ylim(1) temp_ylim(2)], '--', 'Color',rgb('Green'), 'LineWidth',1.5 );
        
end