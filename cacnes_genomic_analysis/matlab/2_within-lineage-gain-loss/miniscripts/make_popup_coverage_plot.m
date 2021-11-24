function make_popup_coverage_plot( subplot_num, subplot_index, ...
    coverage_matrix, matrix_type, ...
    SampleNames, strains, starts, ends, isdeletion, ind, ...
    chr_starts, color_amp, color_del )

% Absolute coverage plot
subplot(subplot_num,1,subplot_index)
hold on
% Plot parameters
extraonsides=2500;
plotstart=max(1,starts(ind)-extraonsides);
plotend=min(ends(ind)+extraonsides,size(coverage_matrix,2));
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
for j=1:numel(SampleNames)
    plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',.1, 'Color', [ 0 0 0 max(1/numel(SampleNames),0.05) ] ) 
%    plot(plotstart:plotend,coverage_matrix(j,plotstart:plotend), 'LineWidth',.1,'Color', (.5*j/numel(SampleNames))*[1 1 1]) 
end
% This plots a the average coverage of all samples
plot(plotstart:plotend,mean(coverage_matrix(:,plotstart:plotend)),'LineWidth',5, 'Color', rgb('Black'))
% This highlights the sample you clicked on
if(isdeletion(ind)>0)
    plot(plotstart:plotend,coverage_matrix(strains(ind),plotstart:plotend), 'LineWidth',5, 'Color', color_del)
else
    plot(plotstart:plotend,coverage_matrix(strains(ind),plotstart:plotend), 'LineWidth',5, 'Color', color_amp)
end
% This adds vertical lines on the region boundaries
temp=ylim;
plot([starts(ind) starts(ind)], [temp(1) temp(2)], '--', 'Color',rgb('Green'), 'LineWidth',1.5 );
plot([ends(ind) ends(ind)],[temp(1) temp(2)], '--', 'Color',rgb('Green'), 'LineWidth',1.5 );
        
end