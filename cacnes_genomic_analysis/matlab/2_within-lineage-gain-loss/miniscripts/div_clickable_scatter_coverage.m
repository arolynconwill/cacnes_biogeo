function div_clickable_scatter_coverage( strains, starts, ends, isdeletion, ...
    SampleNames, coverage_matrix, copynumber_matrix, twice_normalized, ...
    dataset_name, genome_length, chr_starts, genome_contig_seqs, CDS, ...
    dir_save_figs )

color_del = rgb('Red');
color_amp = rgb('Blue');

% Initialize clickable figure
figure(1); clf(1); hold on; box on;

% Plot lines in background for contig boundaries
xlim( [ 0 genome_length+1 ] )
for i=1:numel( chr_starts )
    line( [chr_starts(i) chr_starts(i)], [0 numel(SampleNames)+1], 'Color', [ 0 0 0 0.2 ] )
end

% Plot line for each amp/del event
lw=10;
for i=1:numel(strains)
    if isdeletion(i)>0 % line segments for deletions
        plot([starts(i) ends(i)], [strains(i) strains(i)], 'Color', color_del, 'LineWidth',lw)
    else % line segments for amplifications
        plot([starts(i) ends(i)], [strains(i) strains(i)], 'Color', color_amp, 'LineWidth',lw)
    end
    
end

% Plot clickable dot at beginning of each amp/del event
dot_colors = ~isdeletion'.*color_amp + isdeletion'.*color_del;
p1=scatter(starts,strains,50,dot_colors,'filled');
set( p1, 'ButtonDownFcn', @clicked ); % make dots clickable
ylim([0 numel(SampleNames)+1])

% General plot appearance
fs=20;
set(gca,'FontSize',fs)
set(gca,'TickLabelInterpreter','none')
yticks(1:1:numel(SampleNames))
yticklabels( SampleNames );
ax = gca;
ax.YAxis.FontSize = 5+15*(200-numel(SampleNames))/200;
ylabel('sample name','FontSize',fs)
xlabel('position on genome')
title([ dataset_name ': candidate amplifications and deletions' ])
% Shift whitespace on sides of plot
% pos=get(gca,'position');  % retrieve the current values
% pos(1)=1.25*pos(1);
% set(gca,'position',pos);
% Legend 
text( 0.025*(max(xlim)-min(xlim)), 0.95*max(ylim), 'gain', 'Color', color_amp, 'FontSize', 20 )
text( 0.025*(max(xlim)-min(xlim)), 0.975*max(ylim), 'loss', 'Color', color_del, 'FontSize', 20 )

% Save
%print([dir_save_figs '/' dataset_name '_main.png'],'-dpng')


% Prepare and create pop-up figure upon click

x=starts;
y=strains;

xrange=max(ends)-min(starts);
yrange=numel(SampleNames);

    % Defines what happens upon clicking
    function clicked(src,event)
        
        % Find data point clicked
        ac=get(gca,'position') ;
        fc=get(gcf,'position') ;
        pc=get(gcf,'CurrentPoint') ;
        xl=xlim ;
        yl=ylim ;
        ax = [fc(3) * [ac(1), ac(1)+ac(3)]]  ;
        ay = [fc(4) * [ac(2), ac(2)+ac(4)]]  ;
        x0 = xl(1) + (xl(2)-xl(1))*(pc(1)-ax(1))/(ax(2)-ax(1)) ;
        y0 = yl(1) + (yl(2)-yl(1))*(pc(2)-ay(1))/(ay(2)-ay(1)) ;
        
        [~,ind] = min(((x-x0)/xrange).^2+((y-y0)/yrange).^2) ;
        
        % Get region positions
        start_index = starts(ind);
        start_contig_pos = p2chrpos(start_index,chr_starts);
        end_index = ends(ind);
        end_contig_pos = p2chrpos(end_index,chr_starts);
        % Get genes in region
        [ genes_annotations, genes_translations ] = get_gene_annnotations( start_index, end_index, chr_starts, CDS );
        
        % Print info to console
        fprintf(1,'\n\n') % add some headspace
        if(isdeletion(ind)>0)
            fprintf(1,[ 'Loss:' '\n' ])
        else
            fprintf(1,[ 'Gain:' '\n' ])
        end
        fprintf(1,[ '  Index: ' num2str(ind) '\n' ])
        fprintf(1,[ '  Sample: ' SampleNames{strains(ind)} '\n' ])
        fprintf(1,[ '  Start (index): ' num2str(starts(ind)) '\n' ])
        fprintf(1,[ '  Start (contig num): ' num2str(start_contig_pos(1)) '\n' ])
        fprintf(1,[ '  Start (contig pos): ' num2str(start_contig_pos(2)) '\n' ])
        fprintf(1,[ '  End (index): ' num2str(ends(ind)) '\n' ])
        fprintf(1,[ '  End (contig num): ' num2str(end_contig_pos(1)) '\n' ])
        fprintf(1,[ '  End (contig pos): ' num2str(end_contig_pos(2)) '\n' ])
        fprintf(1,[ '  Length: ' num2str(ends(ind)-starts(ind)+1) '\n' ])
        if start_contig_pos(1)==end_contig_pos(1) % start and end on same contig
            contig_seq = genome_contig_seqs{ start_contig_pos(1) };
            contig_seq_region = contig_seq( start_contig_pos(2):1:end_contig_pos(2) );
            fprintf(1,[ '  Sequence: ' contig_seq_region '\n' ])
        end
        if isempty( genes_annotations)
            fprintf(1,[ '  No gene annotations. ' '\n' ])
        else
            fprintf(1,[ '  Gene annotation(s): ' '\n' ])
            for g=1:numel(genes_annotations)
                fprintf(1,[ '    ' genes_annotations{g} '\n' ])
                fprintf(1,[ '      translation: ' genes_translations{g} '\n' ])
            end
        end
        
        
        % Coverage plot
        fig_num = 10;
        figure(fig_num); clf(fig_num); hold on;
        st=suptitle([ dataset_name ': candidate gain/loss region ' num2str(ind)]);
        st.FontSize=20;
        subplot_num = 3;
        % 1. Absolute coverage plot
        subplot_index = 1;
        make_popup_coverage_plot( subplot_num, subplot_index, coverage_matrix, 'raw coverage', SampleNames, strains, starts, ends, isdeletion, ind, chr_starts, color_amp, color_del )
        % 2. Copy number coverage plot
        subplot_index = 2;
        make_popup_coverage_plot( subplot_num, subplot_index, copynumber_matrix, 'copy number', SampleNames, strains, starts, ends, isdeletion, ind, chr_starts, color_amp, color_del )
        % 3. Normalized coverage plot
        subplot_index = 3;
        make_popup_coverage_plot( subplot_num, subplot_index, twice_normalized, 'normalized coverage', SampleNames, strains, starts, ends, isdeletion, ind, chr_starts, color_amp, color_del )

    end

end
