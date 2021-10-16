function div_bar_charts(c, names, position, avg_cov, showall)
%showall is an optional argument to display more of counts matrix

%% Changes
% % Arolyn, 2019.01.17: Displays position on genome on bar chart; position
% added as a new input

% Tami, 2019.01.23 ylim on bar chart (not always the best on idea)
% Velina, 2020.11.30: Displays 0.5*median of non-zero coverage for each
% sample along with bars; avg_cov summary values passed in as 4th argument
%%

LABEL_SIZE= 14; 


figure(660); clf; hold on;

if nargin < 5
    showall=0;
    num_call = subplot(1,1,1);
else
    num_call = subplot(4,1,1);
end

% fix this argument mess later -- velina
if nargin < 4
    show_avg_cov = 0;
else
    show_avg_cov = 1;
end


%Plot number of calls
a=squeeze(c(1:8,:));
bar(reshape([a; nan(4,size(a,2))],4,[])','stacked')
legend('A','T','C','G', 'Location', 'BestOutside')
ylabel('Number of reads', 'FontSize', 14)
xlabel('Samples', 'FontSize', 14) % GTV added this label
text(size(a,2)*0.1, cast(max(max(a))*0.95,'double'), ['Position on genome: ' num2str(position)])
set(num_call,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
set(num_call,'Xtick',2:3:(3*numel(names)-1))
xlim([0 (3*numel(names)+3)])
%ylim([0 100]) #TDL commented this out March 2020

if show_avg_cov > 0
    hold on
    xt = xticks;
    plot(xt,0.5*avg_cov, 'k.', 'MarkerSize', 12) %take half to account for forward/reverse being plotted separately
    legend('A','T','C','G', '"avg" cov', 'Location', 'BestOutside')
    hold off
end


%Plot other options
if showall > 0
    %Plot call quality
    call_qual = subplot(4,1,2); hold on;
    title('Average call quality');
    a=squeeze(c(9:16,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_bq params.min_bq], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
    legend('Aq','Tq','Cq','Gq', 'Location', 'BestOutside')
    ylabel('Average Base Quality')
    set(call_qual,'Xtick',2:3:(3*numel(names)-1))
    set(call_qual,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(call_qual, 'XTick', []);
    
    
    %Plot mapping quality
    map_qual = subplot(4,1,3); hold on;
    title('Average mapping quality');
    a=squeeze(c(17:24,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_mq params.min_mq], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped','LineStyle', 'none')
    legend('Am','Tm','Cm','Gm', 'Location', 'BestOutside')
    ylabel('Average Mapping Quality')
    set(map_qual,'Xtick',2:3:(3*numel(names)-1))
    set(map_qual,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(map_qual, 'XTick', []);
    
    
    %Plot tail distance f
    tail_dist = subplot(4,1,4); hold on;
    title('Average tail distance');
    a=squeeze(c(25:32,:));
    if nargin>3
        plot([0 3*size(a,2)], [params.min_td params.min_td], 'k:')
        plot([0 3*size(a,2)], [params.max_td params.max_td], 'k:')
    end
    bar(reshape([a; nan(4,size(a,2))],4,[])','grouped', 'LineStyle', 'none')
    legend('Atd','Ttd','Ctd','Gtd', 'Location', 'BestOutside')
    ylabel('Average Tail Distance')
    set(tail_dist,'Xtick',2:3:(3*numel(names)-1))
    set(tail_dist,'Xticklabel',names, 'FontSize', LABEL_SIZE, 'XTickLabelRotation', 90)
    xlim([0 (3*numel(names)+3)])
    set(tail_dist, 'XTick', []);
    
end

end


