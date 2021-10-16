function [ slope, slope_interval ] = fit_clock( times_list, dmrca_list, dir_lineage_spec, next_lineage, dataset_name, fig_num )

% Molecular clock
% % Linear fit
% % Makes a plot
% % Saves fit in text file


%% Perform fit

% Use "fit" to perform a best fit
[ fitobject, gof ] = fit( times_list, dmrca_list, 'Poly1' );
coeff_values = coeffvalues(fitobject);
coeff_intervals = confint(fitobject);
rsquared = gof.rsquare;
% Fitted values
fittedX = linspace(min(times_list), max(times_list), 100);
fittedY = polyval(coeff_values, fittedX);    

% Values to return
slope = 12*coeff_values(1); % snvs/mb/yr
slope_interval = 12*(coeff_intervals(2,1)-coeff_intervals(1,1))/2; % snvs/mb/yr

%% Make a plot showing data and best fit line

% Figure
figure(fig_num)
clf(fig_num)
hold on
% Scatter plot of raw data
color_dots = 0.5*[ 1 1 1 ];
jitter = 0.67*rand( numel(times_list),1 )-0.335;
scatter( times_list+jitter, dmrca_list , 100, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerEdgeAlpha',1,... % no edges
    'MarkerFaceColor', color_dots, ...
    'MarkerFaceAlpha',0.25,...
    'LineWidth', 1.0 ...
    )
% Appearance
box on
set(gca, 'FontSize', 20);
x_buffer = 2;
xlim( [ -x_buffer max(times_list)+x_buffer ] )
y_buffer = x_buffer;
ylim( [ 0 max(dmrca_list)+y_buffer ] )
xlabel('sampling time (months)')
if contains( dataset_name, '_subclades' )
    if contains( dataset_name, '-norm' )
        ylabel('SNVs per megabase from subclade ancestor')
    else
        ylabel('SNVs from subclade ancestor')
    end
else
    if contains( dataset_name, '-norm' )
        ylabel('SNVs per megabase')
    else
        ylabel('SNVs from lineage ancestor')
    end
end
% Add best fit line
plot(fittedX, fittedY, 'k-', 'LineWidth', 2);
yaxis_max = max(ylim);
hold off

print([ dir_lineage_spec '/' 'MolecClocl_' next_lineage '_' dataset_name '.png'],'-dpng')


%% Save fit info in text file

fid = fopen( [ dir_lineage_spec '/' 'MolecClocl_' next_lineage '_' dataset_name '.txt' ], 'w' );
if contains( dataset_name, '-norm' )
    fprintf(fid, ['dMRCA = ' num2str(12*coeff_values(1)) '*t +' num2str(coeff_values(2)) '\n']);
    fprintf(fid, ['c1 range (95%% CI): ' num2str(12*coeff_values(1)) ' +/- ' num2str(12*(coeff_intervals(2,1)-coeff_intervals(1,1))/2) ' SNVs/mb/yr' '\n']);
    fprintf(fid, ['c1 range (95%% CI): ' num2str(12*coeff_intervals(1,1)) ' to ' num2str(12*coeff_intervals(2,1)) ' SNVs/mb/yr' '\n']);
    fprintf(fid, ['dMRCA = ' num2str(coeff_values(1)) '*t +' num2str(coeff_values(2)) '\n']);
    fprintf(fid, ['c1 range (95%% CI): ' num2str(coeff_intervals(1,1)) ' to ' num2str(coeff_intervals(2,1)) ' SNVs/mb/mo' '\n']);
    fprintf(fid, ['c2 range (95%% CI): ' num2str(coeff_intervals(1,2)) ' to ' num2str(coeff_intervals(2,2)) ' SNVs/mb \n']);
    fprintf(fid, ['r2 : ' num2str(rsquared) '\n']);
else
    fprintf(fid, ['dMRCA = ' num2str(coeff_values(1)) '*t +' num2str(coeff_values(2)) '\n']);
    fprintf(fid, ['c1 range (95%% CI): ' num2str(coeff_intervals(1,1)) ' to ' num2str(coeff_intervals(2,1)) ' SNVs/mo' '\n']);
    fprintf(fid, ['c1 range (95%% CI): ' num2str(12*coeff_intervals(1,1)) ' to ' num2str(12*coeff_intervals(2,1)) ' SNVs/yr' '\n']);
    fprintf(fid, ['c2 range (95%% CI): ' num2str(coeff_intervals(1,2)) ' to ' num2str(coeff_intervals(2,2)) ' SNVs/mb \n']);
    fprintf(fid, ['r2 : ' num2str(rsquared) '\n']);
end
fclose(fid);

end