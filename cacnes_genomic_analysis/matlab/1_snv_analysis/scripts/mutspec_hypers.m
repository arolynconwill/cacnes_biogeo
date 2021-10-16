%% Purpose

% Count the number of transition and transversion mutations

% Definitions:
% Transitions: ring structure stays the same
% A<->G or C<->T
% Transversions: change in ring structure
% everything else

% Want to evaluate: 
% % Transition/transversion ratio for hypermutators in Clade 1
% % Transition/transversion ratio for other in Clade 1
% % Transition/transversion ratio for all in Clade 1


%% Directory setup

% Directory for Lieberman Lab scripts:
dir_scripts_liebermanlab = '../lab_scripts';
path(dir_scripts_liebermanlab,path);

% Add my scripts
dir_my_scripts = [ pwd '/' 'scripts/myscripts_evo'];
path(dir_my_scripts,path);

% Where to find cluster data
dir_lineages = '2_snvs';

% Names of lineages
load( [ dir_lineages '/' 'cluster_names' ] )
num_lineages = numel(cluster_names);

% Where to save findings
dir_hyper = '5_parallel_evo/mutspec_hypers';
if ~exist(dir_hyper,'dir')
    mkdir(dir_hyper)
end


%% Options

save_figs = true;


%% Get data: whole clade

clade_of_interest = cluster_names{1};
mut_observed = load([ '5_parallel_evo/mutspec/' clade_of_interest '_mutation-spectrum.mat' ], 'mut_observed' );
mut_observed = mut_observed.mut_observed;

% mut_observed
%AT, TA % transversion
%AC, TG % transversion
%AG, TC % transition
%GC, CG % transversion
%GT, CA % transversion
%GA, CT % transition
% from mutation_spectrum_module via div_matrix2_6types 

transition_bool = [ 0 0 1 0 0 1 ];
transversion_bool = ~transition_bool;

transitions_all = sum( mut_observed'.*transition_bool );
transversions_all = sum( mut_observed'.*transversion_bool );
tt_ratio_all = transitions_all/transversions_all


%% Get data: calls for individual samples

load([dir_lineages '/' clade_of_interest '/data_' clade_of_interest '.mat'])

Calls_for_clade = Calls_for_analysis( :,~outgroup_isolates );

specimens_clade = specimen_numbers(~outgroup_isolates);

Calls_anc = anc_nti_goodpos;


%% Find positions corresponding to hypermutator SNPs

% Hypermutators are all samples belonging to specimen 31 or 41
hypermuts = ( specimens_clade == 31 | specimens_clade == 41 ); % bool
num_hypermuts = sum(hypermuts)

% Calls for just the hypermutators
Calls_hypermuts = Calls_for_clade( :,hypermuts );
if sum(sum( Calls_hypermuts == 0 )) > 0
    fprintf(1,'Warning! Ns exist in hypermutator calls. \n')
else
    fprintf(1,'No Ns in hypermutator calls. \n')
end

% Find positions where there are is variation within hypermutator set
nonvariablep_hypermuts = (sum(Calls_hypermuts==1,2)==num_hypermuts | sum(Calls_hypermuts==2,2)==num_hypermuts | sum(Calls_hypermuts==3,2)==num_hypermuts | sum(Calls_hypermuts==4,2)==num_hypermuts );
fprintf(1,['Number of hypermutator SNPs: ' num2str(sum(~nonvariablep_hypermuts)) '/' num2str(numel(nonvariablep_hypermuts)) '\n'])


%%

% Initialize counters
transition_tally_regular = 0;
transversion_tally_regular = 0;
transition_tally_hypermuts = 0;
transversion_tally_hypermuts = 0;

% Define transitions
transition_list = { 'AG', 'GA', 'CT', 'TC' };
% FMI: https://www.mun.ca/biology/scarr/Transitions_vs_Transversions.html

% Go through mutation positions and count...
for i=1:numel(annotation_full) % loop through annotation positions
    next_hypermut = ~nonvariablep_hypermuts(i);
    next_nts = annotation_full(i).nts;
    % Check to make sure there are only two NTs...
    if numel(next_nts) > 2
        fprintf(1,['More than two NTs at index ' num2str(i) ': ' next_nts '.\n'])
        continue
    end
    % Is this mutation a transition?
    next_transition = ismember( next_nts, transition_list );
    % Record in tally
    if next_hypermut % if a hypermutator SNP
        if next_transition
            transition_tally_hypermuts = transition_tally_hypermuts+1;
        else
            transversion_tally_hypermuts = transversion_tally_hypermuts+1;
        end
    else % if a regular SNP
        if next_transition
            transition_tally_regular = transition_tally_regular+1;
        else
            transversion_tally_regular = transversion_tally_regular+1;
        end
    end
end

% Check number of SNPs
tally_total = transition_tally_hypermuts+transversion_tally_hypermuts+transition_tally_regular+transversion_tally_regular;
tally_total_hypermuts = transition_tally_hypermuts+transversion_tally_hypermuts;
tally_total_regular = transition_tally_regular+transversion_tally_regular;
fprintf(1,['Total number of SNPs: ' num2str(tally_total) '\n'] )
fprintf(1,['Total number of SNPs (hypermutators only): ' num2str(tally_total_hypermuts) '\n'] )
fprintf(1,['Total number of SNPs (regular only): ' num2str(tally_total_regular) '\n'] )

% Compute ratios
ratio_hypermuts = transition_tally_hypermuts/transversion_tally_hypermuts;
ratio_regular = transition_tally_regular/transversion_tally_regular;
fprintf(1,['Hypermutator transition/transversion ratio: ' num2str(transition_tally_hypermuts) '/' num2str(transversion_tally_hypermuts) '=' num2str(ratio_hypermuts) '\n'])
fprintf(1,['Regular transition/transversion ratio: ' num2str(transition_tally_regular) '/' num2str(transversion_tally_regular) '=' num2str(ratio_regular) '\n'])
% And for overall
transition_tally_all = transition_tally_hypermuts+transition_tally_regular;
transversion_tally_all = transversion_tally_hypermuts+transversion_tally_regular;
ratio_all = transition_tally_all/transversion_tally_all;


%% Fisher exact test

fisher_table = table([3;1],[6;7],'VariableNames',{'Flu','NoFlu'},'RowNames',{'NoShot','Shot'})
fisher_table = table([transition_tally_hypermuts;transition_tally_regular],[transversion_tally_hypermuts;transversion_tally_regular],'VariableNames',{'Transition','Transversion'},'RowNames',{'Hypermutators','Others'})

%                      Transition    Transversion
%                      __________    ____________
% 
%     Hypermutators       119              4     
%     Others              360             62     

[h,p,stats] = fishertest(fisher_table,'Tail','right','Alpha',0.01)

% p =
% 
%    1.5254e-04
%    .00015
% 
% 
% stats = 
% 
%   struct with fields:
% 
%              OddsRatio: 5.1236
%     ConfidenceInterval: [1.3197 19.8916]
%     

%% Calculate uncertainty on ratios

% x = a/b

% hypermutators
x = ratio_hypermuts
a = transition_tally_hypermuts;
delta_a = sqrt(a);
b = transversion_tally_hypermuts;
delta_b = sqrt(b);
delta_x_hypermuts = x * sqrt( (delta_a/a)^2 + (delta_b/b)^2 )

% regular
x = ratio_regular
a = transition_tally_regular;
delta_a = sqrt(a);
b = transversion_tally_regular;
delta_b = sqrt(b);
delta_x_regular = x * sqrt( (delta_a/a)^2 + (delta_b/b)^2 )

% all
x = ratio_all
a = transition_tally_all;
delta_a = sqrt(a);
b = transversion_tally_all;
delta_b = sqrt(b);
delta_x_all = x * sqrt( (delta_a/a)^2 + (delta_b/b)^2 )


%% Plot ratios (all three)

% Data to plot
bar_data = [ ratio_hypermuts, ratio_regular, ratio_all ];
bar_names = { 'hypermutators', 'other', 'all' };
bar_error_high = [ delta_x_hypermuts, delta_x_regular, delta_x_all ];
bar_error_low = [ delta_x_hypermuts, delta_x_regular, delta_x_all ];

% Plot bar chart
figure(2)
clf(2)
hold on
box on
bar(1:numel(bar_data),bar_data,'BaseValue',1,'FaceColor',.9*[1 1 1],'LineWidth',2);
xticks(1:numel(bar_names))
xticklabels(bar_names)
yticks([1 10 100])
ylabel('transitions/transversions')
set(gca,'FontSize',25,'LineWidth',1.5)
set(gca,'YScale','log')
axis([ 0.5 3.5 1 100 ])
% Add error bars
er=errorbar(1:numel(bar_data),bar_data,bar_error_low,bar_error_high,'.');
er.Color = [0 0 0]; 
er.LineWidth = 2;
hold off

% Save figure
if save_figs
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([pwd '/' dir_hyper '/' 'Plot_Hyper_Hyper+Other+All.png'],'-dpng')
end


%% Plot ratios (first two only)

% Data to plot
bar_data = [ ratio_hypermuts, ratio_regular ];
bar_names = { 'hypermutator SNVs', 'other SNVs' };
bar_error_high = [ delta_x_hypermuts, delta_x_regular ];
bar_error_low = [ delta_x_hypermuts, delta_x_regular ];

% Plot bar chart
figure(1)
clf(1)
hold on
box on
bar(1:numel(bar_data),bar_data,'BaseValue',1,'FaceColor',.9*[1 1 1],'LineWidth',2);
xticks(1:numel(bar_names))
xticklabels(bar_names)
ylabel('transitions/transversions')
set(gca,'FontSize',25,'LineWidth',1.5)
set(gca,'YScale','log')
axis([ 0.5 2.5 1 100 ])
% Add error bars
er=errorbar(1:numel(bar_data),bar_data,bar_error_low,bar_error_high,'.');
er.Color = [0 0 0]; 
er.LineWidth = 2;
hold off

% Save figure
if save_figs
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [8 6]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 8 6]);
    print([pwd '/' dir_hyper '/' 'Plot_Hyper_Hyper+Other.png'],'-dpng')
end



