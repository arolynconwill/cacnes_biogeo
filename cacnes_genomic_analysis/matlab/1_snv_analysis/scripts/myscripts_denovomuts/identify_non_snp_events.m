function p_involved_in_non_snp_event = identify_non_snp_events( ...
    recombination_block_size, correlation_cutoff, goodpos, p, mutant_AF, ...
    GenomeLength, plot_bool, dir_diagnostic, this_cluster_name )

%% Notes

% Inputs:
% Mutant allele frequency (mutant_AF):
% % mutant_AF = mutant allele frequency
% % mutant_AF should ONLY include ingroup samples (remove outgroups!)
% % mutant_AF dimensions = (positions in goodpos) x (ingroup samples)
% Note: first give inputs are necessary; remaining inputs are optional and
% are used to make a plot to visualize where suspected recombination
% positions are on the genome.

% Filters: recommended starting points for C acnes lineages
% recombination_block_size=500;
% correlation_cutoff=0.75;

% Changes from Tami's version:
% Uses mutant_AF instead of maf (major allele frequency); not sure why
% original version used maf???
% Allow region to have 2 SNPs (previously required at least 3)
% Using correlation instead of covariance (better normalization and also
% symmetric)
% Makes a plot of positions suspected of recombination


%% Function

% Initialize 
involved_in_non_snp_event=[];

% Consider SNPs one at a time
for j=1:numel(goodpos)
    i=goodpos(j); % index in p of SNP
    % Find other indices within recombination_block_size bp on genome
    region=find(p(goodpos)>p(i)-recombination_block_size,1):find(p(goodpos)<p(i)+recombination_block_size,1,'last');
    %p(goodpos(region));
    % Consider pairs of SNPs within this set
    if numel(region)>=2
        % Correlation matrix: dimensions = # SNPs in region x # SNPs in region
        correlationmatrix = corrcoef(mutant_AF(region,:)');
        [a,b]= find( correlationmatrix > correlation_cutoff ); 
        involved_in_non_snp_event=[involved_in_non_snp_event region(a(a~=b)) ]; % remove self and require pairs
    end
end

% Remove duplicates
involved_in_non_snp_event=unique(involved_in_non_snp_event);

% Get positions on genome
p_involved_in_non_snp_event = p(goodpos(involved_in_non_snp_event));


%% Make a plot

if nargin > 5

    if plot_bool && numel(involved_in_non_snp_event)>0
    
        figure(101);
        clf(101);
        subplot(2,1,1)
        hold on
        histogram(p(goodpos(involved_in_non_snp_event)),0:10000:GenomeLength,'FaceColor',.25*[1 1 1]); xlabel('pos on genome'); 
        ylabel('num recombo SNPs')
        temp_p_plasmid = p( goodpos( p(goodpos) >= 1398569 & p(goodpos) <= 1418224 ) );
        temp_num_plasmid = sum( p(goodpos) >= 1398569 & p(goodpos) <= 1418224 );
        temp_p_plasmid_found = p( goodpos( involved_in_non_snp_event( p(goodpos(involved_in_non_snp_event)) >= 1398569 & p(goodpos(involved_in_non_snp_event)) <= 1418224 ) ) );
        temp_num_plasmid_found = sum( p(goodpos(involved_in_non_snp_event)) >= 1398569 & p(goodpos(involved_in_non_snp_event)) <= 1418224 );
        setdiff( temp_p_plasmid, temp_p_plasmid_found ) 
        temp_num_L11 = sum( p(goodpos) >= 1793873 & p(goodpos) <= 1793881 );
        temp_num_L11_found = sum( p(goodpos(involved_in_non_snp_event)) >= 1793873 & p(goodpos(involved_in_non_snp_event)) <= 1793881 );
        line([1415000,1415000,],[0,numel(involved_in_non_snp_event)],'Color','b'); text(1410000+1000,0.9*numel(involved_in_non_snp_event),['plasmid ' num2str(temp_num_plasmid_found) '/' num2str(temp_num_plasmid)],'Color','b')
        line([1793879,1793879,],[0,numel(involved_in_non_snp_event)],'Color','b'); text(1793879+1000,0.8*numel(involved_in_non_snp_event),['L11 ' num2str(temp_num_L11_found) '/' num2str(temp_num_L11)],'Color','b')
        ylim([0 numel(involved_in_non_snp_event)])
        set(gca,'FontSize',14)
        hold off
        subplot(2,1,2)
        hold on
        histogram(setdiff(p(goodpos),p(goodpos(involved_in_non_snp_event))),0:10000:GenomeLength,'FaceColor',.25*[1 1 1]); 
        xlabel('pos on genome'); 
        ylabel('num unfiltered SNPs')
        ylim([0 numel(involved_in_non_snp_event)])
        set(gca,'FontSize',14)
        hold off

        % Save
        print([ dir_diagnostic '/' 'Figure_RecombFilter_' this_cluster_name '.png'],'-dpng')

    end

end

end