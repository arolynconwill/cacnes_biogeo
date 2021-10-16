function [ p_nonsyn, CI_nonsyn, num_muts_N, num_muts_S ] = compute_observed_dnds( annotation_full, gene_nums_of_interest )


%% Summary

% Computes the observed probability of a nonsynonymous mutation on a
% protein coding region based on all the mutations provided in
% annotation_full.

% Optional third input allows p_nonsyn to be calculated over a subset of
% genes.


%% Evaluate mutations in annotation_full

% Number of mutations
num_muts_total = length( annotation_full );

% Which mutations to consider
muts_genenums = [annotation_full.gene_num];
if nargin == 2 % consider specific genes if provided
    muts_filter = ismember( muts_genenums, gene_nums_of_interest );
else % otherwise consider all genes
    muts_filter = muts_genenums>0;
end

% Tally number of nonsynonymous vs synonymous mutations
muts_types = [annotation_full.type];
num_muts_N = sum( muts_types(muts_filter)=='N' );
num_muts_S = sum( muts_types(muts_filter)=='S' );


%% Compute probability of a nonsynonymous mutation with 95% CIs using binofit
% [ p, CI ] = binofit( x, n, alpha )
% p = max likelihood prob of succcess
% CI = confidence interval (100(1-alpha))
% x = number of successes
% n = number of trials
% alpha = to specify confidence interval (optional)

x = num_muts_N;
n = num_muts_N + num_muts_S;
alpha = 0.05;
[ p_nonsyn, CI_nonsyn ] = binofit( x, n, alpha );

    
end