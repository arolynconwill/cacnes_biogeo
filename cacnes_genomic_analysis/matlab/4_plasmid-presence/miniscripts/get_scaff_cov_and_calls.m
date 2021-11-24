function [ Calls, coverage, p ] = get_scaff_cov_and_calls( filename, min_qual_for_call, min_maf_for_call, min_cov_each_strand_for_call )

% Load data
load( filename )

% Get coverage from counts
coverage = squeeze( sum( counts ) );
cov_fwd_strand=squeeze(sum(counts(1:4,:,:)));
cov_rev_strand=squeeze(sum(counts(5:8,:,:)));

% Compute major allele and major allele frequency
[maf, maNT, ~, ~] = div_major_allele_freq( counts );

% Switch sign of quals
Quals = -1*Quals;
% For manually setting Quals if VCF truncated and before cluster step re-runs
% if isequal( filename, 'case_step/case_A_440S-Ch-3_scaffold_34/candidate_mutation_table.mat' )
%     temp_index = find(ismember( SampleNames, 'A:415-Scr-trBa-4_t35' ));
%     Quals(:,temp_index) = 100;
%     fprintf(1,'WARNING! Quals for sample A:415-Scr-trBa-4_t35 scaffold A_440S-Ch-3_scaffold_34 all zero. Fixed manually. \n')
% end

% Get loose calls for each position
Calls = maNT;
Calls( Quals < min_qual_for_call ...
    | maf < min_maf_for_call ...
    | cov_fwd_strand < min_cov_each_strand_for_call ...
    | cov_rev_strand < min_cov_each_strand_for_call ...
    ) = 0; % set to N

end