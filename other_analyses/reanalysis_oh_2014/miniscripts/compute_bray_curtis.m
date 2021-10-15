function bray_curtis_val = compute_bray_cutis( read_counts_1, read_counts_2 )

% Use relative abundances
fracs_1 = read_counts_1/sum(read_counts_1);
fracs_2 = read_counts_2/sum(read_counts_2);
% Bray-Curtis
bray_curtis_val = sum(abs(fracs_1-fracs_2))/(sum(fracs_1)+sum(fracs_2));


% Absolute abundance (number of reads)
% counts_common = sum( min( read_counts_1, read_counts_2 ) );
% counts_total = sum( read_counts_1 ) + sum( read_counts_2 );
% bray_curtis_val=1-2*counts_common/counts_total;
% Alternative equation (same result)
% counts_diff = sum( abs( read_counts_1 - read_counts_2 ) );
% bray_curtis_val_2 = counts_diff/counts_total

end