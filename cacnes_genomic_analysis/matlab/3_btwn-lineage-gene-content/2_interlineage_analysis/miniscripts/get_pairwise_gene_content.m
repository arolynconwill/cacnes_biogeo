function [ pair_num_same, pair_num_diff, pair_frac_rel, pair_frac_tot ] = get_pairwise_gene_content( vec_1, vec_2)

pair_num_same = sum( vec_1 & vec_2 );
pair_num_diff = sum( vec_1 | vec_2 ) - pair_num_same;

pair_frac_rel = sum( vec_1 & vec_2 )/sum( vec_1 | vec_2 );
pair_frac_tot = pair_num_same/numel(vec_1);

end