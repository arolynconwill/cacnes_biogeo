function chi2_p = chi_squared_test( contingency_table )


%% Chi squared test

% Size of contingency table
table_size = size(contingency_table);
n_rows = table_size(1);
n_cols = table_size(2);

% Expected values
matrix_expected = zeros( n_rows,n_cols );
for r=1:n_rows
    for c=1:n_cols
        matrix_expected(r,c) = ...
            sum(contingency_table(r,:)) * sum(contingency_table(:,c)) / sum(sum(contingency_table));
    end
end

% Chi squared
chi2_val = sum( sum( ...
    ( contingency_table - matrix_expected ).^2 ./ matrix_expected ...
    ) );

% P value
chi2_dof = (n_rows-1)*(n_cols-1);
chi2_p = chi2cdf(chi2_val,chi2_dof,'upper');

end
