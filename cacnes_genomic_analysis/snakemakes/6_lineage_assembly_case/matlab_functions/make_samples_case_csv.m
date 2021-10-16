%% OBJECTIVE

% Make samples_case_clade_N.csv for each clade

% Import samples.csv
% Provide directories on cluster here
% Output samples_case_clade_N.csv's


%% Import original samples.csv from assembly mapping step

% Get original samples.csv
csv_original = table2cell(readtable(['../5_lineage_assembly_mapping/samples.csv']));
original_paths = csv_original(:,1);
original_samples = csv_original(:,2);
original_refgenomes = csv_original(:,3);
original_providernames = csv_original(:,4);
original_subjects = csv_original(:,5);


%% Provide mapping directory

path_mapping = '/scratch/mit_lieberman/projects/aro_cacnes_biogeo/5_lineage_assembly_mapping';


%% Make table for samples_case.csv

% List of clades
list_clades = unique( cell2mat( original_subjects ));

% Make directory
mkdir( 'samples_case_csvs' )

% Loop through clades
for c=1:numel(list_clades)
    
    % Which clade
    this_clade = list_clades(c);

    % Make table to keep track of info for new samples.csv
    new_table = {};
    new_table{1,1} = 'Path';
    new_table{1,2} = 'Sample';
    new_table{1,3} = 'ReferenceGenome';
    new_table{1,4} = 'Outgroup';
    new_table_row = 2;

    % Pull info from samples.csv
    this_clade_size = sum( cell2mat(original_subjects) == this_clade );
    clade_samples = original_samples( cell2mat(original_subjects) == this_clade );
    clade_refgenomes = cell2mat(original_subjects( cell2mat(original_subjects) == this_clade ));

    % Put it into new_table
    for i=1:this_clade_size
        new_table{new_table_row,1} = path_mapping;
        new_table{new_table_row,2} = clade_samples{i};
        new_table{new_table_row,3} = clade_refgenomes(i);
        new_table{new_table_row,4} = 0;
        new_table_row = new_table_row+1;
    end

    % Write file
    fid = fopen(['samples_case_csvs/samples_case_clade_' num2str(this_clade) '.csv'],'w');
    for i=1:size(new_table,1)
        fprintf(fid,[new_table{i,1} ',' ...
            new_table{i,2} ',' ...
            num2str(new_table{i,3}) ',' ...
            num2str(new_table{i,4}) ',' '\n' ]);
    end
    fclose(fid);

end