function annotation_full_everything = concatenate_annotation_fulls( paths_to_files, file_tags )

% Generate annotation_full_everything combining annotation_full across all lineages
for file_index=1:numel(paths_to_files)

    % Next file and file tag
    next_file = paths_to_files{file_index};
    next_file_tag = file_tags{file_index}; 
    fprintf(1,['Loading file ' next_file '...\n'])
    
    % Check if clade has enough SNPs for annotations; if so, load and append
    load( next_file, 'goodpos' );
    if numel(goodpos) > 0 % if within lineage SNPs exist
        
        % Load annotation_full
        load(next_file, 'annotation_full' );
        % Get other information about this set of mutations
        load(next_file, 'outgroup_isolates');
        load(next_file, 'subjects' );
        next_subject = char( 64+ unique(subjects(~outgroup_isolates)) );
        load(next_file, 'slst' );
        super_slst = cellfun(@(x) x(1), slst(~outgroup_isolates));
        next_superslst = unique(super_slst);
     
        % Add rows to annotation_full_everything
        if file_index==1
            % Intialize annotation_full_everything
            annotation_full_everything = struct;
            % Add rows to annotation_full_everything
            num_snps = 0;
            fields_next = fieldnames(annotation_full); % get fields here because they could be a subset...
            for snp_index=1:length(annotation_full) % loop through eacn new SNP
                for field_index=1:length(fields_next) % update each field
                    annotation_full_everything(num_snps+snp_index).(fields_next{field_index}) = annotation_full(snp_index).(fields_next{field_index});
                end
                % Add new fields
                annotation_full_everything(num_snps+snp_index).data_origin = next_file;
                annotation_full_everything(num_snps+snp_index).data_index = file_index;
                annotation_full_everything(num_snps+snp_index).data_tag = next_file_tag;
                annotation_full_everything(num_snps+snp_index).subject = next_subject;
                annotation_full_everything(num_snps+snp_index).superSLST = next_superslst;
            end
        else
            % Add rows to annotation_full_everything
            num_snps = length(annotation_full_everything);
            fields_next = fieldnames(annotation_full); % get fields here because they could be a subset...
            for snp_index=1:length(annotation_full) % loop through eacn new SNP
                for field_index=1:length(fields_next) % update each field
                    annotation_full_everything(num_snps+snp_index).(fields_next{field_index}) = annotation_full(snp_index).(fields_next{field_index});
                end
                % Add new fields
                annotation_full_everything(num_snps+snp_index).data_origin = next_file;
                annotation_full_everything(num_snps+snp_index).data_index = file_index;
                annotation_full_everything(num_snps+snp_index).data_tag = next_file_tag;
                annotation_full_everything(num_snps+snp_index).subject = next_subject;
                annotation_full_everything(num_snps+snp_index).superSLST = next_superslst;
            end
        end
        
    end

end