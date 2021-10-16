function [ slst_all, clusters_all_slst ] = get_slst_from_assemblies( clusters_all, slst_all_original, slst_all )

% Gets SLST assignments from clade assemblies. 
% Compares to initial SLST assignments from reference genome alignment.


%% Load data from assemblies

% Get SLST info from assemblies
load('data_assemblies/slsts_from_assemblies')
% get_slst_cladenum
% get_slst_types


%% Compare SLST from assembly to SLST from referenge genome alignment
% Update any SLSTs that were initially missing

% Loop through clades and make sure SLSTs from identify_clusters agree with
% SLSTs from assemblies
clusters_all_slst = {};
for this_clade = 1:numel(clusters_all)
    fprintf(1,['Clade: ' num2str(this_clade) '\n'])
    % Get most abundant SLST
    this_clade_slsts = slst_all_original(clusters_all{this_clade});
    this_clade_slsts_unique = unique( this_clade_slsts );
    if numel(this_clade_slsts_unique) == 1
        this_clade_slst = this_clade_slsts_unique{1};
    else
        this_clade_slsts_counts = cellfun(@(x) sum(cellfun(@(y) isequal(x,y), this_clade_slsts)), this_clade_slsts_unique);
        [~,this_clade_slst_unique_index] = max( this_clade_slsts_counts );
        this_clade_slst = this_clade_slsts_unique{ this_clade_slst_unique_index };
    end
    fprintf(1,['SLST (ref genome): ' this_clade_slst '\n'])
    % Compare to SLST from assembly
    this_clade_slst_assembly = get_slst_types{this_clade};
    if numel( this_clade_slst_assembly ) == 1
        this_clade_slst_assembly = [ this_clade_slst_assembly 'x' ];
    end
    fprintf(1,['SLST (assembly): ' this_clade_slst_assembly '\n'])
    clusters_all_slst{end+1} = this_clade_slst_assembly;
    % Update slst_all in two cases:
    % % 1) if clade ref genome SLST is XX, update all according to assembly
    % (typical cause = SLST sequence with different length than reference
    % genome)
    % % 2) if only some samples are XX, update those specifically according
    % to assembly (typical cause = N in SLST sequence)
    this_clade_indices = clusters_all{this_clade};
    if this_clade_slst == 'XX'
        for i=1:numel(this_clade_indices)
            slst_all{this_clade_indices(i)} = this_clade_slst_assembly;
        end
    else
        for i=1:numel(this_clade_indices)
            if slst_all{this_clade_indices(i)} == 'XX'
                slst_all{this_clade_indices(i)} = this_clade_slst_assembly;
            end
        end        
    end
end


%% Report which SLSTs were found

fprintf(1,['All SLSTs found: '])
all_clades_slst_unique = unique( slst_all );
for i=1:numel(all_clades_slst_unique)
    fprintf(1,[ all_clades_slst_unique{i} ', ' ])
end
fprintf(1,'\n')

end