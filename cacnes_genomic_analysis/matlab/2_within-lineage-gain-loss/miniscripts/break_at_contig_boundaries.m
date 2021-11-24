function [ reg_starts, reg_ends ] = break_at_contig_boundaries( reg_starts_0, reg_ends_0, ...
    contig_num_by_index, contig_start_indices, contig_end_indices, start_pos )

% Convert region starts/ends to genome indices
reg_starts_0_genome = reg_starts_0 + start_pos - 1;
reg_ends_0_genome = reg_ends_0 + start_pos - 1;

% Find indices of candidates covering multiple contigs
reg_has_multiple_contigs_bool = ( contig_num_by_index(reg_starts_0_genome) ~= contig_num_by_index(reg_ends_0_genome) );
reg_has_multiple_contigs = find( reg_has_multiple_contigs_bool );
num_reg_spanning_mult_contigs = sum( reg_has_multiple_contigs_bool );

% Initialize
num_original = numel( reg_starts_0_genome );
num_final = sum( contig_num_by_index(reg_ends_0_genome)-contig_num_by_index(reg_starts_0_genome)+1 );
% num_final = sum(~reg_has_multiple_contigs_bool) + sum( contig_num_by_index(reg_ends_0_genome(reg_has_multiple_contigs_bool))-contig_num_by_index(reg_starts_0_genome(reg_has_multiple_contigs_bool))+1 ); % alternative computation method
% num_final = num_original - num_reg_spanning_mult_contigs + sum( contig_num_by_index(reg_ends_0_genome(reg_has_multiple_contigs_bool))-contig_num_by_index(reg_starts_0_genome(reg_has_multiple_contigs_bool))+1 ); % alternative computation method
reg_starts_genome = zeros( 1,num_final ); % initialize
reg_ends_genome = zeros( 1,num_final ); % initialize
next_reg_index = 1;

if num_reg_spanning_mult_contigs > 0
    
    % First add all candidate regions that only lie on one contig
    num_to_add = sum( ~reg_has_multiple_contigs_bool );
    reg_starts_genome(1:num_to_add) = reg_starts_0_genome( ~reg_has_multiple_contigs_bool );
    reg_ends_genome(1:num_to_add) = reg_ends_0_genome( ~reg_has_multiple_contigs_bool );
    next_reg_index = next_reg_index + num_to_add;
%     num_to_add = del_has_multiple_contigs(1)-1;
%     reg_starts_genome(next_reg_index:next_reg_index+num_to_add-1) = reg_starts_0_genome( 1:1:num_to_add ); % initialize
%     reg_ends_genome(next_reg_index:next_reg_index+num_to_add-1) = reg_ends_0_genome( 1:1:num_to_add); % initialize
%     next_reg_index = next_reg_index + num_to_add;
    
    % Loop through candidate regions that span multiple contigs and split
    % them at contig boundaries
    for n=1:num_reg_spanning_mult_contigs
        first_contig = contig_num_by_index(reg_starts_0_genome(reg_has_multiple_contigs(n)));
        last_contig = contig_num_by_index(reg_ends_0_genome(reg_has_multiple_contigs(n)));
        for k=first_contig:1:last_contig
            if k==first_contig
                reg_starts_genome(next_reg_index) = reg_starts_0_genome(reg_has_multiple_contigs(n));
                reg_ends_genome(next_reg_index) = contig_end_indices(k);
                next_reg_index = next_reg_index+1;
            elseif k==last_contig
                reg_starts_genome(next_reg_index) = contig_start_indices(k);
                reg_ends_genome(next_reg_index) = reg_ends_0_genome(reg_has_multiple_contigs(n));
                next_reg_index = next_reg_index+1;
            else
                reg_starts_genome(next_reg_index) = contig_start_indices(k);
                reg_ends_genome(next_reg_index) = contig_end_indices(k);
                next_reg_index = next_reg_index+1;
            end
        end
%         % Fill in more regions until we hit another multi-contig one
%         if del_has_multiple_contigs(n)~=numel(reg_starts_0_genome) % not last candidate region
%             if n~=numel( del_has_multiple_contigs ) % not last candidate region with mutliple contigs
%                 if del_has_multiple_contigs(n+1) ~= del_has_multiple_contigs(n)+1 % next region does not contain multiple contigs
%                     num_to_add = del_has_multiple_contigs(n+1)-del_has_multiple_contigs(n)-1;
%                     reg_starts_genome(next_reg_index:next_reg_index+num_to_add-1) = reg_starts_0_genome( del_has_multiple_contigs(n)+1:1:del_has_multiple_contigs(n+1)-1 );
%                     reg_ends_genome(next_reg_index:next_reg_index+num_to_add-1) = reg_ends_0_genome( del_has_multiple_contigs(n)+1:1:del_has_multiple_contigs(n+1)-1 );
%                     next_reg_index = next_reg_index + num_to_add;
%                 end
%             else % no more regions with multiple contigs
%                 reg_starts_genome(next_reg_index:end) = reg_starts_0_genome( del_has_multiple_contigs(n):end );
%                 reg_ends_genome(next_reg_index:end) = reg_ends_0_genome( del_has_multiple_contigs(n):end );
%             end
%         end
    end
    
else % no candidate positions spanning multiple contigs
    
    reg_starts_genome = reg_starts_0_genome;
    reg_ends_genome = reg_ends_0_genome;
    
end

% Convert region stards/ends back to matrix indices
reg_starts = reg_starts_genome - start_pos + 1;
reg_ends = reg_ends_genome - start_pos + 1;

% Check that no regions start/end on different contigs
%sum( contig_num_by_index(reg_starts_genome) ~= contig_num_by_index(reg_ends_genome) ) % should be zero

end