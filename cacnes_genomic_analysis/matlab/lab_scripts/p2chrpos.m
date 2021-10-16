function pos = p2chrpos(p, cstarts)


%% 
% pos = p2chrpos(p, cstarts)
% Given the candidate SNP positions (p) in the chromosome and the positions
% of the chromosome beginnings (cstarts), returns a table with chromosome
% number and position for each entry in p: pos = [chr_num p]
%
% Edited 2021_03_24 by GTV added description
% Edited 2018_08_13 by TDL to automatically transpose p if in wrong direction

%%

if size(p,1)<size(p,2)
    p=p';
end

chr_num=ones(size(p));
if numel(cstarts) > 1
    for i=2:numel(cstarts)
        chr_num=chr_num+(p>cstarts(i));
    end
    positions= p - cstarts(chr_num)'; %Position on Chromosome
    pos=[chr_num positions];
else
    pos=[chr_num p];
end
end