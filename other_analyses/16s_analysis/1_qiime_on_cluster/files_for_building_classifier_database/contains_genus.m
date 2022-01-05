function contain_logical = contains_genus(specstring, genus_cats)
contain_logical_temp = zeros(size(genus_cats));
for ii = 1:length(genus_cats)
    contain_logical_temp(ii) = contains(specstring, genus_cats{ii});
end
contain_logical = contain_logical_temp>0;
end