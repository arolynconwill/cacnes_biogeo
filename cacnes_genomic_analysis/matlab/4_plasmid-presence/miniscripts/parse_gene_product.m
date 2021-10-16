function product_str = parse_gene_product( product )

% concatenate if multiple cells
if size(product,1)>1
    product_str = []; % initialize
    for i=1:size(product,1)
        product_str = [ product_str, strtrim(product(i,:)) ];
        if i~=size(product,1)
            product_str = [ product_str, ' ' ];
        end
    end
else
    product_str = strtrim( product );
end

% remove commas
product_str = strrep( product_str, ',', '' );

% remove whitespace inside
while contains( product_str, '  ' )
    product_str = strrep( product_str, '  ', ' ' );
end

% avoid strings that start with -
if contains( product_str, '-PAC1_' )
    product_str = [ 'intergenic; ' product_str ];
end

end