function [MAF, majorNT,minorNT, minorAF] = div_major_allele_freq(cnts)
%
% [MAF, majorNT,minorNT, minorAF] = div_major_allele_freq(cnts)
%
% Returns major and minor allele frequencies (MAF & MinorAF) and the
% nucleotide corresponding to those frequencies: majorNT and minorNT,
% respectively.
%
% Takes forward and reverse counts into account.

c=cnts(1:4,:,:)+cnts(5:8,:,:);

[sorted, sortedpositions] = sort(c,1);
maxcount = sorted(end,:,:);
minorcount = sorted(end-1,:,:);

MAF=double(maxcount)./sum(c,1);
minorAF=double(minorcount)./sum(c,1);

% Note: If counts for all bases are zero, then sort won't change the order
% (since there is nothing to sort), so the majorNT will end up as a 4 and
% the minorNT will end up as a 3.
majorNT = squeeze(sortedpositions(end,:,:));
minorNT = squeeze(sortedpositions(end-1,:,:));

MAF=squeeze(MAF);
MAF(isnan(MAF))=0; %set to 0 to indicate no data

minorAF=squeeze(minorAF); %leave no datat as 'nan' because 0 here could mean a pure position


return

