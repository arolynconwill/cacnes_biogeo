function features = parse_location_gb(s)

s(s=='<') = [] ;

if s(1)=='c'       % if 'true', gene is found on complement strand. ('features.strand = false'). 
    loc = s(12:end-1) ;     % used to skip 'complement('
    features.strand=false ;         
else
    loc = s ;
    features.strand=true ;
end

if loc(1)=='j'       % if 'true', gene overlaps origin
    loc = loc(6:end-1) ;     % used to skip 'complement('
end


isnum=loc>=48 & loc<58;
loc1end=find(~isnum,1);
loc2start=find(~isnum,1,'last');

features.loc1 = str2num(loc(1:loc1end-1)) ;
features.loc2 = str2num(loc(loc2start+1:end)) ;

return
end
