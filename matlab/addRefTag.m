function formattedBasicSeqs=addRefTag(basicSeqs)
%
%  takes sequence header info and creates friendly seq names
%      terminates at first whitespace
%      replaces underscores with hyphens
%      appends _Ref tag
%
nSeqs=length(basicSeqs);
formattedBasicSeqs=basicSeqs;  %  I always make a copy
%
for iSeq=1:nSeqs
    thisHeader=basicSeqs(iSeq).Header;
    %
    %  truncate at first whitespace char
    %
    thisHeaderSpaces=isspace(thisHeader);
    if any(thisHeaderSpaces)
        firstSpace=find(thisHeaderSpaces,1,'first');
        thisHeader(firstSpace:end)='';  % remove first space and rest of string
    end
    %
    %  replace underscores with hyphens
    %
    thisHeader=strrep(thisHeader,'_','-');
    %
    %  append reference sequence identifier
    %
    thisHeader=[thisHeader '_Ref'];
    %
    %  update header
    %
    formattedBasicSeqs(iSeq).Header=thisHeader;
end