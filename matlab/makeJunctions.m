function junctionSeqs=makeJunctions(publicSeqs,privateSeqs,readLength)
%
%  creates junctional sequences using readLength-1 long tail to head
%  combinations
%
nPublicSeqs=length(publicSeqs);
nPrivateSeqs=length(privateSeqs);
nTotalSeqs=nPublicSeqs+nPrivateSeqs;
totalSeqs=[publicSeqs(:);privateSeqs(:)];

htLength=readLength-1;
heads=char(zeros(nTotalSeqs,htLength,'uint8'));
tails=char(zeros(nTotalSeqs,htLength,'uint8'));
for iSeq=1:nTotalSeqs
    heads(iSeq,:)=totalSeqs(iSeq).Sequence(1:htLength);
    tails(iSeq,:)=totalSeqs(iSeq).Sequence((end-htLength)+(1:htLength));
end
iJunc=0;
for hSeq=1:nTotalSeqs
    thisHeadSeq=heads(hSeq,:);
    thisHeadName=totalSeqs(hSeq).Header;
    %
    %  strip any _ tag identifier
    %
    underLocs=strfind(thisHeadName,'_');
    if ~isempty(underLocs)
        firstUnder=underLocs(1);
        thisHeadName(firstUnder:end)='';
    end
    for tSeq=1:nTotalSeqs
        thisTailSeq=tails(tSeq,:);
        thisTailName=totalSeqs(tSeq).Header;
        %
        %  strip any _ tag identifier
        %
        underLocs=strfind(thisTailName,'_');
        if ~isempty(underLocs)
            firstUnder=underLocs(1);
            thisTailName(firstUnder:end)='';
        end
        %
        %  create the junctions
        %
        iJunc=iJunc+1;
        junctionSeqs(iJunc).Header=sprintf('Junction_%s_%s',thisTailName,thisHeadName);
        junctionSeqs(iJunc).Sequence=[thisTailSeq,thisHeadSeq];
    end
end


    
    
