function [fragsConverged,snvsConverged,currFrags,currSNVs]=DDiMAPconvergenceTest(fragFile,snvFile,prevFrags,prevSNVs)
%
%   forms histogram of frags by ref seq and snvs by gene
%   compares to results from prior iteration
%

%
%   construct frag histogram
%
if exist(fragFile,'file')
fragInfo=fastaread(fragFile);
else
    fragInfo=[];
end
nFrags=length(fragInfo);
currRefSeq='';
nRefSeqF=0;
kFrags=zeros(nFrags,1); %  allocate counter spece
for iFrag=1:nFrags
    thisFragName=fragInfo(iFrag).Header;
    refSeqEnd=strfind(thisFragName,'_Frag');
    thisRefSeq=thisFragName(1:refSeqEnd);
    if strcmp(thisRefSeq,currRefSeq)  %  found another one
        kFrags(nRefSeqF)=kFrags(nRefSeqF)+1;
    else   % found a different one
        nRefSeqF=nRefSeqF+1;
        kFrags(nRefSeqF)=1;
        currRefSeq=thisRefSeq;
    end
end
currFrags=kFrags(1:nRefSeqF);
%
%   construct snv histogram
%
if exist(snvFile,'file')
    
    fidSnv=fopen(snvFile,'r');
    resultCell=textscan(fidSnv,'%s %*[^\n]','Delimiter',',','HeaderLines',1);  % cell array containing the reg seq names
    fclose(fidSnv);
else
    resultCell={[]};
end
nSnvs=length(resultCell{1});

currRefSeq='';
nRefSeqS=0;
kSnvs=zeros(nSnvs,1);
for iSnv=1:nSnvs
    thisRefSeq=resultCell{1}(iSnv);
    if strcmp(thisRefSeq,currRefSeq)  %  found another one
        kSnvs(nRefSeqS)=kSnvs(nRefSeqS)+1;
    else
        currRefSeq=thisRefSeq;
        nRefSeqS=nRefSeqS+1;
        kSnvs(nRefSeqS)=1;
    end
end
currSNVs=kSnvs(1:nRefSeqS);
%
%   compare histograms
%
if length(prevFrags)==nRefSeqF %  same number of ref seqs with frags found
    fragsConverged=all(currFrags==prevFrags); % all counts must match
else
    fragsConverged=false;
end
if length(prevSNVs)==nRefSeqS  %  same number of ref seqs with SNVs found
    snvsConverged=all(currSNVs==prevSNVs);  %  all counts must match
else
    snvsConverged=false;
end
