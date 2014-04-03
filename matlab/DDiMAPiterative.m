%
%  iterative mapping script
%

%
%  basic parameters (these would come from a calling script or user interface)
%
% specimenID='127';
% inputDir='/home/jspence/Data/Burack_127/'  
% workingDirPattern='/home/jspence/Data/Burack_127/DDiMAP_%s_Gen%d/'  % will be used to put in aligner order and iteration number
% outputDir='/home/jspence/Data/Burack_127/DDiMAP_Results/'
% 
% fastqFile='Burack_127.fastq'         %  the raw read data, in inputDir
% publicRefSeqFile='127_Public_DDiMAP.fa'  %  just the ref sequences, extra newlines between seqs OK for MATLAB, in inputDir
% germlineRefSeqFile='127_Germline_DDiMAP.fa'  %  just the ref sequences that may be used for mapping but never for junction making
% privateRefSeqFile='127_Private_DDiMAP.fa'   %  just the sangered private sequences that are always used for junction making and may be used for mapping
% privateMapped=false                %  if true, then private is mapped and
                                     %     used for junctions - useful for
                                     %     mapping IGH to Sanger
                                     %  if false, germline is mapped and
                                     %     private is used for junctions - if
                                     %     it exists - useful for mapping
                                     %     IGH to germline Vh and J

% concatenatedUsed=true                %  sample prep used concatenation (or not) so junctions will be made
% readType='cs'                        %  read type is 'cs' or 'c' (ABI SOLiD) or 'nt' or 'n' (Illumina, etc.)
% readLength=50;                       %  nominal max read length
% alignerOrder='C'                     %  valid chars are B (BFAST), C (CUSHAW3), S (SHRiMP2)
% maxIters=6;                           %  MAXIMUM number of iters (will be adjusted to accomodate length of alignerOrder
% firstIter=2;                         %  can restart using this
% roaSize=34;                          %  must be even, if odd, increment
                                     %     - overlapSize is roaSize/2
                                     %     - first track start is 1
                                     %     - seond track start is 1+
                                     %                  roaSize/4 (chopped)
                                     %     - may consider having a gene
                                     %     specific choice like first track
                                     %     start at geneLength/overlapSize remainder
                                     %      divided by 2, (chopped) to center
                                     %      the track pair in the ref sequence
                                     %    - second track start offset is roaSize/4, chopped.
                                     %
%
%   environment stuff
%
scriptPath='/home/jspence/DDiMAPscripts';
DDiMAPpath='/home/jspence/DDiMAP/bin/DDiMAP';
%
%  do the work - set up for the iteration
%
iterMin=length(alignerOrder);
iterMax=max(maxIters,iterMin);  %  always do as many iters as are in alignerOrder string
aligners=zeros(1,iterMax,'uint8');  % will strecth the string by rightward extenstion out to iterMax
lastAligner=alignerOrder(end);
aligners(1:length(alignerOrder))=alignerOrder;
if iterMax>length(alignerOrder)
    aligners((length(alignerOrder)+1):end)=lastAligner;
end
%
%  set up the loop
%
thisFastqFile=fullfile(inputDir,fastqFile);
switch lower(readType(1))
    case 'c'
        alignerReadType='CS';
    case 'n'
        alignerReadType='NT';
end
converged=false;
prevFrags=[];
prevSNVs=[];
thisIter=firstIter-1;
if thisIter==0
    prevWorkingDir=inputDir;
else
    prevWorkingDir=sprintf(workingDirPattern,alignerOrder,thisIter);
end
while ~converged && thisIter<iterMax
    thisIter=thisIter+1;
    %
    %   creates working dir if not present
    %
    thisWorkingDir=sprintf(workingDirPattern,alignerOrder,thisIter);
    if ~isdir(thisWorkingDir)
        mkdir(thisWorkingDir)
    end
    %
    %  sequence preparation - reads original fasta files (private and public), creates junctions if
    %  concatenated used, adds privaate seq or alternativer seqs, adds frags if present
    %
    basicSeqs=fastaread(fullfile(inputDir,publicRefSeqFile));
    formattedBasicSeqs=addRefTag(basicSeqs);
    enhancedFastaFile=fullfile(thisWorkingDir,'refSeqEnhanced.fa');
    if exist(enhancedFastaFile,'file')  %  see if one is already here - need to zap it
        delete(enhancedFastaFile)   %  remove if present because fastawrite appends to existing files
    end
    fastawriteNN(enhancedFastaFile,formattedBasicSeqs);  % modified fastawrite to not put in extra newlines
    %
    %  check for existence of private sequences that are to be used for
    %  junction making and optionally added to ref seq file for alignment
    %
    if exist(fullfile(inputDir,privateRefSeqFile),'file')  %  add to sequences used for junction making
        privateSeqs=fastaread(fullfile(inputDir,privateRefSeqFile));
        formattedPrivateSeqs=addRefTag(privateSeqs);
        if privateMapped  % add to mapping targets
            fastawriteNN(enhancedFastaFile,formattedPrivateSeqs)  %  append to file
        end
    else
        formattedPrivateSeqs=[];    % NULL pointer    
    end
    %
    %  check for existence of any Germline sequences to be used for mapping
    %  only if the private sequences are not being used.  Add to ref seq
    %  file in this case.
    %
    if ~privateMapped && exist(fullfile(inputDir,germlineRefSeqFile),'file')  %  add to sequences used for mapping
        germlineSeqs=fastaread(fullfile(inputDir,germlineRefSeqFile));
        formattedGermlineSeqs=addRefTag(germlineSeqs);
        fastawriteNN(enhancedFastaFile,formattedGermlineSeqs)
    else
        formattedGermlineSeqs=[];    % NULL pointer
    end
    %
    %  call to a function that creates junctions if they are needed
    %  then add to ref seq file
    %
    if concatenatedUsed
        junctionSeqs=makeJunctions(formattedBasicSeqs,formattedPrivateSeqs,readLength);
        fastawriteNN(enhancedFastaFile,junctionSeqs)
    end
    %
    %  append frags from previous iteration if any
    %
    prevFragFile=fullfile(prevWorkingDir,'fasta.fa')
    if exist(prevFragFile,'file')
        fragSeqs=fastaread(prevFragFile);
        fastawriteNN(enhancedFastaFile,fragSeqs)
    end
    prevWorkingDir=thisWorkingDir;
    %
    %  prepare aligner environment
    %
    thisAligner=aligners(thisIter);
    thisAligned=sprintf('DDiMAP_%s_%s',specimenID,thisAligner);
    %
    %  set environment variables for script execution
    %
    setenv('myWorkingDir',thisWorkingDir);
    setenv('mySpecimenID',specimenID)
    setenv('myFasta',enhancedFastaFile)
    setenv('myFastq',thisFastqFile)
    setenv('myAligned',thisAligned)
    
    switch thisAligner
        case 'C' % use cushaw3
            alignerProgram='CUSHAW3'
        case 'B' % use BFAST
            alignerProgram='BFAST'
        case 'S' % use SHRiMP2
            alignerProgram='SHRiMP'
    end
    unixCmd=sprintf('%s/run%sforDDiMAP_%s.sh',scriptPath,alignerProgram,alignerReadType)
    %
    %   run the command and capture the output for a log file
    [status,alignerResult]=unix(unixCmd,'-echo');
    logFile=fullfile(thisWorkingDir,sprintf('%s.log',alignerProgram));
    fid=fopen(logFile,'w');
    fwrite(fid,alignerResult,'char');
    fclose(fid);

    %
    %  now run the DDiMAP code
    %
    setenv('LD_LIBRARY_PATH','/usr/local/lib/bamtools/');  % this probably should be changed to adding it to the path if it is not found
    thisAlignedFile=fullfile(thisWorkingDir,sprintf('%s.bam',thisAligned));
    %
    %  this should be changed to some code that constructs the command line
    %  according to the parameters set in the calling script or UI
    %
    DDiMAPcmd=sprintf('%s -r %d -f %s -b %s -n %f -o %s\n',DDiMAPpath,roaSize,enhancedFastaFile,thisAlignedFile,fragThresh,thisWorkingDir)
    [status,DDiMAPresult]=unix(DDiMAPcmd,'-echo');
    logFile=fullfile(thisWorkingDir,'DDiMAP.log');
    fid=fopen(logFile,'w');
    fwrite(fid,DDiMAPresult,'char');
    fclose(fid);

    %
    %  now check for convergence
    %
    fragFile=fullfile(thisWorkingDir,'fasta.fa');
    snvFile=fullfile(thisWorkingDir,'snv.csv');
    %
    %  call to the convergence test code
    %      result history kept in currFrags/prevFrags and currSNVs/prevSNVs
    %
    [fragsConverged,snvsConverged,currFrags,currSNVs]=DDiMAPconvergenceTest(fragFile,snvFile,prevFrags,prevSNVs)
    prevFrags=currFrags;
    prevSNVs=currSNVs;
    %
    %  convergence based on stabilized SNV calls (this may be driven by a
    %  UI/calling script parameter)
    %
    converged=snvsConverged;
end
%
%  put final results into outputDir
%
if ~isdir(outputDir)
    mkdir(outputDir)
end
%
copyfile(fullfile(thisWorkingDir,'fasta.fa'),fullfile(outputDir,sprintf('convergedFrags%s.fa',specimenID)));
copyfile(fullfile(thisWorkingDir,'snv.csv'),fullfile(outputDir,sprintf('convergedSNVs%s.csv',specimenID)));
copyfile(fullfile(thisWorkingDir,'coverage.csv'),fullfile(outputDir,sprintf('convergedCoverage%s.csv',specimenID)));
copyfile(fullfile(thisWorkingDir,'refSeqEnhanced.fa'),fullfile(outputDir,sprintf('convergedEnhancedRefSeqs%s.fa',specimenID)));
logfiles=dir(fullfile(thisWorkingDir,'*.log'));
for iLog=1:length(logfiles)
    copyfile(fullfile(thisWorkingDir,logfiles(iLog).name),fullfile(outputDir,logfiles(iLog).name));
end   