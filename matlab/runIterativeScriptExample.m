%
%  loop over folders to run DDiMAPiterative
%
%specimenIDlist={'126','127','128','129','131','132','133','134','135','136','137','139','140','155','158','293'};
%specimenIDlist={'126','127','128','129','131','132','133','134','135','136','137','139'};
%specimenIDlist={'128','133','134','139','136','132','126','127','129','131','135','137'};
%specimenIDlist={'128','133','134','139','136','132'};
specimenIDlist={'134','136','132','128'};
%specimenIDlist={'RB159','RB160'};
%specimenIDlist={'139'};
nSpecs=length(specimenIDlist);
for iSpec=1:nSpecs
    specimenID=specimenIDlist{iSpec}
    
    inputDir=sprintf('/home/jspence/Data/Burack_%s/',specimenID);
    workingDirPattern=fullfile(inputDir,'DDiMAP32_CUSHAW3_VhJwithHighFrag_%s_Gen%d/');
    outputDir=fullfile(inputDir,'DDiMAP32results_CUSHAW3_VhJwithHighFrag/');
    
    fastqFile=sprintf('Burack_%s.fastq',specimenID);
    
    publicRefSeqFile=sprintf('%s_Public_DDiMAP.fa',specimenID);  %  just the ref sequences, extra newlines between seqs OK for MATLAB, in inputDir
    germlineRefSeqFile=sprintf('%s_Germline_DDiMAP.fa',specimenID);  %  just the ref sequences that may be used for mapping but never for junction making
    privateRefSeqFile=sprintf('%s_Private_DDiMAP.fa',specimenID);   %  just the sangered private sequences that are always used for junction making and may be used for mapping

    privateMapped=false;                  %  if true, use private for mapping, if false use germline for mapping
    concatenatedUsed=true;               %  sample prep used concatenation (or not) so junctions will be made
    readType='cs';                       %  read type is 'cs' or 'c' (ABI SOLiD) or 'nt' or 'n' (Illumina, etc.)
    readLength=50;                       %  nominal max read length
    alignerOrder='C';                    %  valid chars are B (BFAST), C (CUSHAW3), S (SHRiMP2)
    maxIters=7;                          %  MAXIMUM number of iters (will be adjusted to accomodate length of alignerOrder
    firstIter=1;                         %  can restart using this
    roaSize=32;                          %  must be even, if odd, increment
    fragThresh=0.1;                      %  default value is 0.1, 1.0 turns off unverified frag cores
    %
    %  run the iterative script
    %
    DDiMAPiterative
end