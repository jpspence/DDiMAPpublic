#!/usr/bin/env python
from sys import exit, argv
from os import path, makedirs, remove, stat
from itertools import repeat
from subprocess import Popen, call, STDOUT
from shutil import copyfile, rmtree

import csv
import ConfigParser
import argparse
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main(argv):
    """
    Main function that reads in configuration files from config file and command line and performs iterations
    """

    mapperAbbrs = {'C':'cushaw', 'S':'shrimp', 'B':'bfast', 'W':'bwa-mem', 'N':'novoalign'}

    #Dictionary of commands to use for various mappers - configure your mapper commands here
    aligner_dict = {
	'B,CS,S':[
		'bfast fasta2brg -f DDiFasta -A 0',
		'bfast fasta2brg -f DDiFasta -A 1',
		'bfast index -f DDiFasta -m 1111111111111111111111 -w 14 -i 1 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 111110100111110011111111111 -w 14 -i 2 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 10111111011001100011111000111111 -w 14 -i 3 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 1111111100101111000001100011111011 -w 14 -i 4 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 111111110001111110011111111 -w 14 -i 5 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 11111011010011000011000110011111111 -w 14 -i 6 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 1111111111110011101111111 -w 14 -i 7 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 111011000011111111001111011111 -w 14 -i 8 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 1110110001011010011100101111101111 -w 14 -i 9 -A 1 -n DDiProcs',
		'bfast index -f DDiFasta -m 111111001000110001011100110001100011111 -w 14 -i 10 -A 1 -n DDiProcs',
		'bfast match -f DDiFasta -A 1 -i 1-10 -k 18 -K 100000 -w 0 -t -n DDiProcs -Q 100000 -l -r DDiFastq1 > DDiBMF',
		'bfast localalign -f DDiFasta -m DDiBMF -A 1 -n DDiProcs -U -q 20 -Q 100000 -t > DDiBAF',
		'rm DDiBMF',
		'bfast postprocess -f DDiFasta -i DDiBAF -o DDiAligned -O 1 -a 3 -z -n DDiProcs -q 20 -Q 100000 -t > DDiSAM',
		'rm DDiBAF'
	    ],
        'C,CS,S':[
            'cushaw3 index DDiFasta -c -p bwtindex',
            'cushaw3 calign -r bwtindex -f DDiFastq1 -t DDiProcs -multi 1 CushawOpts -o DDiSAM'
            ],
        'C,NT,S':[
            'cushaw3 index DDiFasta -p bwtindex',
            'cushaw3 align -r bwtindex -f DDiFastq1 -t DDiProcs -multi 1 CushawOpts -o DDiSAM'
            ],
        'C,NT,P':[
            'cushaw3 index DDiFasta -p bwtindex',
            'cushaw3 align -r bwtindex -q DDiFastq1 DDiFastq2 -t DDiProcs -multi 1 CushawOpts -o DDiSAM'
            ],
        'S,CS,S':[
            'gmapper-cs -N DDiProcs -Q -o 1 --strata --all-contigs ShrimpOpts DDiFastq1 DDiFasta > DDiSAM'
            ],
        'S,NT,S':[
            'gmapper-ls -N DDiProcs -Q -o 1 --strata --all-contigs ShrimpOpts DDiFastq1 DDiFasta > DDiSAM'
            ],
        'S,NT,P':[
            'gmapper-ls -N DDiProcs -Q -o 1 --strata --all-contigs ShrimpOpts -1 DDiFastq1 -2 DDiFastq2 DDiFasta  > DDiSAM'
            ],
	'W,NT,S':[
            'bwa index DDiFasta',
	    'bwa mem -t DDiProcs BwaMemOpts DDiFasta DDiFastq1 > DDiSAM'
            ],
	'W,NT,P':[
            'bwa index DDiFasta',
	    'bwa mem -t DDiProcs BwaMemOpts DDiFasta DDiFastq1 DDiFastq2 > DDiSAM'
            ],
	'N,NT,S':[
	    'novoindex DDiNIX DDiFasta',
            'novoalign -r Random -n 100 -o SAM -d DDiNIX -f DDiFastq1 > DDiSAM'
            ],
	'N,NT,P':[
	    'novoindex DDiNIX DDiFasta',
            'novoalign -r Random -n 100 -o SAM NovoOpts -d DDiNIX -f DDiFastq1 DDiFastq2 > DDiSAM'
            ]
        }

    #Arguments that are required
    required = ['fastqFiles', 'mappingRefSeqFiles', 'outputDir']

    parser = argparse.ArgumentParser(description='Iteratively calls 3rd party mappers and DDiMap executable')

    #Argument options
    parser.add_argument('-q', type=str, metavar='file', nargs='+', help='list of fastq files', dest='fastqFiles')
    parser.add_argument('-r', type=str, metavar='file', nargs='+', help='list of files to use for reference sequences', dest='mappingRefSeqFiles')
    parser.add_argument('-j', type=str, metavar='file', nargs='+', help='list of files to use for junctions', dest='junctionRefSeqFiles')
    parser.add_argument('-o', type=str, metavar='directory', help='output directory', dest='outputDir')
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', '--paired', action='store_true', help='fastq files have paired ends', dest='pairedEnds')
    group.add_argument('-s', '--single', action='store_false', help='fastq files have single ends', dest='pairedEnds')
    parser.add_argument('-n', type=int, metavar='cpus', help='number of processors to use', dest='nProcs')
    parser.add_argument('-c', type=str, metavar='config_file', help='location of config file', dest='configFile')
    parser.add_argument('-v', action='store_true', help='turns on verbosity', dest='verbose')

    parser.add_argument('--aligner_order', type=str, metavar='{'+','.join(mapperAbbrs.keys())+'}', help='mapper sequence as a string. ie CSC', dest='alignerOrder')
    parser.add_argument('--first_iter', metavar='n', type=int, help='first iteration', dest='firstIter')
    parser.add_argument('--max_iters', metavar='n', type=int, help='maximum iterations', dest='maxIters')
    parser.add_argument('--read_length', metavar='n', type=int, help='read length', dest='readLength')
    parser.add_argument('--read_type', type=str, help='read type', choices=['CS','NT'], dest='readType')
    parser.add_argument('--req_frag_conv', help='require frags to converge as well as SNVs', action='store_true', dest='reqFragConv')
    parser.add_argument('--no-req_frag_conv', help='does not require frags to converge as well as SNVs', action='store_false', dest='reqFragConv')

    parser.add_argument('--frag_maker_thresh',type=float, metavar='threshold', help='verified frag maker threshold', dest='fragMakerThresh')
    parser.add_argument('--frag_thresh', type=float, metavar='threshold', help='unverified frag maker threshold', dest='fragThresh')
    parser.add_argument('--min_absolute_cover', type=int, metavar='n', help='minimum absolute cover', dest='minAbsoluteCover')
    parser.add_argument('--snv_thresh', type=float, metavar='threshold', help='SNV threshold', dest='SNVthresh')
    parser.add_argument('--snv_type2_thresh', type=float, metavar='threshold', help='SNV type 2 threshold', dest='SNVtype2thresh')
    parser.add_argument('--snv_type3_thresh', type=float, metavar='threshold', help='SNV type 3 threshold', dest='SNVtype3thresh')
    parser.add_argument('--roa_size', type=int, metavar='size', help='Size to use for region of analysis in DDiMAP', dest='roaSize')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--use_DI', action='store_true', help='use reads mapped with deletion and insertion', dest='useDI')
    group.add_argument('--no-use_DI', action='store_false', help='do not use reads mapped with deletion and insertion', dest='useDI')

    parser.add_argument('--cushaw_opts', type=str, metavar="'options'", help='cushaw specific options', dest='cushawOpts')
    parser.add_argument('--shrimp_opts', type=str, metavar="'options'", help='shrimp specific options', dest='shrimpOpts')
    parser.add_argument('--bwamem_opts', type=str, metavar="'options'", help='bwa-mem specific options', dest='bwaMemOpts')
    parser.add_argument('--novo_opts', type=str, metavar="'options'", help='novoalign specific options', dest='novoOpts')


    #Parse args and check for config file
    args = parser.parse_args()
    if args.configFile:
        configFile = args.configFile
        if not path.isfile(configFile):
            print 'config file specified, but not found'
            exit(1)
    else:
        configFile = 'DDiMap.cfg'

    #Read in settings from config file
    Settings = read_config(configFile)

    #  Loop over each section and replace values with those passed in on command line.  
    #  Also create a local variable that matches the keys in the settings dictionary.

    for section in Settings.keys():
        for key in Settings[section].keys():
            if getattr(args, key):
                Settings[section][key] = getattr(args, key)
            exec '%s = Settings[section][key]' % key
            if key in required and not Settings[section][key]:
                print '%s not specified on command line or in config file.  Aborting...' % key
                print Settings[section][key]
                parser.print_help()
                exit(1)
            if (type(Settings[section][key]) == list):
                Settings[section][key] = ', '.join(Settings[section][key])

    if useDI:         #  reads with CIGARs containing both I and D are processed
        kFlag='-k'
    else:             #  reads with CIGARs containing both I and D are not processed
        kFlag=''

    if pairedEnds:
        pair_str='P'
    else:
        pair_str='S'

    #   do the work - set up for the iteration
    aligners = list(alignerOrder)
    iterMin = len(aligners)
    iterMax = max(maxIters, iterMin); #  always do as many iters as are in alignerOrder string
    aligners = aligners + list(repeat(aligners[-1], iterMax - iterMin)) #   define the aligner ID sequence to be used over the iterations


    # Make paths absolute
    fastqFiles = [path.abspath(x) for x in fastqFiles]
    mappingRefSeqFiles = [path.abspath(x) for x in mappingRefSeqFiles]
    junctionRefSeqFiles = [path.abspath(x) for x in junctionRefSeqFiles]
    outputDir = path.abspath(outputDir) + '/'

    #  Make sure the output directory exists

    if not path.isdir(outputDir):
        makedirs(outputDir)

    # Write configuration file in outputDir
    write_config(outputDir, Settings)

    #  INITIAL VALUES OF LOOP CONTROL PARAMETERS
    converged = False
    prevFragList = []   #  this will be replaced by counts of fragments created for each baseline refernce sequence
    prevSNVList = []    #  this will be replaced by counts of SNV candidates found for each baseline reference sequence

    thisIter = firstIter


    for RefSeqFile in fastqFiles:
        if not path.isfile(RefSeqFile):
            print 'Unable to find fastqFile at ' + RefSeqFile
            exit(1)

    #  Delete old enhanced fast file if present.  It should never be...

    enhancedFastaFile = outputDir + 'refSeqEnhanced.fa'
    if path.isfile(enhancedFastaFile): #  see if one is already here - need to zap it
        remove(enhancedFastaFile)      #  remove if present because fastawrite appends to existing files
    output_handle = open(enhancedFastaFile, 'a')

    #  Add reference sequences to file with _Ref tag
    RefSeqs=[]
    for RefSeqFile in mappingRefSeqFiles:
	print 'ref seq file = ' + RefSeqFile
        if not path.isfile(RefSeqFile):
            print 'Unable to find RefSeqFile at ' + RefSeqFile
            exit(1)
        RefSeqs = RefSeqs + list(SeqIO.parse(RefSeqFile, 'fasta'))
    if (RefSeqs):
        formattedRefSeqs = add_ref_tag(RefSeqs)
        SeqIO.write(formattedRefSeqs, output_handle, 'fasta')  # modified MATLAB fastawrite to not put in extra newlines

    
    #  Create junctions if they are needed and then add to ref seq file as mapping targets for chimeric reads
    RefSeqs=[]
    for RefSeqFile in junctionRefSeqFiles:
        if not path.isfile(RefSeqFile):
            print 'Unable to find RefSeqFile at ' + RefSeqFile
            exit(1)
        RefSeqs = RefSeqs + list(SeqIO.parse(RefSeqFile, 'fasta'))
    if (RefSeqs):
        formattedRefSeqs = add_ref_tag(RefSeqs)
        junctionSeqs = make_junctions(formattedRefSeqs,readLength);
        SeqIO.write(junctionSeqs, output_handle, 'fasta')  # modified MATLAB fastawrite to not put in extra newlines

    output_handle.close()            


    #  allows restarts
    if thisIter > 1:  #  there is no previous iteration, so start fresh
        prevWorkingDir = outputDir + ('Gen%d/' % (thisIter-1))
        for i in range(1, thisIter):
            prevWorkingDir = '%sGen%d/' % (outputDir, i) 
            fragFile = prevWorkingDir + 'fasta.fa'
            snvFile = prevWorkingDir + 'snv.csv'
            ddimap_convergence_test(fragFile, snvFile, prevFragList, prevSNVList, reqFragConv)


    while not converged and thisIter <= iterMax:
        
        print '======= Iteration %d of %d  ========' % (thisIter, iterMax)

        # creates working dir if not present
        thisWorkingDir = outputDir + ('Gen%d/' % thisIter)
        if path.isdir(thisWorkingDir):
            rmtree(thisWorkingDir)
        makedirs(thisWorkingDir)
        
        # Delete old enhanced fast file if present.  It should never be...
        enhancedFastaFile = thisWorkingDir + 'refSeqEnhanced.fa'
        if path.isfile(enhancedFastaFile): 
            remove(enhancedFastaFile)      
        copyfile(outputDir + 'refSeqEnhanced.fa', enhancedFastaFile)

        output_handle = open(enhancedFastaFile, 'a')
            
        # Append frags from previous iteration if any (these sequences are tagged as fragments when the file is written by DDiMAP)
        if (thisIter > 1):
            prevFragFile=prevWorkingDir + '/fasta.fa'
            if path.isfile(prevFragFile) and path.getsize(prevFragFile) > 0:
                fragSeqs=list(SeqIO.parse(prevFragFile, 'fasta'))
                SeqIO.write(fragSeqs, output_handle, 'fasta')  # modified MATLAB fastawrite to not put in extra newlines

        output_handle.close()    

        # Setup variables for aligner
        thisAligner=aligners[thisIter-1]
        thisAligned='DDiMAP_%s' % thisAligner
 
        if path.isfile(thisWorkingDir + 'mapper.log'):
            remove(thisWorkingDir + 'mapper.log')

        if not ','.join([thisAligner,readType,pair_str]) in aligner_dict.keys():
            print mapperAbbrs[thisAligner] + ' does not support ' + readType + ' read type with ' + ('paired ends' if pairedEnds else 'non paired ends')
            exit(1)


        # execute commands for aligner

        open(thisWorkingDir + 'mapper.log', 'w').close()
        if verbose:
            b=Popen(['tail', '-F', thisWorkingDir + 'mapper.log'])

        # set substitutions for aligner commands
        commandsubs={'DDiFastq1':fastqFiles[0], 
                     'DDiProcs':nProcs, 
                     'DDiFasta':enhancedFastaFile, 
                     'DDiBMF':thisAligned + '.bmf', 
                     'DDiBAF':thisAligned + '.baf', 
                     'DDiSAM':thisAligned + '.sam',
                     'DDiNIX':thisAligned + '.nix', 
                     'DDiAligned':thisAligned, 
                     'CushawOpts':cushawOpts, 
                     'ShrimpOpts':shrimpOpts, 
                     'BwaMemOpts':bwaMemOpts, 
                     'NovoOpts':novoOpts}

        if (len(fastqFiles) > 1):
            commandsubs['DDiFastq2']=fastqFiles[1]

        for command in aligner_dict[','.join([thisAligner,readType,pair_str])]:
            cmdlist=re.split('\s*',command)
            #remove empty arguments and subsitute in values from commandsubs            
            args=filter(None,[str(commandsubs[x]) if x in commandsubs.keys() else x for x in cmdlist])
            args=re.split('\s*',' '.join(args)) 
            print ' '.join(args)   # output actual command
            if 'DDiFastq2' in args: #This hasn't been substituted because one wasn't provided
                print mapperAbbrs[thisAligner] + ' expects 2 fastq files for use with ' + readType + ' read type with ' + ('paired ends' if pairedEnds else 'non paired ends')
                exit(1)

            # Now we need to detect stdout redirection and do it properly using pOpen
            if '>' in args: 
                i = args.index('>')
                outfile = args[i+1]
                del args[i:i+2]
            else:
                outfile = None
   
            log_file = open(thisWorkingDir + 'mapper.log', 'a')
            
            if (outfile):
                with open(thisWorkingDir + outfile, 'w') as output_file:
                    a=Popen(args, cwd=thisWorkingDir, stdout=output_file, stderr=log_file)
            else:
                a=Popen(args, cwd=thisWorkingDir, stderr=log_file, stdout=log_file)

            success=a.wait()
            log_file.close()
            if not success == 0:
                print '*** mapper exited with error', success
                print 'See ' + thisWorkingDir + 'mapper.log' + ' for more details'
                exit(success)

        if verbose:
            b.terminate()
        # Perform sam to bam conversion for DDiMap
        args=['samtools', 'view', '-b', '-S', '-o', thisAligned + '.bam', thisAligned + '.sam']
        print ' '.join(args)   

        open(thisWorkingDir + 'samtools.log', 'w').close()
        if verbose:
            b=Popen(['tail', '-F', thisWorkingDir + 'samtools.log'])
        log_file = open(thisWorkingDir + 'samtools.log', 'w')
        a=Popen(args, cwd=thisWorkingDir, stderr=log_file, stdout=log_file)
        success=a.wait()
        log_file.close()
        if verbose:
            b.terminate()
        if not success == 0:
            print '*** samtools exited with error', success
            print 'See ' + thisWorkingDir + 'samtools.log' + ' for more details'            
            exit(success)
        # remove the uncompressed sam file
        args=['rm', thisAligned + '.sam'];
        a=Popen(args, cwd=thisWorkingDir)

        # now run the DDiMAP code
        thisAlignedFile = thisWorkingDir + thisAligned + '.bam'
        args = (['DDiMAP', kFlag, '-r', roaSize, '-f', enhancedFastaFile, '-b', 
                 thisAlignedFile, '-c', minAbsoluteCover, '-n', fragThresh, '-a', 
                 fragMakerThresh, '-p', SNVthresh, '-s', SNVtype2thresh, '-l', 
                 SNVtype3thresh, '-o', thisWorkingDir])
        args = [str(x) for x in args]
        print ' '.join(args)
        open(thisWorkingDir + 'DDiMap.log', 'w').close()
        if verbose:
            b=Popen(['tail', '-F', thisWorkingDir + 'DDiMap.log'])
        log_file = open(thisWorkingDir + 'DDiMap.log', 'a')
        a = Popen(args, cwd=thisWorkingDir, stdout=log_file, stderr=log_file)
        success=a.wait()
        if verbose:
            b.terminate()
        log_file.close()
        if not success == 0:
            print '*** DDiMap exited with error', success
            print 'See ' + thisWorkingDir + 'DDiMap.log' + ' for more details'
            exit(success)
        
        # now check for convergence
    
        fragFile = thisWorkingDir + 'fasta.fa'
        snvFile = thisWorkingDir + 'snv.csv'
    
        #  call to the convergence test matlab function
        #  result history kept in currFrags/prevFrags and currSNVs/prevSNVs
    
        if ddimap_convergence_test(fragFile, snvFile, prevFragList, prevSNVList, reqFragConv):
            print 'Convergence found.  Stopping...'
            break

        prevWorkingDir = thisWorkingDir;  #  all done with the previous, this will be the next iteration previous directory
        thisIter = thisIter+1
    else:
        print 'Failed to converge'

    print '%10s %10s %10s' % ('Iteration', 'nFrags', 'nSNVs')
    for i,  (frags, snvs) in enumerate(zip(prevFragList, prevSNVList)):
        print '%10d %10d %10d' % (i+1, sum(frags), sum(snvs))

    #  put final results into outputDir
    #  make renamed copies of the final iteration result files, naming them using
    copyfile(thisWorkingDir+'fasta.fa',outputDir+'convergedFrags.fa')
    copyfile(thisWorkingDir+'dictionary.csv',outputDir+'convergedDictionary.csv')
    copyfile(thisWorkingDir+'snv.csv',outputDir+'convergedSNVs.csv')
    copyfile(thisWorkingDir+'coverage.csv',outputDir+'convergedCoverage.csv')
    copyfile(thisWorkingDir+'refSeqEnhanced.fa',outputDir+'convergedEnhancedRefSeqs.fa')
    
    
def read_config(configFile):
    Config = ConfigParser.ConfigParser() 
    Config.optionxform = str 

    #Set default values for vars
    Settings={
        'Files':{
            'junctionRefSeqFiles':'',
            'mappingRefSeqFiles':'',
            'fastqFiles':'',
            'outputDir':'',
            'pairedEnds':False},
        'Settings':{
            'readLength':50,
            'nProcs':'8',
            'firstIter':1,
            'maxIters':20,
            'readType':'NT',
            'alignerOrder':'C',
            'reqFragConv':True,
            'verbose':False},
        'DDiMap':{
            'fragMakerThresh':'0.01',
            'SNVtype3thresh':'0.1',
            'SNVtype2thresh':'0.003000',
            'minAbsoluteCover':'2',
            'SNVthresh':'0.000750',
            'roaSize':'32',
            'useDI':False,
            'fragThresh':'0.1'},
        'MapperOptions':{
            'cushawOpts':'',           
            'bwaMemOpts':'',
            'novoOpts':'',
            'shrimpOpts':''}
        }
    
    if path.isfile(configFile):
        Config.read(configFile)
        for section, settings in Settings.items():
            if Config.has_section(section):
                for key, value in Config.items(section):
                    if type(Settings[section][key]) == bool: #casting string as boolean always returns true
                        Settings[section][key]=(value in ['True','true','yes','y']);
                    else:
                        Settings[section][key]=type(Settings[section][key])(value)
                    
    #Loop over input files to parse comma, whitespace, or new line delimited lists (with ';' comments) to make up for parser shortcomings
    for key in ['mappingRefSeqFiles','junctionRefSeqFiles','fastqFiles']:
        Settings['Files'][key]=[x for x in re.split(r'[\s,]+',re.sub(r'\s*;.*','', Settings['Files'][key])) if x]
        
    return Settings


def write_config(outputDir, Settings):
    Config = ConfigParser.ConfigParser() 
    Config.optionxform = str 
    for section, settings in Settings.items():
        Config.add_section(section)
        for key, value in settings.items():
            Config.set(section, key, value)        
    with open(outputDir+'DDiMap.cfg', 'wb') as configfile:
        Config.write(configfile)


def add_ref_tag(basicSeqs):
    """
    takes sequence header info and creates friendly seq names
    terminates at first white space & replace underscore with hyphens
    appends _Ref tag
    """

    formattedBasicSeqs=list(basicSeqs)   
    for record in formattedBasicSeqs:
        record.id=record.id+'_Ref'
        record.name=record.name+'_Ref'
        record.description=record.description+'_Ref'
    return formattedBasicSeqs


def make_junctions(RefSeqs,readLength):
    """
    publicSeqs and privateSeqs are Sequences lists 
    create junctional seqences using readLength-1 long tail to head 
    """

    nRefSeqs = len(RefSeqs)
    htLength = readLength - 1 #in matlab readLength-1, but index starts from 0 for python 

    junctionSeqs = []
    heads = []
    tails = []
    names = []

    for iSeq in range(nRefSeqs):
        heads.append(RefSeqs[iSeq].seq[0:htLength])
        tails.append(RefSeqs[iSeq].seq[-htLength:])
        names.append(list(RefSeqs[iSeq].id))
        if '_' in names[iSeq]:
            names[iSeq][names[iSeq].index('_'):]='' #Removes everything from the name after the '_'
            
    for hSeq in range(nRefSeqs): 
        for tSeq in range(nRefSeqs): 
            # create the junctions
            junctionNm = 'Junction_'+''.join(names[tSeq])+'_'+''.join(names[hSeq])
            junctionSeq = tails[tSeq] + heads[hSeq]
            junctionSeqs.append(SeqRecord(id=junctionNm,name=junctionNm, description=junctionNm, seq=junctionSeq)) #miss description

    return junctionSeqs 


def ddimap_convergence_test(fragFile,snvFile, prevFragList, prevSNVList, reqFragConv):
    """
    forms histogram of frags by ref seq and snvs by gene 
    compares to results from prior iteration 
    """

    # construct frag histogram

    if (path.exists(fragFile) and stat(fragFile).st_size != 0):
        fragHandle = open(fragFile, 'rU')
        fragInfo = list(SeqIO.parse(fragHandle, 'fasta'))
        fragids=[x.id for x in fragInfo]
    else:
        fragids=[]

    kFrags=[fragids.count(x) for x in set(fragids)]
    currFrags = kFrags 


    # construct snv histogram

    fragids = []
    if (path.exists(snvFile) and stat(snvFile).st_size != 0):
        fidSnv = open(snvFile, 'rb')
        reader = csv.reader(fidSnv) 
        fragids=[x[0] for x in list(reader)[1:]]
        fidSnv.close()

    kSnvs=[fragids.count(x) for x in set(fragids)]
    currSNVs = kSnvs


    #	compare histograms

    fragsConverged = False
    snvsConverged = False


    for i, (prevFrags,prevSNVs) in enumerate(zip(prevFragList, prevSNVList)):
        fragsConverged = len(prevFrags) == len(currFrags) and currFrags == prevFrags
        snvsConverged = len(prevSNVs) == len(currSNVs) and currSNVs == prevSNVs
        if (snvsConverged and (not reqFragConv or fragsConverged)):
            print 'Gen %d and Gen %d match' % (len(prevFragList)+1, i+1)
            break
        
    prevFragList.append(currFrags)
    prevSNVList.append(currSNVs)      
    return snvsConverged and (not reqFragConv or fragsConverged)
    


if __name__=='__main__':
    main(argv)
