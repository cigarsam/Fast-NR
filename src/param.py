import sys
import os
import argparse
from math import log
#import logging


def defOpt () :
    """Define the required and optional parameters.   
    """
    parser = argparse.ArgumentParser(description='Find negative regulatory element in STARR-seq related method.', 
    prog='Fast-NR', 
    usage='''%(prog)s [options] -t <treatmentFile> -c <controlFile> -g <genomeFile> -o <outputName>
    Example: %(prog)s -t treatment.bed -c control.bed -g dm3.chrom.sizes -o dm3_silencer''',
    epilog='''The output file is tab splitted and contain 9 column:
    chromosome  start   end summit  tCount  cCount  FoldChange  pValue  correctedP
    ''')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument("-t", '--treatment', type=str, required=True, help="Control file path and name, usually the STARR-seq plasmid library. Only .bam and .bed file accpected.")
    required.add_argument("-c", '--control', type=str, required=True, help="Treatment file path and name, usually the STARR-seq cDNA library. Only .bam and 3 column or 6 column .bed file accpected. But this program will not consider the strand information.")
    required.add_argument("-g", '--genome', type=str, required=True, help="Genome chromosome size file, which can download from the UCSC. File include two column, chromosome, length of chromosome.")
    required.add_argument("-o", '--output', type=str, required=True, help="Output path and prefix name of the final result peaks.")
    optional.add_argument("-p", '--pvalue', type=float, default=1e-05, help="The cut-off of p value. Default: %(default)s.")
    optional.add_argument("-cp",'--correctPvalue', type=str, choices=['Bonferroni', 'BH'], default='Bonferroni', help="The correct method of p value. Include 'BH', 'Bonferroni'. Default: %(default)s.")
    optional.add_argument("-cm",'--curveMethod', type=str, choices=['Cosine', 'Pearson', 'Euclidean', 'Gradient'], default='Cosine', help="The method used to calculate the curve similarity. Include 'Cosine', 'Pearson', 'Euclidean', and 'Gradient'. Default: %(default)s.")
    optional.add_argument("-ct", '--curveThreshold', type=float, default=0.25, help="The percent cut-off of similarity distance. From 0 to 1. 1 means accept all regions, do not compute the similarity. Default: %(default)s.")
    optional.add_argument("-ws", '--windowSize', type=int, default=600, help="The size of window, used to find differential coverage region. Default: %(default)s bp.")
    optional.add_argument("-l", '--length', type=int, default=0, help="Extend the fragment length to this fixed bp length from 5'. Default: %(default)s.")
    optional.add_argument('--remove', action="store_true",help="Remove the slop or not. Set to 'True' if want to remove the results in slop regions. Default: false.")
    optional.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    #parser.print_help()
    return parser


def option (parser) :
    args = parser.parse_args()
    
    #if input contain the required parameter
    if not args.treatment or not args.control or not args.genome or not args.output:       
        parser.print_help()
        sys.exit(1)
    #if input format is correct and accessible.
    if os.path.isfile(args.treatment) and os.access(args.treatment, os.R_OK):
        basename = os.path.basename(args.treatment)
        Format = basename.split(".")[-1]  
        if Format == "bed" or Format == "bam":
            args.tFormat = Format
        else:
            print("Error: Format of treatment file is not correct!")
            sys.exit(1)
    else: 
        print("Error: Treatment file %s missing or not accessible." % args.treatment)
        sys.exit(1)
    if os.path.isfile(args.control) and os.access(args.control, os.R_OK):
        basename = os.path.basename(args.control)
        Format = basename.split(".")[-1]  
        if Format == "bed" or Format == "bam":
            args.cFormat = Format
        else:
            print("Error: Format of control file is not correct!")
            sys.exit(1)    
    else:   
        print("Control file %s is missing or not accessible." % args.control)
        sys.exit(1)

    if not os.path.isfile(args.genome) or not os.access(args.genome, os.R_OK):
        print("Genome size file %s is missing or not accessible." % args.genome)
        sys.exit(1)
    
    if (args.curveThreshold <=0) or (args.curveThreshold>1):
        print("The input threshold of curve similarity [%0.2e] is not correct." % args.curveThreshold)
        sys.exit(1)        
    
    outPath = os.path.split(args.output)[0]      
    if len(outPath) >0 and (not os.path.exists(outPath) or not os.access(outPath, os.R_OK)):
        print("OutputPath %s is missing or not accessible." % outPath)
        sys.exit(1)
    elif len(outPath) == 0: #only the name, no path.
        outPath = "./"
 
    args.log_pvalue = log(args.pvalue,10)*-10
   
    # output filename
    args.outfile = args.output+"_peaks.bed"
    if os.path.isfile(args.outfile) and os.access(args.outfile, os.R_OK):
        print("The output file %s already exist!" % args.outfile)
        sys.exit(1)    
    
    #Show the command line
    print("Command: Fast-NR -t %s -c %s -g %s -o %s -p %0.2e -ws %d -cp %s -cm %s -ct %0.2e -l %d -rm %s." % (args.treatment, args.control, args.genome, args.output, args.pvalue, args.windowSize, args.correctPvalue, args.curveMethod, args.curveThreshold, args.length, args.remove))

    return args
  
  
  
