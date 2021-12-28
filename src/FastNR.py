import os
import pysam
import sys

import math
import numpy as np
import multiprocessing
import timeit
from functools import partial 
import gzip
from .param import *
from .statis import *

global controlTotal,treatTotal
     

def main ():
    """This is the main function of Fast-NR"""
    #Get the defined parameter.
    defPar = defOpt()
    args = option(defPar)
    #Start to read in inputs.
    print("Step 1: Read in input files.")
    
    #The chromosome size file.
    print("Read in chromosome size file.")   
    
    genome = GenomeIn(args.genome)
    chrSizeList = genome.genomeInf() #all chromosome size store in a dict
    pThres = args.pvalue
    correctM = args.correctPvalue
    similarity = args.curveMethod
    bw = args.windowSize
    
    print("Read in control file.")
    
    if args.cFormat == "bed":               
        controlPro = BEDProcess(args.control,args.length)
        control=controlPro.Genome_Frag_all_pos(chrSizeList)
    elif args.cFormat == "bam":
        controlPro = BAMProcess(args.control,args.length)
        control=controlPro.Genome_Frag_all_pos(chrSizeList)

    print("Read in treatmant file.")
    if args.tFormat == "bed":
        treatmentPro = BEDProcess(args.treatment,args.length) 
        treat=treatmentPro.Genome_Frag_all_pos(chrSizeList)
    elif args.tFormat == "bam":
        treatmentPro = BAMProcess(args.treatment,args.length)
        treat=treatmentPro.Genome_Frag_all_pos(chrSizeList)
    
    #get the scale value of treatment and control.
    scale = (float)(treat.total)/(control.total)
    tMeanLen = treat.length/treat.total
    cMeanLen = control.length/control.total
    controlTotal = control.total*scale
    treatTotal = treat.total
    print("Total fragment in treatment is: [%d], with mean length: [%0.2e]" % (treat.total,tMeanLen))
    print("Total fragment in control is: [%d], with mean length: [%0.2e]" % (control.total,cMeanLen))
    print("The scale of treat to control is [%0.2e]", scale) 
    #get the background count in 
    bgTreat = bw*treat.total/genome.sumlen
    bgControl = bw*control.total*scale/genome.sumlen
    
    #The min control is the 
    #minControl = nb_cdf_inv_i(args.pvalue,treat.total,bgTreat)
    #The min treatment is 0
    if args.pvalue > 1e-6 :
        minControl = math.ceil(nb_cdf_inv_i(args.pvalue,treat.total,0))
    else: #if p value less than 1e-6, then use 1e-6 as input to get minControl. This will not affect final result, only filter less sites when compute p value.
        minControl = math.ceil(nb_cdf_inv_i(1e-6,treat.total,0))
     
    #maxTreat = nb_cdf_inv_k(args.pvalue,treat.total,bgControl,control.total*scale)
    #the max control is max in control pos
    #maxTreat = math.floor(nb_cdf_inv_k(args.pvalue,treat.total,max(control.pos.values()[0][0]),control.total*scale))-1 #This is for python 2.x. max(control.pos.values()) for python 3.x
     
    #calculate the p value

    candiNRregion = {}
     
    if args.curveMethod == "Euclidean" or args.curveMethod =="Cosine":
        percentile = 100*(1-args.curveThreshold)
    elif args.curveMethod == "Pearson" or args.curveMethod == "Gradient": #in the correlation part, smaller value means lower similarity
        percentile = 100*(args.curveThreshold)
        
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)       
    
    ###################################################################################################################    
    #####Split the numpy array according to the value of treat and control.
    #first, sort the pos by treat and control count, and return the index of pos
    par_NB_pvalue = partial(cal_NB_pvalue, treatTotal,controlTotal)
    pvalueNum = 0
    for chromosome in genome.get_chr_names():
        filterPArr = []
        #filter the position according to the min control and max treatment cover.
        #The max treatment much smaller than the largest control generated count according to threshold p value.
        maxChrTreat = math.floor(nb_cdf_inv_k(args.pvalue,treat.total,np.amax(treat.pos[chromosome][0]),control.total*scale))-1        
        #The min control much larger than the smallest treatment (0 in all chromosome) generated count according to threshold p value. 
        minChrControl = minControl        
        #filter
        print("Filter the chromosome [%s] with treatment <= [%d] and control >= [%d] and control-treatment >= [%d]. " % (chromosome, maxChrTreat, minChrControl,minChrControl))
        # pos = np.where((control.pos[chromosome]>=minChrControl)&(treat.pos[chromosome]<=maxChrTreat)) #this only remove little sites, use following command instead. For the minus of control and treat also need to larger than the smallest control count.
        chrPos = np.where((control.pos[chromosome]>=minChrControl)&(treat.pos[chromosome]<=maxChrTreat)&(control.pos[chromosome][0]-treat.pos[chromosome][0]>minChrControl))
        #get unique pos
        if len(chrPos[1]) <=0 :
            print("No region find in this chromosome [%s]" %(chromosome))
            continue    
        else:        
            tC = np.take(treat.pos[chromosome][0],chrPos[1])
            cC = np.take(control.pos[chromosome][0],chrPos[1])
            items = np.column_stack((tC,cC*scale,chrPos[1])) 
            sortItems = items[np.lexsort(np.transpose(items)[::-1])]
            # #or 
            # sortItems = items[items[:, 2].argsort(kind='mergesort')]  # sort by pos
            # sortItems = sortItems[sortItems[:, 1].argsort(kind='mergesort')]  # sort by control
            # sortItems = sortItems[sortItems[:, 0].argsort(kind='mergesort')]  # sort by treat  
            #second, the split this numpy by treat & control pair value.
            diff = np.vstack((np.diff(sortItems[:,0]),np.diff(sortItems[:,1])))
            ItemSplit = np.split(sortItems, np.where((diff[0]!=0) | (diff[1]!=0))[0]+1)
            #Third, get the first value (that is the unique treat & control pair) and store it in a numpy array.
            uniqPair = np.array([item[0] for item in ItemSplit])
            #Fourth, get the pvalue of uniqPair.
            #par_NB_pvalue = partial(cal_NB_pvalue, treatTotal,controlTotal)
            print("Calculate the p value.")
            pvalueArray = np.array(pool.map(par_NB_pvalue,uniqPair))
            #Fifth, filter the pvalue
            print("Filter the p value.")
            filterInd = np.where(pvalueArray[:,0]<pThres)
            filterPArr = pvalueArray[filterInd]
            #filter the pos of all filtered value.
            if len(filterInd[0]) == 0:
                print("No region find in this chromosome [%s]" %(chromosome))
                continue
            else: 
                filterPos = [ItemSplit[i] for i in filterInd[0]]
                filterPosList = np.array(list(np.concatenate([item[:,2] for item in filterPos]).flat)) #flatten the pos list
                filterPosList = np.sort(filterPosList)
                #Sixth, extend the pos in each filterPos to bw, and calculate the region similarity.   
                #Put the all filtered positions in a shared memory.
                size = genome.get_size_of_chr(chromosome)[1] 
                #store all pvalue in chromosome each position, less than threshold mark as zero, others mark as 1
                chrPvalue=np.ones((len(treat.pos[chromosome][0]),1))
                np.put(chrPvalue,filterPosList.astype(int),0)

                RegionR=[]
                print("Extend the filter pos to fixed window size.")
                for i in range(0,len(filterPos)):
                    #extend the sites. Extend the pos in each filterPos to bw.  
                    exStart = filterPos[i][:,2]-bw/2
                    exStart[np.where(exStart<0)] = 0
                    exEnd = filterPos[i][:,2]+bw/2
                    exEnd[np.where(exEnd>size)] = size
                    for j in range(0,len(exEnd)):
                        regionS,regionE = int(exStart[j]),int(exEnd[j])
                        if len(chrPvalue[regionS:regionE+1,0][np.where(chrPvalue[regionS:regionE+1,0]<pThres)]) >= (regionE-regionS)*3/4:
                            RegionR.append([regionS,regionE,filterPos[i][j][2],filterPos[i][j][0],filterPos[i][j][1],filterPArr[i][0]])
                #Then remove the overlap.
                if len(RegionR) == 0:
                    print("No region find in this chromosome [%s]" %(chromosome))
                    continue     
                else:           
                    RegionR = np.array(RegionR)
                    sortRegionR = RegionR[np.argsort(RegionR[:,0])]
                    keepInterval = remove_intervals(sortRegionR)
                    similarArr = np.ones((1,len(keepInterval))) 
                    #Then get the similar.
                    if percentile == 100 or percentile == 0:       #if use all regions, then, do not use similarity compute.
                        candiNRregion[chromosome] = np.column_stack((keepInterval[:,0],keepInterval[:,1],keepInterval[:,2],keepInterval[:,3],keepInterval[:,4],(keepInterval[:,3]/keepInterval[:,4]),keepInterval[:,5]))
                        pvalueNum += len(candiNRregion[chromosome])                    
                    else:
                        print("Calculate the curve similarity.")
                        i=0
                        for tmpR in keepInterval: 
                            regionS,regionE = int(tmpR[0]),int(tmpR[1])
                            treatPoint = treat.get_interval_cover(chromosome,regionS,regionE,chrSizeList)
                            controlPoint = control.get_interval_cover(chromosome,regionS,regionE,chrSizeList)
                            controlPoint = controlPoint*scale
                            if args.curveMethod == "Pearson": #when there is little count in a bin, this similarity is also low, so, maybe this is not suitable for less count.
                                similar = PearsonCoef(treatPoint,controlPoint)[0,1]
                            elif args.curveMethod == "Euclidean":
                                #similar = EuclideanDis(len(treatPoint),treatPoint,controlPoint)
                                similar = euclidean(treatPoint,controlPoint)
                            elif args.curveMethod == "Cosine": #new, 2020.10.3
                                treatPoint = treatPoint + 0.01  #add, 2021,4,12, in case of zero vector 
                                controlPoint = controlPoint + 0.01  #add, 2021,4,12, in case of zero vector 
                                similar = cosine(treatPoint,controlPoint)
                            elif args.curveMethod == "Gradient":
                                if len(range((int)(regionS),(int)(regionE)+1)) > len(treatPoint): #for the region at end of genome
                                   tT = np.row_stack((range((int)(regionS),(int)(regionE)),treatPoint))
                                   tC = np.row_stack((range((int)(regionS),(int)(regionE)),controlPoint))
                                else:
                                   tT = np.row_stack((range((int)(regionS),(int)(regionE)+1),treatPoint))
                                   tC = np.row_stack((range((int)(regionS),(int)(regionE)+1),controlPoint))                           
                                slope_treatment = SlopeLine(tT)
                                slope_control = SlopeLine(tC)
                                similar = PearsonCoef(slope_treatment,slope_control)[0,1]
                            similarArr[0][i] = similar  
                            i += 1
                        keepInterval = np.column_stack((keepInterval,similarArr[0]))
                        #seventh, get the percentile threshold of similarity.
                        #filter by similarity
                        percentVec = np.percentile(keepInterval[:,6],(percentile))
                        if args.curveMethod == "Euclidean" or args.curveMethod == "Cosine":
                            filterR = keepInterval[np.where(keepInterval[:,6] > percentVec)]
                        elif args.curveMethod == "Pearson" or args.curveMethod == "Gradient":
                            filterR = keepInterval[np.where(keepInterval[:,6] < percentVec)]   
                        #sort the regions, then remove the overlapped region, only keep the smallest p value.
                        #candiNRregion[chromosome] = filterR
                        pvalueNum += len(filterR)
                        candiNRregion[chromosome] = np.column_stack((filterR[:,0],filterR[:,1],filterR[:,2],filterR[:,3],filterR[:,4],(filterR[:,3]/filterR[:,4]),filterR[:,5]))

    pool.close()
    pool.join() 

    #step 5. correct the filtered p value with signed method.   
    if pvalueNum <=0 :
        print("No region find!")
        sys.exit()
    else:    
        correctRegion = {}
        totalP = 0
        totalRP = 0
        if args.correctPvalue == "Bonferroni":
            for chromosome in candiNRregion.keys():           
                totalP += len(candiNRregion[chromosome]) #use filtered p value to get the length.
            for chromosome in candiNRregion.keys():
                correctRegion[chromosome] = np.column_stack((candiNRregion[chromosome],np.empty((len(candiNRregion[chromosome]),1))))
                correctRegion[chromosome][:,7] = candiNRregion[chromosome][:,6]*totalP
                #larger than one will go to one
                correctRegion[chromosome][:,7][np.where(correctRegion[chromosome][:,7] > 1)] = 1 
        elif args.correctPvalue == "BH":
            #use the same method that used in R p.adjust(,method="BH")
            #first, build a sorted p value with index of they original position.
            pD = {}
            i = 0
            for chromosome in candiNRregion.keys():
                tmpLen = len(candiNRregion[chromosome])
                tmpD = dict(list(zip(range(i,i+tmpLen), candiNRregion[chromosome][:,6]))) #test from here.
                pD.update(tmpD) #it may takes few times.
                i += tmpLen
            #sort by value,if sort by key sorted(pD.items(), key=lambda d: d[0])
            sortPD = sorted(pD.items(), key=lambda d: d[1])
            #from the biggest to smallest, calculate the BH correct p value
            i = len(sortPD)
            n = len(sortPD)
            BHpvalueD = {}
            prevP = (sortPD[i-1][1]*n)/i
            while i> 0:
                tmpCP = min((sortPD[i-1][1]*n)/i,prevP)
                prevP = tmpCP
                BHpvalueD[sortPD[i-1][0]] = tmpCP #give orignal index to the corrected p value
                i -= 1
            #return the corrected p value to all p value.
            i = 0
            for chromosome in candiNRregion.keys():
                correctRegion[chromosome] = np.column_stack((candiNRregion[chromosome],np.empty((len(candiNRregion[chromosome]),1))))
                for tn in range(0,len(candiNRregion[chromosome][:,6])):
                    correctRegion[chromosome][tn,7] = BHpvalueD[i]
                    i += 1
        #filter the corrected p value.
        totalR = 0 
        for chromosome in correctRegion.keys():
            correctRegion[chromosome] = correctRegion[chromosome][np.where(correctRegion[chromosome][:,7] < pThres)]   #second, filter the sites with fixed p value
            
        finalNR={}
        #step 6. if set remove slope region, the last step. add in 2020.10.14. not sure if it necessary.
        TotalNR=0
        if args.remove:
            for chromosome in correctRegion.keys():
                filterFR = []
                for region in correctRegion[chromosome]:
                    regionS,regionE = int(region[0]),int(region[1])
                    treatPoint = treat.get_interval_cover(chromosome,regionS,regionE,chrSizeList)
                    controlPoint = control.get_interval_cover(chromosome,regionS,regionE,chrSizeList)
                    #get eleven point, all point seems not usfull
                    pointV = np.arange(0,1.1,0.1)*(regionE-regionS)
                    tmpCTE = np.take(treatPoint,pointV.astype(int))
                    tmpCCE = np.take(controlPoint,pointV.astype(int))
                    minusTR = tmpCTE[1:11,]-tmpCTE[0:10,]
                    minusCR = tmpCCE[1:11,]-tmpCCE[0:10,]
                    #if all TR and CR have similar slope like positive or negative tend, then remove it.
                    lenTP = len(np.where(minusTR>0)[0])
                    lenCP = len(np.where(minusCR>0)[0])
                    lenTN = len(np.where(minusTR<0)[0])
                    lenCN =len(np.where(minusCR<0)[0])
                    if (lenTP>=9 and lenCP>=9) or (lenTN>=9 and lenCN>=9): 
                        if np.corrcoef(minusTR,minusCR)[0][1] >= 0.5: #and the correlation is high
                            continue
                        else:
                            filterFR.append(region)
                    else:
                        filterFR.append(region)
                finalNR[chromosome]=np.array(filterFR)  
                TotalNR += len(filterFR)
                print("Find %d regions in chromosome %s." % (len(filterFR),chromosome))
        else:
            finalNR = correctRegion
            for chromosome in finalNR.keys():
                TotalNR += len(finalNR[chromosome])
                print("Find %d regions in chromosome %s." % (len(finalNR[chromosome]),chromosome))        
        print("Total %d regions in this sample." % TotalNR)
        #store the final region in a file.
        # negRe = negRegion()
        chrs = sorted(finalNR.keys())
        text = ""

        for chromosome in chrs:
            for line in finalNR[chromosome]:
                chrom = chromosome
                (start,end,summit,treatCount,controlCount,foldChange,pValue,correctpValue)=line
                #negRe.__Region_line(chromosome,start,end,treatCount,controlCount,foldChange,pValue,correctpValue)
                text+= "%s\t%d\t%d\t%d\t%d\t%d\t%f\t%e\t%e\n" % (chromosome,start,end,summit,treatCount,controlCount,foldChange,pValue,correctpValue)

        print("Output result in file %s." % args.outfile)
        outBED = open(args.outfile,"w")
        outBED.write(text)
        outBED.close()   

class BAMProcess():
    """Read in bam, and build a list.
    """    
    def __init__ (self,bfile,extend):
        self.bf = pysam.AlignmentFile(bfile, 'rb')
        self.extend=extend        
    def Genome_Frag_all_pos (self,chrSizeList):
        """Build genome cover from all lines, return a Cover object.
        """
        cover = Cover()
        #the process of paired or single end data.
        #first, judge if it is paired or single end.
        for tmpline in self.bf:
            
            fChr=""
            fStart,fEnd=0,0    

            #for single end reads        
            if not tmpline.is_paired and not tmpline.is_unmapped and not tmpline.is_duplicate: # single mapped unique read.
                fChr = tmpline.reference_name
                fStart = tmpline.reference_start+1
                fEnd = tmpline.reference_end+1   
                fStrand = "-" if tmpline.is_reverse else "+"
                #if used extend parameter, then, extend from 5' to 3' of fragment.
                if self.extend >0 :
                    if tmpline.is_reverse:
                        fStart = 1 if fStart <1 else fEnd-self.extend 
                    else:
                        fEnd = fStart+self.extend              
                        if fEnd > chrSizeList[tmpline.reference_name][1]:
                            fEnd = chrSizeList[tmpline.reference_name][1]   
                cover.add_frag_all_pos(fChr,fStart,fEnd,fStrand,chrSizeList) 
            #for paired end reads. Only continue when read1 and read2 all obtained.
            # if not tmpline.is_paired or not tmpline.is_proper_pair or tmpline.mate_is_unmapped or tmpline.is_duplicate:
                # continue      
            elif tmpline.is_paired and tmpline.is_proper_pair and not tmpline.mate_is_unmapped and not tmpline.is_duplicate: # paired mapped unique read. 
                #obtain both read1 and read2, the go to judge
                if tmpline.is_read2:
                    read2 = tmpline
                else:
                    read1 = tmpline
                    read2 = None
                    continue
                if not read1 is None and not read2 is None: 
                    if read1.query_name == read2.query_name and read1.reference_name == read1.reference_name : # find a pair. 
                        fChr = read1.reference_name
                        if(read1.reference_start>=read2.reference_start):
                            fStart,fEnd=read2.reference_start,read1.reference_end
                        else:
                            fStart,fEnd=read1.reference_start,read2.reference_end                        
                        #if used extend parameter, then, extend from 5' to 3' of fragment. judge, if the extend frag out of chromosome size.
                        if read2.is_reverse and not read1.is_reverse:
                                fStrand = "+"
                                if self.extend >0 :
                                    fEnd = fStart + self.extend
                                    if fEnd > chrSizeList[read1.reference_name][1]:
                                        fEnd = chrSizeList[read1.reference_name][1]
                        elif not read2.is_reverse and read1.is_reverse:
                                fStrand = "-"
                                if self.extend >0 :
                                    fStart = 1 if fStart <1 else fEnd-self.extend 
                        else:
                                print("Warning: read 1 and read 2 not in a forward-reverse type! Discard it [%s]." % read1.query_name)
                                continue
                        #then call function to obtain the coverage of this fragment   
                        cover.add_frag_all_pos(fChr,fStart,fEnd,fStrand,chrSizeList) 
                    elif read1.query_name != read2.query_name: 
                        print("Error: not sorted by read names! You can use \"samtools sort -n\" to sort it.")
                        sys.exit()
                    elif read1.reference_name != read2.reference_name: 
                        print("Error: read 1 and read 2 not mapped in the same chromosome!")
                        sys.exit()
                
        self.bf.close()         
#       self.fhd.seek(0)
        return cover
        
  
class BEDProcess():
    """Read in BED, and build a list.
    """
    def __init__ (self,path,extend):  
        # try gzip first
        f = gzip.open(path, 'r')
        try:
            f.read(10)
        except IOError:
            # not a gzipped file
            f.close()
            f = open(path, 'r', -1)
        else:
            f.seek(0)
        self.fhd=f
        self.extend=extend
    def __Read_bed_line (self, tmpline):
        tmpline = tmpline.rstrip()
        spt = tmpline.split('\t')
        if len(spt) < 6 : # if no strand information include 6 column is ".", then default +
            tmpStrand = 0 
        else:
            if spt[5] == "+":
                tmpStrand = 0 
            elif spt[5] == "-":
                tmpStrand = 1
            elif spt[5] == ".":
                tmpStrand = 0    
            else:
                print("Error: Strand information not correct in fragment: \"%s\",\"%s\"" % (tmpline,spt[5]))
                sys.exit(1)
        return (spt[0],int(spt[1]),int(spt[2]),tmpStrand)
    def Genome_Frag_all_pos (self,chrSizeList):
        """Build genome cover from all lines, return a Cover object.
        """
        cover = Cover()
        for tmpline in self.fhd:
            (chromosome,start,end,strand) = self.__Read_bed_line(tmpline)
            if not chromosome or start is "" or end is "":   #except the possible that when it is equal to zero as start pos. 
                continue
            #If extend, then extend the start, and end.
            if self.extend >0 :
                if strand == "-":
                    start = 1 if start <1 else end-self.extend 
                else:
                    end = start+self.extend              
                    if end > chrSizeList[chromosome][1]:
                        end = chrSizeList[chromosome][1] 
            cover.add_frag_all_pos(chromosome,start,end,strand,chrSizeList)
        self.fhd.seek(0)
        return cover
   
class Cover:
    """the coverage of each bp in genome by fragment,stored in a dict.
    This dic is stored and organized by sequence names. 
    And can be sorted.
    """
    def __init__ (self):
        """       
        """    
        self.pos = {}
        self.total = 0 
        self.length = 0
    def add_frag_all_pos(self,chromosome,start,end,strand,chrSizeList):
        """Add fragment information to each chromosome position   
        size -- the size of this chromosome, which will help to build the zero numpy array
        """
        size = chrSizeList[chromosome][1] 
        #if not self.pos.has_key(chromosome):
        if chromosome not in self.pos:
            self.pos[chromosome] = np.zeros((1,size))
        self.pos[chromosome][0][start:(end+1)] +=1
        self.total += 1 #add total fragment for scale
        self.length += (end-start+1)
    def get_interval_all_count(self,chromosome,start,end,chrSizeList):
        """Return sum count and mean count of position in chromosome intervals
        """
        size = chrSizeList[chromosome][1] 
        #if self.pos.has_key(chromosome):
        if chromosome in self.pos:
            if start > size or end > size:
                raise Exception("Beyond the range of this chromosome (%s) in Cover object!\n" % (chromosome)) 
            else:
                tmpCount = np.sum(self.pos[chromosome][0][start:(end+1)])
                tmpMean = np.mean(self.pos[chromosome][0][start:(end+1)])
                return (tmpCount,tmpMean)          
        else:
            raise Exception("No such chromosome name (%s) in Cover object!\n" % (chromosome)) 
    def get_interval_cover(self,chromosome,start,end,chrSizeList):
        """Return all cover in chromosome intervals
         """
        size = chrSizeList[chromosome][1] 
        #if self.pos.has_key(chromosome):
        if chromosome in self.pos:
            if start > size or end > size:
                raise Exception("Beyond the range of this chromosome (%s) in Cover object!\n" % (chromosome)) 
            else:
                return self.pos[chromosome][0][start:(end+1)]           
        else:
            raise Exception("No such chromosome name (%s) in Cover object!\n" % (chromosome))
    def get_pos_cover(self,chromosome,position,chrSizeList):
        """Return cover in one position in chromosome
        """
        size = chrSizeList[chromosome][1] 
        #if self.pos.has_key(chromosome):
        if chromosome in self.pos:
            if pos > size:
                raise Exception("Beyond the range of this chromosome (%s) in Cover object!\n" % (chromosome)) 
            else:
                return self.pos[chromosome][0][position]           
        else:
            raise Exception("No such chromosome name (%s) in Cover object!\n" % (chromosome))            
    def get_chr_cover(self,chromosome):
        """Return coverage to each chromosome 
        """
        #if self.pos.has_key(chromosome):
        if chromosome in self.pos:
            return self.pos[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in Cover object!\n" % (chromosome))


class GenomeIn():
    """Read in genome chromosome size. 
    """
    def __init__ (self,path):
        f = gzip.open(path, 'r')
        try:
            f.read(10)
        except IOError:
            # not a gzipped file
            f.close()
            f = open(path, 'r', -1)
        else:
            f.seek(0)
        self.fhd=f  
        self.chrSize = {} 
        self.sumlen = 0         
    def __Read_genome_line (self, tmpline):
        tmpline = tmpline.rstrip()
        spt = tmpline.split('\t')
        if len(spt) != 2:
            print("Error: The chromosome size file not in a correct format with two column: chromosome length.")
            sys.exit(1)
        #all the start of genome is 1
        else:    
            return (spt[0],1,int(spt[1])) 

    def genomeInf(self):
        for tmpline in self.fhd:
            (chromosome,Start,End) = self.__Read_genome_line(tmpline)
            if not End or not chromosome:
                continue
            self.__add_chr(chromosome,Start,End)
            self.sumlen += End
        self.fhd.seek(0)
        return self.chrSize
        
    def __add_chr (self, chromosome, Start, End):
        """Add a chromosome size to the list.
        """
        #if not self.chrSize.has_key(chromosome):
        if chromosome not in self.chrSize:
            self.chrSize[chromosome] = []
        else:
            print("Error: Have two line of chromosome [%s], please check your chromosome size file!" % chromosome)
            sys.exit(1)
        self.chrSize[chromosome].append(Start)
        self.chrSize[chromosome].append(End)

    def get_size_of_chr (self, chromosome):
        """Return the size of certain chromosome.
        """
        #if self.chrSize.has_key(chromosome):
        if chromosome in self.chrSize:
            return self.chrSize[chromosome]
        else:
            raise Exception("No such chromosome name (%s) in this genome!\n" % (chromosome))
           
    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        l = list(self.chrSize.keys())
        l.sort()
        return l

if __name__ == "__main__":
    sys.exit(main())   


