from scipy.stats import nbinom
import numpy as np
from scipy.spatial.distance import euclidean
from scipy.spatial.distance import cosine 
from scipy.special import nbdtrik
from scipy.special import nbdtri


def PearsonCoef(vA,vB):
    """Pearson correlation coefficient between two vector
    """
    pcc = np.corrcoef(vA, vB)
    return pcc

def SlopePoint(pt1,pt2):
    """The slope between two point (x1,y1),(x2,y2)
    (y2-y1)/(x2-x1)
    """
    tmp_slope=0
    tmp_slope = (pt2[1]-pt1[1])/(pt2[0]-pt1[0])   
    return tmp_slope

def SlopeLine(curveA):
    """The slope of each point in a curve
    """
    slope = np.zeros(len(curveA[0])-1)
    for i in range(0,len(curveA[0])-1):
        slope[i] = SlopePoint(curveA[:,i],curveA[:,i+1])   
    return slope

def negativeBinomail(sampleCount,sampleTotal,ControlCount,ControlTotal):
    """Negative binomail CDF     
    """
    #we need the sampleCount as failure, so, use the following one.
    #need float before the calculate the ratio, otherwise it will return 0.
    pvalue = nbinom.cdf(sampleCount,(sampleTotal-sampleCount),(float)(ControlTotal-ControlCount)/ControlTotal)
    return pvalue

def nb_cdf_inv_k(pvalue,sampleTotal,ControlCount,ControlTotal):
    """the inverse of Negative binomail CDF vs k.
    The maximum number of allowed failures.
    y - is the pvalue
    n - is the target number of successes (here is total sample count, 
    in deed the (sampleTotal-sampleCount), but the sample count is what we need k,
    and total sample count minus the sampleCount will not affect too much).
    p - is probability of success in a single event (here is probability of control,
    that (ControlTotal-ControlCount)/ControlTotal).
    it will return the max k, that smaller than k will give a smaller pvalue.
    it means that sampleCount need smaller than max k, at given ControlCount       
    """
    maxK = nbdtrik(pvalue,sampleTotal,(float)(ControlTotal-ControlCount)/ControlTotal)
    return maxK  

def nb_cdf_inv_i(pvalue,sampleTotal,sampleCount):
    """the inverse of Negative binomial CDF vs p. 
    Probability of success in a single event (float) such that nbdtr(k, n, p) = y.
    k - is the The maximum number of allowed failures (nonnegative int).
    n - is the target number of successes (here is sampleTotal-sampleCount).
    y - is the p value.
    it will return the min count, that larger than m will give a smaller pvalue.
    it means that controlCount need larger than min C, at given sampleCount
    """
    p = nbdtri(sampleCount,sampleTotal-sampleCount,pvalue)
    minC = (1-p)*sampleTotal
    return minC

def cal_NB_pvalue (treatTotal,controlTotal,items):
    """calculate the pvalue in pos of chromosome.    
    """
    pvalue = 1
    (treatCount,controlCount,pos)=items  
    pvalue = negativeBinomail(treatCount,treatTotal,controlCount,controlTotal)
    return (pvalue,treatCount,controlCount,pos)
   
def remove_intervals(intervals): 
    """For overlapped intervals, return the smallest pvalue intervals.
    intervals -- the numpy array, contain all position information. pV,revP,tcount,ccount,pos,start,end
    """
    keep = np.empty((1,len(intervals[0])))
    keep[0] = intervals[0]
    for current in intervals:
        previous = keep[-1]
        if current[0] <= previous[1]:## if overlap 
            if current[5] < previous[5]: ##if p value is smaller, then move to this, else, keep the original one
                keep[-1] = current
            else:   ## if overlap, but p value is larger, then skip this interval continue,keep the original one.
               continue 
        else: ## if no overlap, then add to end of keep
            keep = np.vstack ((keep, current) )
    return keep
    
