__author__ = 'yqb7'
import math, Tools
def heuristicLikelihoodCal(sample1Name, sample2Name, sample1FreqHeuristic, sample2FreqHeuristic,
                           bayesianFreq, lociname2level):

    epsilon = 0.3072
    nloci = len(lociname2level)
    delta = 0
    delta_holder = {}
    loci1Ploidy = Tools.findlocil2Ploidy(lociname2level)
    for lociIndex in range(nloci):
        # ---find necessary info for the likelihood calculation-------
        lociName = lociname2level[lociIndex]
        delta_holder[lociName] = []
        # find the hap types for one loci for both samples and sum the logValue
        hapsValueDict1 = Tools.findHapsValueDictHeuris2(sample1FreqHeuristic, lociName)
        hapsValueDict2 = Tools.findHapsValueDictHeuris2(sample2FreqHeuristic, lociName)
        hapsSet1 = Tools.findHapCount(hapsValueDict1)
        hapsSet2 = Tools.findHapCount(hapsValueDict2)
        hapsCount1 = len(hapsSet1)
        hapsCount2 = len(hapsSet2)
        minimumHapCount = Tools.findMinimumValue(hapsCount1, hapsCount2)
        hapsSumNumber = hapsCount1 + hapsCount2    #x
        sharedalleles = Tools.findNsharedallelesHeuris(hapsValueDict1, hapsValueDict2)
        nsharedalleles = len(sharedalleles)
        # H-nu requires Bayesian Frequency information to do the calculation
        bayesianFreqList = Tools.findHeurH_nu2(bayesianFreq, lociName)[1]
        # sum of log of probability of sample1
        # n is the minimum number of the hap,
        # x is the sum hap number,
        # m is the constant 1,
        # y is the shared locus number between v1 and v2
        # w
        # k
        w = 0
        delta_nu_raw = 0
        delta_nu = 0
        z = 0
        jj = 0
        P_nu = 1
        H_nu = Tools.findHeurH_nu2(bayesianFreq, lociName)[0]
        k = 0
        delta_ex_raw = 0
        delta_ex = 0
        P_ex = 1
        m = 1
        if hapsCount1 > 0 and hapsCount2 > 0:
            if lociName not in loci1Ploidy:  # ploid = 2
                if minimumHapCount > 1:
                    w = hapsSumNumber
                if minimumHapCount == 1 and hapsSumNumber == 2:
                    w = 4
                if minimumHapCount == 1 and hapsSumNumber > 2:
                    w = 1 + hapsSumNumber
                if m == 1:
                    jj = 2
                if m > 1:
                    jj = m
                if minimumHapCount == 1 and nsharedalleles == 1 and hapsSumNumber > 2:
                    z = z + 3
                if minimumHapCount > 1 and jj >= nsharedalleles: #(n>1)*(jj>=y)
                    z = z + 2 * nsharedalleles
                if minimumHapCount > 1 and nsharedalleles > jj:  #(n > 1)* (y > jj)
                    z = z + 2 * jj
                if minimumHapCount == 1 and hapsSumNumber == 2 and nsharedalleles == 1: #(n==1) * (x==2) * (y==1))
                    z = z + 2 * jj
                if nsharedalleles == 0:
                    delta_nu_raw = w + conditionForDeltaNuRaw(jj, z)
                    P_nu = 1
                if nsharedalleles > 0:
                    delta_nu_raw = 2 * jj + conditionForDeltaNuRaw(jj, z)
                    for asharedallele in sharedalleles:
                        P_nu = (hapsValueDict1[asharedallele]) ** 2
                k = P_nu
                if delta_nu_raw > 0:
                    delta_nu = H_nu * delta_nu_raw * k
                if delta_nu_raw == 0:
                    delta_nu = H_nu * P_nu * k
                delta = delta_nu
            if lociName in loci1Ploidy:  # ploid = 1
                if nsharedalleles == 0:  # y =  nsharedalleles: the shared locus number between v1 and v2
                    delta_ex_raw = 2 * hapsSumNumber
                    P_ex = 1
                    k = 1
                if nsharedalleles > 0:
                    for asharedallele in sharedalleles:
                        P_ex = (hapsValueDict1[asharedallele]) ** 2
                        k = P_ex
                if delta_ex_raw > 0:
                    delta_ex = H_nu * delta_ex_raw * k
                if delta_ex_raw == 0:
                    delta_ex = H_nu * P_ex * k
                delta = delta_ex
        else:
            delta = 0
        delta_holder[lociName].append(delta)

    return delta_holder

def conditionForDeltaNuRaw(jj, z):
    sum = 0
    for avalue in range (2*jj):
        avalue = avalue + 1
        if avalue == z:
            sum = sum - avalue
    return sum