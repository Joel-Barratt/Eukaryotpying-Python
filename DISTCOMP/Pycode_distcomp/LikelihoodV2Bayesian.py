__author__ = 'yqb7'
import math, Tools

def likelihoodCal(sample1Name, sample2Name, sample1Freq, sample2Freq,lociname2level):
    nloci = len(lociname2level)
    likelihoodinfor = []
    likelihoodhold = {}
    epsilon = 0.3072
    loci2Ploidy = Tools.findlocil2Ploidy(lociname2level)

    for lociIndex in range(nloci):
        # ---find necessary info for the likelihood calculation-------
        likelihood0 = 0
        likelihood1 = 0
        likelihood2 = 0
        lociname = lociname2level[lociIndex]
        # find the hap types for one loci for both samples and sum the logValue
        likelihoodTmp = []
        valueLH1 = 0
        valueLH2 = 0
        # sum of log of probability of sample1
        hapsCount1, hapsValueList1, hapsValueDict1 = findValueList(sample1Freq, lociname2level[lociIndex])
        sumhaps1 = 0
        # ploid is 2
        if len(hapsValueList1) == 1 and lociname not in loci2Ploidy:
            for haps1index in range(len(hapsValueList1)):
                sumhaps1 = sumhaps1 + math.log(hapsValueList1[haps1index])
            sumhaps1 = sumhaps1 * 2
        if len(hapsValueList1) > 1 and lociname not in loci2Ploidy:
            for haps1index in range(len(hapsValueList1)):
                sumhaps1 = sumhaps1 + math.log(hapsValueList1[haps1index])
        if lociname in loci2Ploidy:
            for haps1index in range(len(hapsValueList1)):
                sumhaps1 = sumhaps1 + math.log(hapsValueList1[haps1index])
        if len(hapsValueList1) == 0:
            sumhaps1 = 0
        # sum of log of probability of sample2
        hapsCount2, hapsValueList2, hapsValueDict2 = findValueList(sample2Freq, lociname2level[lociIndex])
        sumhaps2 = 0
        if len(hapsValueList2) == 1 and lociname not in loci2Ploidy:
            for haps2index in range(len(hapsValueList2)):
                sumhaps2 = sumhaps2 + math.log(hapsValueList2[haps2index])
            sumhaps2 = sumhaps2 * 2
        if len(hapsValueList2) > 1 and lociname not in loci2Ploidy:
            for haps2index in range(len(hapsValueList2)):
                sumhaps2 = sumhaps2 + math.log(hapsValueList2[haps2index])
        if lociname in loci2Ploidy:
            for haps2index in range(len(hapsValueList2)):
                sumhaps2 = sumhaps2 + math.log(hapsValueList2[haps2index])
        if len(hapsValueList2) == 0:
            sumhaps2 = 0
        hapsValueList1.sort()  # get a ascending list
        hapsValueList2.sort()  # get a ascending list
        listtmp = hapsValueList1 + hapsValueList2
        listtmp.sort()
        # find the shared alleles and unshared values
        nsharedalleles = Tools.compareDict(hapsValueDict1, hapsValueDict2)
        if hapsCount1 != 0 and hapsCount2 != 0:
            valueLH1 = Tools.compareDictUnsharedLH1V2(hapsValueDict1, hapsValueDict2, lociname, loci2Ploidy)
            valueLH2 = Tools.compareDictUnsharedLH2V2(hapsValueDict1, hapsValueDict2, lociname, loci2Ploidy)
        # ---end of info search-----------------
        # calculate likelihood0, likelihood1, likelihood2
        # this case include as 1). sample1 = sample2 & 2). haps number is equal & 3). ploid is 1 or 2
        if hapsCount1 == 0 or hapsCount2 == 0:
            likelihood0 = 0
            likelihood1 = 0
            likelihood2 = 0
        if lociIndex == 11:
            a = 0
        if hapsCount1 > 0 and hapsCount2 > 0:
            likelihood0 = sumhaps1 + sumhaps2
            if nsharedalleles == 0 and lociname not in loci2Ploidy:
                likelihood1 = math.log(epsilon * listtmp[0])
                likelihood2 = math.log((epsilon * listtmp[0]) ** (2 - nsharedalleles))
            if nsharedalleles == 0 and lociname in loci2Ploidy:
                likelihood1 = math.log(epsilon * listtmp[0])
                likelihood2 = likelihood1
            if nsharedalleles != 0:
                likelihood1 = valueLH1 + sumhaps2
                likelihood2 = valueLH2 + sumhaps2
                if valueLH2 == 1 and lociname not in loci2Ploidy:
                    likelihood2 = math.log((epsilon * listtmp[0]) ** (2 - nsharedalleles))
                if lociname in loci2Ploidy:
                    likelihood2 = likelihood1

        # this case include as 1). sample1 != sample2 & 2). ploid is 1 or 2
        # case1 sample1 has one hap , sample2 has more than 1 haps
        # case1 sample2 has one hap , sample1 has more than 1 haps
        # case1 sample1 has more than one haps, sample2 has more than 1 haps
        likelihoodTmp.append(likelihood0)
        likelihoodTmp.append(likelihood1)
        likelihoodTmp.append(likelihood2)
        likelihoodhold[lociname] = likelihoodTmp
    likeallloci = calLikelihood(likelihoodhold)
    likelihoodinfor.append(likeallloci)
    likelihoodinfor.append(likelihoodhold)
    return likelihoodinfor
def findValueList(sample1Freq, lociName):
    hapsCount = 0
    hapsValueList = []
    hapsValueDict = {}
    for haps in sample1Freq.keys():
        if haps.__contains__(lociName):
            hapsCount = hapsCount + 1
            hapsValueList.append(sample1Freq.get(haps))
            hapsValueDict[haps] = sample1Freq.get(haps)
    return hapsCount, hapsValueList, hapsValueDict
def calLikelihood(likelihoodhold):
    sumCol0 = 0
    sumCol1 = 0
    sumCol2 = 0
    for likelihoodlist in likelihoodhold.values():
        sumCol0 = sumCol0 + likelihoodlist[0]
        sumCol1 = sumCol1 + likelihoodlist[1]
        sumCol2 = sumCol2 + likelihoodlist[2]
    sumlikallloci = math.exp(sumCol0) + math.exp(sumCol1) + math.exp(sumCol2)
    exp0 = math.exp(sumCol0)/sumlikallloci
    exp1 = math.exp(sumCol1)/sumlikallloci
    exp2 = math.exp(sumCol2)/sumlikallloci
    likeallloci = exp0 * 0 + exp1 * 1 + exp2 * 2
    return likeallloci