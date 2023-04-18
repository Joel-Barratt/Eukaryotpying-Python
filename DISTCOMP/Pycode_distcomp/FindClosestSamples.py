__author__ = 'yqb7'
import random, Tools
from orderedset import OrderedSet
# This function is to find the closest samples for the sample with missing loci
def findClosestSample(cleanSampleHFreq):
    sharedHap = {}
    sharedComHap = {}
    lociname2level = random.choice(list(cleanSampleHFreq.values())).lociname2level
    samplesMissLoci, samplesCompLociList = findSampleMisLoci(cleanSampleHFreq, lociname2level)
    for IDindex1 in range(len(samplesMissLoci)):
        pattern1_double = doubleSingleHap(cleanSampleHFreq[samplesMissLoci[IDindex1]].codedPattern, lociname2level)
        pattern1 = OrderedSet(pattern1_double)
        for IDindex2 in range(IDindex1 + 1, len(samplesMissLoci)):
            if IDindex2 < len(samplesMissLoci):
                pattern2 = OrderedSet(cleanSampleHFreq[samplesMissLoci[IDindex2]].codedPattern)
                samplePair = str(samplesMissLoci[IDindex1]) + '__' + str(samplesMissLoci[IDindex2])
                sharedhaptmp = pattern1.intersection(pattern2)
                sharedHap[samplePair] = sharedhaptmp

    # The result will be used to fill in the "0"s in the pairwise distance.
    for IDindex1 in range(len(samplesMissLoci)):
        pattern1 = OrderedSet(cleanSampleHFreq[samplesMissLoci[IDindex1]].codedPattern)
        for IDindex2 in range(len(samplesCompLociList)):
            pattern2 = OrderedSet(cleanSampleHFreq[samplesCompLociList[IDindex2]].codedPattern)
            samplePair = str(samplesMissLoci[IDindex1]) + '__' + str(samplesCompLociList[IDindex2])
            sharedComhaptmp = pattern1.intersection(pattern2)
            if len(sharedComhaptmp) == len(pattern1):
                sharedComHap[samplePair] = sharedComhaptmp
    realsharedComHap = checkRealInclude(sharedComHap, cleanSampleHFreq, lociname2level)
    return realsharedComHap, samplesMissLoci, samplesCompLociList

def findSampleMisLoci(cleansampleHFreq, lociname2level):
    lociIDlist = []
    sampleMissLociList = []
    sampleCompLociList = []
    for lociIndex in range(len(lociname2level)):
        lociIndex += 1
        lociIDlist.append(lociIndex)
    for sampleID, acontent in cleansampleHFreq.items():
        samplelociList = OrderedSet()
        for acodedPattern in acontent.codedPattern:
            samplelociList.add(acodedPattern.split('_')[0])
        if len(samplelociList) < len(lociname2level):
            sampleMissLociList.append(sampleID)
        if len(samplelociList) == len(lociname2level):
            sampleCompLociList.append(sampleID)
    return sampleMissLociList, sampleCompLociList

def checkRealInclude(sharedComHap, cleanSampleHFreq, lociname2level):
    realSharedComp = {}
    for asharedPair in sharedComHap.keys():
        flag = 'true'
        for lociindex in range(len(lociname2level)):
            lociindex += 1
            hap1Count = 0
            hap2Count = 0

            pattern1 = cleanSampleHFreq[asharedPair.split('__')[0]].codedPattern
            for aHapP1 in pattern1:
                ahap1 = aHapP1.split('_')[0]
                if ahap1 == str(lociindex):
                    hap1Count += 1
            pattern2 = cleanSampleHFreq[asharedPair.split('__')[1]].codedPattern
            for aHapP2 in pattern2:
                ahap2 = aHapP2.split('_')[0]
                if ahap2 == str(lociindex):
                    hap2Count += 1
            if hap1Count > 0 and hap2Count > 0 and hap1Count != hap2Count:
                flag = 'false'
        if flag == 'true':
            realSharedComp[asharedPair] = sharedComHap[asharedPair]
    return realSharedComp

def doubleSingleHap(pattern, lociname2level):
    loci2Ploidy = Tools.findlocil2Ploidy(lociname2level)
    doublePattern = []
    for lociIndex in range(len(lociname2level)):
        lociname = lociname2level[lociIndex]
        hapHolder = []
        for hapIndex in range(len(pattern)):
            if pattern[hapIndex].split('_')[0] == str(lociIndex):
                hapHolder.append(pattern[hapIndex])
        if len(hapHolder) == 1 and lociname not in loci2Ploidy:
            doublePattern = doublePattern + hapHolder + hapHolder
        if len(hapHolder) > 1 and lociname not in loci2Ploidy:
            doublePattern = doublePattern + hapHolder
        if lociname in loci2Ploidy:
            doublePattern = doublePattern + hapHolder
    return doublePattern