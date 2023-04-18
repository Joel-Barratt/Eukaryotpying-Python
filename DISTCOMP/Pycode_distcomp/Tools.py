__author__ = 'Yueli Zheng'
import math, pandas as pd
from itertools import combinations
def compareDict(dict1, dict2):
    nSharedKey = 0
    for key1 in dict1.keys():
        for key2 in dict2.keys():
            if key1 == key2:
                nSharedKey = nSharedKey + 1
    return nSharedKey

def compareDictUnsharedLH1V2(dict1, dict2, lociname, loci2Ploidy):
    # shared case1. hap1, hap1;
    # case2. hap1 hap2, hap1 hap2;
    # case3. hap1 hap2 hap3, hap1 hap2 hap3;
    unshared = {}
    shared = {}
    logUnshared = 0
    sumMax = 0
    dict1valuelist = list(dict1.values())
    dict1valuelist.sort()
    for i in range(len(dict1valuelist)):
        sumMax = sumMax + math.log(dict1valuelist[i])
    for key1 in dict1.keys():
        if key1 not in dict2.keys():
            unshared[key1] = dict1[key1]
        if key1 in dict2.keys():
            shared[key1] = dict1[key1]
    if lociname not in loci2Ploidy and len(list(unshared.keys())) == 0 and len(list(dict1.keys())) == 1:
        #print("case1")
        logUnshared = math.log(dict1[list(dict1.keys())[0]])
        return logUnshared
    # if 2 haps in the loci, the max number is needed to calculate LH1
    if len(list(unshared.keys())) >= 0 and len(list(shared.keys())) > 0:
        sharedvaluelist = list(shared.values())
        sharedvaluelist.sort()
        logUnshared = sumMax - math.log(sharedvaluelist[0])
        return logUnshared

def compareDictUnsharedLH2V2(dict1, dict2, lociname, loci2Ploidy):
    # shared case1. hap1, hap1;
    # case2. hap1 hap2, hap1 hap2;
    # case3. hap1 hap2 hap3, hap1 hap2 hap3;
    epsilon = 0.3072
    unshared1 = {}
    shared1 = {}
    unshared2 = {}
    shared2 = {}
    pairflag = 0
    for key1 in dict1.keys():
        if key1 not in dict2.keys():
            unshared1[key1] = dict1[key1]
        if key1 in dict2.keys():
            shared1[key1] = dict1[key1]
    for key2 in dict2.keys():
        if key2 not in dict1.keys():
            unshared2[key2] = dict2[key2]
        if key2 in dict1.keys():
            shared2[key2] = dict2[key2]
    #case1, loci < 12, one hap in the locus, one hap forms one pair
    if len(list(shared1.keys())) == len(list(shared2.keys())) == 1 \
            and len(list(unshared1.keys())) == len(list(unshared2.keys())) == 0:
        pairflag = len(list(shared1.keys()))
    #case2, more than 2 haps in the locus. 2 shared haps form one pair
    if len(list(shared1.keys())) >= 2:
        pairflag = len(list(shared1.keys()))
    # assign value to logUnshared
    if pairflag == 0:
        logUnshared = 1
        return logUnshared
    if lociname not in loci2Ploidy and pairflag == 1:
        logUnshared = math.log(dict1[list(dict1.keys())[0]])
        return logUnshared
    if lociname in loci2Ploidy and pairflag >= 1:
        logUnshared = 0
        return logUnshared
    if lociname not in loci2Ploidy and pairflag >= 2:
        sumMax = 0
        dict1valuelist = list(dict1.values())
        dict1valuelist.sort()
        shared1valuelist = list(shared1.values())
        shared1valuelist.sort()
        for i in range(len(dict1valuelist)):
            sumMax = sumMax + math.log(dict1valuelist[i])
        logUnshared = sumMax - math.log(shared1valuelist[0]) - math.log(shared1valuelist[1])
        return logUnshared

def findSharedMaxLH2(dict1, dict2):
    shared = {}
    sharedMax = 0
    sharedMaxPair = 0
    for key1 in dict1.keys():
        if key1 in dict2.keys():
            shared[key1] = dict1[key1]
    sharedList = list(shared.values())
    sharedList.sort()
    if len(sharedList) >= 1:
        sharedMax = math.log(sharedList[len(sharedList)-1])
        return sharedMax
    else:
        sharedMax = 0
        return sharedMax
def findSharedMaxPairLH1(dict1, dict2):
    #shared case1. hap1, hap1; hap1 hap2, hap1 hap2; hap1 hap2 hap3;
    shared = {}
    sharedMaxPair = 0
    for key1 in dict1.keys():
        if key1 in dict2.keys():
            shared[key1] = dict1[key1]
    sharedList = list(shared.values())
    sharedList.sort()
    if len(sharedList) == len(list(dict1.keys())) == 2:
        sharedMaxPair = 0
        return sharedMaxPair
    if len(sharedList) > 1 and len(list(dict1.keys())) > 2:
        sharedMaxPair = math.log(sharedList[len(sharedList) - 1]) + math.log(sharedList[len(sharedList) - 2])
        return sharedMaxPair
    if len(sharedList) == 1 == len(list(dict1.keys())):
        sharedMaxPair = math.log(sharedList[len(sharedList) - 1])
        return sharedMaxPair
    else:
        sharedMaxPair = 0
        return sharedMaxPair
def findPositionOfPattern(lociPattern):
    position = []
    for index in range(len(lociPattern)):
        if lociPattern[index] == 'X':
            position.append(index)
    return position
#sampattern form different comparation
def pairList(positionList):
    combinationList = []
    for index1 in range(len(positionList)):
        length = index1 + 1
        combi = combinations(positionList,length)
        combinationList.append(list(combi))
    return combinationList
#find the frequency according to the shared pattern.
def findSameLociPattern(valueSetLoci, listValueLoci):
    # listValueLoci holds all loci patterns of all samples in the same loci
    # {pattern:{position tuple : count}}
    lociPatternHold = {}
    # value1 is a pattern in the pattern set
    for value1 in valueSetLoci:
        if str(value1) == str('0000XX0000'):
            a = 0
        value1Pos = findPositionOfPattern(value1)
        keyTupleHolder = []
        combinationList = pairList(value1Pos) # To do combination for all possible patterns
        for atuple in combinationList: # combinationList 1st level: key
            for tupleIndex in range(len(atuple)):
                keyTupleHolder.append(atuple[tupleIndex])
        tupleCountDict = {}
        for keyList in keyTupleHolder:
            count_tuple = 0
            value1List_key = {}
            value1List_key[keyList] = []
            for keyIndex in range(len(keyList)):
                value1List_key[keyList].append(value1[keyList[keyIndex]])
            for value2 in listValueLoci:
                value2List_key = {}
                value2List_key[keyList] = []
                for keyIndex in range(len(keyList)):
                    value2List_key[keyList].append(value2[keyList[keyIndex]])
                if value1List_key == value2List_key:
                    count_tuple += 1
            tupleCountDict[keyList] = count_tuple
        lociPatternHold[value1] = tupleCountDict
    return lociPatternHold

def findHapsValueDictHeuris(sampleFreqHeuristic, lociIndex):
    hapsValueDict = {}
    for lociname, pfreq in sampleFreqHeuristic.items():
        if lociname.split('CDC')[1] == str(lociIndex):
            for pattern, freqs in pfreq.items():
                for position, freq in freqs.items():
                    hapsValueDict[position] = freq
    return hapsValueDict

def findHapsValueDictHeuris2(sampleFreqHeuristic, lociName):
    hapsValueDict = {}
    for lociname, pfreq in sampleFreqHeuristic.items():
        if lociname == str(lociName):
            for pattern, freqs in pfreq.items():
                for position, freq in freqs.items():
                    hapsValueDict[position] = freq
    return hapsValueDict

def findNsharedallelesHeuris(hapsValueDict1, hapsValueDict2):
    shared = 0
    sharedAllele = []
    for key1 in hapsValueDict1.keys():
        for key2 in hapsValueDict2.keys():
            if key1 == key2:
                shared = shared + 1
                sharedAllele.append(key1)
    return sharedAllele

def findHapCount(hapsValueDict):
    hapList = set()
    for key in hapsValueDict.keys():
        for each in list(key):
            hapList.add(each)
    return hapList

def findMinimumValue(hapsCount1, hapsCount2):
    tmpCount = 0
    if hapsCount1 >= hapsCount2:
        tmpCount = hapsCount2
    else: tmpCount = hapsCount1
    return tmpCount
def findHeurH_nu(bayesianFreq, lociindex):
    bayesianFreqList = []
    H_nu = 0
    for key, value in bayesianFreq.items():
        keyIndex = int(key.split("_")[1].split('CDC')[1])
        if keyIndex == lociindex and value > 0:
            bayesianFreqList.append(value)
    for avalue in bayesianFreqList:
        H_nu = H_nu + avalue * math.log(avalue,2)
    return -H_nu, bayesianFreqList
def findHeurH_nu2(bayesianFreq, lociName):
    bayesianFreqList = []
    H_nu = 0
    for key, value in bayesianFreq.items():
        if key.__contains__(lociName) and value > 0:
            bayesianFreqList.append(value)
    for avalue in bayesianFreqList:
        H_nu = H_nu + avalue * math.log(avalue,2)
    return -H_nu, bayesianFreqList
def findLevel3LociName(lociname3level, level2Loci, hapID):
    for alevel3Lociname in lociname3level:
        hapMarker = alevel3Lociname.split('_')[len(alevel3Lociname.split('_')) - 1]
        if alevel3Lociname.__contains__(level2Loci) and int(hapMarker) == int(hapID):
            return alevel3Lociname

def findlocil2Ploidy(lociname2level):
    loci1Ploidy = list()
    for eachname in lociname2level:
        if eachname.__contains__('MSR') or eachname.__contains__('Junction') or 'Mt' in eachname:
            loci1Ploidy.append(eachname)
    return loci1Ploidy

def output(sampleID, matrix, outputname):
    sampleNum = sampleID.__len__()
    df = pd.DataFrame(matrix.reshape(sampleNum, sampleNum), columns = sampleID, index = sampleID)
    df.to_csv(outputname)


