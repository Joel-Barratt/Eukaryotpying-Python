__author__ = 'Yueli Zheng'

import math, pandas as pd
class SampleContent2:
    def __init__(self,sampleID, aSample, lociNameSet1level,lociNameSet2level, lociNameSet3level):
        self.sampleID = sampleID
        self.aSample = aSample
        self.lociname1level = lociNameSet1level
        self.lociname2level = lociNameSet2level
        self.lociname3level = lociNameSet3level
        self.lociCount_level1 = self.calLociV2()
        self.totalloci = self.calTotalLoci()
        self.lociCount, self.codedPattern, self.lociPattern = self.caleachLociTotal()
        self.lociPatternWCount = {}
        self.sampleFreq = {}
        self.sampleFreqHeuristic = {}
        self.sampleFreqHeurisPattern = {}
    # calculate the overall loci info
    def calTotalLoci(self):
        totallociTmp = 0
        for aname in self.lociname1level:
            loci = 'false'
            # print(self.aSample.items())
            for key, value in self.aSample.items():
                if key.__contains__(aname) and value == "X":
                    loci = 'true'
                    # print(key,value)
                    break
            if loci == 'true':
                totallociTmp = totallociTmp + 1
        return totallociTmp
    # calculate the loci info
    def calLociV2(self):
        # initialize lociCount_l1
        lociCount_l1 = dict()
        for aname in self.lociname1level:
            lociCount_l1[aname] = 0
        for key, value in self.aSample.items():
            for lociname in self.lociname1level:
                if key.__contains__(lociname) and value == "X":
                    lociCount_l1[lociname] = lociCount_l1[lociname] + 1
        return lociCount_l1
    def caleachLociTotal(self):
        lociCount = dict()
        codedPattern = list()
        lociPattern = dict()
        #initialize lociCount
        for aname in self.lociname2level:
            lociCount[aname] = 0
            lociPattern[aname] = ''
        anameIndex = 0
        for aname in self.lociname2level:
            anameIndex = anameIndex + 1
            for lociname, lociValue in self.aSample.items():
                if lociname.__contains__(aname) and lociValue == "X":
                    lociCount[aname] = lociCount[aname] + 1
                    codedPattern.append(str(anameIndex) + '_' + lociname.split('Hap_')[1])
                    lociPattern[aname]= lociPattern[aname] + lociValue
                if lociname.__contains__(aname) and pd.isna(lociValue):
                    lociPattern[aname] = lociPattern[aname] + "0"
        return lociCount, codedPattern, lociPattern

def main():
    ob = SampleContent2()
#main
#main()