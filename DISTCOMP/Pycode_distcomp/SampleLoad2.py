import sys, re, random
import collections
from orderedset import OrderedSet
import pandas as pd
import numpy as np
from SampleContent2 import SampleContent2

__author__ = 'yqb7'
"To do " \
"1. read in original data " \
"2. clean data to trim off the samples with less haps detected."
class SampleLoad2:
    def __init__(self, inputfile):
        # hold the cleandata in a dict type variable
        #self.cleanSamples = self.exportFile(inputfile)
        a = 0
    def importFile(self, inputfile):
        cleanSamples = self.subimport(inputfile)
        return cleanSamples
    def subimport(self, inputfile):
        cleanSamples = {}
        data = pd.read_csv(inputfile, sep='\t')
        lociNameList = list(data.columns)
        locicolumn0 = self.columncount(data, lociNameList)
        # remove the columns without data
        data = self.dropcolumn(data,locicolumn0)
        lociNameList = list(data.columns)
        if "CDC1" in lociNameList[1]:
            lociNameSet1level, lociNameSet2level, lociNameSet3level = self.findLociName(lociNameList)
        else:
            lociNameSet1level, lociNameSet2level, lociNameSet3level = self.findLociNamev2(lociNameList)
        #### remove the loci with only 1 hap = 1column.
        locil2dict1hap = self.countl3hpsl2(lociNameSet2level, lociNameSet3level)
        for adictvalue in locil2dict1hap.values():
            data = self.dropcolumn(data, adictvalue)
        cleanSamples = self.getCleanSample(data)
        # Remove loci with one haplotype in the cleanSamples
        lociNamel2 = list(random.choice(list(cleanSamples.values())).lociname2level)
        lociNamel3 = list(random.choice(list(cleanSamples.values())).lociname3level)
        locisinglehap = self.checkhapdiversity(lociNamel2, lociNamel3, cleanSamples)
        for aloci in locisinglehap.values():
            for onehap in aloci:
                data = data.drop(onehap, axis=1)
        cleanSamples = self.getCleanSample(data)
        return cleanSamples
    def getCleanSample(self, data):
        cleanSamples = {}
        lociNameList = list(data.columns)
        # Refresh the locinames
        if "CDC1" in lociNameList[1]:
            lociNameSet1level, lociNameSet2level, lociNameSet3level = self.findLociName(lociNameList)
        else:
            lociNameSet1level, lociNameSet2level, lociNameSet3level = self.findLociNamev2(lociNameList)
        for eachrow in data.values.tolist():
            asample = {}
            sampleID = eachrow[0].strip('_').strip()
            for rowelementCount in range(len(eachrow)):
                rowelementCount = rowelementCount + 1
                if rowelementCount < len(eachrow):
                    asample[lociNameList[rowelementCount]] = eachrow[rowelementCount]
            aSampleContent = SampleContent2(sampleID, asample, lociNameSet1level, lociNameSet2level,
                                            lociNameSet3level)
            if aSampleContent.totalloci >= 5:
                cleanSamples[sampleID] = aSampleContent
            else:
                l1_3loci_count = 0
                for locinamel1, countl1 in aSampleContent.lociCount_level1.items():
                    if 'CDS' not in locinamel1 and 'CDC' not in locinamel1 and countl1 > 0:
                        l1_3loci_count += 1
                if l1_3loci_count >= 3 and aSampleContent.totalloci >= 4:
                    cleanSamples[sampleID] = aSampleContent
        return cleanSamples
    def checkhapdiversity(self, locinamel2, lociNamel3, cleansamples):
        singleHaploci = dict()
        multihaploci = list()
        for eachloci in locinamel2:
            patternhold = set()
            for eachsample in cleansamples:
                if cleansamples[eachsample].lociPattern[eachloci] != len(cleansamples[eachsample].lociPattern[eachloci]) * "0":
                    patternhold.add(cleansamples[eachsample].lociPattern[eachloci])
                if len(patternhold) > 1:
                    multihaploci.append(eachloci)
                    break
        for eachloci in locinamel2:
            if eachloci not in multihaploci:
                singleHaploci[eachloci] = []
                for aloci in lociNamel3:
                    if eachloci in aloci:
                        singleHaploci[eachloci].append(aloci)
        return singleHaploci

    def dropcolumn(self, data, collist):
        for acol in collist:
            data = data.drop(acol, axis=1)
        return data

    def findLociName(self, locinameList):
        locinameLevel3 = OrderedSet()
        locinameLevel2 = OrderedSet()
        locinameLevel1 = OrderedSet()
        for index in range(len(locinameList)):
            index = index + 1
            if index < len(locinameList):
                # full name CDC1_Hap_1
                locinameLevel3.add(locinameList[index].strip())
                # CDC1_
                locinameLevel2.add(locinameList[index].strip().split('Hap')[0])
                # CDC1
                locinameLevel1.add(locinameList[index].strip().split('_')[0])
        return locinameLevel1, locinameLevel2, locinameLevel3
    #update the method findLociName
    def findLociNamev2(self, locinameList):
        locinameLevel3 = OrderedSet()
        locinameLevel2 = OrderedSet()
        locinameLevel1 = OrderedSet()
        for index in range(len(locinameList)):
            if index < len(locinameList) and locinameList[index] != "Seq_ID":
                locinameLevel3.add(locinameList[index].strip())
                if '.' in locinameList[index]:
                    locinameLevel2.add(re.sub('\d','', locinameList[index].replace('.', '_PART_').split('_PART')[0]))
                    locinameLevel1.add(re.sub('\d','', locinameList[index].replace('.', '_PART_').split('_PART')[0]))
                else:
                    locinameLevel2.add(locinameList[index].split('Hap')[0])
                    locinameLevel1.add(locinameList[index].split('_PART')[0])
        return locinameLevel1, locinameLevel2, locinameLevel3
    def reorderlociname(self, lociname):
        newlociname = ''
        locinamelist = lociname.split('_')
        locinamelist[1:3] = reversed(locinamelist[1:3])
        for index in range(len(locinamelist)):
            newlociname = newlociname + locinamelist[index] + '_'
        newlociname = newlociname.strip('_')
        return newlociname
    #count haps of level3 per loci level2 and return the loci list of level2 with only one hap inside.
    def countl3hpsl2(self, locil2, locil3):
        locil2dict = dict()
        locil2dict1hap = dict()
        for alocil2 in locil2:
            locil2dict[alocil2] = []
            for alocil3 in locil3:
                if alocil2 in alocil3:
                    locil2dict[alocil2].append(alocil3)
        for aloci,value in locil2dict.items():
            if len(value) == 1:
                locil2dict1hap[aloci] = locil2dict[aloci]
        return locil2dict1hap
    # count each column and get the lociname if the sum of column is 0
    def columncount(self, data, locinamelist):
        locinamelistnon0 = dict()
        locicolumn0 = list()
        for key, value in data.iteritems():
            a = 0
        for aloci in locinamelist:
            if aloci != 'Seq_ID':
                #print (aloci)
                if len(list(data[aloci].value_counts())) == 0:
                    locicolumn0.append(aloci)
        return locicolumn0
def main():
    input = sys.argv[1]
    ob = SampleLoad2(input)
#call main
#main()
