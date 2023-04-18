__author__ = 'yqb7'

import math
import os, sys, random
from FrequencyCalBayesian2 import FrequencyCalBayesian2
import Tools, LikelihoodV2Bayesian
from SampleLoad3 import SampleLoad3
# SampleContent,
# SampleLoad,
# FrequencyCalPloid,
# LikelihoodV1Bayesian,
# LikelihoodV1Bayesian
# Bayesian use fixed frequecy values even with different number of shrared alleles,
# Heuristic use different frequency values when
class LikelihoodCalBayesian2:
    def __init__(self, input, expectlociNum, out, expectLociFile, arg):
        if arg == 1:
            ob1 = SampleLoad3(input, expectlociNum, out, expectLociFile)
            cleanSample = ob1.importFile(input, expectlociNum, expectLociFile)
            ob = FrequencyCalBayesian2(cleanSample)
            self.freqDict = ob.freqencyCal(cleanSample)
            # assign the values for sampleFreq of class sampleContent
            cleanSampleWFreq = self.findFreqEachSample(cleanSample, self.freqDict)
            self.bayesianLH = self.pairdistanceCal(cleanSampleWFreq)
            print("finish Bayesian calculation")
        if arg == 2:
            ob1 = SampleLoad3(input, expectlociNum, out, expectLociFile)
            cleanSample = ob1.importFile(input, expectlociNum, expectLociFile)
            ob = FrequencyCalBayesian2(cleanSample)
            self.freqDict = ob.freqencyCal(cleanSample)
    def findFreqEachSample(self, cleanSample, freqDict):
        lociname3level = random.choice(list(cleanSample.values())).lociname3level
        lociname2level = random.choice(list(cleanSample.values())).lociname2level
        for sampleID, value in cleanSample.items():
            #print("clean items value assign", 3)
            for apattern in value.codedPattern:
                hapID = apattern.split('_')[1]
                hapIndex = int(apattern.split('_')[0])
                hapname = Tools.findLevel3LociName(lociname3level, lociname2level[hapIndex-1], hapID)
                value.sampleFreq[hapname] = freqDict[hapname]
        return cleanSample

    def pairdistanceCal(self, cleanSampleWFreq):
        #print('pairdistance baye: ', 1)
        keys = list(cleanSampleWFreq)
        lociname2level = random.choice(list(cleanSampleWFreq.values())).lociname2level
        likeliHood012 = {}
        likeliallloci = {}
        for i in range(len(keys)):
            sample1 = cleanSampleWFreq.get(keys[i])
            for j in range(i, len(keys)):
                #print ("j are ", j)
                sample2 = cleanSampleWFreq.get(keys[j])
                pairFileKey = keys[i] + "__" + keys[j]
                #print("LikelihoodCalBayesian pairFileKey: ", pairFileKey)
                likelihoodinfo = LikelihoodV2Bayesian.likelihoodCal\
                    (keys[i], keys[j], sample1.sampleFreq, sample2.sampleFreq, lociname2level)
                likeliHood012[pairFileKey] = likelihoodinfo[1]
                likeliallloci[pairFileKey] = likelihoodinfo[0]
        return likeliallloci
    '''likeliHood012 holds 3 types of numbers,
    likelihood0: sum of log probabilites of sample1 + sum of log probability of sample2,
    likelihood1: sum of 2 max log probabilites of sample1 + sum of log probability of sample2,
    likelihood2: 1 max log probabilites of sample1 + sum of log probability of sample2.'''

def main():
    input = sys.argv[1]
    ob = LikelihoodCalBayesian2(input)
# call main
#main()

