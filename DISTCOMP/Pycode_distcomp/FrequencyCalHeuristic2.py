__author__ = 'yqb7'
import sys, random
#from SampleLoad2 import SampleLoad2
import FindSameLociSamples2
from FrequencyCalBayesian2 import FrequencyCalBayesian2
#This class will calculate the 14 loci's frequencies
class FrequencyCalHeuristic2:
    def __init__(self, cleanSample):
        #ob = SampleLoad2(input)
        #FindSameLociSamples2.findSameLociSamples2(ob.cleanSamples)
        #self.findCDCTotalHapHeuris(ob.cleanSamples)
        #self.freqencyCalHeuristic(ob.cleanSamples)
        #calculate the total hap for each loci and prepare the number to calculate
        #the frequency according to different pattern.
        #print("FrequencyCalHeuristic", 1)
        a = 0
    # To do: find out the total samples with confirmed sequence for each loci with heuristic rule
    def findCDCTotalHapHeuris(self, cleanSamples):
        lociname2level = random.choice(list(cleanSamples.values())).lociname2level
        CDCTotalHapHeuris = {}
        for eachLoci in lociname2level:
            CDCTotalHapHeuris[eachLoci] = 0
        for onesample in cleanSamples.values():
            for lociName, lociValue in onesample.lociCount.items():
                if int(lociValue) > 0:
                    CDCTotalHapHeuris[lociName] += 1
        return CDCTotalHapHeuris

    def freqencyCalHeuristic(self,cleanSamples):
        #print("FrequencyCalHeuristic", 2)
        CDCTotalHapHeuris = self.findCDCTotalHapHeuris(cleanSamples)
        # Calculate the frequency according to the combination list
        for onesample in cleanSamples.values():
            lociDictFreqHeuriPattern = {}
            for lociName, patternWCount in onesample.lociPatternWCount.items():
                for pattern, posWcount in patternWCount.items():
                    if pattern.__contains__('X'):
                        value2SameLevel = {}
                        value3SameLevel = {}
                        for lociName1, lociCount in CDCTotalHapHeuris.items():
                            for pos, count in posWcount.items():
                                if lociName1 == lociName:
                                    tmpFreq = count/lociCount
                                    value3SameLevel[pos] = tmpFreq
                            value2SameLevel[pattern] = value3SameLevel
                        lociDictFreqHeuriPattern[lociName] = value2SameLevel
            onesample.sampleFreqHeurisPattern = lociDictFreqHeuriPattern
        #print("heuristic frequency calculation finish")
        return cleanSamples
def main():
    input = sys.argv[1]
    ob = FrequencyCalHeuristic2(input)
#call main
#main()
