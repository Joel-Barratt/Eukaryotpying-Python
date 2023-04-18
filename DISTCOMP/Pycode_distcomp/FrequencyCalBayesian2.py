__author__ = 'yqb7'
import sys
#from SampleLoad2 import SampleLoad2
import random, collections, Tools

#This class will calculate the 14 loci's frequencies for each sample according to the clean samples
class FrequencyCalBayesian2:
    def __init__(self,cleanSample):
        #ob = SampleLoad2(input)
        #self.cleanSample = ob.exportFile(input)
        #self.freqDict = self.freqencyCal(self.cleanSample)
        #self.freqencyCal(self.cleanSample)
        #print("freqencyCalPloid()", 1)
        cleanSample = cleanSample

    def freqencyCal(self, cleanSample):
        #print("freqencyCal()", "inside")
        lociname3level = random.choice(list(cleanSample.values())).lociname3level
        lociname2level = random.choice(list(cleanSample.values())).lociname2level
        #auto detect the loci names with ploidy = 2
        loci2Ploidy = Tools.findlocil2Ploidy(lociname2level)
        freqDict = {}
        for lociIndex in range(len(lociname2level)):
            singleLociPattern = {}
            for oneSampleID, sampleContent in cleanSample.items():
                singleLociPattern[oneSampleID] = []
                for apattern in sampleContent.codedPattern:
                    if int(apattern.split("_")[0]) == lociIndex + 1:
                        singleLociPattern[oneSampleID].append(apattern)
            #separate the sample patterns according to the ploid conditions
            singlePloid, multiPloid = self.separatePattern(singleLociPattern)
            sum_sp = 0
            sum_mp = 0
            resultsp = collections.Counter(singlePloid)
            resultmp = collections.Counter(multiPloid)
            hap_sm = set() #find out all the haps including both single and multiple ploid haps
            for haps, frequencys in collections.Counter(singlePloid).items():
                hap_sm.add(haps)
            for hapm, frequencym in collections.Counter(multiPloid).items():
                hap_sm.add(hapm)
            #for multiploid sample, frequency does not need *2
            for haps, frequencys in collections.Counter(multiPloid).items():
                sum_mp = sum_mp + frequencys
            lociname = lociname2level[lociIndex]
            #single ploid samples: frequency *2 for non MRT and non Junction loci
            if lociname not in loci2Ploidy:
                for haps, frequencys in collections.Counter(singlePloid).items():
                    sum_sp = sum_sp + frequencys * 2
                for eachHap in hap_sm:
                    hapname = Tools.findLevel3LociName(lociname3level, lociname2level[lociIndex], eachHap)
                    ratio = (resultsp[eachHap] * 2 + resultmp[eachHap])/(sum_sp + sum_mp)
                    freqDict[hapname] = ratio
            # loci 12, 13, 14 are calculated as single ploid, frequency does not need *2
            if lociname in loci2Ploidy:
                for haps, frequencys in collections.Counter(singlePloid).items():
                    sum_sp = sum_sp + frequencys
                for eachHap in hap_sm:
                    hapname = Tools.findLevel3LociName(lociname3level, lociname2level[lociIndex], eachHap)
                    ratio = (resultsp[eachHap] + resultmp[eachHap])/(sum_sp + sum_mp)
                    freqDict[hapname] = ratio
        #print("finish bayesian frequency calculation")
        return freqDict
    def separatePattern(self, singleLociPattern):
        singlePloid = []
        multiplePloid = []
        for value in singleLociPattern.values():
            if len(value) == 1:
                singlePloid.append(int(value[0].split('_')[1]))
            if len(value) > 1:
                for subValue in value:
                    multiplePloid.append(int(subValue.split('_')[1]))
        return singlePloid, multiplePloid

def main():
    input = sys.argv[1]
    ob = FrequencyCalBayesian2(input)
    print("freqencyCalPloid()", 3)
#call main
#main()
