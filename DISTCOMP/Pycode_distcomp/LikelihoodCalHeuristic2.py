__author__ = 'yqb7'

import copy
import math, Tools, time
import os, sys, random, collections
from orderedset import OrderedSet
import string, numpy as np
from FrequencyCalBayesian2 import FrequencyCalBayesian2
from FrequencyCalHeuristic2 import FrequencyCalHeuristic2
import Tools, LikelihoodV2Heuristic
import LikelihoodCalBayesian2, re
import FindSameLociSamples2, FindClosestSamples
from LikelihoodCalBayesian2 import LikelihoodCalBayesian2
from SampleLoad3 import SampleLoad3

class LikelihoodCalHeuristic2:
    def __init__(self, input, expectlociNum, out, expectLociFile):
        # assign the values for sampleFreq of class sampleContent
        # self.cleanSampleWFreq = LikelihoodCalBayesian(input).cleanSampleWFreq
        ob = SampleLoad3(input,expectlociNum, out, expectLociFile)
        cleanSample = ob.importFile(input, expectlociNum, expectLociFile)
        # print out the correct clean samples 1 time here
        ob.out(cleanSample, out)
        cleanSampleName = list(cleanSample.keys())
        FindSameLociSamples2.findSameLociSamples2(cleanSample)
        obfch = FrequencyCalHeuristic2(cleanSample)
        self.CDCTotalHapHeuris = obfch.findCDCTotalHapHeuris(cleanSample)
        self.cleanSampleHFreq = obfch.freqencyCalHeuristic(cleanSample)
        sampleIDlist = list(self.cleanSampleHFreq.keys()) #complete cleansample IDs
        obBay = LikelihoodCalBayesian2(input, expectlociNum, out, expectLociFile, 2)
        #obBay.findFreqEachSample(self.cleanSampleHFreq, obBay.freqDict)
        self.bayesianFreq = obBay.freqDict
        # keep the original likelihood in pairdistanceHeuristic
        pairdistanceHeuristic, lociname2level = self.pairdistanceCalHeuristic(self.cleanSampleHFreq)
        likelihoodheuristic = dict()
        counter = 0
        likelihoodheuristic[str(counter)] = copy.deepcopy(pairdistanceHeuristic)
        samplelh0dict = dict()
        #check0 = self.check0valuenp(likelihoodheuristic.get(str(counter)))
        print("start empty value fill")
        while self.check0valuenp(likelihoodheuristic[str(counter)], counter) == 'true':
            #print("cycle-fill  counter ", counter)
            if counter == 0:
                samplelh0 = self.check0valueSamp(cleanSample)
                samplelh0dict[str(counter)] = samplelh0
            if counter > 0 :
                samplelh0 = self.check0value(likelihoodheuristic[str(counter)],samplelh0dict['0'], counter)
                samplelh0dict[str(counter)] = samplelh0
            counter += 1
            matchingSampleList, samplesMissLoci, samplesCompLociList = \
                FindClosestSamples.findClosestSample(self.cleanSampleHFreq)
            #print("finish find matching samples")
            likelihoodheuristic[str(counter)] = \
                self.fillMissedValue(likelihoodheuristic[str(counter-1)], matchingSampleList,
                samplelh0dict[str(counter-1)], lociname2level, counter - 1, cleanSampleName)
            aaa = 0
        print('heauristc matrix finish updating ')
        self.heurisLH = self.normheurislh(likelihoodheuristic[str(counter)])
    def pairdistanceCalHeuristic(self, cleanSampleHFreq):
        start = time.time()
        likeliHood = {}
        keys = list(cleanSampleHFreq)
        lociname2level = random.choice(list(cleanSampleHFreq.values())).lociname2level
        for i in range(len(keys)):
            #print("pairdistanheuristic", i)
            sample1 = self.cleanSampleHFreq.get(keys[i])
            for j in range(i, len(keys)):
                sample2 = self.cleanSampleHFreq.get(keys[j])
                pairFileKey = keys[i] + "__" + keys[j]
                #print("heuristic pairFileKey: ", pairFileKey)
                likelihoodinfo = LikelihoodV2Heuristic.\
                    heuristicLikelihoodCal(keys[i], keys[j], sample1.sampleFreqHeurisPattern,
                                           sample2.sampleFreqHeurisPattern, self.bayesianFreq,lociname2level)
                likeliHood[pairFileKey] = likelihoodinfo
        end = time.time()
        print("heuristic pairdistance calculation complete, time used: ", end - start)
        return likeliHood, lociname2level

    # check the samples having missing loci with the likelihood
    def check0valuenp(self, likelihood, counter):
        #print('84 check0valuenp starts. counter ', counter)
        check0 = 'false'
        for pairkey, locivaluelist in likelihood.items():
            cleanvalues = self.convertDictV2ListV(pairkey, locivaluelist.values())
            #print('pairkey 88 ', pairkey, ' counter ', counter)
            if 0 in cleanvalues:
                check0 = 'true'
                return check0
        return check0
    # obtain the orignal samples with missing locis from the clean samples
    def check0valueSamp(self, cleanSample):
        start = time.time()
        mssamploci = dict()
        for sampleID, content in cleanSample.items():
            msloci = list()
            for loci, count in content.lociCount.items():
                if count == 0:
                    msloci.append(loci)
            if len(msloci) > 0:
                mssamploci[sampleID] = msloci
        end = time.time()
        #print("check0valueSamp ends: ", end - start)
        return mssamploci

    #Find out the samples having missing loci based on searching the likelihood
    def check0value(self, likelihood, sampleidlh0list, counter):
        #print("check0value starts")
        start = time.time()
        samplelist = list()
        for pairkey, value in likelihood.items():
            cleanValues = self.convertDictV2ListV(pairkey, value.values())
            if 0 in cleanValues:
                samplelist.append(pairkey)
                a = 0
        newmissloci = set()
        newmissdict = dict()
        for apair in samplelist:
            pair = apair.split('__')
            newmissloci.add(pair[0])
            newmissloci.add(pair[1])
        for asam in newmissloci:
            newmissdict[asam] = sampleidlh0list[asam]
        #to solve the problem that matrix outputs at different runs will have slightly value difference
        #fix the order of the miss loci samples
        newmissdict = dict(sorted(newmissdict.items(), key = lambda x: list(sampleidlh0list.keys()).index(x[0])))
        end = time.time()
        print("check0value ends: ", end - start)
        return newmissdict

    def findmisloci(self, mslocisamlist, mslocisamples):
        start = time.time()
        mssamploci = dict()
        for sample in mslocisamlist:
            mslocilist = []
            for loci, samplelist in mslocisamples.items():
                if sample in samplelist:
                    mslocilist.append(loci)
            mssamploci[sample] = mslocilist
        end = time.time()
        print("findmisloci ends: ", end - start)
        return mssamploci
    # to reduce the size of the likelihoodHeuris and save time in filling values with matching samples.
    def shrinkMatchingSample(self, matchingSample, likelihoodheur):
        machingset = set(a.split('__')[1] for a in matchingSample)
        sublhheur = dict()
        for a in machingset:
            sublhheur[a] = list()
        sublhheur_diag = dict()
        for a in machingset:
            tmp = dict()
            for key, value in likelihoodheur.items():
                if key.__contains__(a) and key.split('__')[0] != key.split('__')[1]:
                    tmp[key] = value
                if key.__contains__(a) and key.split('__')[0] == key.split('__')[1]:
                    sublhheur_diag[key] = value
            sublhheur[a].append(tmp.copy())
        return sublhheur, sublhheur_diag
    # updating on 01/27/2021 with updating matching samples search method
    def fillMissedValue2(self,likelihood, matchingSample, samplelh0, lociname2level, counter, sublhheur, sublhheur_diag):
        #print("fillMissedValue, l83")
        start = time.time()
        # case1, find match sample and assign the value
        # case2, can not find the match sample, average the corresponding loci and assign the value
        meanN0Loci = self.meanNon0Loci(likelihood, lociname2level)
        meanN0Loci_diag = self.meanNon0_diag_loci(likelihood, lociname2level)
        filledlikelihood = copy.deepcopy(likelihood)
        #self.checkmatchedsample(samplelh0, matchingSample)
        if counter == 1:
            a = 0
        for misslociS, missloci in samplelh0.items():
            msample = self.checkmatchedsample(misslociS, matchingSample)
            if len(msample) >= 1:
                filledlikelihood = self.fillWMatch2(filledlikelihood, likelihood, msample, misslociS, missloci, sublhheur, sublhheur_diag )
            else:
                filledlikelihood = self.fillmeanN0loci(filledlikelihood, meanN0Loci, misslociS, missloci)
                filledlikelihood = self.fillmeanN0loci_diag(filledlikelihood, meanN0Loci_diag, misslociS, missloci)
        end = time.time()
        print("fillMissedValue ends: ", end - start)
        return filledlikelihood

    def fillWMatch2(self,filledlikelihood, likelihood, smatchsample, misslociS, missloci, sublhheur, sublhheur_diag):
        # to update the likelihoods(diagonal and non-diagonol )of the sample with missing loci with the matching sample's likelihood value.
        #To do: average the matchind sample loci values
        #case1. 1 match sample
        if len(smatchsample) == 1:
            for pairkey, value in filledlikelihood.items():
                if pairkey.__contains__(misslociS):
                    matchpairkey1 = pairkey.replace(misslociS, smatchsample[0])
                    anotherPK = self.new2posibPK(matchpairkey1)
                    matchpairkey = self.findMkey(matchpairkey1, anotherPK, likelihood.keys())
                    for amsloci in missloci:
                        value[amsloci] = likelihood[matchpairkey].get(amsloci)
        #case2. multiple match samples
        if len(smatchsample) > 1:
            msampleHold = dict()
            locivalue = dict()
            for asamp in smatchsample:
                msampleHold[asamp] = sublhheur[asamp]
            for amsloci in missloci:
                locivalue[amsloci] = list()
        return filledlikelihood

    #one sample with miss loci could have one or multiple matching samples.
    def checkmatchedsample(self, misslociS, matchingSample):
        start = time.time()
        msample = list()
        for smatchsample in matchingSample.keys():
            if misslociS == smatchsample.split('__')[0]:
                msample.append(smatchsample.split('__')[1])
        end = time.time()
        #print("checkmatchedsample ends: ", end - start)
        return msample

    # fill in the positions with likelihood = 0
    # samplelh0 holds the samples with missing loci.
    def fillMissedValue(self,likelihood, matchingSample, samplelh0, lociname2level, counter, cleanSampleName):
        #print("fillMissedValue, l220")
        start = time.time()
        # case1, find match sample and assign the value
        # case2, can not find the match sample, average the corresponding loci and assign the value
        meanN0Loci = self.meanNon0Loci(likelihood, lociname2level, counter)
        meanN0Loci_diag = self.meanNon0_diag_loci(likelihood, lociname2level)
        filledlikelihood = copy.deepcopy(likelihood)
        #collect the all filledvalues.
        filledvaluehold = dict()
        #self.checkmatchedsample(samplelh0, matchingSample)
        if counter == 0:
            a = 0
        if counter == 1:
            a = 0
        for misslociS, missloci in samplelh0.items():
            msample = self.checkmatchedsample(misslociS, matchingSample)
            if len(msample) >= 1:
                filledlikelihood = \
                    self.fillWMatch(filledlikelihood, likelihood, msample, misslociS, missloci, cleanSampleName, filledvaluehold)
                a = 0
                #value0hold.update(value0)
            else:
                filledlikelihood = self.fillmeanN0loci(filledlikelihood, meanN0Loci, misslociS, missloci, cleanSampleName, filledvaluehold)
                filledlikelihood = self.fillmeanN0loci_diag(filledlikelihood, meanN0Loci_diag, misslociS, missloci, filledvaluehold)
        #newmissloci = self.find0values(filledvaluehold, samplelh0, filledlikelihood)
        end = time.time()
        #print("fillMissedValue ends: ", end - start)
        return filledlikelihood
    #01/27/2021 update multiple matching samples
    def fillWMatch(self,filledlikelihood, likelihood, smatchsample, misslociS, missloci, cleanSampleName, filledvaluehold):
        # to update the likelihoods(diagonal and non-diagonol )of the sample with missing loci with the matching sample's likelihood value.
        #To do: average the matchind sample loci values
        #case1. 1 match sample
        # value0 hold all pairkey values with likelihood = 0
        value0 = dict()
        if len(smatchsample) == 1:
            pairkeyList = list()
            for asname in cleanSampleName:
                # build pairkeys single matching samples
                # matchpairkey1 is the pairkey with missing loci sample
                # matchpairkey2 is the pairkey with maching sample
                tmppair1 = asname + '__' + misslociS
                anotherPK1 = self.new2posibPK(tmppair1)
                matchpairkey1 = self.findMkey(tmppair1, anotherPK1, likelihood.keys())
                #if asname != smatchsample[0]:
                tmppair2 = asname + '__' + smatchsample[0]
                anotherPK2 = self.new2posibPK(tmppair2)
                matchpairkey2 = self.findMkey(tmppair2, anotherPK2, likelihood.keys())
                pairkeyList.append(matchpairkey1 + '|||' + matchpairkey2)
            for eachpairkey in pairkeyList:
                twopairkey = eachpairkey.split('|||')
                tmpfilledvalue = dict()
                tmp = dict()
                for amsloci in missloci:
                    filledlikelihood[twopairkey[0]][amsloci] = likelihood[twopairkey[1]][amsloci]
                    tmpfilledvalue[amsloci] = likelihood[twopairkey[1]][amsloci]
                    '''if likelihood[twopairkey[1]][amsloci] == 0:
                        tmp[amsloci] = 0
                        value0[twopairkey[0]] = tmp'''
                filledvaluehold[twopairkey[0]] = tmpfilledvalue
                #filledvaluehold[twopairkey[0]] = []
                #filledvaluehold[twopairkey[0]].append(tmpfilledvalue)
                a = 0
            #updating diagnal value
            tmpfilledvalue = dict()
            tmp = dict()
            for amsloci in missloci:
                filledlikelihood[misslociS + '__' + misslociS][amsloci] = \
                    likelihood[smatchsample[0] + '__' + smatchsample[0]][amsloci]
                tmpfilledvalue[amsloci] = likelihood[smatchsample[0] + '__' + smatchsample[0]][amsloci]
            filledvaluehold[misslociS + '__' + misslociS] = tmpfilledvalue
            #filledvaluehold[misslociS + '__' + misslociS] = []
            #filledvaluehold[misslociS + '__' + misslociS].append(tmpfilledvalue)
            a = 0
        #case2. multiple match samples
        if len(smatchsample) > 1:
            #build the pairkeys for the multiple matching samples.
            valuehold = dict()
            for asampleName in cleanSampleName:
                tmppair1 = asampleName + '__' + misslociS
                anotherPK1 = self.new2posibPK(tmppair1)
                matchpairkey1 = self.findMkey(tmppair1, anotherPK1, likelihood.keys())
                valuehold[matchpairkey1] = list()
                # find out the matching samples for the misslociS pairkeys
                for amatchsam in smatchsample:
                    tmppair2 = asampleName + '__' + amatchsam
                    anotherPK2 = self.new2posibPK(tmppair2)
                    matchpairkey2 = self.findMkey(tmppair2, anotherPK2, likelihood.keys())
                    valuehold[matchpairkey1].append(matchpairkey2)
                # find out the corresponding values to fill the missing values
                tmp = dict()
                tmpfilledvalue = dict()
                for amisloci in missloci:
                    filledlikelihood[matchpairkey1][amisloci] = \
                        self.meanLoci(amisloci, valuehold[matchpairkey1], likelihood)
                    tmpfilledvalue[amisloci] = self.meanLoci(amisloci, valuehold[matchpairkey1], likelihood)
                    '''if self.meanLoci(amisloci, valuehold[matchpairkey1], likelihood) == 0:
                        tmp[amisloci] = 0
                        value0[matchpairkey1] = tmp'''
                filledvaluehold[matchpairkey1] = tmpfilledvalue
                #filledvaluehold[matchpairkey1] = []
                #filledvaluehold[matchpairkey1].append(tmpfilledvalue)
            #build pairkeys of diagonal values for matching samples
            matchSpairkeys = list()
            for amatchsam in smatchsample:
                matchSpairkeys.append(amatchsam + '__' + amatchsam)
            tmpfilledvalue = dict()
            tmp = dict()
            for amisloci in missloci:
                filledlikelihood[misslociS + '__' + misslociS][amisloci] = \
                    self.meanloci_diag(amisloci, matchSpairkeys, likelihood)
                tmpfilledvalue[amisloci] = self.meanloci_diag(amisloci, matchSpairkeys, likelihood)
            filledvaluehold[misslociS + '__' + misslociS] = tmpfilledvalue
            #filledvaluehold[misslociS + '__' + misslociS] = []
            #filledvaluehold[misslociS + '__' + misslociS].append(tmpfilledvalue)
            a = 0
        return filledlikelihood
    # this function is for missinglociS with multiple matching samples
    def meanLoci(self, loci, valuehold, likelihood):
        locisum = 0
        for apairkey in valuehold:
            locisum = locisum + likelihood[apairkey][loci]
        meanv = locisum/len(valuehold)
        return meanv

    # this function is for missinglociS with multiple matching samples
    def meanloci_diag(self, loci, pairkeylist, likelihood):
        locisum_diag = 0
        for apairkey in pairkeylist:
            locisum_diag = locisum_diag + likelihood[apairkey][loci]
        mean_diag = locisum_diag/len(pairkeylist)
        return mean_diag
    def fillmeanN0loci(self, filledlikelihood, meanN0Loci, misslociS, missloci, cleanSampleName, filledvaluehold):
        #build pairkey for the misslociS
        pairkeyList = list()
        for asname in cleanSampleName:
            tmppair = asname + '__' + misslociS
            anotherPK = self.new2posibPK(tmppair)
            matchpairkey = self.findMkey(tmppair, anotherPK, filledlikelihood.keys())
            pairkeyList.append(matchpairkey)
        for pairkey in pairkeyList:
            #print("pairkey fillmeanN0loci ", pairkey)
            tmpfilledValue = dict()
            for amsloci in missloci:
                filledlikelihood[pairkey][amsloci] = meanN0Loci[amsloci]
                tmpfilledValue[amsloci] = meanN0Loci[amsloci]
            filledvaluehold[pairkey] = tmpfilledValue
        return filledlikelihood
    # updating in the 0 values in diagonal lines of the matrix
    def fillmeanN0loci_diag(self, filledlikelihood, meanN0Loci_diag, misslociS, missloci, filledvaluehold):
        # build pairkey for the misslociS
        pairkey = misslociS + '__' + misslociS
        tmpfilledValue = dict()
        for amsloci in missloci:
            filledlikelihood[pairkey][amsloci] = meanN0Loci_diag[amsloci]
            tmpfilledValue[amsloci] = meanN0Loci_diag[amsloci]
        filledvaluehold[pairkey] = tmpfilledValue
        return filledlikelihood
    #This function is used to find the missloci samples according to the filled in values.
    def find0values(self, filledvaluehold, originalmisslociS, lhfilled):
        samplelist = list()
        for pairkey, value in filledvaluehold.items():
            #print('pairkey ', pairkey)
            cleanValues = self.convertDictV2ListV(pairkey, value.values())
            if 0 in cleanValues:
                samplelist.append(pairkey)
                a = 0
        newmissloci = set()
        for apair in samplelist:
            pair = apair.split('__')
            newmissloci.add(pair[0])
            newmissloci.add(pair[1])
        # obtain the different samples by comparing the missloci samples got from filledvalue with the original missloci samples
        difsamp = [ele for ele in list(originalmisslociS.keys()) if ele not in newmissloci]
        difsampdict = dict()
        for asam in difsamp:
            difsampdict[asam] = originalmisslociS[asam]
        #Check if the difsamp has complete non 0 values in the filled likelihood if it still has 0 then combine it
        #with newmissloci. To do:
        #self.checkdiffsam0(difsampdict)

        return newmissloci
    def convertDictV2ListV(self, pairkey, dictvalues):
        cleanvalues = list()
        for avalue in list(dictvalues):
            #print('type ', type(avalue))
            if type(avalue) is list:
                cleanvalues.append(avalue[0])
            else:
                cleanvalues.append(avalue)
        return cleanvalues

    def meanNon0Loci(self, likelihood, lociname2level, counter):
        meanNon0loci = dict()
        for lociname in lociname2level:
            #print("l508 counter + loci ", counter, ' | ', lociname)
            likelihood_loci = []
            likelihood_diag_loci = []
            for akeypair, avalue in likelihood.items():
                #after filling value types are mixed
                avalue[lociname] = self.cleanlocivalue(avalue[lociname])
                if (avalue[lociname] > 0) and akeypair.split('__')[0] != akeypair.split('__')[1]:
                    likelihood_loci.append(avalue[lociname])
                if (avalue[lociname] > 0) and akeypair.split('__')[0] == akeypair.split('__')[1]:
                    likelihood_diag_loci.append(avalue[lociname])
            mean_loci = (sum(np.array(likelihood_loci)) * 2 + sum(np.array(likelihood_diag_loci)))/\
                        (len(likelihood_loci) * 2 + len(likelihood_diag_loci))
            meanNon0loci[lociname] = mean_loci
        return meanNon0loci
    def cleanlocivalue(self, avalue):
        if type(avalue) is list:
            cleanvalue = avalue[0]
            return cleanvalue
        else: cleanvalue = avalue
        return cleanvalue
    def meanNon0_diag_loci(self, likelihood, lociname2level):
        meanNon0_diag_loci = dict()
        for lociname in lociname2level:
            likelihood_diag_loci = []
            for akeypair, avalue in likelihood.items():
                avalue[lociname] = self.cleanlocivalue(avalue[lociname])
                if (avalue[lociname] > 0) and akeypair.split('__')[0] == akeypair.split('__')[1]:
                    likelihood_diag_loci.append(avalue[lociname])
            mean_diag_loci = sum(np.array(likelihood_diag_loci))/len(likelihood_diag_loci)
            meanNon0_diag_loci[lociname] = mean_diag_loci
        return meanNon0_diag_loci

    def normheurislh(self, likelihood):
        normheurlh = dict()
        for pairkey, value in likelihood.items():
            sumloci = sum(np.array(list(value.values())))
            normheurlh[pairkey] = sumloci
        return normheurlh

    def newPKreorder(self, newpairkey):
        newpairkey1 = ''
        key1 = newpairkey.split('__')[0].split('_')[1]
        key1num = re.sub('[A-Za-z]+', '', key1).lstrip('0')
        key2 = newpairkey.split('__')[1].split('_')[1]
        key2num = re.sub('[A-Za-z]+', '', key2).lstrip('0')
        if int(key1num) < int(key2num):
            newpairkey1 = newpairkey
        if int(key1num) >= int(key2num):
            newpairkey1 = newpairkey.split('__')[1] + '__' + newpairkey.split('__')[0]
        a = 0
        return newpairkey1
    #This function is to find another possible pairkey
    def new2posibPK(self, newpairkey):
        twosampID = newpairkey.split('__')
        twosampID = twosampID[::-1]
        anotherPPK = twosampID[0] + '__' + twosampID[1]
        return anotherPPK

    def findMkey(self, key1, key2, likelihoodkeys):
        if key1 in likelihoodkeys:
            return key1
        if key2 in likelihoodkeys:
            return key2

def main():
    input = sys.argv[1]
    ob = LikelihoodCalHeuristic2(input)
# call main
# main()