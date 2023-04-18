__author__ = 'Yueli Zheng'
import os, sys, numpy as np, collections, pandas as pd, time, math, argparse
from argparse import ArgumentParser
from LikelihoodCalBayesian2 import LikelihoodCalBayesian2
from LikelihoodCalHeuristic2 import LikelihoodCalHeuristic2
class BayesianHeuristic:
    def __init__(self, input, outputHH, outputHB, outputB, outputH, expectlociNum,
                 out, expectLociFile):
        start = time.time()
        expectlociNum = expectlociNum
        # assign the values for sampleFreq of class sampleContent
        bayesianLH = LikelihoodCalBayesian2(input, expectlociNum, out, expectLociFile, 1).bayesianLH
        bayesianLHnorm = dict()
        for pairkey in bayesianLH.keys():
            bayesianLHnorm[pairkey] = 2 - bayesianLH[pairkey]
        minBtmp = min(np.array(list(bayesianLHnorm.values())))
        for pairkey in bayesianLHnorm.keys():
            bayesianLHnorm[pairkey] = bayesianLHnorm[pairkey] - minBtmp
        maxBtmp = max(np.array(list(bayesianLHnorm.values())))
        for pairkey in bayesianLHnorm.keys():
            bayesianLHnorm[pairkey] = bayesianLHnorm[pairkey]/maxBtmp
        np.array(list(bayesianLHnorm.values()))
        end1 = time.time()
        print("bayesian Run finish. Time used to complete: ", end1 - start)
        heurisLHnorm = dict()
        heurisLH = LikelihoodCalHeuristic2(input, expectlociNum,
                 out, expectLociFile).heurisLH
        minHeuris = min(np.array(list(heurisLH.values())))
        for pairkey in heurisLH.keys():
            heurisLHnorm[pairkey] = heurisLH[pairkey] - minHeuris
        maxHeuris = max(np.array(list(heurisLHnorm.values())))
        for pairkey in heurisLHnorm.keys():
            heurisLHnorm[pairkey] = heurisLHnorm[pairkey]/maxHeuris
        np.array(list(heurisLHnorm.values()))
        # hbmatrix1 is the average with H-norm + sorted H with rank order of B-norm
        # hbmatrix2 is the average with H-norm + sorted B with rank order of H-norm
        sampleID, hbmatrix1, hbmatrix2 = self.averagebayheu(bayesianLHnorm, heurisLHnorm, outputB, outputH)
        #hbmatrix2 = hbmatrix2.reshape(-1,1)
        self.output(sampleID, hbmatrix1, outputHH)
        self.output(sampleID, hbmatrix2, outputHB)
        end = time.time()
        print("heuristic Run finish. Time used to complete: ", end - end1)
        print("Run finish. Time used for the whole process to complete: ", end - start)

    def completeMatrix(self, halfMatrix):
        compmatrix = dict()
        for samplePair, value in halfMatrix.items():
            if samplePair.split('__')[0] != samplePair.split('__')[1]:
                compmatrix[samplePair] = value
                compmatrix[samplePair.split('__')[1] + '__' + samplePair.split('__')[0]] = value
            if samplePair.split('__')[0] == samplePair.split('__')[1]:
                compmatrix[samplePair] = value
        # sort the commatrix to match the original sample lists
        commatrix = collections.OrderedDict(sorted(compmatrix.items()))
        return commatrix

    # sort bayesianLHnorm and average it with the heurisLHnorm. The sample pairs of bayesianHnorm is not considered.
    def averagebayheu(self, bayesianLHnorm, heurisLHnorm, output1, output2):
        bayesianLHnormcomp = self.completeMatrix(bayesianLHnorm)
        heurisLHnormcomp = self.completeMatrix(heurisLHnorm)
        sampleID, average1 = self.averageMethod1(bayesianLHnormcomp, heurisLHnormcomp)
        average2 = self.averageMethod2(bayesianLHnormcomp, heurisLHnormcomp)
        bayesianLHnormcompnp = np.array(list(bayesianLHnormcomp.values()))
        heurisLHnormcompnp = np.array(list(heurisLHnormcomp.values()))
        self.output(sampleID, bayesianLHnormcompnp, output1)
        self.output(sampleID, heurisLHnormcompnp, output2)
        return sampleID, average1, average2

    # average1: heuris norm values and heuris sort values, do not use bayesian values. HH
    def averageMethod1(self, bayesianLHnorm, heurisLHnorm):
        # rank the bayesian matrix
        bayeRanks = pd.Series(list(bayesianLHnorm.values()))
        bayeRank = bayeRanks.sample(frac=1).rank(method = 'first').reindex_like(bayeRanks)

        heuriSort = sorted(heurisLHnorm.items(), key=lambda kv: kv[1])
        bayRank = list(bayeRank.array.astype(int))
        heurisSort_D = list()
        # sort the heuristic according to the bayRank
        for element in bayRank:
            heurisSort_D.append(heuriSort[element -1])
        heurisSort_D = dict(heurisSort_D)
        # list holds two matrixs ,heurisLHnormcomp and heurisSort_D
        listD_heuris = list()
        listD_heuris.append(list(heurisLHnorm.values()))
        listD_heuris.append(list(heurisSort_D.values()))
        tmp1 = np.array(list(heurisLHnorm.values()))
        tmp2 = np.array(list(heurisSort_D.values()))
        # average heuris norm and heuris sort, do not use bayesian values
        averageHB = np.add(list(heurisLHnorm.values()), list(heurisSort_D.values()))/2
        maxlh = max(averageHB)
        averageHBratio = averageHB / maxlh
        sampleID = list(bayesianLHnorm.keys())[0:int(math.sqrt(len(bayesianLHnorm.keys())))]
        sampleIDnew = [i.split('__', 1)[1] for i in sampleID]
        return sampleIDnew, averageHBratio

    # average2 the bayesian sort value and heurisLHnorm results, HB.
    def averageMethod2(self, bayesianLHnorm, heurisLHnorm):
        # rank the heuristic matrix
        heurisRanks = pd.Series(list(heurisLHnorm.values()))
        heurisRank = heurisRanks.sample(frac=1).rank(method='first').reindex_like(heurisRanks)
        # Rank heuristic matrix, sort b-matrix according to ranked h-matrix
        bayeSort = sorted(bayesianLHnorm.items(), key=lambda kv: kv[1])
        heurisRank = list(heurisRank.array.astype(int))
        bayeSort_D = list()
        # sort the bayesian according to the bayRank
        for element in heurisRank:
            bayeSort_D.append(bayeSort[element - 1])
        bayeSort_D = dict(bayeSort_D)
        # list holds two matrixs ,heurisLHnormcomp and bayeSort_D
        listD_ensemble = list()
        listD_ensemble.append(list(heurisLHnorm.values()))
        listD_ensemble.append(list(bayeSort_D.values()))
        # average the bayesian sort value and heurisLHnorm results.
        averageHB = np.add(list(heurisLHnorm.values()), list(bayeSort_D.values()))/2
        maxlh = max(averageHB)
        averageHBratio = averageHB / maxlh
        return averageHBratio

    def output(self, sampleID, matrix, outputname):
        sampleNum = sampleID.__len__()
        df = pd.DataFrame(matrix.reshape(sampleNum, sampleNum), columns = sampleID, index = sampleID)
        df.to_csv(outputname)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-infile', type=argparse.FileType('r', encoding='UTF-8'), required=True)
    parser.add_argument('-outhh', type=argparse.FileType('w', encoding='UTF-8'), required=True)
    parser.add_argument('-outbh', type=argparse.FileType('w', encoding='UTF-8'), required=True)
    parser.add_argument('-outB', type=argparse.FileType('w', encoding='UTF-8'), required=True)
    parser.add_argument('-outH', type=argparse.FileType('w', encoding='UTF-8'), required=True)
    parser.add_argument('-expectlocinumber', type=int, required=True)
    parser.add_argument('-sampmeetlocirequire', type=argparse.FileType('w', encoding='UTF-8'), required=True)
    parser.add_argument('-expectlocifile', type=argparse.FileType('r', encoding='UTF-8'), required=True)
    args = parser.parse_args()
    ob = BayesianHeuristic(args.infile.name, args.outhh.name, args.outbh.name, args.outB.name, args.outH.name,
                           args.expectlocinumber, args.sampmeetlocirequire.name, args.expectlocifile.name)

# call main
main()