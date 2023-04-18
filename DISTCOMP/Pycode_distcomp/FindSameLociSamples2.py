__author__ = 'Yueli Zheng'

import random
import collections
from orderedset import OrderedSet

import Tools
#Assign values to locipatternWcount
def findSameLociSamples2(cleanSamples):
    lociPatternList = {} # hold all patern of all samples for each loci
    lociPatternSet = {} # hold all pattern set for each loci
    lociPatternHold = {}
    lociname2level = random.choice(list(cleanSamples.values())).lociname2level
    for eachLociName in lociname2level:
        lociPatternList[eachLociName] = []
        lociPatternSet[eachLociName] = OrderedSet()
    for key, avalue in cleanSamples.items():
        for hapname, alociPattern in avalue.lociPattern.items():
            if alociPattern.__contains__('X'):
                lociPatternList[hapname].append(alociPattern)
                lociPatternSet[hapname].add(alociPattern)
    for eachLociName in lociname2level:
        lociPatternHold[eachLociName] = Tools.findSameLociPattern(lociPatternSet[eachLociName], lociPatternList[eachLociName])
    # assign the calculated frequency for the loci of each clean samples
    for key, avalue in cleanSamples.items():
        for hapname in avalue.lociPattern.keys():
            avalue.lociPatternWCount[hapname] = {}
            for ahapPattern in lociPatternHold[hapname].keys():
                if avalue.lociPattern[hapname] == ahapPattern:
                    avalue.lociPatternWCount[hapname][ahapPattern] = lociPatternHold[hapname][ahapPattern]
