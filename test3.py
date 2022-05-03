#! /usr/bin/env python3

import random

import QaryPolarEncoderDecoder
from ScalarDistributions import QaryMemorylessDistribution


def make_xVectorDistribution_fromQaryMemorylessDistribution(q, xyDistribution, length):
    def make_xVectorDistribution():
        xDistribution = QaryMemorylessDistribution.QaryMemorylessDistribution(q)
        xDistribution.probs = [xyDistribution.calcXMarginals()]
        xVectorDistribution = xDistribution.makeQaryMemorylessVectorDistribution(length, None)
        return xVectorDistribution
    return make_xVectorDistribution


def make_codeword_noprocessing(encodedVector):
    return encodedVector

def simulateChannel_fromQaryMemorylessDistribution(xyDistribution):
    def simulateChannel(codeword):
        receivedWord = []
        length = len(codeword)

        for j in range(length):
            x = codeword[j]

            rand = random.random()
            probSum = 0.0

            for y in range(len(xyDistribution.probs)):
                if probSum + xyDistribution.probXGivenY(x,y) >= rand:
                    receivedWord.append(y)
                    # print("x = ", x, ", y = ", y, " probXGivenY(x,y) = ", xyDistribution.probXGivenY(x,y), ", rand = ", rand)
                    break
                else:
                    probSum += xyDistribution.probXGivenY(x,y)

        return receivedWord

    return simulateChannel

def make_xyVectorDistribution_fromQaryMemorylessDistribution(xyDistribution):
    def make_xyVectrorDistribution(receivedWord):
        length = len(receivedWord)
        useTrellis = False

        if useTrellis:
            xyVectorDistribution = xyDistribution.makeQaryTrellisDistribution(length, receivedWord)
        else:
            xyVectorDistribution = xyDistribution.makeQaryMemorylessVectorDistribution(length, receivedWord)

        return xyVectorDistribution
    return make_xyVectrorDistribution

q = 2
p = 0.11
L = 100
n = 7
N = 2 ** n

upperBoundOnErrorProbability = 0.1

xDistribution = None
xyDistribution = QaryMemorylessDistribution.makeQSC(q, p)

frozenSet = QaryMemorylessDistribution.calcFrozenSet_degradingUpgrading(n, L, upperBoundOnErrorProbability, xDistribution, xyDistribution)

# print("Rate = ", N - len(frozenSet), "/", N, " = ", (N - len(frozenSet)) / N)

numberOfTrials = 4000

make_xVectorDistribution = make_xVectorDistribution_fromQaryMemorylessDistribution(q, xyDistribution, N)
make_codeword = make_codeword_noprocessing
simulateChannel = simulateChannel_fromQaryMemorylessDistribution(xyDistribution)
make_xyVectorDistribution = make_xyVectorDistribution_fromQaryMemorylessDistribution(xyDistribution)

QaryPolarEncoderDecoder.encodeDecodeSimulation(q, N, make_xVectorDistribution, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfTrials, frozenSet, verbosity=0)

# # trustXYProbs = False
# trustXYProbs = True
# PolarEncoderDecoder.genieEncodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfTrials, upperBoundOnErrorProbability, trustXYProbs)