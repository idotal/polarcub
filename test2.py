#! /usr/bin/env python3

from ScalarDistributions import BinaryMemorylessDistribution
from VectorDistributions import BinaryMemorylessVectorDistribution as bmvd
from VectorDistributions import BinaryTrellis
import BinaryPolarEncoderDecoder
import random


# TODO: move these into BinaryMemorylessDistribution, and perhaps make them class methods

def make_xVectorDistribuiton_fromBinaryMemorylessDistribution(xyDistribution, length):
    def make_xVectorDistribuiton():
        xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

        xDistribution.probs.append( [-1.0,-1.0] )
        for x in range(2):
            xDistribution.probs[0][x] = xyDistribution.calcXMarginal(x)

        xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)
        return xVectorDistribution
    return make_xVectorDistribuiton

def make_codeword_noprocessing(encodedVector):
    return encodedVector

def simulateChannel_fromBinaryMemorylessDistribution(xyDistribution):
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

def make_xyVectorDistribution_fromBinaryMemorylessDistribution(xyDistribution):
    def make_xyVectrorDistribution(receivedWord):
        length = len(receivedWord)
        useTrellis = False

        if useTrellis:
            xyVectorDistribution =  xyDistribution.makeBinaryTrellisDistribution(length, receivedWord)
        else:
            xyVectorDistribution = xyDistribution.makeBinaryMemorylessVectorDistribution(length, receivedWord)

        return xyVectorDistribution
    return make_xyVectrorDistribution

p = 0.11
L = 100
n = 7
N = 2 ** n

upperBoundOnErrorProbability = 0.1

xDistribution = None
xyDistribution = BinaryMemorylessDistribution.makeBSC(p)

frozenSet = BinaryMemorylessDistribution.calcFrozenSet_degradingUpgrading(n, L, upperBoundOnErrorProbability, xDistribution, xyDistribution)

# print("Rate = ", N - len(frozenSet), "/", N, " = ", (N - len(frozenSet)) / N)

numberOfTrials = 4000

make_xVectorDistribuiton = make_xVectorDistribuiton_fromBinaryMemorylessDistribution(xyDistribution, N)
make_codeword = make_codeword_noprocessing
simulateChannel = simulateChannel_fromBinaryMemorylessDistribution(xyDistribution)
make_xyVectorDistribution = make_xyVectorDistribution_fromBinaryMemorylessDistribution(xyDistribution)

BinaryPolarEncoderDecoder.encodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfTrials, frozenSet)

# # trustXYProbs = False
# trustXYProbs = True
# PolarEncoderDecoder.genieEncodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfTrials, upperBoundOnErrorProbability, trustXYProbs)
