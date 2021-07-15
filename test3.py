#! /usr/bin/env python3

import BinaryTrellis
import CollectionOfBinaryTrellises
import random
import Guardbands
import PolarEncoderDecoder
import BinaryMemorylessDistribution

# def deletionChannelSimulation(codeword, p, seed, trimmedZerosAtEdges=False):
#     N = len(codeword)
#     receivedWord = []
#
#     random.seed(seed)
#
#
#     for i in range(N):
#         r = random.random()
#
#         if r < p:
#             pass
#         else:
#             receivedWord.append(codeword[i])
#
#
#     if trimmedZerosAtEdges == True:
#         trimmedReceivedWord = Guardbands.trimZerosAtEdges(receivedWord)
#
#         return trimmedReceivedWord
#
#     else:
#         return receivedWord
#
# def temp():
#     bt = BinaryTrellis.BinaryTrellis(2)
#     
#     # print(bt)
#     
#     fromVertex_stateId = 3
#     fromVertex_verticalPosInLayer = 17
#     fromVertex_layer = 0
#     
#     toVertex_stateId = 4
#     toVertex_verticalPosInLayer = 14
#     toVertex_layer = 1
#     
#     edgeLabel = 0
#     edgeProb = 0.5
#     
#     bt.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, edgeProb)
#     
#     print(bt)
#
# codeword = [0,0,1,0]
#
# xi = 0.1
# n = 2
# n0 = 2
#
# withGuardBand = Guardbands.addDeletionGuardBands(codeword, n, n0, xi)
#
# # print(codeword)
# # print(withGuardBand)
# # print(Guardbands.removeDeletionGuardBands(withGuardBand, n, n0))
#
# deletionProb = 0.3
# seed = 0
# trimmedZerosAtEdges=False
#
# receivedWord = deletionChannelSimulation(codeword, deletionProb, seed, trimmedZerosAtEdges)
# trellis = BinaryTrellis.buildTrellis_uniformInput_deletion(receivedWord, len(codeword), deletionProb, trimmedZerosAtEdges)
#
# print(codeword)
# print(receivedWord)
# print(trellis)
# # print(trellis.minusTransform())
# print(trellis.plusTransform([0,0]))
# print(trellis.plusTransform([0,0]).plusTransform([1]))

def make_xVectorDistribuiton_deletion_uniform(length):
    def make_xVectorDistribuiton():
        xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

        xDistribution.probs.append( [0.5,0.5] )
        
        xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)
        return xVectorDistribution
    return make_xVectorDistribuiton

def make_codeword_addDeletionGuardBands(xi, n, n0):
    def make_codeword(encodedVector):
        return Guardbands.addDeletionGuardBands(encodedVector, n, n0, xi)

    return make_codeword

def make_simulateChannel_deletion(p, seed=None):
    def simulateChannel(codeword):
        return BinaryTrellis.deletionChannelSimulation(codeword, p, seed)

    return simulateChannel

def make_xyVectorDistribution_deletion(deletionProbability, xi, n, n0):
    codewordLength = 2 ** n
    def make_xyVectorDistribution(receivedWord):
        return CollectionOfBinaryTrellises.buildCollectionOfBinaryTrellises_uniformInput_deletion(receivedWord, deletionProbability, xi, n, n0)

    return make_xyVectorDistribution

deletionProbability = 0.001
numberOfGenieTrials = 1
numberOfEncodingDecodingTrials = 2000
n = 2
N = 2 ** n

# guardband parameters
xi = 0.1
n0 = 1

upperBoundOnErrorProbability = 0.1

make_xVectorDistribuiton = make_xVectorDistribuiton_deletion_uniform(N)
make_codeword = make_codeword_addDeletionGuardBands(xi, n, n0)

simulateChannel = make_simulateChannel_deletion(deletionProbability)
make_xyVectorDistribution = make_xyVectorDistribution_deletion(deletionProbability, xi, n, n0)

frozenSet = PolarEncoderDecoder.genieEncodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfGenieTrials, upperBoundOnErrorProbability)
# PolarEncoderDecoder.encodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfEncodingDecodingTrials, frozenSet)
