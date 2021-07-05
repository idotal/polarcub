#! /usr/bin/env python3

import BinaryTrellis
import random

def deletionChannelSimulation(codeword, p, seed):
    N = len(codeword)
    receivedWord = []

    random.seed(seed)


    for i in range(N):
        r = random.random()

        if r < p:
            pass
        else:
            receivedWord.append(codeword[i])

    return receivedWord

def temp():
    bt = BinaryTrellis.BinaryTrellis(2)
    
    # print(bt)
    
    fromVertex_stateId = 3
    fromVertex_verticalPosInLayer = 17
    fromVertex_layer = 0
    
    toVertex_stateId = 4
    toVertex_verticalPosInLayer = 14
    toVertex_layer = 1
    
    edgeLabel = 0
    edgeProb = 0.5
    
    bt.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, edgeProb)
    
    print(bt)

codeword = [1,0,1,1,1,0]
p = 0.9
seed = 0

receivedWord = deletionChannelSimulation(codeword, p, seed)
print(codeword)
print(receivedWord)
