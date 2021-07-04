#! /usr/bin/env python3

import BinaryTrellis

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
