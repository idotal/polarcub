#! /usr/bin/env python3

import BinaryMemorylessDistribution
import LinkedListHeap

# llh = LinkedListHeap.LinkedListHeap()
# llh.insertAtTail(8, ["a"])
# llh.insertAtTail(5, ["b"])
# llh.insertAtTail(13, ["c"])
# llh.insertAtTail(9, ["d"])
# llh.insertAtTail(4, ["e"])
# llh.insertAtTail(100, ["f"])
# print( llh )
#
# llh.extractHeapMin()
# print( llh )
#
# newMin = llh.getHeapMin()
#
# llh.updateKey(newMin.leftElementInList, 103)
# print( llh )
#
# print( llh.returnData() ) 
#
# data = [[0.1,0.2], [0.4,0.3], [0.1,0.7], [0.9,0.8]]
# key  = [   3     ,     9    ,     2    ,    4     ]
#
# llh = LinkedListHeap.LinkedListHeap(key, data)
# print( llh )
#
# dataLeft = [0.01,0.02]
# dataCenter = [0.11,0.25]
# dataRight = None #[0.3,0.65]

p = 0.1

bsc = BinaryMemorylessDistribution.makeBSC(p)
minus = bsc.minusTransform()
plus = bsc.plusTransform()
plusplus = plus.plusTransform()
plusplusplus = plusplus.plusTransform()

# print( minus.probs )
# print( plus.probs )
# print( bsc.conditionalEntropy() )
# print( minus.conditionalEntropy() )
# print( plus.conditionalEntropy() )

d = plusplusplus.degrade(4)
print( plusplusplus.conditionalEntropy() )
print( d.conditionalEntropy() )




