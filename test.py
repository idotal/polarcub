#! /usr/bin/env python3

import BinaryMemorylessDistribution
import LinkedListHeap

llh = LinkedListHeap.LinkedListHeap()
llh.insertAtTail(8, ["a"])
llh.insertAtTail(5, ["b"])
llh.insertAtTail(13, ["c"])
llh.insertAtTail(9, ["d"])
llh.insertAtTail(4, ["e"])
llh.insertAtTail(100, ["f"])
print ( llh )

llh.extractHeapMin()
print ( llh )


# bmd = BinaryMemorylessDistribution.BinaryMemorylessDistribution()
#
# p = 0.1
#
# bsc = BinaryMemorylessDistribution.makeBSC(p)
# minus = bsc.minusTransform()
# plus = bsc.plusTransform()
#
# # print( minus.probs )
# # print( plus.probs )
# print( bsc.conditionalEntropy() )
# print( minus.conditionalEntropy() )
# print( plus.conditionalEntropy() )
#



