#! /usr/bin/env python3

import BinaryMemorylessDistribution
import QaryMemorylessDistribution
import LinkedListHeap

q = 3

p = 0.11
qchannel = QaryMemorylessDistribution.makeQEC(q, p) 
print(qchannel)

# p = 0.11
#
# bsc = BinaryMemorylessDistribution.makeBSC(p)
#
# print( "base capacity = ", 1.0 - bsc.conditionalEntropy() )
#
# n = 8
# L = 100
#
# channels = []
# channels.append([])
# channels[0].append(bsc)
#
# for m in range(1,n):
#     channels.append([])
#     for channel in channels[m-1]:
#         # channels[m].append(channel.minusTransform().degrade(L))
#         # channels[m].append(channel.plusTransform().degrade(L))
#         channels[m].append(channel.minusTransform().upgrade(L))
#         channels[m].append(channel.plusTransform().upgrade(L))
#
# entropySum = 0.0
#
# for channel in channels[m]:
#     print( channel.conditionalEntropy() )
#     entropySum += channel.conditionalEntropy()
#
# print( "average capacity = ", 1.0 - entropySum / 2**(n-1) )


