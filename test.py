#! /usr/bin/env python3

import BinaryMemorylessDistribution
import QaryMemorylessDistribution
import LinkedListHeap
from math import log2

# q = 3
#
# p = 0.11
# qchannel = QaryMemorylessDistribution.makeQEC(q, p) 
# print( "* original" )
# print(qchannel)
#
# qminus = qchannel.minusTransform()
# degraded = qminus.degrade(9)
# print( "* minus transform" )
# print( qminus )
# print( "* degraded" )
# print( degraded )
#
# qplus = qchannel.plusTransform()
# degraded = qplus.degrade(9)
# print( "* plus transform" )
# print( qplus )
# print( "* degraded" )
# print( degraded )

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
#         channels[m].append(channel.minusTransform().degrade(L))
#         channels[m].append(channel.plusTransform().degrade(L))
#         # channels[m].append(channel.minusTransform().upgrade(L))
#         # channels[m].append(channel.plusTransform().upgrade(L))
#
# entropySum = 0.0
#
# for channel in channels[m]:
#     print( channel.conditionalEntropy() )
#     entropySum += channel.conditionalEntropy()
#
# print( "average capacity = ", 1.0 - entropySum / 2**(n-1) )

q = 3
p = 0.11

qsc = QaryMemorylessDistribution.makeQSC(q, p)

print( "base capacity = ", log2(q) - qsc.conditionalEntropy() )

n = 6
L = 200

channels = []
channels.append([])
channels[0].append(qsc)

for m in range(1,n):
    channels.append([])
    for channel in channels[m-1]:
        channels[m].append(channel.minusTransform().degrade(L))
        channels[m].append(channel.plusTransform().degrade(L))
        # channels[m].append(channel.minusTransform().upgrade(L))
        # channels[m].append(channel.plusTransform().upgrade(L))

entropySum = 0.0

for channel in channels[m]:
    print( channel.conditionalEntropy() )
    entropySum += channel.conditionalEntropy()

print( "average capacity = ", log2(3.0) - entropySum / 2**(n-1) )

