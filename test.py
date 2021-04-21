#! /usr/bin/env python3

import BinaryMemorylessDistribution
import QaryMemorylessDistribution
import LinkedListHeap
from math import log2
import cProfile

def bdegrade():
    p = 0.11
    
    bsc = BinaryMemorylessDistribution.makeBSC(p)
    
    print( "base capacity = ", 1.0 - bsc.conditionalEntropy() )
    
    n = 8
    L = 100
    
    channels = []
    channels.append([])
    channels[0].append(bsc)
    
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
    
    print( "average capacity = ", 1.0 - entropySum / 2**(n-1) )

def qdegrade():
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

qdegrade()
# cProfile.run('qdegrade()')
