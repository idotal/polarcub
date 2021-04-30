#! /usr/bin/env python3

import BinaryMemorylessDistribution
import QaryMemorylessDistribution
import LinkedListHeap
from math import log2
import cProfile

def bupgrade():
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
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", 1.0 - entropySum / 2**(n-1) )

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
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(3.0) - entropySum / 2**(n-1) )

def qupgrade():
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
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(3.0) - entropySum / 2**(n-1) )

def upgradeSimple():

    p = 0.11
    bsc = BinaryMemorylessDistribution.makeBSC(p)
    transformed = bsc

    transformed = transformed.plusTransform()
    transformed = transformed.minusTransform()

    print(transformed)
    upgraded = transformed.upgrade(3)
    print(upgraded)

def qupgradeSimple():

    q = 3
    p = 0.11

    # qsc = QaryMemorylessDistribution.makeQSC(q, p)
    qec = QaryMemorylessDistribution.makeQEC(q, p)
    transformed = qec

    transformed = transformed.plusTransform()
    # transformed = transformed.minusTransform()

    # print(transformed)
    oneHot = transformed.oneHotBinaryMemorylessDistributions()

    upgraded = transformed.upgrade(9)

    # print( oneHot[0] )
    # upgraded = oneHot[0].upgrade(3)

    print( upgraded )
    # upgraded = transformed.upgrade(3)
    # print(upgraded)


# bdegrade()
# bupgrade()
qupgradeSimple()

# cProfile.run('qdegrade()')
