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
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", 1.0 - entropySum / 2**n )

def bdegrade():
    p = 0.11
    
    bsc = BinaryMemorylessDistribution.makeBSC(p)
    
    print( "base capacity = ", 1.0 - bsc.conditionalEntropy() )
    
    n = 8
    L = 100
    
    channels = []
    channels.append([])
    channels[0].append(bsc)
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().degrade(L))
            channels[m].append(channel.plusTransform().degrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", 1.0 - entropySum / 2**n )

def qdegrade_static():
    q = 3
    p = 0.11
    
    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    
    print( "base capacity = ", log2(q) - qsc.conditionalEntropy() )
    
    n = 5
    # n = 2
    L = 400
    
    channels = []
    channels.append([])
    channels[0].append(qsc)
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().degrade_static(L))
            channels[m].append(channel.plusTransform().degrade_static(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( log2(q) - channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(q) - entropySum / 2**n )

def qdegrade():
    q = 3
    p = 0.11
    
    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    
    print( "base capacity = ", log2(q) - qsc.conditionalEntropy() )
    
    n = 5
    L = 400
    
    channels = []
    channels.append([])
    channels[0].append(qsc)
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().degrade(L))
            channels[m].append(channel.plusTransform().degrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( log2(q) - channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(q) - entropySum / 2**n )

def qupgrade():
    q = 3
    p = 0.11
    
    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    
    print( "base capacity = ", log2(q) - qsc.conditionalEntropy() )
    
    n = 5
    L = 400
    
    channels = []
    channels.append([])
    channels[0].append(qsc)
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( log2(q) - channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(q) - entropySum / 2**n )

def qupgrade_static():
    q = 3
    p = 0.11
    
    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    
    print( "base capacity = ", log2(q) - qsc.conditionalEntropy() )
    
    # n = 5
    n = 2
    L = 400
    
    channels = []
    channels.append([])
    channels[0].append(qsc)
    
    for m in range(1,n+1):
        channels.append([])
        for channel in channels[m-1]:
            channels[m].append(channel.minusTransform().upgrade_static(L))
            channels[m].append(channel.plusTransform().upgrade_static(L))
    
    entropySum = 0.0
    
    for channel in channels[m]:
        print( log2(q) - channel.conditionalEntropy() )
        entropySum += channel.conditionalEntropy()
    
    print( "average capacity = ", log2(q) - entropySum / 2**n )

def upgradeSimple():

    p = 0.11
    bsc = BinaryMemorylessDistribution.makeBSC(p)
    transformed = bsc

    transformed = transformed.plusTransform()
    transformed = transformed.minusTransform()

    print(transformed)
    upgraded = transformed.upgrade(3)
    print(upgraded)

def qdegradeSimple():

    q = 3
    p = 0.11
    L = 16

    # myset = {1, 5, 7}
    # print(myset)
    # for i in myset:
    #     print(i)

    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    transformed = qsc

    # qec = QaryMemorylessDistribution.makeQEC(q, p)
    # transformed = qec

    transformed = transformed.plusTransform()
    # transformed = transformed.plusTransform()
    print(transformed)
    transformed = transformed.degrade_static(100)
    print(transformed)
    # transformed = transformed.minusTransform()
    #
    # print("original")
    # print(transformed)
    # oneHot = transformed.oneHotBinaryMemorylessDistributions()
    #
    # upgraded = transformed.upgrade(L)
    #
    # # print( oneHot[0] )
    # # upgraded = oneHot[0].upgrade(3)
    #
    # print("upgraded")
    # print( upgraded )
    # # upgraded = transformed.upgrade(3)
    # # print(upgraded)
    #
    # oneHotUpgraded = upgraded.oneHotBinaryMemorylessDistributions()
    # print("transformed one-hot")
    # print(oneHotUpgraded[0])
    # print(oneHotUpgraded[1])

def qupgradeSimple():

    q = 3
    p = 0.11
    L = 16

    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    transformed = qsc

    # qec = QaryMemorylessDistribution.makeQEC(q, p)
    # transformed = qec

    transformed = transformed.plusTransform()
    transformed = transformed.degrade(100)
    # transformed = transformed.minusTransform()

    print("original")
    print(transformed)
    oneHot = transformed.oneHotBinaryMemorylessDistributions()

    upgraded = transformed.upgrade(L)

    # print( oneHot[0] )
    # upgraded = oneHot[0].upgrade(3)

    print("upgraded")
    print( upgraded )
    # upgraded = transformed.upgrade(3)
    # print(upgraded)

    oneHotUpgraded = upgraded.oneHotBinaryMemorylessDistributions()
    print("transformed one-hot")
    print(oneHotUpgraded[0])
    print(oneHotUpgraded[1])

# bdegrade()
# bupgrade()
# qupgradeSimple()
# qdegrade()
# qdegrade_static()
# qupgrade()
qupgrade_static()
# qdegradeSimple()

# cProfile.run('qupgrade()')
