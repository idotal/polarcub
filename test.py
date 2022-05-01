#! /usr/bin/env python3

from math import log2

from ScalarDistributions import BinaryMemorylessDistribution
from ScalarDistributions import QaryMemorylessDistribution


def bupgrade():
    p = 0.11

    bsc = BinaryMemorylessDistribution.makeBSC(p)

    print("base capacity = ", 1.0 - bsc.conditionalEntropy())

    n = 8
    L = 100

    channels = []
    channels.append([])
    channels[0].append(bsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))

    entropySum = 0.0

    for channel in channels[m]:
        print(channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", 1.0 - entropySum / 2 ** n)


def bdegrade():
    p = 0.11

    bsc = BinaryMemorylessDistribution.makeBSC(p)

    print("base capacity = ", 1.0 - bsc.conditionalEntropy())

    n = 8
    L = 100

    channels = []
    channels.append([])
    channels[0].append(bsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().degrade(L))
            channels[m].append(channel.plusTransform().degrade(L))

    entropySum = 0.0

    for channel in channels[m]:
        print(channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", 1.0 - entropySum / 2 ** n)


def qdegrade_static():
    q = 3
    p = 0.11

    qsc = QaryMemorylessDistribution.makeQSC(q, p)

    print("base capacity = ", log2(q) - qsc.conditionalEntropy())

    n = 5
    L = 400

    binningToUse = QaryMemorylessDistribution.Binning.TalSharovVardy  # standard for degrading
    # binningToUse = QaryMemorylessDistribution.Binning.PeregTal # non-standard

    channels = []
    channels.append([])
    channels[0].append(qsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().degrade_static(L, binningToUse))
            channels[m].append(channel.plusTransform().degrade_static(L, binningToUse))

    entropySum = 0.0

    for channel in channels[m]:
        print(log2(q) - channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", log2(q) - entropySum / 2 ** n)


def qdegrade():
    q = 3
    p = 0.11

    qsc = QaryMemorylessDistribution.makeQSC(q, p)

    print("base capacity = ", log2(q) - qsc.conditionalEntropy())

    n = 5
    L = 400

    channels = []
    channels.append([])
    channels[0].append(qsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().degrade(L))
            channels[m].append(channel.plusTransform().degrade(L))

    entropySum = 0.0

    for channel in channels[m]:
        print(log2(q) - channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", log2(q) - entropySum / 2 ** n)


def qupgrade():
    q = 3
    p = 0.11

    qsc = QaryMemorylessDistribution.makeQSC(q, p)

    print("base capacity = ", log2(q) - qsc.conditionalEntropy())

    n = 5
    L = 400

    channels = []
    channels.append([])
    channels[0].append(qsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))

    entropySum = 0.0

    for channel in channels[m]:
        print(log2(q) - channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", log2(q) - entropySum / 2 ** n)


def qupgrade_static():
    q = 3
    p = 0.11

    qsc = QaryMemorylessDistribution.makeQSC(q, p)

    print("base capacity = ", log2(q) - qsc.conditionalEntropy())

    n = 5
    L = 400

    binningToUse = QaryMemorylessDistribution.Binning.PeregTal  # standard for upgrading
    # binningToUse = QaryMemorylessDistribution.Binning.TalSharovVardy # non-standard

    channels = []
    channels.append([])
    channels[0].append(qsc)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().upgrade_static(L, binningToUse))
            channels[m].append(channel.plusTransform().upgrade_static(L, binningToUse))

    entropySum = 0.0

    for channel in channels[m]:
        print(log2(q) - channel.conditionalEntropy())
        entropySum += channel.conditionalEntropy()

    print("average capacity = ", log2(q) - entropySum / 2 ** n)


def upgradeSimple():
    p = 0.11
    bsc = BinaryMemorylessDistribution.makeBSC(p)
    transformed = bsc

    transformed = transformed.plusTransform()
    transformed = transformed.minusTransform()

    print(transformed)
    upgraded = transformed.upgrade_static(400)
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
    L = 400
    # binningToUse = QaryMemorylessDistribution.Binning.TalSharovVardy 
    binningToUse = QaryMemorylessDistribution.Binning.PeregTal

    # probs = [0.8,0.1,0.1]
    # qsc = QaryMemorylessDistribution.makeInputDistribution(probs)
    qsc = QaryMemorylessDistribution.makeQSC(q, p)
    transformed = qsc

    # qec = QaryMemorylessDistribution.makeQEC(q, p)
    # transformed = qec

    transformed = transformed.plusTransform()
    # transformed = transformed.degrade(L)
    transformed = transformed.upgrade_static(L)
    transformed = transformed.minusTransform()

    print("original")
    print(transformed)

    M = transformed.calcMFromL(L)
    mu = transformed.calcMuForPeregTal(M)

    print("mu =", mu)

    upgraded = transformed.upgrade_static(L)

    print("upgraded")
    print(upgraded)

    # degraded = transformed.degrade_static(L, binningToUse)
    #
    # print("degraded")
    # print( degraded )


def qupgradeInputDistribution():
    probs = [0.34, 0.33, 0.33]
    inputDist = QaryMemorylessDistribution.makeInputDistribution(probs)

    print("base entropy = ", inputDist.conditionalEntropy())
    print("probs = ", probs)

    n = 8
    L = 100
    tvLimit = 0.001

    channels = []
    channels.append([])
    channels[0].append(inputDist)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().upgrade(L))
            channels[m].append(channel.plusTransform().upgrade(L))

    entropySum = 0.0
    goodCount = 0

    for channel in channels[m]:
        # print(channel.probs)
        print(channel.conditionalEntropy(), channel.totalVariation())
        entropySum += channel.conditionalEntropy()
        if channel.totalVariation() < tvLimit:
            goodCount += 1

    print("average entropy = ", entropySum / 2 ** n)
    print("pass tvLimit of ", tvLimit, " = ", goodCount)


def qupgradeInputDistribution_static():
    probs = [0.34, 0.33, 0.33]
    inputDist = QaryMemorylessDistribution.makeInputDistribution(probs)

    print("base entropy = ", inputDist.conditionalEntropy())
    print("probs = ", probs)

    n = 8
    L = 100
    tvLimit = 0.001

    channels = []
    channels.append([])
    channels[0].append(inputDist)

    for m in range(1, n + 1):
        channels.append([])
        for channel in channels[m - 1]:
            channels[m].append(channel.minusTransform().upgrade_static(L))
            channels[m].append(channel.plusTransform().upgrade_static(L))

    entropySum = 0.0
    goodCount = 0

    for channel in channels[m]:
        # print(channel.probs)
        print(channel.conditionalEntropy(), channel.totalVariation())
        entropySum += channel.conditionalEntropy()
        if channel.totalVariation() < tvLimit:
            goodCount += 1

    print("average entropy = ", entropySum / 2 ** n)
    print("pass tvLimit of ", tvLimit, " = ", goodCount)


def qupgradeMultivariateNormal():
    sigmaSquare = 0.1
    gridMax = 2.0
    gridMin = -2.0
    gridPoints = 100

    q = 3
    MM = 20
    original = QaryMemorylessDistribution.makeMultivariateNormal(sigmaSquare, gridMax, gridMin, gridPoints)

    originalEntropy = original.conditionalEntropy()
    print("*", originalEntropy)

    # print("* upper bound on difference, equation (13) in Ordentlich and Tal")
    # for L in range(10,730,20):
    #     print(L, QaryMemorylessDistribution.upgrade_dynamic_upper_bound(q, L))
    #
    # print("* lower bound on upgrading cost, equation (43) in Kartowsky and Tal")
    # for L in range(10,730,20):
    #     print(L, QaryMemorylessDistribution.upgrade_cost_lower_bound(q, L))

    print("* static upgrade")
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.upgrade_static(L)
        print(transformed.calcOutputAlphabetSize(), originalEntropy - transformed.conditionalEntropy())

    print("* dynamic upgrade")

    MM = 27
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.upgrade(L)
        print(transformed.calcOutputAlphabetSize(), originalEntropy - transformed.conditionalEntropy())


# def qdegradeMultivariateNormal():
#     sigmaSquare = 0.05
#     gridMax = 1.2
#     gridMin = -1.2
#     gridPoints = 400
#
#     q = 3
#     MM = 20
#     original = QaryMemorylessDistribution.makeMultivariateNormal(sigmaSquare, gridMax, gridMin, gridPoints)
#
#     originalEntropy = original.conditionalEntropy()
#     print("*",  originalEntropy)
#
#     # print("* upper bound on difference, equation (13) in Ordentlich and Tal")
#     # for L in range(10,730,20):
#     #     print(L, QaryMemorylessDistribution.upgrade_dynamic_upper_bound(q, L))
#     #
#     # print("* lower bound on upgrading cost, equation (43) in Kartowsky and Tal")
#     # for L in range(10,730,20):
#     #     print(L, QaryMemorylessDistribution.upgrade_cost_lower_bound(q, L))
#
#     print("* static degrade")
#     for M in range(2,MM+1):
#         L = M ** (q-1)
#         transformed = original.degrade_static(L)
#         print( transformed.calcOutputAlphabetSize(), transformed.conditionalEntropy() - originalEntropy )
#
#     print("* dynamic degrade")
#
#     MM = 27
#     for M in range(2,MM+1):
#         L = M ** (q-1)
#         transformed = original.degrade(L)
#         print( transformed.calcOutputAlphabetSize(), transformed.conditionalEntropy() - originalEntropy )
#

def qupgradeUniform():
    q = 3
    T = 400
    MM = 20
    original = QaryMemorylessDistribution.makeQuantizedUniform(q, T)

    originalEntropy = original.conditionalEntropy()
    print("*", originalEntropy)

    print("* upper bound on difference, equation (13) in Ordentlich and Tal")
    for L in range(10, 730, 20):
        print(L, QaryMemorylessDistribution.upgrade_dynamic_upper_bound(q, L))

    print("* lower bound on upgrading cost, equation (43) in Kartowsky and Tal")
    for L in range(10, 730, 20):
        print(L, QaryMemorylessDistribution.upgrade_cost_lower_bound(q, L))

    print("* static upgrade")
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.upgrade_static(L)
        print(transformed.calcOutputAlphabetSize(), originalEntropy - transformed.conditionalEntropy())

    print("* dynamic upgrade")

    MM = 27
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.upgrade(L)
        print(transformed.calcOutputAlphabetSize(), originalEntropy - transformed.conditionalEntropy())


def qdegradeUniform():
    q = 3
    T = 400
    MM = 20
    original = QaryMemorylessDistribution.makeQuantizedUniform(q, T)

    originalEntropy = original.conditionalEntropy()
    print("*", originalEntropy)

    print("* upper bound on difference, equation (27) in Ordentlich and Tal")
    for L in range(10, 730, 20):
        print(L, QaryMemorylessDistribution.degrade_dynamic_upper_bound(q, L))

    print(
        "* lower bound on degrading cost, equation (3) in Tal, On the Construction of Polar Codes for Channelswith Moderate Input Alphabet Sizes")
    for L in range(10, 730, 20):
        print(L, QaryMemorylessDistribution.degrade_cost_lower_bound(q, L))

    # Strange! There is not a big difference!
    print("* static degrade")
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.degrade_static(L)
        print(transformed.calcOutputAlphabetSize(), transformed.conditionalEntropy() - originalEntropy)

    print("* dynamic degrade")

    MM = 27
    for M in range(2, MM + 1):
        L = M ** (q - 1)
        transformed = original.degrade(L)
        print(transformed.calcOutputAlphabetSize(), transformed.conditionalEntropy() - originalEntropy)


# bdegrade()
# bupgrade()
# qdegrade()
# qdegrade_static()
# qupgrade_static()
# qupgrade()
# qdegradeSimple()
# qupgradeSimple()
# qupgradeUniform()
# qdegradeUniform()
# qupgradeInputDistribution()
# qupgradeInputDistribution_static()
# cProfile.run('qupgrade()')
# qupgradeMultivariateNormal()
# qdegradeMultivariateNormal()

p = 0.11
channel = BinaryMemorylessDistribution.makeBEC(p)
print(channel.mmse())
transformed = channel.plusTransform()
print(transformed.mmse())
transformed = channel.minusTransform()
print(transformed.mmse())
