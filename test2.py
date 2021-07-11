#! /usr/bin/env python3

import BinaryMemorylessDistribution
import BinaryMemorylessVectorDistribution as bmvd
import BinaryTrellis
import PolarEncoderDecoder
import random

def calcFrozenSet_degradingUpgrading(n, L, upperBoundOnErrorProbability, xDistribution, xyDistribution):
    """Calculate the frozen set by degrading and upgrading the a-priori and joint distributions, respectively

    Args:
        n (int): number of polar transforms

        L (int): number of quantization levels, for both encoding and decoding

        upperBoundOnErrorProbability (float): select an index i to be unfrozen if the total-variation K(U_i|U^{i-1}) \leq epsilon/(2N) and the probatility of error Pe(U_i | U^{i-1}, Y^{N-1}) \leq epsilon/(2N), where epsilon is short for upperBoundOnErrorProbability 

        xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1). If this is set to None, then we assume a uniform input distribution, in which case the above criterion for i being unfrozen is simplified to Pe(U_i | U^{i-1}, Y^{N-1}) \leq epsilon/N.

        xyVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-posteriori entries for P(X=0) and P(X=1). That is, entry i contains P(X=0,Y=y_i) and P(X=1,Y=y_i).
    Returns:
        frozenSet (set): the set of frozen indices
    """

    assert( n >= 0 )
    assert( L > 0 )
    assert( upperBoundOnErrorProbability > 0 )
    assert( xyDistribution is not None )

    if xDistribution is not None:
        xDists = []
        xDists.append([])
        xDists[0].append(xDistribution)

        for m in range(1,n+1):
            xDists.append([])
            for dist in xDists[m-1]:
                xDists[m].append(dist.minusTransform().upgrade(L))
                xDists[m].append(dist.plusTransform().upgrade(L))

    xyDists = []
    xyDists.append([])
    xyDists[0].append(xyDistribution)

    for m in range(1,n+1):
        xyDists.append([])
        for dist in xyDists[m-1]:
            xyDists[m].append(dist.minusTransform().degrade(L))
            xyDists[m].append(dist.plusTransform().degrade(L))

    frozenSet = set()

    N = 1 << n

    if xDistribution is not None:
        delta =  upperBoundOnErrorProbability / (2 * N)
        for i in 2 ** n:
            if xyDists[m][i].errorProb() > delta or xDists[m][i].totalVariation() > delta:
                frozenSet.add(i)

    else:
        delta =  upperBoundOnErrorProbability / N
        for i in range(2 ** n):
            if xyDists[m][i].errorProb() > delta:
                frozenSet.add(i)

    return frozenSet

def genieEncode():
    uLen = 8

    frozenSet = set()

    vecDist = bmvd.BinaryMemorylessVectorDistribution(uLen)
    for i in range(uLen):
        vecDist.probs[i][0] = 0.11
        vecDist.probs[i][1] = 0.89

    encoder = PolarEncoderDecoder.PolarEncoderDecoder(uLen, frozenSet, 0)

    TVvec = None
    Hvec = None
    numTrials = 100

    for rngSeed in range(numTrials):
        (encodedVector, TVvecTemp, HvecTemp) = encoder.genieSingleEncodeSimulatioan(vecDist, rngSeed)
        if  TVvec is None:
            TVvec = TVvecTemp
            Hvec = HvecTemp
        else:
            assert( len(TVvec) == len(TVvecTemp) )
            for i in range(len(TVvec)):
                TVvec[i] += TVvecTemp[i]
                Hvec[i] += HvecTemp[i]

    Hsum = 0.0
    for i in range(len(TVvec)):
        TVvec[i] /= numTrials
        Hvec[i] /= numTrials
        Hsum += Hvec[i]

    print( TVvec )
    print( Hvec )
    print( Hsum /len(Hvec) )

def genieEncodeDecodeSimulation(length, xyDistribution, numberOfTrials): 

    rngSeed = 0

    useTrellis = True
    # useTrellis = False

    xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

    xDistribution.probs.append( [-1.0,-1.0] )
    for x in range(2):
        xDistribution.probs[0][x] = xyDistribution.calcXMarginal(x)

    xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)

    frozenSet = set()
    TVvec = None
    HEncvec = None
    HDecvec = None

    encDec = PolarEncoderDecoder.PolarEncoderDecoder(length, frozenSet, rngSeed)

    # DO *NOT* SET SEED TO 0, AS THIS TAKES THE SYSTEM CLOCK AS A SEED!!!
    for rngSeed in range(1, numberOfTrials+1):
        (encodedVector, TVvecTemp, HencvecTemp) = encDec.genieSingleEncodeSimulatioan(xVectorDistribution, rngSeed)

        receivedWord = []
        # TODO: move this somewhere else (this code is duplicated, which is bad!)
        for j in range(length):
            x = encodedVector[j]

            rand = random.random()
            probSum = 0.0

            for y in range(len(xyDistribution.probs)):
                if probSum + xyDistribution.probXGivenY(x,y) >= rand:
                    receivedWord.append(y)
                    # print("x = ", x, ", y = ", y, " probXGivenY(x,y) = ", xyDistribution.probXGivenY(x,y), ", rand = ", rand)
                    break
                else:
                    probSum += xyDistribution.probXGivenY(x,y)

        if useTrellis:
            xyVectorDistribution =  xyDistribution.makeBinaryTrellisDistribution(length, receivedWord)
        else:
            xyVectorDistribution = xyDistribution.makeBinaryMemorylessVectorDistribution(length, receivedWord)
        
        (decodedVector, PevecTemp, HdecvecTemp) = encDec.genieSingleDecodeSimulatioan(xVectorDistribution, xyVectorDistribution, rngSeed)

        if  TVvec is None:
            TVvec = TVvecTemp
            Pevec = PevecTemp
            HEncvec = HencvecTemp
            HDecvec = HdecvecTemp
        else:
            assert( len(TVvec) == len(TVvecTemp) )
            for i in range(len(TVvec)):
                TVvec[i] += TVvecTemp[i]
                Pevec[i] += PevecTemp[i]
                HEncvec[i] += HencvecTemp[i]
                HDecvec[i] += HdecvecTemp[i]

    HEncsum = 0.0
    HDecsum = 0.0
    for i in range(len(TVvec)):
        TVvec[i] /= numberOfTrials
        Pevec[i] /= numberOfTrials
        HEncvec[i] /= numberOfTrials
        HDecvec[i] /= numberOfTrials
        HEncsum += HEncvec[i]
        HDecsum += HDecvec[i]

    print( "TVVec = ", TVvec )
    print( "pevec = ", Pevec )
    print( "HEncvec = ", HEncvec )
    print( "HDecvec = ", HDecvec )
    print( "Normalized HEncsum = ",  HEncsum /len(HEncvec) )
    print( "Normalized HDecsum = ", HDecsum /len(HDecvec) )

    TVPlusPeVec = []

    for i in range(len(TVvec)):
        TVPlusPeVec.append(TVvec[i] + Pevec[i])

    errorUpperBound = 0.1

    sortedIndices = sorted(range(len(TVPlusPeVec)), key=lambda k: TVPlusPeVec[k]) 

    print( sortedIndices )

    errorSum = 0.0
    indexInSortedIndicesArray = -1
    frozenSet = set()

    while errorSum < errorUpperBound and indexInSortedIndicesArray + 1 < len(TVPlusPeVec):
        i = sortedIndices[indexInSortedIndicesArray + 1]
        if TVPlusPeVec[i] + errorSum <= errorUpperBound:
            errorSum += TVPlusPeVec[i]
            indexInSortedIndicesArray += 1
        else:
            break

    for j in range(indexInSortedIndicesArray+1,  len(TVPlusPeVec)):
        i = sortedIndices[j]
        frozenSet.add(i)

    print( frozenSet )
    print( 1.0 - len(frozenSet) / len(HEncvec) )

    

def testEncode():
    uLen = 8

    vecDist = bmvd.BinaryMemorylessVectorDistribution(uLen)
    frozenSet = set()

    for i in range(uLen):
        vecDist.probs[i][0] = 0.1
        vecDist.probs[i][1] = 0.9

    for i in range(uLen):
        frozenSet.add(i)

    rngSeed = 1
    encoder = encdec.PolarEncoderDecoder(uLen, frozenSet, rngSeed)

    information = []


    encodedVector = encoder.encode(vecDist, information)

    print( "got here" )
    print( encodedVector )

def encodeDecodeSimulation(length, frozenSet, xyDistribution, numberOfTrials): 

    rngSeed = 0
    misdecodedWords = 0

    # useTrellis = True
    useTrellis = False

    xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

    xDistribution.probs.append( [-1.0,-1.0] )
    for x in range(2):
        xDistribution.probs[0][x] = xyDistribution.calcXMarginal(x)

    xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)

    encDec = PolarEncoderDecoder.PolarEncoderDecoder(length, frozenSet, rngSeed)

    random.seed(1)

    for t in range(numberOfTrials):
        information = []
        for i in range( encDec.k ):
            inf = 0 if random.random() < 0.5 else 1
            information.append(inf)

        codeword = encDec.encode(xVectorDistribution, information)
        receivedWord = []

        # TODO: move this somewhere else (this code is duplicated, which is bad!)
        for j in range(length):
            x = codeword[j]

            rand = random.random()
            probSum = 0.0

            for y in range(len(xyDistribution.probs)):
                if probSum + xyDistribution.probXGivenY(x,y) >= rand:
                    receivedWord.append(y)
                    # print("x = ", x, ", y = ", y, " probXGivenY(x,y) = ", xyDistribution.probXGivenY(x,y), ", rand = ", rand)
                    break
                else:
                    probSum += xyDistribution.probXGivenY(x,y)

        if useTrellis:
            xyVectorDistribution =  xyDistribution.makeBinaryTrellisDistribution(length, receivedWord)
        else:
            xyVectorDistribution = xyDistribution.makeBinaryMemorylessVectorDistribution(length, receivedWord)

        (decodedVector, decodedInformation) = encDec.decode(xVectorDistribution, xyVectorDistribution)

        for i in range( encDec.k ):
            if information[i] != decodedInformation[i]:
                misdecodedWords += 1
                print( t, ") error, transmitted information: ", information, ", decoded information: ", decodedInformation, ", transmitted codeword: ", codeword, ", received word: ", receivedWord )
                break

    print( "Error probability = ", misdecodedWords, "/", numberOfTrials, " = ", misdecodedWords/numberOfTrials )

# testEncode()

p = 0.11
L = 1000
n = 7
N = 2 ** n
# upperBoundOnErrorProbability = 1.0

xDistribution = None
xyDistribution = BinaryMemorylessDistribution.makeBSC(p)

# frozenSet = calcFrozenSet_degradingUpgrading(n, L, upperBoundOnErrorProbability, xDistribution, xyDistribution)
#
# print("Rate = ", N - len(frozenSet), "/", N, " = ", (N - len(frozenSet)) / N)

# numberOfTrials = 2000

genieEncodeDecodeSimulation(N, xyDistribution, L)
