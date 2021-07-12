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

    # TODO: This is sub-optimal. Change this to a calculation similar to what we do with the genie.
    # TODO: Actually, factor common code out into a new function
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

def genieEncodeDecodeSimulation(length, make_xVectorDistribution, make_codeword, simulateChannel, make_xyVectrorDistribution, numberOfTrials, errorUpperBoundForFrozenSet):
    """Run a genie encoder and corresponding decoder, and return frozen set

    Args:
       length (int): the number of indices in the polar transformed vector

       make_xVectorDistribution (function): return xVectorDistribution, and takes no arguments

       make_codeword (function): make a codeword out of the encodedVector (for example, by doing nothing, or by adding guard bands)

       simulateChannel (function): transfroms a codeword to a recieved word, using the current state of the random number generator

       make_xyVectrorDistribution (function): return xyVectorDistribution, as a function of the received word

       numberOfTrials (int): number of Monte-Carlo simulations
    """

    rngSeed = 0

    xVectorDistribution = make_xVectorDistribution()

    frozenSet = set()
    TVvec = None
    HEncvec = None
    HDecvec = None

    encDec = PolarEncoderDecoder.PolarEncoderDecoder(length, frozenSet, rngSeed)

    # DO *NOT* SET SEED TO 0, AS THIS TAKES THE SYSTEM CLOCK AS A SEED!!!
    for rngSeed in range(1, numberOfTrials+1):
        (encodedVector, TVvecTemp, HencvecTemp) = encDec.genieSingleEncodeSimulatioan(xVectorDistribution, rngSeed)

        codeword = make_codeword(encodedVector)

        receivedWord = simulateChannel(codeword)

        xyVectorDistribution = make_xyVectrorDistribution(receivedWord)

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

    sortedIndices = sorted(range(len(TVPlusPeVec)), key=lambda k: TVPlusPeVec[k]) 

    print( sortedIndices )

    errorSum = 0.0
    indexInSortedIndicesArray = -1
    frozenSet = set()

    while errorSum < errorUpperBound and indexInSortedIndicesArray + 1 < len(TVPlusPeVec):
        i = sortedIndices[indexInSortedIndicesArray + 1]
        if TVPlusPeVec[i] + errorSum <= errorUpperBoundForFrozenSet:
            errorSum += TVPlusPeVec[i]
            indexInSortedIndicesArray += 1
        else:
            break

    for j in range(indexInSortedIndicesArray+1,  len(TVPlusPeVec)):
        i = sortedIndices[j]
        frozenSet.add(i)

    print( frozenSet )
    print( 1.0 - len(frozenSet) / len(HEncvec) )

    return frozenSet

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

def encodeDecodeSimulation(length, make_xVectorDistribution, make_codeword, simulateChannel, make_xyVectrorDistribution, numberOfTrials, frozenSet, rngSeed=1):
    """Run a polar encoder and a corresponding decoder (SC, not SCL)

    Args:
       length (int): the number of indices in the polar transformed vector

       make_xVectorDistribution (function): return xVectorDistribution, and takes no arguments

       make_codeword (function): make a codeword out of the encodedVector (for example, by doing nothing, or by adding guard bands)

       simulateChannel (function): transfroms a codeword to a recieved word, using the current state of the random number generator

       make_xyVectrorDistribution (function): return xyVectorDistribution, as a function of the received word

       frozenSet (set): the set of (dynamically) frozen indices
    """

    misdecodedWords = 0

    xVectorDistribution = make_xVectorDistribution()

    encDec = PolarEncoderDecoder.PolarEncoderDecoder(length, frozenSet, rngSeed)

    random.seed(1)

    for t in range(numberOfTrials):
        information = []
        for i in range( encDec.k ):
            inf = 0 if random.random() < 0.5 else 1
            information.append(inf)

        encodedVector = encDec.encode(xVectorDistribution, information)

        codeword = make_codeword(encodedVector)

        receivedWord = simulateChannel(codeword)

        xyVectorDistribution = make_xyVectrorDistribution(receivedWord)

        (decodedVector, decodedInformation) = encDec.decode(xVectorDistribution, xyVectorDistribution)

        for i in range( encDec.k ):
            if information[i] != decodedInformation[i]:
                misdecodedWords += 1
                print( t, ") error, transmitted information: ", information, ", decoded information: ", decodedInformation, ", transmitted codeword: ", codeword, ", received word: ", receivedWord )
                break

    print( "Error probability = ", misdecodedWords, "/", numberOfTrials, " = ", misdecodedWords/numberOfTrials )

# TODO: move these into BinaryMemorylessDistribution, and perhaps make them class methods

def make_xVectorDistribuiton_fromBinaryMemorylessDistribution(xyDistribution, length):
    def make_xVectorDistribuiton():
        xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

        xDistribution.probs.append( [-1.0,-1.0] )
        for x in range(2):
            xDistribution.probs[0][x] = xyDistribution.calcXMarginal(x)
        
        xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)
        return xVectorDistribution
    return make_xVectorDistribuiton

def make_codeword_noprocessing(encodedVector):
    return encodedVector

def simulateChannel_fromBinaryMemorylessDistribution(xyDistribution):
    def simulateChannel(codeword):
        receivedWord = []
        length = len(codeword)

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
        
        return receivedWord

    return simulateChannel

def make_xyVectorDistribution_fromBinaryMemorylessDistribution(xyDistribution):
    def make_xyVectrorDistribution(receivedWord):
        length = len(receivedWord)
        useTrellis = False

        if useTrellis:
            xyVectorDistribution =  xyDistribution.makeBinaryTrellisDistribution(length, receivedWord)
        else:
            xyVectorDistribution = xyDistribution.makeBinaryMemorylessVectorDistribution(length, receivedWord)

        return xyVectorDistribution
    return make_xyVectrorDistribution

p = 0.11
L = 1000
n = 6
N = 2 ** n

upperBoundOnErrorProbability = 0.1

xDistribution = None
xyDistribution = BinaryMemorylessDistribution.makeBSC(p)

frozenSet = calcFrozenSet_degradingUpgrading(n, L, upperBoundOnErrorProbability, xDistribution, xyDistribution)

print("Rate = ", N - len(frozenSet), "/", N, " = ", (N - len(frozenSet)) / N)

numberOfTrials = 2000

make_xVectorDistribuiton = make_xVectorDistribuiton_fromBinaryMemorylessDistribution(xyDistribution, N)
make_codeword = make_codeword_noprocessing
simulateChannel = simulateChannel_fromBinaryMemorylessDistribution(xyDistribution)
make_xyVectorDistribution = make_xyVectorDistribution_fromBinaryMemorylessDistribution(xyDistribution)

encodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, numberOfTrials, frozenSet)

# genieEncodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel, make_xyVectorDistribution, L)
