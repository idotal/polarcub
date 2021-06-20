#! /usr/bin/env python3

import BinaryMemorylessVectorDistribution as bmvd
import PolarEncoderDecoder as encdec

def testEncode():
    uLen = 4

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

testEncode()

    


