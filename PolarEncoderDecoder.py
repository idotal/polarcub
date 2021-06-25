import numpy as np
import random as rand
from enum import Enum

class uIndexType(Enum):
    frozen = 0
    information = 1

class PolarEncoderDecoder():
    def __init__(self, length, frozenSet, rngSeed ): # length is the length of the U vector, if rngSeed is set to 0, then we freeze all frozen bits to zero
        self.rngSeed = rngSeed
        self.frozenSet = frozenSet
        self.length = length
        
        self.frozenOrInformation = np.empty(length, uIndexType)
        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()


    def initializeFrozenOrInformationAndRandomlyGeneratedNumbers(self):
        self.k = 0
        for i in range(self.length):
            if i in self.frozenSet:
                self.frozenOrInformation[i] = uIndexType.frozen
            else:
                self.frozenOrInformation[i] = uIndexType.information
                self.k += 1

        self.randomlyGeneratedNumbers = np.empty(self.length)
        self.randomlyGeneratedNumbers[:] = np.nan

        if self.rngSeed != 0:
            rand.seed(self.rngSeed)

            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = rand.random()
        else:
            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = 1.0

    # returns encodedVector
    def encode(self, xVectorDistribution, information):
        uIndex = 0
        informationVectorIndex = 0
        assert( len(xVectorDistribution) == self.length )


        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(xVectorDistribution, information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers)

        assert( next_uIndex == len(encodedVector) == len(xVectorDistribution) )
        assert( next_informationVectorIndex == len(information) )

        return encodedVector

    # returns (encodedVector, information)
    def decode(self, xVectorDistribution, xyVectorDistriubiton):
        uIndex = 0
        informationVectorIndex = 0

        information = np.empty(self.k, np.int64)
        information[:] = -1

        assert( len(xVectorDistribution) == len(xyVectorDistribution) == self.length )


        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(xVectorDistribution, information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers, xyVectorDistribution)

        assert( next_uIndex == len(encodedVector) == len(vectorDistribution) )
        assert( next_informationVectorIndex == len(information) )

        return (encodedVector, information)

    # returns (encodedVector, next_uIndex, next_informationVectorIndex)
    def recursiveEncode(self, xVectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xyVectorDistribution=None):
        encodedVector = np.empty(len(xVectorDistribution), np.int64)
        encodedVector[:] = -1

        if len(xVectorDistribution) == 1:
            if self.frozenOrInformation[uIndex] == uIndexType.information:

                # For decoding
                if xyVectorVectorDistribution != None:
                    marginalizedVector = xyVectorDistribution.calcMarginalizedProbabilities()
                    encodedVector[0] = 0 if marginalizedVector[0] >= marginalizedVector[1] else 1

                encodedVector[0] = information[informationVectorIndex]
                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex + 1
                return (encodedVector, next_uIndex, next_informationVectorIndex)
            else:
                marginalizedVector = xVectorDistribution.calcMarginalizedProbabilities()
                print( marginalizedVector[0] )
                if marginalizedVector[0] >= randomlyGeneratedNumbers[uIndex]:
                    encodedVector[0] = 0
                else:
                    encodedVector[0] = 1

                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex
                return (encodedVector, next_uIndex, next_informationVectorIndex)
        else:
            xMinusVectorDistribution = xVectorDistribution.minusTransform()
            normalization = xMinusVectorDistribution.calcNormalizationVector()
            xMinusVectorDistribution.normalize(normalization)

            # For decoding
            if xyVectorDistribution != None:
                xyMinusVectorDistribution = xyVectorDistribution.minusTransform()
                normalization = xyMinusVectorDistribution.calcNormalizationVector()
                xyMinusVectorDistribution.normalize(normalization)
            else:
                xyMinusVectorDistribution = None

            (minusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(xMinusVectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xyMinusVectorDistribution)

            xPlusVectorDistribution = xVectorDistribution.plusTransform(minusEncodedVector)
            normalization = xPlusVectorDistribution.calcNormalizationVector()
            xPlusVectorDistribution.normalize(normalization)

            # For decoding
            if xyVectorDistribution != None:
                xyPlusVectorDistribution = xyVectorDistribution.plusTransform(minusEncodedVector)
                normalization = xyPlusVectorDistribution.calcNormalizationVector()
                xyPlusVectorDistribution.normalize(normalization)
            else:
                xyPlusVectorDistribution = None

            uIndex = next_uIndex
            informationVectorIndex = next_informationVectorIndex
            (plusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(xPlusVectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xyPlusVectorDistribution)

            halfLength = len(xVectorDistribution) // 2

            for halfi in range(halfLength):
                encodedVector[2*halfi] = (minusEncodedVector[halfi] + plusEncodedVector[halfi]) % 2
                encodedVector[2*halfi + 1] = plusEncodedVector[halfi]

            return (encodedVector, next_uIndex, next_informationVectorIndex)

