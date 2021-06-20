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


    def initializeFrozenOrInformationAndRandomlyGeneratedNumbers(self):
        for i in range(self.length):
            if i in self.frozenSet:
                self.frozenOrInformation[i] = uIndexType.frozen
            else:
                self.frozenOrInformation[i] = uIndexType.information

        self.randomlyGeneratedNumbers = np.empty(self.length)
        self.randomlyGeneratedNumbers[:] = np.nan

        if self.rngSeed != 0:
            rand.seed(self.rngSeed)

            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = rand.random()
        else:
            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = 1.0

    def encode(self, vectorDistribution, information):
        uIndex = 0
        informationVectorIndex = 0
        assert( len(vectorDistribution) == self.length )

        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()

        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(vectorDistribution, information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers)

        assert( next_uIndex == len(encodedVector) == len(vectorDistribution) )
        assert( next_informationVectorIndex == len(information) )

        return encodedVector

    def decode(self, vectorDistribution):
        pass # TODO

    # returns (encodedVector, next_uIndex, next_informationVectorIndex)
    def recursiveEncode(self, vectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers):
        encodedVector = np.empty(len(vectorDistribution), np.int64)
        encodedVector[:] = -1

        if len(vectorDistribution) == 1:
            if self.frozenOrInformation[uIndex] == uIndexType.information:
                encodedVector[0] = information[informationVectorIndex]
                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex + 1
                return (encodedVector, next_uIndex, next_informationVectorIndex)
            else:
                marginalizedVector = vectorDistribution.calcMarginalizedProbabilities()
                print( marginalizedVector[0] )
                if marginalizedVector[0] >= randomlyGeneratedNumbers[uIndex]:
                    encodedVector[0] = 0
                else:
                    encodedVector[0] = 1

                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex
                return (encodedVector, next_uIndex, next_informationVectorIndex)
        else:
            minusVectorDistribution = vectorDistribution.minusTransform()
            normalization = minusVectorDistribution.calcNormalizationVector()
            minusVectorDistribution.normalize(normalization)

            (minusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(minusVectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers)

            plusVectorDistribution = vectorDistribution.plusTransform(minusEncodedVector)
            normalization = plusVectorDistribution.calcNormalizationVector()
            plusVectorDistribution.normalize(normalization)

            uIndex = next_uIndex
            informationVectorIndex = next_informationVectorIndex
            (plusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncode(plusVectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers)

            halfLength = len(vectorDistribution) // 2

            for halfi in range(halfLength):
                encodedVector[2*halfi] = (minusEncodedVector[halfi] + plusEncodedVector[halfi]) % 2
                encodedVector[2*halfi + 1] = plusEncodedVector[halfi]

            return (encodedVector, next_uIndex, next_informationVectorIndex)

