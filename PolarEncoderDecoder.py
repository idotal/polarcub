def class PolarEncoderDecoder():
    def init(self, frozenSet, rngSeed ):
        self.rngSeed = rngSeed
        self.frozenSet = frozenSet

    def encode(self, vectorDistribution, information):

        pass # TODO

    def decode(self, vectorDistribution):
        pass # TODO

    # returns (encodedVector, next_uIndex, next_informationVectorIndex)
    def recursiveEncode(self, vectorDistribution, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers):
        if len(vectorDistribution) == 1:
            encodedVector = np.vector(1)
            if frozenSet[uIndex] == informationIndex:
                encodedVector[0] = information[informationVectorIndex]
                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex + 1
                return (encodedVector, next_uIndex, next_informationVectorIndex)
            else:
                marginalizedVector = vectorDistribution.marginalizedVector()
                if marginalizedVector[0] >= randomlyGeneratedNumbers[uIndex]:
                    encodedVector[0] = 0
                else:
                    encodedVector[0] = 1
                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex
                return (encodedVector, next_uIndex, next_informationVectorIndex)
        else:
            # TODO: stopped here










