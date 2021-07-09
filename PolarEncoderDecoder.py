import numpy as np
import random as rand
from enum import Enum

class uIndexType(Enum):
    frozen = 0
    information = 1

class PolarEncoderDecoder():
    def __init__(self, length, frozenSet, rngSeed): # length is the length of the U vector, if rngSeed is set to 0, then we freeze all frozen bits to zero
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
        """Encode k information bits according to a-priori input distribution

        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            information (numpy array of Int64): the k information bits to encode

        Returns:
            The encoded vector (reverse polar transform of the U)
        """

        uIndex = 0
        informationVectorIndex = 0
        assert( len(xVectorDistribution) == self.length )


        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers, xVectorDistribution)

        assert( next_uIndex == len(encodedVector) == len(xVectorDistribution) )
        assert( next_informationVectorIndex == len(information) )

        return encodedVector

    # returns (encodedVector, information)
    def decode(self, xVectorDistribution, xyVectorDistribution):
        """Decode k information bits according to a-priori input distribution and a-posteriori input distribution

        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            xyVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-posteriori entries for P(X=0) and P(X=1). That is, entry i contains P(X=0,Y=y_i) and P(X=1,Y=y_i).

        Returns:
            (encodedVector, information): The encoded vector (reverse polar transform of the U) and the corresponding information bits.
        """

        uIndex = 0
        informationVectorIndex = 0

        information = np.empty(self.k, np.int64)
        information[:] = -1

        assert( len(xVectorDistribution) == len(xyVectorDistribution) == self.length )


        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers, xVectorDistribution, xyVectorDistribution)

        assert( next_uIndex == len(encodedVector) == self.length )
        assert( next_informationVectorIndex == len(information) )

        return (encodedVector, information)

    def geniePreSteps(self, simulationSeed):
        self.backupFrozenSet = self.frozenSet
        self.backupSeed = self.rngSeed

        self.frozenSet = { i for i in range(self.length) }
        self.rngSeed = simulationSeed
        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()

    def geniePostSteps(self):
        self.frozenSet = self.backupFrozenSet 
        self.rngSeed = self.backupSeed 
        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()

    def genieSingleDecodeSimulatioan(self, xVectorDistribution, xyVectorDistribution, simulationSeed):
        """Pick up statistics of a single decoding run
        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            xyVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-posteriori entries for P(X=0) and P(X=1). That is, entry i contains P(X=0,Y=y_i) and P(X=1,Y=y_i).
            simulationSeed (int): The seed used to randomly pick the value of the u_i

        Returns:
            (decodedVector, Pe): a pair of arrays. The first array is the codeword we have produced. Entry i of the second array is the probability of error, min{P(U_i=0|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1}), P(U_i=1|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1})}.
        """

        #TODO: get the encoding simulation running, and then do decoding
        pass

    def genieSingleEncodeSimulatioan(self, xVectorDistribution, xyVectorDistribution, simulationSeed):
        """Pick up statistics of a single encoding run
        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            simulationSeed (int): The seed used to randomly pick the value of the u_i

        Returns:
            (encodedVector, K): a pair of arrays. The first array is the codeword we have produced. Entry i of the second array is the the total variation |P(U_i=0|U_0^{i-1} = u_0^{i-1})-P(U_i=1|U_0^{i-1} = u_0^{i-1})|.
        """

        marginalizedUProbs = []
        uIndex = 0
        informationVectorIndex = 0
        infromation = []

        self.geniePreSteps()

        assert( len(xVectorDistribution) == self.length )

        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, genieRandomlyGeneratedNumbers, vectorDistribution, None, marginalizedUProbs)

        assert( next_uIndex == len(encodedVector) == len(xVectorDistribution) )
        assert( next_informationVectorIndex == len(information) == 0 )
        assert( len(marginalizedUProbs) ==  self.length )

        Pevec = []

        for probPair in marginalizedUProbs:
            Pevec.append( min( probPair[0], probPair[1] )

        # return things to the way they were
        self.geniePostSteps()

        return (encodedVector, Pevec)

    def recursiveEncodeDecode(self, information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xVectorDistribution, xyVectorDistribution=None, marginalizedUProbs=None):
        """Encode/decode according to supplied vector distributions

        Args:
            information (numpy array of Int64): an array of inforamation bits to either read from when encoding or write to when decoding

            uIndex (int): the first releavant index in the polar transformed U vector of the *whole* codeword (non-recursive)

            informationVectorIndex (int): the first relevant index in the information vector associated with the *whole* codeword (non-recursive)

            randomlyGeneratedNumbers (numpy array of floats): the random numbers generated by rngSeed, corresponding to the *whole* codeword (non-recursive)

            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector (whose length is a function of the recursion depth) with a-priori entries for P(X=0) and P(X=1)

            xyVectorDistribution (VectorDistribution): in a memorylyess setting, this is essentially a vector (whose length is a function of the recursion depth) with a-posteriori entries for P(x=0) and P(x=1). A None value means we are encoding.

            marginalizedUProbs (empty array, or None): If not None, we populate (return) this array so that if xyVectorDistribution is None (encoding), then marginalizedUProbs[i][x] = P(U_i=x|U_0^{i-1} = \hat{u}_0^{i-1}). Otherwise, marginalizedUProbs[i][x] = P(U_i=x|U_0^{i-1} = \hat{u}_0^{i-1}, Y_0^{N-1} = y_0^{N-1}). For genie decoding, we will have \hat{u}_i = u_i, as the frozen set contains all indices.

        Returns:
            (encodedVector, next_uIndex, next_informationVectorIndex): the recursive encoding of the relevant part of the information vector, as well as updated values for the parameters uIndex and informationVectorIndex
        """
    
        # By defualt, we assume encoding, and add small corrections for decoding.

        encodedVector = np.empty(len(xVectorDistribution), np.int64)
        encodedVector[:] = -1

        if len(xVectorDistribution) == 1:
            if self.frozenOrInformation[uIndex] == uIndexType.information:

                # For decoding
                if xyVectorDistribution is not None:
                    marginalizedVector = xyVectorDistribution.calcMarginalizedProbabilities()
                    # print( marginalizedVector )
                    information[informationVectorIndex] = 0 if marginalizedVector[0] >= marginalizedVector[1] else 1

                encodedVector[0] = information[informationVectorIndex]
                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex + 1
            else:
                marginalizedVector = xVectorDistribution.calcMarginalizedProbabilities()
                if marginalizedVector[0] >= randomlyGeneratedNumbers[uIndex]:
                    encodedVector[0] = 0
                else:
                    encodedVector[0] = 1


                next_uIndex = uIndex + 1
                next_informationVectorIndex = informationVectorIndex

            # should we return the marginalized probabilities?
            if marginalizedUProbs is not None:
                if xyVectorDistribution is not None:
                    marginalizedVector = xyVectorDistribution.calcMarginalizedProbabilities()
                else:
                    marginalizedVector = xVectorDistribution.calcMarginalizedProbabilities()
                marginalizedUProbs.append( [marginalizedVector[0], marginalizedVector[1]] )

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

            (minusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xMinusVectorDistribution, xyMinusVectorDistribution)

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
            (plusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xPlusVectorDistribution, xyPlusVectorDistribution)

            halfLength = len(xVectorDistribution) // 2

            for halfi in range(halfLength):
                encodedVector[2*halfi] = (minusEncodedVector[halfi] + plusEncodedVector[halfi]) % 2
                encodedVector[2*halfi + 1] = plusEncodedVector[halfi]

            return (encodedVector, next_uIndex, next_informationVectorIndex)

