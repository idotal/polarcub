import numpy as np
from ScalarDistributions import BinaryMemorylessDistribution
import random
from enum import Enum

import sys

class uIndexType(Enum):
    frozen = 0
    information = 1

class BinaryPolarEncoderDecoder():
    def __init__(self, length, frozenSet, commonRandomnessSeed): # length is the length of the U vector, if rngSeed is set to -1, then we freeze all frozen bits to zero
        self.commonRandomnessSeed = commonRandomnessSeed
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

        commonRandomnessRNG = random.Random()
        if self.commonRandomnessSeed != -1:
            commonRandomnessRNG.seed(self.commonRandomnessSeed)

            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = commonRandomnessRNG.random()
        else:
            for i in range(self.length):
                self.randomlyGeneratedNumbers[i] = 1.0

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

    def geniePreSteps(self, genieSingleRunSeed):
        self.backupFrozenSet = self.frozenSet
        self.backupCommonRandomnessSeed = self.commonRandomnessSeed

        self.frozenSet = { i for i in range(self.length) }
        self.commonRandomnessSeed = genieSingleRunSeed
        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()

    def geniePostSteps(self):
        self.frozenSet = self.backupFrozenSet
        self.commonRandomnessSeed = self.backupCommonRandomnessSeed
        self.initializeFrozenOrInformationAndRandomlyGeneratedNumbers()

    def genieSingleDecodeSimulatioan(self, xVectorDistribution, xyVectorDistribution, genieSingleRunSeed, trustXYProbs):
        """Pick up statistics of a single decoding run
        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            xyVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-posteriori entries for P(X=0) and P(X=1). That is, entry i contains P(X=0,Y=y_i) and P(X=1,Y=y_i).
            genieSingleRunSeed (int): The seed used to randomly pick the common randomness for this genie run

            trustXYProbs (bool): Do we trust the probabilities of U_i=0 and U_i = 1 given past U and all Y (we usually should), or don't we (in case we have guard bands, which can be parsed wrong, and then result in garbage probs).

        Returns:
            (decodedVector, Pe, H): a triplet of arrays. The first array is the codeword we have produced. If we trust the probabilities, entry i of the second array is the probability of error, min{P(U_i=0|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1}), P(U_i=1|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1})}. Entry i of the third array is the entropy, -P(U_i=0|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1}) * log_2(P(U_i=0|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1})) - P(U_i=1|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1}) * log_2(P(U_i=1|U_0^{i-1} = u_0^{i-1}, Y_0^{N-1} = y_0^{N-1})). If we don't trust the probabilities, entry i of the second array is 1.0, and the third array is empty.
        """

        marginalizedUProbs = []
        uIndex = 0
        informationVectorIndex = 0
        information = []

        self.geniePreSteps(genieSingleRunSeed)

        assert( len(xVectorDistribution) == self.length )

        # print(xyVectorDistribution)

        (decodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers, xVectorDistribution, xyVectorDistribution, marginalizedUProbs)

        # print( decodedVector, marginalizedUProbs )
        assert( next_uIndex == len(decodedVector) == len(xVectorDistribution) )
        assert( next_informationVectorIndex == len(information) == 0 )
        assert( len(marginalizedUProbs) ==  self.length )

        Pevec = []
        Hvec = []

        if trustXYProbs == True:
            for probPair in marginalizedUProbs:
                Pevec.append( min( probPair[0], probPair[1]) )
                Hvec.append( BinaryMemorylessDistribution.eta(probPair[0]) + BinaryMemorylessDistribution.eta(probPair[1]) )
        else:
            Uvec = polarTransformOfBits( decodedVector )
            # print( "decodedVector = ", decodedVector, ", Uvec = ", Uvec )
            i = 0
            for probPair in marginalizedUProbs:
                decision = Uvec[i]
                # print( i, decision, probPair)
                i += 1
                if probPair[decision] > probPair[1 - decision]:
                    Pevec.append( 0.0 )
                elif probPair[decision] == probPair[1 - decision]:
                    Pevec.append( 0.5 )
                else:
                    Pevec.append( 1.0 )

                # Hvec.append( BinaryMemorylessDistribution.eta(probPair[0]) + BinaryMemorylessDistribution.eta(probPair[1]) )


        # return things to the way they were
        self.geniePostSteps()

        return (decodedVector, Pevec, Hvec)

    def genieSingleEncodeSimulatioan(self, xVectorDistribution, genieSingleRunSeed):
        """Pick up statistics of a single encoding run
        Args:
            xVectorDistribution (VectorDistribution): in a memoryless setting, this is essentially a vector with a-priori entries for P(X=0) and P(X=1)

            genieSingleRunSeed (int): The seed used to randomly pick the common randomness for this genie run

        Returns:
            (encodedVector, TVvec, H): a triplet of arrays. The first array is the codeword we have produced. Entry i of the second array is the the total variation |P(U_i=0|U_0^{i-1} = u_0^{i-1})-P(U_i=1|U_0^{i-1} = u_0^{i-1})|. Entry i of the third array is the entropy, -P(U_i=0|U_0^{i-1} = u_0^{i-1} * log_2(P(U_i=0|U_0^{i-1} = u_0^{i-1}) - P(U_i=1|U_0^{i-1} = u_0^{i-1} * log_2(P(U_i=1|U_0^{i-1} = u_0^{i-1})
        """

        marginalizedUProbs = []
        uIndex = 0
        informationVectorIndex = 0
        information = []

        self.geniePreSteps(genieSingleRunSeed)

        assert( len(xVectorDistribution) == self.length )

        (encodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, self.randomlyGeneratedNumbers, xVectorDistribution, None, marginalizedUProbs)

        # print( encodedVector, marginalizedUProbs )
        assert( next_uIndex == len(encodedVector) == len(xVectorDistribution) )
        assert( next_informationVectorIndex == len(information) == 0 )
        assert( len(marginalizedUProbs) ==  self.length )

        TVvec = []
        Hvec = []

        for probPair in marginalizedUProbs:
            TVvec.append( abs( probPair[0] - probPair[1]) )
            Hvec.append( BinaryMemorylessDistribution.eta(probPair[0]) + BinaryMemorylessDistribution.eta(probPair[1]) )

        # return things to the way they were
        self.geniePostSteps()

        return (encodedVector, TVvec, Hvec)

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

        # By default, we assume encoding, and add small corrections for decoding.

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

            (minusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xMinusVectorDistribution, xyMinusVectorDistribution, marginalizedUProbs)

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
            (plusEncodedVector, next_uIndex, next_informationVectorIndex) = self.recursiveEncodeDecode(information, uIndex, informationVectorIndex, randomlyGeneratedNumbers, xPlusVectorDistribution, xyPlusVectorDistribution, marginalizedUProbs)

            halfLength = len(xVectorDistribution) // 2

            for halfi in range(halfLength):
                encodedVector[2*halfi] = (minusEncodedVector[halfi] + plusEncodedVector[halfi]) % 2
                encodedVector[2*halfi + 1] = plusEncodedVector[halfi]

            return (encodedVector, next_uIndex, next_informationVectorIndex)

def encodeDecodeSimulation(length, make_xVectorDistribution, make_codeword, simulateChannel, make_xyVectrorDistribution, numberOfTrials, frozenSet, commonRandomnessSeed=1, randomInformationSeed=1, verbosity=0):
    """Run a polar encoder and a corresponding decoder (SC, not SCL)

    Args:
       length (int): the number of indices in the polar transformed vector

       make_xVectorDistribution (function): return xVectorDistribution, and takes no arguments

       make_codeword (function): make a codeword out of the encodedVector (for example, by doing nothing, or by adding guard bands)

       simulateChannel (function): transfroms a codeword to a recieved word, using the current state of the random number generator

       make_xyVectrorDistribution (function): return xyVectorDistribution, as a function of the received word

       frozenSet (set): the set of (dynamically) frozen indices

       commonRandomnessSeed (int): the seed used for defining the encoder/decoder common randomness

       randomInformationSeed (int): the seed used to create the random information to be encoded

    """

    misdecodedWords = 0

    xVectorDistribution = make_xVectorDistribution()

    encDec = BinaryPolarEncoderDecoder(length, frozenSet, commonRandomnessSeed)

    informationRNG = random.Random()
    informationRNG.seed(randomInformationSeed)

    # Note that we set a random seed, which is in charge of both setting the information bits as well as the channel output.
    for t in range(numberOfTrials):
        information = []
        for i in range( encDec.k ):
            inf = 0 if informationRNG.random() < 0.5 else 1
            information.append(inf)

        encodedVector = encDec.encode(xVectorDistribution, information)

        codeword = make_codeword(encodedVector)

        receivedWord = simulateChannel(codeword)

        # if t == 844:
        #     xyVectorDistribution = make_xyVectrorDistribution(receivedWord, verbosity=1)
        # else:
        #     xyVectorDistribution = make_xyVectrorDistribution(receivedWord)
        xyVectorDistribution = make_xyVectrorDistribution(receivedWord)

        (decodedVector, decodedInformation) = encDec.decode(xVectorDistribution, xyVectorDistribution)

        for i in range( encDec.k ):
            if information[i] != decodedInformation[i]:
                misdecodedWords += 1
                if verbosity > 0:
                    s = str(t) + ") error, transmitted inforamtion:\n" + str(information)
                    s += "\ndecoded information:\n" + str(decodedInformation)
                    s += "\nencoded vector before guard bands added:\n" + str(encodedVector)
                    s += "\ncodeword:\n" + str(codeword)
                    s += "\nreceived word:\n" + str(receivedWord)
                    print(s)
                    # print( t, ") error, transmitted information: ", information", ", decoded information: ", decodedInformation, ", transmitted codeword: ", codeword, ", received word: ", receivedWord )
                break

    print( "Error probability = ", misdecodedWords, "/", numberOfTrials, " = ", misdecodedWords/numberOfTrials )

def genieEncodeDecodeSimulation(length, make_xVectorDistribution, make_codeword, simulateChannel, make_xyVectrorDistribution, numberOfTrials, errorUpperBoundForFrozenSet, genieSeed, trustXYProbs=True, filename=None):
    """Run a genie encoder and corresponding decoder, and return frozen set

    Args:
       length (int): the number of indices in the polar transformed vector

       make_xVectorDistribution (function): return xVectorDistribution, and takes no arguments

       make_codeword (function): make a codeword out of the encodedVector (for example, by doing nothing, or by adding guard bands)

       simulateChannel (function): transfroms a codeword to a recieved word, using the current state of the random number generator

       make_xyVectrorDistribution (function): return xyVectorDistribution, as a function of the received word

       numberOfTrials (int): number of Monte-Carlo simulations

       errorUpperBoundForFrozenSet (float): choose a frozen set that will result in decoding error not more than this varialble

       genieSeed (int): the seed used by the genie to have different encoding/decoding common randomness in each run

       trustXYProbs (bool): Do we trust the probabilities of U_i=0 and U_i = 1 given past U and all Y (we usually should), or don't we (in case we have guard bands, which can be parsed wrong, and then result in garbage probs).
    """

    commonRandomnessSeed = 0 # doesn't matter, as this common randomness will not be used (will be replaced by a different seed by the genie)

    xVectorDistribution = make_xVectorDistribution()

    frozenSet = set()
    TVvec = None
    HEncvec = None
    HDecvec = None
    codewordLength = 0

    encDec = BinaryPolarEncoderDecoder(length, frozenSet, commonRandomnessSeed)
    genieSingleRunSeedRNG = random.Random()
    genieSingleRunSeedRNG.seed(genieSeed)

    for trialNumber in range(numberOfTrials):
        genieSingleRunSeed = genieSingleRunSeedRNG.randint(1, 1000000)

        (encodedVector, TVvecTemp, HencvecTemp) = encDec.genieSingleEncodeSimulatioan(xVectorDistribution, genieSingleRunSeed)

        codeword = make_codeword(encodedVector)

        receivedWord = simulateChannel(codeword)

        # print("codeword = ", codeword, ", receivedWord = ", receivedWord)

        xyVectorDistribution = make_xyVectrorDistribution(receivedWord)

        (decodedVector, PevecTemp, HdecvecTemp) = encDec.genieSingleDecodeSimulatioan(xVectorDistribution, xyVectorDistribution, genieSingleRunSeed, trustXYProbs)

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
                if trustXYProbs:
                    HDecvec[i] += HdecvecTemp[i]

    HEncsum = 0.0
    HDecsum = 0.0
    for i in range(len(TVvec)):
        TVvec[i] /= numberOfTrials
        Pevec[i] /= numberOfTrials
        HEncvec[i] /= numberOfTrials
        HEncsum += HEncvec[i]
        if trustXYProbs:
            HDecvec[i] /= numberOfTrials
            HDecsum += HDecvec[i]


    print( "TVVec = ", TVvec )
    print( "pevec = ", Pevec )
    print( "HEncvec = ", HEncvec )
    if trustXYProbs:
        print( "HDecvec = ", HDecvec )
    print( "Normalized HEncsum = ",  HEncsum /len(HEncvec) )
    if trustXYProbs:
        print( "Normalized HDecsum = ", HDecsum /len(HDecvec) )

    frozenSet = frozenSetFromTVAndPe(TVvec, Pevec, errorUpperBoundForFrozenSet)
    print( "code rate = ", (len(TVvec) - len(frozenSet)) /  len(codeword) )
    print( "codeword length = ", len(codeword) )

    if filename is not None:
        f = open( filename, "w" )
        s = "* " + ' '.join(sys.argv[:]) + "\n"
        f.write(s)

        for i in frozenSet:
            f.write(str(i))
            f.write("\n")

        s = "** number of trials = " + str(numberOfTrials) + "\n"
        f.write(s)
        s = "* (TotalVariation+errorProbability) * (number of trials)" + "\n"
        f.write(s)

        for i in range(len(TVvec)):
            s = "*** " + str(i) + " " + str((TVvec[i] + Pevec[i]) * numberOfTrials) +  "\n"
            f.write(s)

        f.close()

    return frozenSet

def polarTransformOfBits( xvec ):
    # print("xvec =", xvec)
    if len(xvec) == 1:
        return xvec
    else:
        assert( len(xvec) % 2 == 0 )

        vfirst = []
        vsecond = []
        for i in range((len(xvec) // 2)):
            vfirst.append( (xvec[2*i] + xvec[2*i+1]) % 2 )
            vsecond.append( xvec[2*i+1] )

        ufirst = polarTransformOfBits(vfirst)
        usecond = polarTransformOfBits(vsecond)

        transformed = []

        transformed.extend(ufirst)
        transformed.extend(usecond)

        # print( "transformed = ", transformed )
        return transformed

def frozenSetFromTVAndPe(TVvec, Pevec, errorUpperBoundForFrozenSet):

    TVPlusPeVec = []

    for i in range(len(TVvec)):
        TVPlusPeVec.append(TVvec[i] + Pevec[i])

    sortedIndices = sorted(range(len(TVPlusPeVec)), key=lambda k: TVPlusPeVec[k])

    # print( sortedIndices )

    errorSum = 0.0
    indexInSortedIndicesArray = -1
    frozenSet = set()

    while errorSum < errorUpperBoundForFrozenSet and indexInSortedIndicesArray + 1 < len(TVPlusPeVec):
        i = sortedIndices[indexInSortedIndicesArray + 1]
        if TVPlusPeVec[i] + errorSum <= errorUpperBoundForFrozenSet:
            errorSum += TVPlusPeVec[i]
            indexInSortedIndicesArray += 1
        else:
            break

    for j in range(indexInSortedIndicesArray+1,  len(TVPlusPeVec)):
        i = sortedIndices[j]
        frozenSet.add(i)

    print( "frozen set =", frozenSet )
    print( "fraction of non-frozen indices =", 1.0 - len(frozenSet) / len(TVPlusPeVec) )

    return frozenSet
