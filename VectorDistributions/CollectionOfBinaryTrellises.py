import Guardbands
import VectorDistribution
from VectorDistributions import BinaryMemorylessVectorDistribution
from VectorDistributions import BinaryTrellis


class CollectionOfBinaryTrellises(VectorDistribution.VectorDistribution):
    def __init__(self, length, numberOfTrellises):
        """Initialize a collection of empty trellises

        Args:
            length (int): the number of inputs (not including the guard bands). That is, 2**n, where n is the number of polar transforms.

            numberOfTrellises (int): the number of trellises in the collection, which must be a power of 2.

        Returns:
            None
        """
        assert (length > 0)
        assert (numberOfTrellises > 0)

        assert (length % numberOfTrellises == 0)

        self.length = length
        self.numberOfTrellises = numberOfTrellises
        self.trellisLength = length // numberOfTrellises

        self.trellises = []

        for i in range(numberOfTrellises):
            newTrellis = BinaryTrellis.BinaryTrellis(self.trellisLength)
            self.trellises.append(newTrellis)

    def __len__(self):
        return self.length

    def __str__(self):
        s = "This collection of Binary trellis contains " + str(
            self.numberOfTrellises) + " trellises. Each trellis has input length " + str(
            self.trellisLength) + ". These trellises are\n"

        for trellis in self.trellises:
            s += "***\n"
            s += trellis.toString()

        s += "***\n"
        return s

    def minusTransform(self):
        return self.__miusPlusTransform()

    def plusTransform(self, decisionVector):
        return self.__miusPlusTransform(decisionVector)

    def __miusPlusTransform(self, decisionVector=None):
        assert (self.length % 2 == 0)

        if decisionVector is not None:
            decisionVectorSubLength = len(decisionVector) // self.numberOfTrellises

        if self.length // 2 > self.numberOfTrellises:
            newCollectionOfBinaryTrellises = CollectionOfBinaryTrellises(self.length // 2, self.numberOfTrellises)
            for i in range(self.numberOfTrellises):
                newCollectionOfBinaryTrellises.trellises[i] = self.trellises[i].minusTransform() if (
                            decisionVector is None) else self.trellises[i].plusTransform(
                    decisionVector[i * decisionVectorSubLength:(i + 1) * decisionVectorSubLength])
            return newCollectionOfBinaryTrellises
        else:
            assert (self.length // 2 == self.numberOfTrellises)

            newBinaryMemorylessVectorDistribution = BinaryMemorylessVectorDistribution.BinaryMemorylessVectorDistribution(
                self.numberOfTrellises)

            for i in range(self.numberOfTrellises):
                tempTrellis = self.trellises[i].minusTransform() if (decisionVector is None) else self.trellises[
                    i].plusTransform(decisionVector[i * decisionVectorSubLength:(i + 1) * decisionVectorSubLength])
                marginalizedProbs = tempTrellis.calcMarginalizedProbabilities(normalize=False)

                for x in range(2):
                    newBinaryMemorylessVectorDistribution.probs[i][x] = marginalizedProbs[x]

            return newBinaryMemorylessVectorDistribution

    def calcMarginalizedProbabilities(self):
        """Since a collection of trellises is eventually colapsed to a simple joint distribution after enough transforms, we should not get to this point
        """

        assert (False)

    def calcNormalizationVector(self):
        normalization = []

        for i in range(self.numberOfTrellises):
            trellisNormalization = self.trellises[i].calcNormalizationVector()
            normalization.append(trellisNormalization)

        return normalization

    def normalize(self, normalization):
        assert (len(normalization) == self.numberOfTrellises)

        for i in range(self.numberOfTrellises):
            self.trellises[i].normalize(normalization[i])


def buildCollectionOfBinaryTrellises_uniformInput_deletion(receivedWord, deletionProb, xi, n, n0,
                                                           numberOfOnesToAddAtBothEndsOfGuardbands, verbosity=0):
    trimmedSubwords = Guardbands.removeDeletionGuardBands(receivedWord, n, n0)

    trellisLength = 2 ** n0
    numberOfTrellises = 2 ** (n - n0)
    totalLength = 2 ** n

    collection = CollectionOfBinaryTrellises(totalLength, numberOfTrellises)
    collection.trellises = []

    trimmedZerosAtEdges = True
    if verbosity > 0:
        print("trimmed subwords")
    for subword in trimmedSubwords:
        if verbosity > 0:
            print(subword)

        tempTrellis = BinaryTrellis.buildTrellis_uniformInput_deletion(subword, trellisLength, deletionProb,
                                                                       trimmedZerosAtEdges,
                                                                       numberOfOnesToAddAtBothEndsOfGuardbands)
        collection.trellises.append(tempTrellis)

    return collection
