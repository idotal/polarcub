import BinaryTrellis
import BinaryMemorylessVectorDistribution

class CollectionOfBinaryTrellises(VectorDistribution.VectorDistribution):
    def __init__(self, length, numberOfTrellises):
        """Initialize a collection of empty trellises

        Args:
            length (int): the number of inputs (not including the guard bands). That is, 2**n, where n is the number of polar transforms. Thus, the number of layers is length + 1.

            numberOfTrellises (int): the number of trellises in the collection, which must be a power of 2.

        Returns:
            None
        """
        assert( length > 0 )
        assert( numberOfTrellises > 0 )

        assert( length % numberOfTrellises == 0)

        self.length = length
        self.numberOfTrellises = numberOfTrellises
        self.trellisLength = length // numberOfTrellises

        self.trellises = []

    
        for i in range(numberOfTrellises):
            newTrellis = BinaryTrellis(trellisLength)
            self.trellises.append(newTrellis)

    def __len__(self):
        return self.length

    def minusTransform(self):
        return self.__miusPlusTransform()

    def plusTransform(self, decisionVector):
        return self.__miusPlusTransform(decisionVector)

    def __miusPlusTransform(self, decisionVector=None):
        assert( self.length % 2 == 0 )

        if self.length // 2 > self.numberOfTrellises:
            newCollectionOfBinaryTrellises = CollectionOfBinaryTrellises(self.lenth//2, self.numberOfTrellises)
            for i in range(self.numberOfTrellies):
                newCollectionOfBinaryTrellises.trellises[i] = self.trellises[i].minusTransform() if decisionVector == None else self.trellises[i].plusTransform(decisionVector)
            return newCollectionOfBinaryTrellises
        else:
            assert(self.length // 2 == self.numberOfTrellises)

            newBinaryMemorylessVectorDistribution = BinaryMemorylessVectorDistribution.BinaryMemorylessVectorDistribution(self.numberOfTrellises)
            for i in range(self.numberOfTrellies):
                tempTrellis = self.trellises[i].minusTransform() if decisionVector == None else self.trellises[i].plusTransform(decisionVector)
                marginzalizedProbs = tempTrellis.calcMarginalizedProbabilities(normalize=False)

                for x in range(2):
                    newBinaryMemorylessVectorDistribution.probs[i][x] = marginalizedProbs[x]

            return newBinaryMemorylessVectorDistribution

    def calcMarginalizedProbabilities(self):
        """Since a collection of trellises is eventually colapsed to a simple joint distribution after enough transforms, we should not get to this point
        """

        assert(False)

    def calcNormalizationVector(self):
        normalization = []

        for i in range(self.numberOfTrellises):
            trellisNormalization = self.trellises[i].calcNormalizationVector()

        return normalization

            
    def normalize(self, normalization):
        assert( len(normalization) == slef.numberOfTrellies)

        for i in range(self.numberOfTrellises):
            self.trellises[i].normalize(normalization[i])

