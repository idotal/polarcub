class VectorDistribution:

    # The following functions are pure virtual

    # Apply the minus transform, and return a new VectorDistribution of half the length
    def minusTransform(self):
        assert (False)

    # Apply the plus transform, according to the hard decisions made in the minus transform, and return a new VectorDistribution
    def plusTransform(self, uminusDecisions):
        assert (False)

    # Return the vector length. We clarify that in a deletion setting or an insertion/deletion/substitution setting, this is the number of inputs (not outputs)
    def __len__(self):
        assert (False)

    # Valid only for len(self) = 1. Returns a vector of probabilities summing to 1, whose length is the output alphabet size. That is, we return P_X(x).
    def calcMarginalizedProbabilities(self):
        assert (False)

    # Return a numpy 1-D vector of length len(self), which specifies the normalization of probabilities
    # That is, in a memoryless setting, we would divide the probabilities in entry i by the number in entry i of the normalizing vector, if we were doing SC decoding.
    # For SCL decoding, we would, for each path, divide the probabilities in entry i by the maximum of the numbers in entry i of the normalization vectors of each path.
    def calcNormalizationVector(self):
        assert (False)

    # normalize according to the above described vector
    def normalize(self, normalization):
        assert (False)
