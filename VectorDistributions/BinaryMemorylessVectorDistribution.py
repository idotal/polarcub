import numpy as np

import VectorDistribution


class BinaryMemorylessVectorDistribution(VectorDistribution.VectorDistribution):

    def __init__(self, length):
        assert (length > 0)
        self.probs = np.empty(
            (length, 2))  # so, probs[i][x] equals the probability of x transmitted or received at time i
        self.probs[:] = np.nan
        self.length = length

    def minusTransform(self):
        assert (self.length % 2 == 0)
        halfLength = self.length // 2

        newVector = BinaryMemorylessVectorDistribution(halfLength)

        for halfi in range(halfLength):
            i = 2 * halfi
            newVector.probs[halfi][0] = self.probs[i][0] * self.probs[i + 1][0] + self.probs[i][1] * self.probs[i + 1][
                1]
            newVector.probs[halfi][1] = self.probs[i][0] * self.probs[i + 1][1] + self.probs[i][1] * self.probs[i + 1][
                0]
            # print(newVector.probs[halfi][0], newVector.probs[halfi][1]) 

        return newVector

    def plusTransform(self, uminusDecisions):
        assert (self.length % 2 == 0)
        halfLength = self.length // 2

        newVector = BinaryMemorylessVectorDistribution(halfLength)

        for halfi in range(halfLength):
            i = 2 * halfi
            if uminusDecisions[halfi] == 0:
                newVector.probs[halfi][0] = self.probs[i][0] * self.probs[i + 1][0]
                newVector.probs[halfi][1] = self.probs[i][1] * self.probs[i + 1][1]
            else:
                newVector.probs[halfi][0] = self.probs[i][1] * self.probs[i + 1][0]
                newVector.probs[halfi][1] = self.probs[i][0] * self.probs[i + 1][1]
            # print(newVector.probs[halfi][0], newVector.probs[halfi][1]) 

        return newVector

    def __len__(self):
        return self.length

    def calcMarginalizedProbabilities(self):
        assert (len(self) == 1)

        marginalizedProbs = np.empty(2)
        marginalizedProbs[:] = np.nan

        s = 0.0
        for x in range(2):
            s += self.probs[0][x]

        if (s > 0.0):
            for x in range(2):
                marginalizedProbs[x] = self.probs[0][x] / s
        else:
            for x in range(2):
                marginalizedProbs[x] = 0.5

        return marginalizedProbs

    def calcNormalizationVector(self):
        normalization = np.zeros(self.length)

        for i in range(self.length):
            normalization[i] = np.maximum(self.probs[i][0], self.probs[i][1])

        return normalization

    def normalize(self, normalization):
        for i in range(self.length):
            t = normalization[i]
            assert (t >= 0)
            if t == 0:
                t = 1

            for x in range(2):
                self.probs[i][x] /= t
