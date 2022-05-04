import numpy as np
import VectorDistribution

class QaryMemorylessVectorDistribution(VectorDistribution.VectorDistribution):

    def __init__(self, q, length):
        assert (q > 1)
        self.q = q
        assert ( length > 0 )
        self.probs = np.empty((length, q), dtype=np.float)  # so, probs[i][x] equals the probability of x transmitted or received at time i
        self.probs[:] = np.nan
        self.length = length

    def minusTransform(self):
        assert( self.length % 2 == 0 )
        halfLength = self.length // 2

        newVector = QaryMemorylessVectorDistribution(self.q, halfLength)
        newVector.probs[:] = 0.0

        for halfi in range(halfLength):
            i = 2 * halfi
            for x1 in range(self.q):
                for x2 in range(self.q):
                    u1 = (x1 + x2) % self.q
                    newVector.probs[halfi][u1] += self.probs[i][x1] * self.probs[i + 1][x2]
        return newVector

    def plusTransform(self, uminusDecisions):
        assert( self.length % 2 == 0 )
        halfLength = self.length // 2

        newVector = QaryMemorylessVectorDistribution(self.q, halfLength)
        newVector.probs[:] = 0.0

        for halfi in range(halfLength):
            i = 2 * halfi

            u1 = uminusDecisions[halfi]
            for u2 in range(self.q):
                x1 = (u1 - u2 + self.q) % self.q
                x2 = u2
                newVector.probs[halfi][u2] += self.probs[i][x1] * self.probs[i + 1][x2]
            # print(newVector.probs[halfi][0], newVector.probs[halfi][1])

        return newVector

    def __len__(self):
        return self.length

    def calcMarginalizedProbabilities(self):
        assert( len(self) == 1 )

        marginalizedProbs = np.empty(self.q)
        marginalizedProbs[:] = np.nan

        s = 0.0
        for x in range(self.q):
            s += self.probs[0][x]

        if (s > 0.0):
            for x in range(self.q):
                marginalizedProbs[x] = self.probs[0][x] / s
        else:
            for x in range(self.q):
                marginalizedProbs[x] = 1.0 / self.q

        return marginalizedProbs

    def calcNormalizationVector(self):
        normalization = np.zeros(self.length)

        for i in range(self.length):
            normalization[i] = self.probs[i].max(axis=0)

        return normalization

    def normalize(self, normalization=None):
        if normalization is None:
            normalization = self.calcNormalizationVector()
        for i in range(self.length):
            t = normalization[i]
            assert( t >= 0 )
            if t == 0:
                t = 1

            for x in range(self.q):
                self.probs[i][x] /= t