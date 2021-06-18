import numpy as np

class BinaryMemorylessVectorDistribution(VectorDistribution):
    
    def __init__(self, length):
        assert( length > 0 )
        self.probs = np.empty(length, 2) # so, probs[i][x] equals the probability of x transmitted or received at time i
        self.probs[:] = NaN
        self.length = length

    def minusTransform(self):
        assert( length % 2 == 0 )
        halfLength = self.length // 2

        newVector = BinaryMemorylessVectorDistribution(halfLength)

        for halfi in range(halfLength):
            newVector.probs[halfi][0] = self.probs[2*i][0] * self.probs[2*i+1][0] + self.probs[2*i][1] * self.probs[2*i+1][1]
            newVector.probs[halfi][1] = self.probs[2*i][0] * self.probs[2*i+1][1] + self.probs[2*i][1] * self.probs[2*i+1][0]

        return newVector

    def plusTransform(self, uminusDecisions):
        assert( length % 2 == 0 )
        halfLength = self.length // 2

        newVector = BinaryMemorylessVectorDistribution(halfLength)

        for halfi in range(halfLength):
            if uminusDecisions[halfi] == 0:
                newVector.probs[halfi][0] = self.probs[2*i][0] * self.probs[2*i+1][0] 
                newVector.probs[halfi][1] = self.probs[2*i][1] * self.probs[2*i+1][1]
            else:
                newVector.probs[halfi][0] = self.probs[2*i][1] * self.probs[2*i+1][0] 
                newVector.probs[halfi][1] = self.probs[2*i][0] * self.probs[2*i+1][1]

        return newVector

    def __len__(self):
        return self.length

    def calcMarginalizedProbabilities(self):
        assert( len(self) == 1 )

        marginalizedProbs = np.empty(length, 2)
        self.probs[:] = NaN

        s = 0.0
        for x in range(2):
            s += self.probs[0][x]

        for x in range(2):
            marginalizedProbs[x] = self.probs[0][x]/s

        return marginalizedProbs

    def calcNormalizationVector(self):
        normalization = np.zeros(length)

        for i in range(self.length):
            normalization[i] = np.maximum( probs[i] )

        return normalization

    def normalize(self, normalization):
        for i in range(self.length):
            t = normalization[y]
            assert( t >= 0 )
            if t == 0:
                t = 1

            for x in range(2):
                probs[i][x] /= t

    
        

