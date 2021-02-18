import math

class BinaryMemorylessDistribution:
    def __init__(self):
        self.probs = []  # probs[yindex][xindex]

    def errorProb(self):
        errorProbSum = 0.0

        for probPair in self.probs:
            errorProbSum += min(probPair)

        return errorProbSum

    def bhattacharyya(self):
        bhattacharyyaSum = 0.0

        for probPair in self.probs:
            bhattacharyyaSum += math.sqrt(probPair[0]*probPair[1])

        return 2.0 * bhattacharyyaSum

    def totalVariationDistance(self):
        tvSum = 0.0

        for probPair in self.probs:
            tvSum += abs(probPair[0] - probPair[1])

        return tvSum

    def normalize(self):
        probSum = 0.0

        for probPair in self.probs:
            probSum += sum(probPair)

        for probPair in self.probs:
            probPair[:] = [ prob / probSum for prob in probPair ]


