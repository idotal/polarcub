class BinaryMemorylessDistribution:
    def __init__(self):
        self.probs = []  # probs[yindex][xindex]

    def errorProb(self):
        errorProbSum = 0.0

        for probPair in self.probs:
            errorProbSum += min(probPair)

        return errorProbSum


