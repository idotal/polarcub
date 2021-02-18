import math
import sys

class BinaryMemorylessDistribution:
    def __init__(self):
        self.probs = []  # probs[yindex][xindex]

    # polar toolbox
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

    def conditionalEntropy(self):
        entropySum = 0.0

        for probPair in self.probs:
            entropySum += ( eta(probPair[0]) + eta(probPair[1]) - eta(probPair[0] + probPair[1]) )

        return entropySum

    # housekeeping

    def normalize(self):
        probSum = 0.0

        for probPair in self.probs:
            probSum += sum(probPair)

        for probPair in self.probs:
            probPair[:] = [ prob / probSum for prob in probPair ]

    def removeZeroProbOutput(self):
        newProbs = []

        for probPair in self.probs:
            if sum(probPair) > 0.0:
                newProbs.append(probPair)

        self.probs = newProbs

    def mergeEquivalentSymbols(self):
        self.removeZeroProbOutput()

        # sort probs according to p(x=0|y)
        self.probs.sort(key = lambda probPair: probPair[0] / (probPair[0] + probPair[1]) )

        # insert the first output letter
        newProbs = [self.probs[0]]

        # loop over all other output letters, and check if we need to append as a new letter, or add to the last letter
        for probPair in self.probs[1:]:
            prevProbPair = newProbs[-1]
            if not math.isclose(probPair[0] / (probPair[0] + probPair[1]), prevProbPair[0] / (prevProbPair[0] + prevProbPair[1])):
                newProbs.append( probPair )
            else:
                newProbs[-1][0] += probPair[0]
                newProbs[-1][1] += probPair[1]

        self.probs = newProbs

        self.normalize() # for good measure

    # polar transforms
    def minusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.probs.append( [y1[0] * y2[0] + y1[1] * y2[1], y1[0] * y2[1] + y1[1] * y2[0]])

        newDistribution.mergeEquivalentSymbols()

        return newDistribution

    def plusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.probs.append( [y1[0] * y2[0], y1[1] * y2[1]])
                newDistribution.probs.append( [y1[1] * y2[0], y1[0] * y2[1]])

        newDistribution.mergeEquivalentSymbols()

        return newDistribution

# useful channels

def makeBSC(p):
    bsc = BinaryMemorylessDistribution()
    bsc.probs.append( [0.5 * p, 0.5 * (1.0-p)] )
    bsc.probs.append( [0.5 * (1.0-p), 0.5 * p] )
    
    return bsc

# useful functions

def eta(p):
    assert 0 <= p <= 1

    if p == 0 or p == 1:
        return p
    else:
        return -p * math.log2(p)

