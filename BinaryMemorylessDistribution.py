import math
import sys
import LinkedListHeap

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

    # public functions for degrading/upgrading
    def degrade(self, L):
        dataList = []
        keyList = []

        for probPair in self.probs:
            datum = [probPair[0], probPair[1]] # make a copy
            dataList.append(datum)
                               
        for i in range(len(dataList)):
            key = _calcKey_degrade(_listIndexingHelper(dataList,i-1), _listIndexingHelper(dataList,i))
            keyList.append(key)

        llh = LinkedListHeap.LinkedListHeap(keyList, dataList)

        while llh.numberOfElements() > L:
            topOfHeap = llh.extractHeapMin()
            leftElement = topOfHeap.leftElementInList
            rightElement = topOfHeap.rightElementInList
            
            # move probs to left element
            for b in range(2):
                leftElement.data[b] += topOfHeap.data[b]

            # recalculate key of left element
            leftLeftElement = leftElement.leftElementInList

            if leftLeftElement != None:
                key = _calcKey_degrade(leftLeftElement.data, leftElement.data)
                llh.updateKey(leftElement, key)

            # recalculate key of right element
            if rightElement != None:
                key = _calcKey_degrade(leftElement.data, rightElement.data)
                llh.updateKey(rightElement, key)

        newDistribution = BinaryMemorylessDistribution()
        newDistribution.probs = llh.returnData()

        return newDistribution

# functions for degrade/upgrade/merge
def eta(p):
    assert 0 <= p <= 1

    if p == 0 or p == 1:
        return p
    else:
        return -p * math.log2(p)

def hxgiveny(data):
    py = data[0] + data[1]
    return py * ( eta(data[0]/py) + eta(data[1]/py) )

# useful channels

def makeBSC(p):
    bsc = BinaryMemorylessDistribution()
    bsc.probs.append( [0.5 * p, 0.5 * (1.0-p)] )
    bsc.probs.append( [0.5 * (1.0-p), 0.5 * p] )
    
    return bsc

def makeBEC(p):
    bec = BinaryMemorylessDistribution()
    bec.probs.append( [0.5 * (1.0-p), 0] )
    bec.probs.append( [0, 0.5 * (1.0-p)] )
    bec.probs.append( [0.5 * p, 0.5 * p] )
    
    return bec

# private functions for degrade
def _calcKey_degrade(dataLeft, dataCenter): # how much would it cost to merge dataLeft and dataCenter
    if dataLeft == None:
        return float("inf")

    assert len(dataLeft) == len(dataCenter) == 2

    dataMerge = [ dataLeft[0] + dataCenter[0], dataLeft[1] + dataCenter[1] ]

    return hxgiveny(dataMerge) - hxgiveny(dataLeft) - hxgiveny(dataCenter)

def _listIndexingHelper(l, i):
    return l[i] if (0 <= i < len(l)) else None

