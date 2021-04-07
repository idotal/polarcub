import math
import sys
import LinkedListHeap

class BinaryMemorylessDistribution:
    def __init__(self):
        self.probs = []  # probs[yindex][xindex]

    def __str__(self):
        s = "Binary memoryless channel with " + str(len(self.probs)) + " symbols and error probability " + str(self.errorProb()) + ". [p(y,x=0), p(y,x=1)]: "

        for probPair in self.probs:
            s += "[" + str(probPair[0]) + ", " + str(probPair[1]) + "], "

        s = s[0:-2] # remove the last ", "

        return s

    def append(self, item):
        self.probs.append(item)

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

        tempProbs = []
        for probPair in self.probs:
            tempProbs.append(probPair[0])
            tempProbs.append(probPair[1])

        tempProbs.sort()

        for tp in tempProbs:
            probSum += tp

        # for probPair in self.probs:
        #     probSum += sum(probPair)

        for probPair in self.probs:
            probPair[:] = [ prob / probSum for prob in probPair ]

    def removeZeroProbOutput(self):
        newProbs = []

        for probPair in self.probs:
            if sum(probPair) > 0.0:
                newProbs.append(probPair)

        self.probs = newProbs

    # sort according to p(x=0|y), in *ascending* order
    # that way, recalling the definition of Delta in upgradeLeftRightProbs, we have deltaLeftMinusRight >= 0
    def sortProbs(self):
    
        # for numerical stability, split into list with p(x=0|y) > 0.5 and p(x=0|y) <= 0.5
        zeroMoreProbable = []
        oneMoreProbable = []
        
        for probPair in self.probs:
            if probPair[0] / sum(probPair) > 0.5:
                zeroMoreProbable.append(probPair)
            else:
                oneMoreProbable.append(probPair)
        
        # sort probs according to p(x=0|y) 
        zeroMoreProbable.sort(key = lambda probPair: probPair[1] / sum(probPair) )
        oneMoreProbable.sort(key = lambda probPair: -probPair[0] / sum(probPair) )
    
        self.probs = zeroMoreProbable + oneMoreProbable 

    def mergeEquivalentSymbols(self):
        self.removeZeroProbOutput()

        self.sortProbs()

        # insert the first output letter
        newProbs = [self.probs[0]]

        # loop over all other output letters, and check if we need to append as a new letter, or add to the last letter
        for probPair in self.probs[1:]:
            prevProbPair = newProbs[-1]

            isclose = True
            for b in range(2):
                if not math.isclose(probPair[b] / sum(probPair), prevProbPair[b] / sum(prevProbPair)):  
                    isclose = False

            if not isclose:
                newProbs.append( [probPair[0], probPair[1]] )
            else:
                for b in range(2):
                    newProbs[-1][b] += probPair[b]

        self.probs = newProbs

        self.normalize() # for good measure

    # polar transforms
    def minusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.append( [y1[0] * y2[0] + y1[1] * y2[1], y1[0] * y2[1] + y1[1] * y2[0]])

        newDistribution.mergeEquivalentSymbols()

        return newDistribution

    def plusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.append( [y1[0] * y2[0], y1[1] * y2[1]])
                newDistribution.append( [y1[1] * y2[0], y1[0] * y2[1]])

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

    def upgrade(self, L):
        dataList = []
        keyList = []

        for probPair in self.probs:
            datum = [probPair[0], probPair[1]] # make a copy
            dataList.append(datum)
                               
        for i in range(len(dataList)):
            key = _calcKey_upgrade(_listIndexingHelper(dataList,i-1), _listIndexingHelper(dataList,i), _listIndexingHelper(dataList,i+1))
            keyList.append(key)

        llh = LinkedListHeap.LinkedListHeap(keyList, dataList)

        while llh.numberOfElements() > L:
            # print("llh.numberOfElements = ", llh.numberOfElements())
            # print(llh)

            topOfHeap = llh.extractHeapMin()
            leftElement = topOfHeap.leftElementInList
            rightElement = topOfHeap.rightElementInList
            
            dataMergeLeft, dataMergeRight = upgradedLeftRightProbs(leftElement.data, topOfHeap.data, rightElement.data)

            # move probs to left and right elements
            for b in range(2):
                leftElement.data[b] += dataMergeLeft[b]
                rightElement.data[b] += dataMergeRight[b]

            # recalculate key of left element
            leftLeftElement = leftElement.leftElementInList

            if leftLeftElement != None:
                key = _calcKey_upgrade(leftLeftElement.data, leftElement.data, rightElement.data)
                llh.updateKey(leftElement, key)

            # recalculate key of right element
            rightRightElement = rightElement.rightElementInList

            if rightRightElement != None:
                key = _calcKey_upgrade(leftElement.data, rightElement.data, rightRightElement.data)
                llh.updateKey(rightElement, key)

        newDistribution = BinaryMemorylessDistribution()
        newDistribution.probs = llh.returnData()

        return newDistribution

# functions for degrade/upgrade/merge
def eta(p):
    assert 0 <= p <= 1

    if p == 0 or p == 1:
        return 0
    else:
        return -p * math.log2(p)

def hxgiveny(data):
    py = data[0] + data[1]
    return py * ( eta(data[0]/py) + eta(data[1]/py) )

# def myisclose(a, b):
#     return True if abs(a-b) < 100.0 * sys.float_info.epsilon else False

# useful channels

def makeBSC(p):
    bsc = BinaryMemorylessDistribution()
    bsc.append( [0.5 * p, 0.5 * (1.0-p)] )
    bsc.append( [0.5 * (1.0-p), 0.5 * p] )
    
    return bsc

def makeBEC(p):
    bec = BinaryMemorylessDistribution()
    bec.append( [0.5 * (1.0-p), 0] )
    bec.append( [0, 0.5 * (1.0-p)] )
    bec.append( [0.5 * p, 0.5 * p] )
    
    return bec

# private functions for degrade/upgrade
def _calcKey_degrade(dataLeft, dataCenter): # how much would it cost to merge dataLeft and dataCenter
    if dataLeft == None:
        return float("inf")

    assert len(dataLeft) == len(dataCenter) == 2

    dataMerge = [ dataLeft[0] + dataCenter[0], dataLeft[1] + dataCenter[1] ]

    return hxgiveny(dataMerge) - hxgiveny(dataLeft) - hxgiveny(dataCenter)

def _calcKey_upgrade(dataLeft, dataCenter, dataRight): # how much would it cost to split dataCenter  into dataLeft and dataRight
    if dataLeft == None or dataRight == None:
        return float("inf")

    assert len(dataLeft) == len(dataCenter) == len(dataRight) == 2

    dataMergeLeft, dataMergeRight = upgradedLeftRightProbs(dataLeft, dataCenter, dataRight)

    return hxgiveny(dataCenter) - hxgiveny(dataMergeLeft) - hxgiveny(dataMergeRight)  

def _listIndexingHelper(l, i):
    return l[i] if (0 <= i < len(l)) else None

def upgradedLeftRightProbs(dataLeft, dataCenter, dataRight):
    # pi = p(y,x=0) + p(y,x=1)
    # Delta = (p(y,x=0) - p(y,x=1))/pi
    # piCenter is split into thetaLeft * piCenter and thetaRight * piCenter

    piLeft = sum(dataLeft)
    piCenter = sum(dataCenter)
    piRight = sum(dataRight)

    normalizedLeft = []
    normalizedCenter = []
    normalizedRight = []

    for b in range(2):
        normalizedLeft.append(dataLeft[b]/piLeft)
        normalizedCenter.append(dataCenter[b]/piCenter)
        normalizedRight.append(dataRight[b]/piRight)

    # calculate Delta_left - Delta_right
    # split into cases, for numerical stability
    if normalizedLeft[0] < 0.5 and normalizedRight[0] < 0.5:
        deltaLeftMinusRight = 2.0 * (normalizedLeft[0] - normalizedRight[0])
    elif normalizedLeft[1] < 0.5 and normalizedRight[1] < 0.5:
        deltaLeftMinusRight = 2.0 * (normalizedRight[1] - normalizedLeft[1])
    else:
        deltaLeft = normalizedLeft[0] - normalizedLeft[1]
        deltaRight = normalizedRight[0] - normalizedRight[1]
        deltaLeftMinusRight = deltaLeft - deltaRight

    assert(deltaLeftMinusRight > 0.0) # should have been resolved by merging equivalent symbols

    # calculate thetaLeft and thetaRight
    # split into cases, for numerical stability
    if normalizedLeft[0] < 0.5 and normalizedRight[0] < 0.5:
        if normalizedLeft[0] - normalizedCenter[0] < normalizedCenter[0] - normalizedRight[0]: # center is closer to left
            # so, the smaller theta is thetaRight
            deltaLeftMinusCenter = 2.0 * ( normalizedLeft[0] - normalizedCenter[0] )
            thetaRight = deltaLeftMinusCenter/deltaLeftMinusRight
            thetaLeft = 1.0 - thetaRight
        else: # the smaller theta is thetaLeft
            deltaCenterMinusRight = 2.0 * ( normalizedCenter[0] - normalizedRight[0] )
            thetaLeft = deltaCenterMinusRight/deltaLeftMinusRight
            thetaRight = 1.0 - thetaLeft
    elif normalizedLeft[1] < 0.5 and normalizedRight[1] < 0.5:
        if normalizedCenter[1] - normalizedLeft[1] <  normalizedRight[1] - normalizedCenter[1]: # center is closer to left
            # so, the smaller theta is thetaRight
            deltaLeftMinusCenter = 2.0 * ( normalizedCenter[1] - normalizedLeft[1] )
            thetaRight = deltaLeftMinusCenter/deltaLeftMinusRight
            thetaLeft = 1.0 - thetaRight
        else: # the smaller theta is thetaLeft
            deltaCenterMinusRight = 2.0 * ( normalizedRight[1] - normalizedCenter[1] )
            thetaLeft = deltaCenterMinusRight/deltaLeftMinusRight
            thetaRight = 1.0 - thetaLeft
    else:
        deltaLeft = normalizedLeft[0] - normalizedLeft[1]
        deltaCenter = normalizedCenter[0] - normalizedCenter[1]
        deltaLeftMinusCenter = deltaLeft - deltaCenter
        thetaRight = deltaLeftMinusCenter/deltaLeftMinusRight
        thetaLeft = 1.0 - thetaRight

    assert(0.0 < thetaLeft < 1.0 and 0.0 < thetaRight < 1.0)

    dataMergeLeft = []
    dataMergeRight = []

    for b in range(2):
        dataMergeLeft.append(thetaLeft * piCenter * normalizedLeft[b])
        dataMergeRight.append(thetaRight * piCenter * normalizedRight[b])

    return dataMergeLeft, dataMergeRight 

