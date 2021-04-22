import math
import sys
import LinkedListHeap
import itertools
from cython.cython_BinaryMemorylessDistribution import eta as fast_eta
from cython.cython_BinaryMemorylessDistribution import hxgiveny as fast_hxgiveny

use_fast = True
# use_fast = False

class BinaryMemorylessDistribution:
    def __init__(self):
        self.probs = []  # probs[yindex][xindex]
        self.auxiliary = None # for keeping track of things needed by degrading or upgrading procedures

    def __str__(self):
        s = "Binary memoryless channel with " + str(len(self.probs)) + " symbols and error probability " + str(self.errorProb()) + ". [p(y,x=0), p(y,x=1)]: "

        for probPair in self.probs:
            s += "[" + str(probPair[0]) + ", " + str(probPair[1]) + "], "

        s = s[0:-2] # remove the last ", "

        if self.auxiliary != None:
            s += "\n the auxiliary data is " + str(self.auxiliary)

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
        """remove output symbols with probability zero.

        If the auxiliary list is not None, we simply drop entries corresponding to removed letters.
        """
        newProbs = []
        newAuxiliary = []

        tempAux = [] if self.auxiliary == None else self.auxiliary
        for probPair, auxDatum in itertools.zip_longest(self.probs, tempAux ):
            if sum(probPair) > 0.0:
                newProbs.append(probPair)
                newAuxiliary.append(auxDatum)

        self.probs = newProbs
        if self.auxiliary != None:
            self.auxiliary = newAuxiliary

    def sortProbs(self):
        """Sort according to p(x=0|y), in *ascending* order.

        This way, recalling the definition of Delta in upgradeLeftRightProbs, we have deltaLeftMinusRight >= 0.
        """
    
        # for numerical stability, split into list with p(x=0|y) > 0.5 and p(x=0|y) <= 0.5
        zeroMoreProbable = []
        oneMoreProbable = []
        
        # for probPair in self.probs:
        #     if probPair[0] / sum(probPair) > 0.5:
        #         zeroMoreProbable.append(probPair)
        #     else:
        #         oneMoreProbable.append(probPair)
        #
        # # sort probs according to p(x=0|y) 
        # zeroMoreProbable.sort(key = lambda probPair: probPair[1] / sum(probPair) )
        # oneMoreProbable.sort(key = lambda probPair: -probPair[0] / sum(probPair) )
        #
        # self.probs = zeroMoreProbable + oneMoreProbable 

        tempAux = [] if self.auxiliary == None else self.auxiliary

        for probPair, auxDatum in itertools.zip_longest(self.probs, tempAux):
            if probPair[0] / sum(probPair) > 0.5:
                zeroMoreProbable.append((probPair, auxDatum))
            else:
                oneMoreProbable.append((probPair, auxDatum))

        def zeroMoreProbableKey(tup):
            probPair = tup[0]
            return probPair[1] / sum(probPair)

        def oneMoreProbableKey(tup):
            probPair = tup[0]
            return -probPair[0] / sum(probPair)

        zeroMoreProbable.sort(key = zeroMoreProbableKey)
        oneMoreProbable.sort(key = oneMoreProbableKey)

        # self.probs = zeroMoreProbable + oneMoreProbable

        self.probs = []
        if self.auxiliary != None:
            self.auxiliary = []

        for probPair, auxDatum in zeroMoreProbable:
            self.probs.append(probPair)
            if self.auxiliary != None:
                self.auxiliary.append(auxDatum)

        for probPair, auxDatum in oneMoreProbable:
            self.probs.append(probPair)
            if self.auxiliary != None:
                self.auxiliary.append(auxDatum)

    def mergeEquivalentSymbols(self):
        """Merge output symbols with very close LLRs into a single output letter.

        If the auxiliary list is not None, it is assumed to contain sets (a set for each output letter). Merging results in the union of these sets.
        The probabilities are normalized at the end, for good measure.
        """
        self.removeZeroProbOutput()

        self.sortProbs()

        # insert the first output letter
        newProbs = [self.probs[0]]
        if self.auxiliary != None:
            newAuxiliary = [self.auxiliary[0]]

        # loop over all other output letters, and check if we need to append as a new letter, or add to the last letter
        tempAux = [] if self.auxiliary == None else self.auxiliary
        for probPair, auxDatum in itertools.zip_longest(self.probs[1:], tempAux[1:]):
            prevProbPair = newProbs[-1]

            isclose = True
            for b in range(2):
                if not math.isclose(probPair[b] / sum(probPair), prevProbPair[b] / sum(prevProbPair)):  
                    isclose = False

            if not isclose:
                newProbs.append( [probPair[0], probPair[1]] )
                if self.auxiliary != None:
                    newAuxiliary.append(auxDatum)
            else:
                for b in range(2):
                    newProbs[-1][b] += probPair[b]

                    if self.auxiliary != None:
                        newAuxiliary[-1] |= auxDatum 

        self.probs = newProbs

        if self.auxiliary != None:
            self.auxiliary = newAuxiliary

        self.normalize() # for good measure

    # polar transforms
    def minusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.append( [y1[0] * y2[0] + y1[1] * y2[1], y1[0] * y2[1] + y1[1] * y2[0]])

        # newDistribution.mergeEquivalentSymbols()

        return newDistribution

    def plusTransform(self):
        
        newDistribution = BinaryMemorylessDistribution()

        for y1 in self.probs:
            for y2 in self.probs:
                newDistribution.append( [y1[0] * y2[0], y1[1] * y2[1]])
                newDistribution.append( [y1[1] * y2[0], y1[0] * y2[1]])

        # newDistribution.mergeEquivalentSymbols()

        return newDistribution

    # public functions for degrading/upgrading
    def degrade(self, L):
        """Degrade into a new channel, containing at most L output letters.

        If auxiliary is not None, it is assumed to contain a list of sets, a set for each output symbol. The new channel contains a corresponding auxiliary list, also containing sets, one for each output symbol. Each such set is the union of the sets corresponding to the output symbols in the original channel which have been merged together into the output symbol of the new channel.
        """
        dataList = []
        keyList = []

        # Merge, and also sort according to LLR
        self.mergeEquivalentSymbols()

        # note that linkedListDatum = [[probPair[0], probPair[1]], auxDatum]

        tempAux = [] if self.auxiliary == None else self.auxiliary
        for probPair, auxDatum in itertools.zip_longest(self.probs, tempAux):
            linkedListDatum = [[probPair[0], probPair[1]], auxDatum] # make a copy
            dataList.append(linkedListDatum)
                               
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
                leftElement.data[0][b] += topOfHeap.data[0][b]

            # update aux, if needed
            if self.auxiliary != None:
                leftElement.data[1] |= topOfHeap.data[1]

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

        if self.auxiliary != None:
            newDistribution.auxiliary = []

        for probPair, auxDatum in llh.returnData():
            newDistribution.probs.append(probPair) 
            if self.auxiliary != None:
                newDistribution.auxiliary.append(auxDatum) 

        return newDistribution

    def upgrade(self, L):
        """Upgrade into a new channel, containing at most L output letters.

        We first remove output letters with zero probability, sort according to ascending order of p(x=0|y), merge symbols with the same LLR, and then upgrade. That is, if the original order to the symbols was
        y0,y1,y2,...
        Then the new symbols z are also sorted according to p(x=0|z). Denote these as
        z0,z1,z2
        Since each zi has the same LLR as some yj, we can merge the above two figures into, say
        z0=y0,y1,y2,z1={y3,y4},y5,y6,z3=y8
        where y3 and y4 have been merged together by the preliminary merge operation, since they had the same LLR. Note that in the above "=" means "same LLR".

        
        If auxiliary is not None, it is assumed to contain a list of sets, a set for each output symbol. The new channel contains a corresponding auxiliary list, each element is a list of three sets, [sybmolsLeft, symbolsCenter, symbolsRight]. The set symbolsCenter contains the indices of the symbols in the original channel having the same LLR as the corresponding symbol in the new channel. That is, for the above example, symbolsCenter(z1) = set(y3) \cup set(y4). Continuing the above example, symbolsLeft(z1) = set(y1) \cup set(y2) and sybmolsRight(z1) = set(y5) \cup set(y6). Note that the symbolsCenter set is always non-empty, the sets symbolsLeft and symbolsRight might be empty, and are always empty for the first and last output letters, respectively.
        """
        dataList = []
        keyList = []

        # Merge, and also sort according to LLR
        self.mergeEquivalentSymbols()

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
if use_fast == True:
    eta = fast_eta
    hxgiveny = fast_hxgiveny
else:
    def eta(p):
        assert 0 <= p <= 1
    
        if p == 0:
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

    probLeft = dataLeft[0]
    probCenter = dataCenter[0]

    assert len(probLeft) == len(probCenter) == 2

    probMerge = [ probLeft[0] + probCenter[0], probLeft[1] + probCenter[1] ]

    return hxgiveny(probMerge) - hxgiveny(probLeft) - hxgiveny(probCenter)

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

