import BinaryMemorylessDistribution
from BinaryMemorylessDistribution import eta
from math import floor
import numpy as np

# constants
lcrLeft = 0
lcrCenter = 1
lcrRight = 2

class QaryMemorylessDistribution:
    def __init__(self, q):
        self.probs = []  # probs[yindex][xindex]
        self.q = q

    def __str__(self):
        s = "Qry memoryless channel with q = " + str(self.q) + " and " + str(len(self.probs)) + " output symbols. The error probability is " + str(self.errorProb()) + ". The conditional entropy is " + str(self.conditionalEntropy()) +  ". [p(y,x=0), p(y,x=1), ..., p(y,x=q-1)]: "

        for probTuple in self.probs:
            s += "[" 
            for x in range(self.q):
                s += str(probTuple[x]) + ", "
            s = s[0:-2] # remove the last ", "
            s += "], "
        s = s[0:-2] # remove the last ", "
        return s

    def append(self, item):
        self.probs.append(item)

    # polar toolbox
    def errorProb(self):
        totalProbSum = 0.0
        noErrorProbSum = 0.0

        for probTuple in self.probs:
            totalProbSum += sum(probTuple)
            noErrorProbSum += max(probTuple)

        return totalProbSum - noErrorProbSum

    def conditionalEntropy(self):
        entropySum = 0.0

        for probTuple in self.probs:
            for p in probTuple:
                entropySum += eta(p)

            entropySum -= eta(sum(probTuple))

        return entropySum

    def oneHotBinaryMemorylessDistributions(self):
        """Return q-1 one-hot binary channels.
        
        The probs list of a channel contains the corresponding probabilities.
        The auxiliary list of a channel contains sets, each set containing a single element: the index of the corresponding output letter in the q-ary channel.
        Note that a binary channel may contain output letters with probability zero.
        """

        binaryMemorylessDistributions = []
        q = self.q

        # create q-1 "empty" channels
        for i in range(q-1):
            binaryMemorylessDistributions.append(BinaryMemorylessDistribution.BinaryMemorylessDistribution())
            binaryMemorylessDistributions[i].auxiliary = []

        # calculate the marginals p(X = j)
        marginals = self.calcXMarginals()

        # calculate p(X > j)
        probGreaterThan = []
        for x in range(q):
            probGreaterThan.append(0)

        for x in range(q-2,-1,-1):
            probGreaterThan[x] = probGreaterThan[x+1] + marginals[x+1]

        # x \equiv (x_0,x_1,...x_{q-2}), where x_j = 1 iff x = j (the vector is all zero iff x = q-1)

        # For each y, first calculate P(X_j = b, Y=y, X^{j-1} = 0) for b in {0,1}, and then use the above calculated marginals to get 
        # P(X_j = b, Y=y | X^{j-1} = 0) = P(X_j = b, Y=y, X^{j-1} = 0) / P(X^{j-1} = 0), where P(X^{j-1} = 0) = P( X > j-1 )

        for yindex, probTuple in enumerate(self.probs):

            prev_pbzero = prev_pbone = None # to catch errors

            for j in range(q-2,-1,-1):
                # Calculate P(X_j = b, Y=y, X^{j-1} = 0), for b = 1. That is, P(X = j, Y=y)
                pbone = probTuple[j]
                # Calculate P(X_j = b, Y=y, X^{j-1} = 0), for b = 0. That is, P(X > j, Y=y)
                if j == q-2:
                    pbzero = probTuple[j+1]
                else:
                    pbzero =  prev_pbone + prev_pbzero # greater than j either means j+1 or greater than j+1

                prev_pbzero, prev_pbone = pbzero, pbone

                if j == 0:
                    probPair = [pbzero, pbone]
                else:
                    probPair = [pbzero/probGreaterThan[j-1], pbone/probGreaterThan[j-1]]

                binaryMemorylessDistributions[j].append(probPair)
                binaryMemorylessDistributions[j].auxiliary.append({yindex})

        return binaryMemorylessDistributions

    def calcXMarginals(self):
        """calculate the marginals p(X = j)"""
        marginals = []

        for x in range(self.q):
            tempSum = 0.0
            for probTuple in self.probs:
                tempSum += probTuple[x]
            marginals.append(tempSum)

        return marginals

    def calcYMarginal(self, y):
        """calculate the marginal p(Y = y)"""
        return sum(self.probs[y])

    # polar transforms
    def minusTransform(self):
        newDistribution = QaryMemorylessDistribution(self.q)

        for y1 in self.probs:
            for y2 in self.probs:
                tempProbs = [0 for x in range(self.q)]

                for x1 in range(self.q):
                    for x2 in range(self.q):
                        u1 = (x1 + x2) % self.q
                        tempProbs[u1] += y1[x1] * y2[x2]

                newDistribution.append(tempProbs)

        return newDistribution

    def plusTransform(self):
        newDistribution = QaryMemorylessDistribution(self.q)

        for y1 in self.probs:
            for y2 in self.probs:

                for u1 in range(self.q):
                    tempProbs = [0 for x in range(self.q)]
                    for u2 in range(self.q):
                        x1 = (u1 - u2 + self.q) % self.q
                        x2 = u2
                        tempProbs[u2] += y1[x1] * y2[x2]
                    newDistribution.append(tempProbs)

        return newDistribution

    # public functions for degrading/upgrading
    def degrade(self, L):
        return self.degrade_dynamic(L)

    def degrade_dynamic(self, L):
        oneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()

        M = floor( L ** (1.0/(self.q-1)) )
        
        degradedOneHotBinaryMemorylessDistributions = []
        for x in range(self.q-1):
            degradedOneHotBinaryMemorylessDistributions.append( oneHotBinaryMemorylessDistributions[x].degrade(M) )

        newDistribution = QaryMemorylessDistribution(self.q)

        newOutputAlphabetSize = self.calcNewOutputAlphabetSize(degradedOneHotBinaryMemorylessDistributions)

        allzeroindexvector = [0 for x in range(self.q -1)]

        yoldMappedTo = []
        for yold in range(len(self.probs)):
            # by default, if yold is not mapped to any symbol in a one-hot channel
            # (because of zero probability), it is mapped to symbol 0 of that channel
            yoldMappedTo.append(allzeroindexvector.copy()) 


        for x in range(self.q - 1):
            for yOneHotIndex, auxDatum in enumerate(degradedOneHotBinaryMemorylessDistributions[x].auxiliary):
                for yold in auxDatum:
                    yoldMappedTo[yold][x] = yOneHotIndex

        allzeroprobvector = [0.0 for x in range(self.q)]

        for ynew in range(newOutputAlphabetSize):
            newDistribution.probs.append(allzeroprobvector.copy())

        conversionToYNewMultipliers = self.calcConversionToYNewMultipliers(degradedOneHotBinaryMemorylessDistributions)

        for yold, yoldprobs in enumerate(self.probs):
            ynew = self.yoldToNew_degrade(yold, yoldMappedTo, conversionToYNewMultipliers)

            for x in range(self.q):
                newDistribution.probs[ynew][x] += yoldprobs[x]

        newDistribution.removeZeroProbOutput()

        return newDistribution

    def degrade_static(self, L):
        M = floor( L ** (1.0/(self.q-1)) )

        dims = [M for i in range(self.q-1)]

        cellsArray = np.array([set() for _ in range(M ** (self.q-1))]).reshape(dims)

        # first, put each output letter into the corresponding cell

        for yold, prob in enumerate(self.probs):
            probsum = sum(prob)
            cell = []
            for x in range(self.q - 1):
                cell.append(floor(M * prob[x]/probsum)) # TODO: stub

            print(yold, prob, cell)
            cellsArray[tuple(cell)] |= {yold}
            # print(cellsArray)

        print(cellsArray)


        # then, merge all the letters in a cell into a single letter

    def upgrade(self, L):
        return self.upgrade_dynamic(L)

    def upgrade_dynamic(self, L):
        oneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()
        M = floor( L ** (1.0/(self.q-1)) )

        upgradedOneHotBinaryMemorylessDistributions = []
        for x in range(self.q-1):
            upgradedOneHotBinaryMemorylessDistributions.append( oneHotBinaryMemorylessDistributions[x].upgrade(M) )

        newDistribution = QaryMemorylessDistribution(self.q)

        newOutputAlphabetSize = self.calcNewOutputAlphabetSize(upgradedOneHotBinaryMemorylessDistributions)

        # yoldMappedTo[yold][x][lcr], where lcr=0,1,2 corresponds to left, center, right
        yoldMappedTo = []
        for yold in range(len(self.probs)):
            yoldMappedTo.append([[None for lcr in range(3)] for x in range(self.q-1)]) 

        for x in range(self.q - 1):
            for yOneHotIndex, auxDatum in enumerate(upgradedOneHotBinaryMemorylessDistributions[x].auxiliary):
                for lcr in range(3): #left-center-right
                    for yold in auxDatum[lcr]:
                        assert( yoldMappedTo[yold][x][lcr] == None )
                        yoldMappedTo[yold][x][lcr] = yOneHotIndex

        allzeroprobvector = [0.0 for x in range(self.q)]

        for ynew in range(newOutputAlphabetSize):
            newDistribution.probs.append(allzeroprobvector.copy())

        conversionToYNewMultipliers = self.calcConversionToYNewMultipliers(upgradedOneHotBinaryMemorylessDistributions)

        # calculate again, since upgrading removes output symbols with probability zero, and then sorts
        originalOneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()

        for yold in range(len(self.probs)):
            self.addToNewDistribution_upgrade(yold, yoldMappedTo, newDistribution, conversionToYNewMultipliers, upgradedOneHotBinaryMemorylessDistributions, originalOneHotBinaryMemorylessDistributions)

        newDistribution.removeZeroProbOutput()

        return newDistribution

    def addToNewDistribution_upgrade(self, yold, yoldMappedTo, newDistribution, conversionToYNewMultipliers, upgradedOneHotBinaryMemorylessDistributions, originalOneHotBinaryMemorylessDistributions):
        yoldMarginal = self.calcYMarginal(yold)
        
        lcrvec = self.initializeLCRVector(yold, yoldMappedTo)

        while True:
            ynew = self.yoldToNew_upgrade(yold, yoldMappedTo, lcrvec, conversionToYNewMultipliers)
            # we've picked yold, and now ynew (through lcrvec). 
            x_ynew_yold_probs = self.calc_probs_of_x_ynew_given_yold(yold, yoldMappedTo, lcrvec, upgradedOneHotBinaryMemorylessDistributions, originalOneHotBinaryMemorylessDistributions)


            marginalFromNowProb = [0.0 for i in range(self.q+1)] # marginalFromNowProb[i] = P(ynew[i],ynew[i+1],...,ynew[q-2]|yold)

            marginalFromNowProb[self.q-1] = 1.0
            for i in range(self.q-2, -1, -1): 
                marginalFromNowProb[i] = marginalFromNowProb[i+1] * (x_ynew_yold_probs[i][0] + x_ynew_yold_probs[i][1])


            zeroTilNowProb = 1.0
            for x in range(self.q): # add to newDistribution[ynew][x]
                if x < self.q - 1:
                    prob = zeroTilNowProb * x_ynew_yold_probs[x][1] * marginalFromNowProb[x+1]
                    zeroTilNowProb *= x_ynew_yold_probs[x][0]
                else:
                    prob = zeroTilNowProb

                prob *= yoldMarginal
                newDistribution.probs[ynew][x] += prob

            if self.iterateLCRVector(yold,yoldMappedTo,lcrvec) == False:
                break

    def calc_probs_of_x_ynew_given_yold(self, yold, yoldMappedTo, lcrvec, upgradedOneHotBinaryMemorylessDistributions, originalOneHotBinaryMemorylessDistributions):
        """For each 0 <= i < q-1, calculate p(x=0, ynew[i]|yold) and p(x=1, ynew[i]|yold), for the one-hot channel
        """
        probs = []

        for i in range(self.q - 1):
            upgradedBinDist = upgradedOneHotBinaryMemorylessDistributions[i]
            originalBinDist = originalOneHotBinaryMemorylessDistributions[i]
            tempProbs = []

            if originalBinDist.calcYMarginal(yold) == 0.0:
                tempProbs = [-1000.0, -1000.0] # to catch errors

            elif lcrvec[i] == lcrCenter:
                for x in range(2):
                    # yold implies ynew[i], with probability 1
                    tempProbs.append(originalBinDist.probXGivenY(x,yold))
            else:
                ynew = yoldMappedTo[yold][i][lcrvec[i]]
                otherynew = yoldMappedTo[yold][i][lcrLeft if lcrvec[i] == lcrRight else lcrRight]
                for x in range(2):
                    fractionMultiplier = upgradedBinDist.probXGivenY(x,ynew)

                    if upgradedBinDist.probXGivenY(x,otherynew) < 0.5:
                        denominator = upgradedBinDist.probXGivenY(x,ynew) - upgradedBinDist.probXGivenY(x,otherynew)
                        numerator = originalBinDist.probXGivenY(x,yold) - upgradedBinDist.probXGivenY(x,otherynew)
                    else:
                        denominator = upgradedBinDist.probXGivenY(1-x,otherynew) - upgradedBinDist.probXGivenY(1-x,ynew)
                        numerator = upgradedBinDist.probXGivenY(1-x,otherynew) - originalBinDist.probXGivenY(1-x,yold) 

                    tempProbs.append(fractionMultiplier * numerator / denominator)

            probs.append(tempProbs)

        return probs

    def initializeLCRVector(self,yold,yoldMappedTo):
        lcrvec = []

        for i in range(self.q-1):
            if yoldMappedTo[yold][i][lcrLeft] == None:
                assert(yoldMappedTo[yold][i][lcrRight] == None)
                lcrvec.append(lcrCenter)
            else:
                assert(yoldMappedTo[yold][i][lcrCenter] == None and yoldMappedTo[yold][i][lcrRight] != None)
                lcrvec.append(lcrLeft)
        return lcrvec

    def iterateLCRVector(self,yold,yoldMappedTo,lcrvec):

        for i in range(self.q-1):
            if lcrvec[i] == lcrLeft:
                assert(yoldMappedTo[yold][i][lcrLeft] != None and yoldMappedTo[yold][i][lcrRight] != None)
                lcrvec[i] = lcrRight
                return True
            elif lcrvec[i] == lcrRight:
                assert(yoldMappedTo[yold][i][lcrLeft] != None and yoldMappedTo[yold][i][lcrRight] != None)
                lcrvec[i] = lcrLeft
                # and don't return (continue to the next i)
            else: # lcrvec[i] == lcrCenter
                assert( lcrvec[i] == lcrCenter )
                assert(yoldMappedTo[yold][i][lcrLeft] == None and yoldMappedTo[yold][i][lcrRight] == None)
                # and don't return (continue to the next i)
        return False

    def calcConversionToYNewMultipliers(self, oneHotBinaryMemorylessDistributions):
        conversionToYNewMultipliers = [1]
        for x in range(1, self.q-1):
            conversionToYNewMultipliers.append(conversionToYNewMultipliers[x-1] * len(oneHotBinaryMemorylessDistributions[x-1].probs))

        return conversionToYNewMultipliers

    def calcNewOutputAlphabetSize(self, oneHotBinaryMemorylessDistributions):
        newOutputAlphabetSize = 1
        for x in range(self.q-1):
             newOutputAlphabetSize *= len(oneHotBinaryMemorylessDistributions[x].probs)

        return newOutputAlphabetSize

    def removeZeroProbOutput(self):
        newProbs = []

        for probTuple in self.probs:
            if sum(probTuple) > 0.0:
                newProbs.append(probTuple)

        self.probs = newProbs

    def yoldToNew_degrade(self, yold, yoldMappedTo, conversionToYNewMultipliers):
        ynew = 0

        for x in range(self.q-1):
            ynew += yoldMappedTo[yold][x] * conversionToYNewMultipliers[x]

        return ynew

    def yoldToNew_upgrade(self, yold, yoldMappedTo, lcrvec, conversionToYNewMultipliers):
        ynew = 0

        for i in range(self.q-1):                            
            # Some output symbols have probability zero, when conditioned on a certain set of
            # possible inputs. These, we map to the first symbol in the corresponding one-hot channel.
            mappedTo = yoldMappedTo[yold][i][lcrvec[i]] if yoldMappedTo[yold][i][lcrvec[i]] != None else 0
            ynew += mappedTo * conversionToYNewMultipliers[i]

        return ynew

# useful channels
def makeQSC(q, p):
    qsc = QaryMemorylessDistribution(q)

    for y in range(q):
        tempProbs = [(1.0 - p)/q if x == y else p/(q*(q-1)) for x in range(q)]
        qsc.append(tempProbs)
    
    return qsc

def makeQEC(q, p):
    qec = QaryMemorylessDistribution(q)

    # first, create the non-erasure symbols
    for y in range(q):
        tempProbs = [(1.0 - p)/q if x == y else 0.0 for x in range(q)]
        qec.append(tempProbs)
    
    tempProbs = [p/q for x in range(q) ]
    qec.append(tempProbs)

    return qec
    
