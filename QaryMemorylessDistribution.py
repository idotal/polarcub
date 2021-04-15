import BinaryMemorylessDistribution
from BinaryMemorylessDistribution import eta
from math import floor

class QaryMemorylessDistribution:
    def __init__(self, q):
        self.probs = []  # probs[yindex][xindex]
        self.q = q

    def __str__(self):
        s = "Qry memoryless channel with q = " + str(self.q) + " and " + str(len(self.probs)) + " output symbols. The error probability is " + str(self.errorProb()) + ". [p(y,x=0), p(y,x=1), ..., p(y,x=q-1)]: "

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
        binaryMemorylessDistributions = []
        q = self.q

        for i in range(q-1):
            binaryMemorylessDistributions.append(BinaryMemorylessDistribution.BinaryMemorylessDistribution())
            binaryMemorylessDistributions[i].auxiliary = []

        # calculate the marginals p(X = j)
        marginals = []

        for x in range(q):
            tempSum = 0.0
            for probTuple in self.probs:
                tempSum += probTuple[x]
            marginals.append(tempSum)

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
        oneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()

        M = floor( L ** (1.0/(self.q-1)) )
        
        degradedOneHotBinaryMemorylessDistributions = []
        for x in range(self.q-1):
            degradedOneHotBinaryMemorylessDistributions.append( oneHotBinaryMemorylessDistributions[x].degrade(M) )

        # for x in range(self.q-1):
        #     print(degradedOneHotBinaryMemorylessDistributions[x])

        newDistribution = QaryMemorylessDistribution(self.q)
        newOutputAlphabetSize = 1
        for x in range(self.q-1):
             newOutputAlphabetSize *= len(degradedOneHotBinaryMemorylessDistributions[x].probs)

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

        yoldToNewBasis = [1]
        for x in range(1, self.q-1):
            yoldToNewBasis.append(yoldToNewBasis[x-1] * len(degradedOneHotBinaryMemorylessDistributions[x-1].probs))

        for yold, yoldprobs in enumerate(self.probs):
            ynew = self.yoldToNew_degrade(yold, yoldMappedTo, yoldToNewBasis)

            for x in range(self.q):
                newDistribution.probs[ynew][x] += yoldprobs[x]

        newDistribution.removeZeroProbOutput()

        return newDistribution

    def upgrade(self, L):
        oneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()
        M = floor( L ** (1.0/(self.q-1)) )

        upgradedOneHotBinaryMemorylessDistributions = []
        for x in range(self.q-1):
            upgradedOneHotBinaryMemorylessDistributions.append( oneHotBinaryMemorylessDistributions[x].upgrade(M) )


    def removeZeroProbOutput(self):
        newProbs = []

        for probTuple in self.probs:
            if sum(probTuple) > 0.0:
                newProbs.append(probTuple)

        self.probs = newProbs

    def yoldToNew_degrade(self, yold, yoldMappedTo, yoldToNewBasis):
        ynew = 0

        for x in range(self.q-1):
            ynew += yoldMappedTo[yold][x] * yoldToNewBasis[x]

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
    
