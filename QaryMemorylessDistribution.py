import BinaryMemorylessDistribution

class QaryMemorylessDistribution:
    def __init__(self, q):
        self.probs = []  # probs[yindex][xindex]
        self.q = q

    def __str__(self):
        s = "Qry memoryless channel with q = " + str(self.q) + " and " + str(len(self.probs)) + " output symbols. The error probability is " + str(self.errorProb()) + ". [p(y,x=0), p(y,x=1), ..., p(y,x=q-1)]: "

        for probPair in self.probs:
            s += "[" 
            for y in range(self.q):
                s += str(probPair[y]) + ", "
            s = s[0:-2] # remove the last ", "
            s += "], "
        s = s[0:-2] # remove the last ", "
        return s

    # polar toolbox
    def errorProb(self):
        totalProbSum = 0.0
        noErrorProbSum = 0.0

        for probPair in self.probs:
            totalProbSum += sum(probPair)
            noErrorProbSum += max(probPair)

        return totalProbSum - noErrorProbSum

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

                newDistribution.probs.append(tempProbs)

        newDistribution.mergeEquivalentSymbols()

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
                    newDistribution.probs.append(tempProbs)

        newDistribution.mergeEquivalentSymbols()

        return newDistribution

    def mergeEquivalentSymbols(self):
        pass #TODO: implement this

# useful channels
def makeQSC(q, p):
    qsc = QaryMemorylessDistribution(q)

    for y in range(q):
        tempProbs = [(1.0 - p)/q if x == y else p/(q*(q-1)) for x in range(q)]
        qsc.probs.append(tempProbs)
    
    return qsc

def makeQEC(q, p):
    qec = QaryMemorylessDistribution(q)

    # first, create the non-erasure symbols
    for y in range(q):
        tempProbs = [(1.0 - p)/q if x == y else 0.0 for x in range(q)]
        qec.probs.append(tempProbs)
    
    tempProbs = [p/q for x in range(q) ]
    qec.probs.append(tempProbs)

    return qec
    
