import BinaryMemorylessDistribution
from BinaryMemorylessDistribution import eta, naturalEta
from math import floor
import numpy as np
import sys
import math
from enum import Enum

# constants
lcrLeft = 0
lcrCenter = 1
lcrRight = 2

class Binning(Enum):
    TalSharovVardy = 1 # standard cell function for static degrade
    PeregTal = 2 # standard cell function for static upgrade

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
        total = 0.0

        for probTuple in self.probs:
            tempProbTuple = probTuple.copy()
            tempProbTuple.sort()
            total += sum( tempProbTuple[:-1] )

        return total

        # totalProbSum = 0.0
        # noErrorProbSum = 0.0
        #
        # for probTuple in self.probs:
        #     totalProbSum += sum(probTuple)
        #     noErrorProbSum += max(probTuple)
        #
        # return totalProbSum - noErrorProbSum

    def conditionalEntropy(self):
        entropySum = 0.0

        for probTuple in self.probs:
            # debugDelta = 0.0
            for p in probTuple:
                entropySum += eta(p)
                # debugDelta += eta(p)

            entropySum -= eta(sum(probTuple))
            # debugDelta -= eta(sum(probTuple))
            # print( debugDelta, probTuple )

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

        M = self.calcMFromL(L)
        
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
        newDistribution.normalize() # for good measure

        return newDistribution

    def putLettersInBins(self, L, binningToUse):
        M = self.calcMFromL(L)

        if binningToUse == Binning.TalSharovVardy:
            mu =  1.0 / (math.e * (M / 2) ) # mu from the paper by Tal, Sharov, and Vardy, but for simplicity, use the natural logarithm (beta = alpha = 1/e)
        elif binningToUse == Binning.PeregTal:
            mu, indexOfBorderCell, maxProbOfBorderCell, gamma =  self.calcMuForPeregTal(M)
        else:
            assert(false)

        dims = [M for i in range(self.q)]

        cellsArray = np.array([set() for _ in range(M ** (self.q))]).reshape(dims)

        # first, put each output letter into the corresponding cell

        for yold, prob in enumerate(self.probs):
            probsum = sum(prob)
            if probsum == 0.0:
                continue

            cell = []
            for x in range(self.q):
                if binningToUse == Binning.TalSharovVardy:
                    cell.append(self.calcCell_TalSharovVardy(prob[x]/probsum, M, mu))
                else:
                    cell.append(self.calcCell_PeregTal(prob[x]/probsum, M, mu, indexOfBorderCell, maxProbOfBorderCell, gamma))

            cellsArray[tuple(cell)] |= {yold}

        return cellsArray

    def degrade_static(self, L, binningToUse = Binning.TalSharovVardy):
        # TODO: cells or bins?
        cellsArray = self.putLettersInBins(L, binningToUse)
        newDistribution = QaryMemorylessDistribution(self.q)

        for setOfY in np.nditer(cellsArray, flags=["refs_ok"]):
            actualSet = setOfY.item()
            if len(actualSet) == 0:
                continue

            ynewProb = [0.0 for i in range(self.q)]

            for yold in actualSet:
                for x in range(self.q):
                    ynewProb[x] += self.probs[yold][x]
            
            newDistribution.probs.append(ynewProb)
        newDistribution.normalize() # for good measure
        return newDistribution

    def calcCell_TalSharovVardy(self, postProb, M, mu):
        if postProb <= 1.0 / math.e:
            cell = floor(naturalEta(postProb) / mu)
        else:
            cell = M - 1 - floor(naturalEta(postProb) / mu)

        assert(cell >= 0 and cell < M)
        return cell

    def upgrade(self, L):
        return self.upgrade_dynamic(L)

    def upgrade_dynamic(self, L):
        oneHotBinaryMemorylessDistributions = self.oneHotBinaryMemorylessDistributions()
        M = self.calcMFromL(L)

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
        newDistribution.normalize() # for good measure

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

    def upgrade_static(self, L, binningToUse = Binning.PeregTal):
        # TODO: cells or bins?
        cellsArray = self.putLettersInBins(L, binningToUse)
        newDistribution = QaryMemorylessDistribution(self.q)

        # The first q symbols are the boost symbols
        allzeroprobvector = [0.0 for i in range(self.q)]
        for x in range(self.q):
            newDistribution.probs.append(allzeroprobvector.copy())

        for setOfY in np.nditer(cellsArray, flags=["refs_ok"]):
            actualSet = setOfY.item()
            if len(actualSet) == 0:
                continue
            self.upgradeCellToSymbolPlusBoosts(actualSet, newDistribution.probs)

        newDistribution.removeZeroProbOutput()
        newDistribution.normalize() # for good measure
        return newDistribution

    def upgradeCellToSymbolPlusBoosts(self, actualSet, newprobs):
        ynewProb = [0.0 for i in range(self.q)]

        if len(actualSet) == 1:
            for yold in actualSet:
                newprobs.append(self.probs[yold])
            return

        # find the leading input symbol
        leadingX = -1
        leadingPostProb = -1.0

        for yold in actualSet:
            probSum = sum(self.probs[yold])
            assert(probSum) > 0.0

            for x in range(self.q):
                postProb = self.probs[yold][x]/probSum
                if postProb > leadingPostProb:
                    leadingX = x
                    leadingPostProb = postProb

        # the leading input symbol is leadingX
        # now, find the posterior of the cell
        cellPosteriorProb = [1.0 for i in range(self.q)]
        cellPosteriorProb[leadingX] = 0.0
        
        for yold in actualSet:
            probSum = sum(self.probs[yold])
            for x in range(self.q):
                if x == leadingX:
                    continue
                postProb = self.probs[yold][x]/probSum
                if postProb < cellPosteriorProb[x]:
                    cellPosteriorProb[x] = postProb

        cellPosteriorProb[leadingX] = 1.0 - sum(cellPosteriorProb)

        # now calculate alpha_x(y), as per (10), and add probabilities 
        # to regular symbols (13a) and boost symbols (13b)
        ynewProb = [0.0 for i in range(self.q)]

        print( actualSet )
        for yold in actualSet:
            debugProb = [0.0 for i in range(self.q)]
            for x in range(self.q):
                if self.probs[yold][x] > 0.0:
                    alphaxy = (cellPosteriorProb[x] / self.probs[yold][x]) * (self.probs[yold][leadingX] / cellPosteriorProb[leadingX] )
                else:
                    alphaxy = 1.0

                alphaxy = min(1.0,alphaxy) # fix floating point rounding errors

                # Now that we've calculate alphaxy, calculate the probability to add
                ynewProb[x] += self.probs[yold][x] * alphaxy
                debugProb[x] = self.probs[yold][x] * alphaxy

                # for the boost symbol
                newprobs[x][x] += (1.0 - alphaxy) * self.probs[yold][x]
            
            debugProbTwo = self.probs[yold].copy()
            debugSum = sum(debugProb)
            debugSumTwo = sum(debugProbTwo)
            for x in range(self.q):
                debugProb[x] /= debugSum
                debugProbTwo[x] /= debugSumTwo

            # print( cellPosteriorProb, debugProb, debugProbTwo )
            print( self.probs[yold], cellPosteriorProb, debugProbTwo )



        print(" * ", ynewProb)

        newprobs.append(ynewProb)

    def calcCell_PeregTal(self, postProb, M, mu, indexOfBorderCell, maxProbOfBorderCell, gamma):
        if postProb > maxProbOfBorderCell:
            cellIndex = indexOfBorderCell + 1 + floor( (postProb - maxProbOfBorderCell) * mu )
        elif postProb > 1.0 / gamma:
            cellIndex = indexOfBorderCell
        else:
            cellIndex = floor( naturalEta(postProb) * mu )

        if cellIndex > M - 1:
            cellIndex = M - 1

        return cellIndex

    def calcMuForPeregTal(self, M):
        """Calculate the value of the optimal (largest) mu, as well as other parameters of interest, for a given M."""

        # Consider Figure 3 in the IEEE-IT paper by Pereg and Tal.
        # For brevity, let gamma = 1/e^2.
        # Since larger mu implies a better bound in terms of the loss in mutual information,
        # we would like mu to be as large as possible, such that the number of cells is M
        # (no point in having the number of cells be less than M, since this means we can
        # enlarge mu). To find this largest possible mu, we will implement the following
        # algorithm. Let mu be given.
        # * Divide the cells to the left of the blue line exactly as in Figure 3. That is
        # all the cells containing only probabilities less than or equal to gamma are ordered such
        # that the difference between naturalEta of the left and right (min and max) probabilities
        # is exactly 1/mu.
        # The number of these cells is floor(2 * gamma * mu), since naturalEta(gamma) = 2 * gamma
        # Next, *unlike* Figure 3, divide the cells to the right of the blue line such that
        # the width of the right most cell (20 in the Figure) is exactly 1/mu, and keep
        # going left as long as you can; that is, as long as the cell does not contain probabilities
        # strictly less than gamma. The number of these cells is exactly floor( (1-gamma) * mu ). Our first
        # observation is that we must have a non-empty region that is not covered by cells. That is, if this
        # were not the case, then both 2 * gamma * mu and (1-gamma) * mu had to be integers, which
        # can't be the case, since gamma = 1/e^2 is irrational. Let us call this region the leftover cell.
        # So, our total number of cells is
        # floor(2 * gamma * mu) + floor( (1-gamma) * mu ) + 1
        # If the total number of cells is larger than M, we reduce mu
        # If the total number of cells is smaller than M, we enlarge mu
        # If the total number of cells is exactly M, then we will shortly describe what we do.
        # In aid of this, first, consider two boundary values of mu: the smallest and largest mu
        # (up to some small epsilon of our choosing) for which the number of cells is exactly M. 
        # Our aim now is to show that some "good" criterion we will shortly define holds for one of
        # these extreme mu and does not hold for the other. Specifically, the good criterion is
        # that the width of the leftover cell is at most 1/mu (the difference between the maximum and
        # minimum probability contained in the cell is at most 1/mu), and also that
        # max_p naturalEta(p) - min_p naturalEta(p) is at most 1/mu, where p ranges over all the probabilities in the
        # cell. For the larger extreme mu (we are about to transition from having M cells into having M+1 cells)
        # the good condition does not hold, since a new cell is about to be born, either because the dotted
        # red line immediately to the right of the blue line is very close to 1/mu, or because the top most
        # horizontal dotted red line to the left of the blue line almost naturalEta(gamma) - 1/mu. This, and the
        # fact that we have some slack on the other side of the blue line (again, from the irrationality of gamma)
        # implies that we are not fulfilling at least one of requirements of the good condition. 
        # Conversely, for the other extreme mu, the good condition is met. To see this, note that in the case,
        # either the dotted red line immediately to the right of the blue
        # line is almost touching it, or the vertical dotted red line immediately to the left of the blue line is
        # almost touching it. For both cases, the good condition holds, since gamma is the point for which the
        # slope equals 1, and it is decreasing in p.
        # So, to sum up, if the number of good cells is M:
        # * if the good condition is met, we enlarge mu
        # * if the good condition is not met, we reduce mu

        assert( M > 1)

        muUpper = 2.0 * M
        muLower = 1.0
        gamma = 1.0/(math.e ** 2)

        while muUpper - muLower > 0.0000001:
            mu = (muUpper + muLower) / 2.0
            cellsToLeftOfGamma = floor(2.0 * gamma * mu) # for a given mu, this is exactly the number of cells to the left of gamma
            cellsToRightOfGamma = floor((1.0 - gamma) * mu)
            totalCells = cellsToLeftOfGamma + cellsToRightOfGamma + 1
            if totalCells < M: # mu too small
                muLower = mu
                continue
            if totalCells > M: # mu too large
                muUpper = mu
                continue
            
            # right side and left side mean left and right of blue line, respectively
            maxProbOfBorderCell = 1.0 - (1.0/mu) * cellsToRightOfGamma
            maxEtaOfBorderCell = naturalEta(maxProbOfBorderCell) if maxProbOfBorderCell < 1.0 / math.e else naturalEta(1.0 / math.e)
            minEtaOnLeftSideOfBorderCell = (1.0 / mu) * cellsToLeftOfGamma

            if maxEtaOfBorderCell - minEtaOnLeftSideOfBorderCell > 1.0 / mu: # mu too large
                muUpper = mu
                continue

            if maxProbOfBorderCell > 1.0 / math.e:
                minEtaOnRightSideOfBorderCell = min(naturalEta(maxProbOfBorderCell), naturalEta(gamma))
                if maxEtaOfBorderCell - minEtaOnRightSideOfBorderCell > 1.0 / mu: # mu too large
                    muUpper = mu
                    continue

            # is the width of the border cell greater than 1/mu?
            if maxProbOfBorderCell > 1.0/mu and naturalEta(maxProbOfBorderCell - 1.0/mu) > minEtaOnLeftSideOfBorderCell: # mu too large
                muUpper = mu
                continue

            # got here, so mu is too small
            muLower = mu

        # to be on the safe side
        mu = muLower

        # re-calculate key parameters
        cellsToLeftOfGamma = floor(2.0 * gamma * mu) # for a given mu, this is exactly the number of cells to the left of gamma
        cellsToRightOfGamma = floor((1.0 - gamma) * mu)

        # print(cellsToLeftOfGamma, cellsToRightOfGamma, mu, gamma)
        assert( cellsToLeftOfGamma + cellsToRightOfGamma + 1 == M )
        indexOfBorderCell  = cellsToLeftOfGamma
        maxProbOfBorderCell = 1.0 - (1.0/mu) * cellsToRightOfGamma

        return mu, indexOfBorderCell, maxProbOfBorderCell, gamma

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

    def normalize(self):
        probSum = 0.0

        tempProbs = []
        for probList in self.probs:
            for prob in probList:
                tempProbs.append(prob)

        tempProbs.sort()

        for tp in tempProbs:
            probSum += tp

        for y in range(len(self.probs)):
            for x in range(self.q):
                self.probs[y][x] /= probSum

    def calcMFromL(self, L):
        M = floor( L ** (1.0/(self.q-1)) + sys.float_info.epsilon )
        return M

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

def makeInputDistribution(probs):
    dist = QaryMemorylessDistribution(len(probs))
    dist.append(probs)
    dist.normalize() # for good measure

    return dist

# My conjecture is that the mutual information of this channel, in the limit is log2(q) - q \cdot \frac{\int_0^1 -p log2(p) (1-p)^{q-1}}{\int_0^1 (1-p)^{q-1}}
def makeQuantizedUniform(q, T):
    """Make a joint distribution such that for an input x and an output (a_0,a_1,...,a_{q-1}), where the a_i are non-negative and sum to T, the joint probability is a_x / M, where M is the number of possible outputs. That is, the number of non-negative vectors  (a_0,a_1,...,a_{q-1}), whose entries sum to L. That is, M = binom{L + q - 1}{q-1}. See the paper by Tal, On the construction of polar codes for channels with moderate input alphabet sizes, and the paper by Kartowsky and Tal, Greedy-Merge Degrading has Optimal Power-Law, in which this channel is proven to be hard to degrade and hard to upgrade, respectively."""
    # M is the number of output letters we will produce
    M = math.comb(T + q - 1, q-1) 

    quantizedUniform = QaryMemorylessDistribution(q)

    outputLetter = []
    level = q

    recursivlyBuildQuantizedUniform(quantizedUniform, outputLetter, level, T, M)

    return quantizedUniform

def recursivlyBuildQuantizedUniform(quantizedUniform, outputLetter, level, T, M):
    if level == 0:
        prob = []
        for x in range(quantizedUniform.q):
            prob.append((1.0 * outputLetter[x])/(M*T))
        quantizedUniform.probs.append(prob)
        return

    if level == 1:
        outputLetter.append(T - sum(outputLetter))
        recursivlyBuildQuantizedUniform(quantizedUniform, outputLetter, level-1, T, M)
        outputLetter.pop()
        return

    for t in range(0, T+1 - sum(outputLetter)):
        outputLetter.append(t)
        recursivlyBuildQuantizedUniform(quantizedUniform, outputLetter, level-1, T, M)
        outputLetter.pop()
