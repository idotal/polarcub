import math

# TODO: change "codeword" to "encodedVector"
def addDeletionGuardBands(codeword, n, n0, xi):
    """Add deletion guard bands, according to equations (68) and (69) in the deletions paper"""

    if n <= n0:
        return codeword
    else:
        assert(len(codeword) % 2 == 0)

        ln = math.floor(2**((1-xi)*(n-1)))

        # print( " codeword = ", codeword, " n = ", n, " ln =", ln )

        leftPartOfCodeword = codeword[0:len(codeword)//2]
        rightPartOfCodeword = codeword[len(codeword)//2:len(codeword)]

        leftAfterGuardBandsAdded = addDeletionGuardBands(leftPartOfCodeword, n-1, n0, xi)
        rightAfterGuardBandsAdded = addDeletionGuardBands(rightPartOfCodeword, n-1, n0, xi)

        guardBand = []

        for i in range(ln):
            guardBand.append(0)

        allTogether = []
        allTogether.extend(leftAfterGuardBandsAdded)
        allTogether.extend(guardBand)
        allTogether.extend(rightAfterGuardBandsAdded)

        return allTogether


def removeDeletionGuardBands(receivedWord, n, n0):
    """Undo the addition of guard bands (plus trim), and return a list of the resulting substrings."""

    trimmedReceivedWord = trimZerosAtEdges(receivedWord)

    if n <= n0:
        return [trimmedReceivedWord]
    else:
        leftHalfOfTrimmedReceivedWord = trimmedReceivedWord[0:len(trimmedReceivedWord)//2]
        rightHalfOfTrimmedReceivedWord = trimmedReceivedWord[len(trimmedReceivedWord)//2:len(trimmedReceivedWord)]

        leftList = removeDeletionGuardBands(leftHalfOfTrimmedReceivedWord, n-1, n0)
        rightList = removeDeletionGuardBands(rightHalfOfTrimmedReceivedWord,n-1,n0)

        return leftList + rightList

def trimZerosAtEdges(receivedWord):
        trimmedReceivedWord = []

        firstOneIndex = -1

        print( receivedWord )
        for i in range(len(receivedWord)):
            if receivedWord[i] == 1:
                firstOneIndex = i
                break
        
        if firstOneIndex == -1:
            return trimmedReceivedWord # which is empty

        lastOneIndex = -1
        for i in range(len(receivedWord)-1,-1,-1):
            if receivedWord[i] == 1:
                lastOneIndex = i
                break

        assert(lastOneIndex != -1)

        for i in range(firstOneIndex, lastOneIndex+1):
            trimmedReceivedWord.append(receivedWord[i])

        return trimmedReceivedWord
