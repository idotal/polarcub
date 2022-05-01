import math


def addDeletionGuardBands(encodedVector, n, n0, xi, numberOfOnesToAddAtBothEndsOfGuardbands=0):
    """Add deletion guard bands, according to equations (68) and (69) in the deletions paper

    Optionally, add ones at both ends (and also at the start and end of the encoded vector)
    """

    if n <= n0:
        if numberOfOnesToAddAtBothEndsOfGuardbands > 0:
            allTogether = []
            allTogether.extend([1 for i in range(numberOfOnesToAddAtBothEndsOfGuardbands)])
            allTogether.extend(encodedVector)
            allTogether.extend([1 for i in range(numberOfOnesToAddAtBothEndsOfGuardbands)])
            return allTogether
        else:
            return encodedVector
    else:
        assert (len(encodedVector) % 2 == 0)

        ln = math.floor(2 ** ((1 - xi) * (n - 1)))

        # print( " encodedVector = ", encodedVector, " n = ", n, " ln =", ln )

        leftPartOfEncodedVector = encodedVector[0:len(encodedVector) // 2]
        rightPartOfEncodedVector = encodedVector[len(encodedVector) // 2:len(encodedVector)]

        leftAfterGuardBandsAdded = addDeletionGuardBands(leftPartOfEncodedVector, n - 1, n0, xi,
                                                         numberOfOnesToAddAtBothEndsOfGuardbands)
        rightAfterGuardBandsAdded = addDeletionGuardBands(rightPartOfEncodedVector, n - 1, n0, xi,
                                                          numberOfOnesToAddAtBothEndsOfGuardbands)

        guardBand = []

        for i in range(ln):
            guardBand.append(0)

        allTogether = []
        allTogether.extend(leftAfterGuardBandsAdded)
        allTogether.extend(guardBand)
        allTogether.extend(rightAfterGuardBandsAdded)

        return allTogether


def removeDeletionGuardBands(receivedWord, n, n0):
    """Undo the addition of guard bands (trim zeros), and return a list of the resulting substrings.

    """

    trimmedReceivedWord = trimZerosAtEdges(receivedWord)

    if n <= n0:
        return [trimmedReceivedWord]
    else:
        leftHalfOfTrimmedReceivedWord = trimmedReceivedWord[0:len(trimmedReceivedWord) // 2]
        rightHalfOfTrimmedReceivedWord = trimmedReceivedWord[len(trimmedReceivedWord) // 2:len(trimmedReceivedWord)]

        leftList = removeDeletionGuardBands(leftHalfOfTrimmedReceivedWord, n - 1, n0)
        rightList = removeDeletionGuardBands(rightHalfOfTrimmedReceivedWord, n - 1, n0)

        return leftList + rightList


def trimZerosAtEdges(receivedWord):
    trimmedReceivedWord = []

    firstOneIndex = -1

    # print("trimZerosAtEdges, receivedWord = ", receivedWord )
    for i in range(len(receivedWord)):
        if receivedWord[i] == 1:
            firstOneIndex = i
            break

    if firstOneIndex == -1:
        return trimmedReceivedWord  # which is empty

    lastOneIndex = -1
    for i in range(len(receivedWord) - 1, -1, -1):
        if receivedWord[i] == 1:
            lastOneIndex = i
            break

    assert (lastOneIndex != -1)

    for i in range(firstOneIndex, lastOneIndex + 1):
        trimmedReceivedWord.append(receivedWord[i])

    # print("trimZerosAtEdges, trimmedReceivedWord = ", trimmedReceivedWord )

    return trimmedReceivedWord
