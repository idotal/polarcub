import numpy as np
import VectorDistribution
import math

class Vertex():
    def __init__(self, stateId = -1, verticalPosInLayer = -1, layer = -1, vertexProb = -1.0):
        """A vertex has three properties that identifies it: statedId, verticalPosInLayer, and layer.

        Args:
            stateId (int): The state in the Markov chain associated with the the input distribution (just before the current input)

            verticalPosInLayer (int): a value of i means that i symbols have already appeared at the output

            layer (int): a value of j means that j symbols have already appeared at the input

            vertexProb (float): applicable to vertices in first and last layer, all other vertices have a defualt of -1.0
        """

        self.stateId = stateId
        self.verticalPosInLayer = verticalPosInLayer
        self.layer = layer
        self.outgoingEdges = {}
        self.incomingEdges = {}
        self.vertexProb = vertexProb

    def sanityCheck(self):
        assert(self.stateId >= 0)
        assert(self.verticalPosInLayer >= 0)
        assert(self.layer >= 0)

    def getKey(self):
        self.sanityCheck()
        return (self.stateId, self.verticalPosInLayer, self.layer)

    def toString(self, printEdges = True):
        s = "* " + str(self.getKey()) + ": stateId = " + str(self.stateId) + ", verticalPosInLayer = " + str(self.verticalPosInLayer) + ", layer = " + str(self.layer) + ", vertexProb = " + str(self.vertexProb) + "\n"

        if printEdges == True:
            s += "  incoming edges (" + str(len(self.incomingEdges)) + "):\n"
            for (edgeKey, edge) in self.incomingEdges.items():
                s += "    " + edge.toString() + "\n"
    
            s += "  outgoing edges (" + str(len(self.outgoingEdges)) + "):\n"
            for (edgeKey, edge) in self.outgoingEdges.items():
                s += "    " + edge.toString() + "\n"

        return s

    def __str__(self):
        return self.toString()

class Edge():
    def __init__(self, fromVertex=None, toVertex=None, edgeLabel=-1, edgeProb=-1.0):
        self.fromVertex = fromVertex
        self.toVertex = toVertex
        self.edgeLabel = edgeLabel
        self.edgeProb = edgeProb

    def sanityCheck(self):
        assert(self.fromVertex is not None)
        assert(self.toVertex is not None)
        assert(self.edgeLabel != -1)
        assert(self.fromVertex.layer + 1 == self.toVertex.layer)
        assert(self.edgeProb >= 0.0  and self.edgeProb <= 1.0)

    def getKey(self):
        self.sanityCheck()
        return (self.fromVertex.stateId, self.fromVertex.verticalPosInLayer, self.fromVertex.layer, self.toVertex.stateId, self.toVertex.verticalPosInLayer, self.edgeLabel)

    def toString(self):
        s = str(self.fromVertex.getKey())

        s += " --[lbl=" + str(self.edgeLabel) + ",p=" + str(self.edgeProb) + "]--> "

        s += str(self.toVertex.getKey())

        return s

        # s = "The edge has label = " + str(self.edgeLabel) + " and edgeProb = " + str(self.edgeProb) + "\n"
        #
        # if printFromVertex == True:
        #     s += "The edge leaves the following vertex:\n"
        #     s += self.fromVertex.toString(printEdges=False)
        #
        # if printToVertex == True:
        #     s += "The edge enters the following vertex:\n"
        #     s += self.toVertex.toString(printEdges=False)
        #
        # return s

    def __str__(self):
        return self.toString()
        
class BinaryTrellis(VectorDistribution.VectorDistribution):
    def __init__(self, length):
        """Initialize an empty trellis

        Args:
            length (int): the number of inputs (not including the guard bands). That is, 2**n, where n is the number of polar transforms. Thus, the number of layers is length + 1.

        Returns:
            None
        """
        assert( length > 0 )
        self.length = length
        self.layers = length + 1
        self.verticesInLayer = [] # vericesInLayer[layer] is a dictionary accessed by vertexKey. So, vericesInLayer[layer][vertexKey].

        for l in range(self.layers):
            self.verticesInLayer.append( {} )

    def setVertexProb(self, vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb):
        vertex = self.__getVertexAndAddIfNeeded(vertex_stateId, vertex_verticalPosInLayer, vertex_layer)
        vertex.vertexProb = vertexProb

    def addToEdgeProb(self, fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, probToAdd):
        fromVertex = self.__getVertexAndAddIfNeeded(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer)
        toVertex =  self.__getVertexAndAddIfNeeded(toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer)
        self.addToEdgeProb_vertexReferences(fromVertex, toVertex, edgeLabel, probToAdd)

    def addToEdgeProb_vertexReferences(self, fromVertex, toVertex, edgeLabel, probToAdd):
        edge = self.__getEdgeAndAddIfNeeded(fromVertex, toVertex, edgeLabel)
        edge.edgeProb += probToAdd

    def getEdgeProb(self, fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel):
        fromVertex_key = Vertex(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer).getKey()
        toVertex_key = Vertex(toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer).getKey()

        fromVertex = self.verticesInLayer[fromVertex_layer][fromVertex_key]
        toVertex = self.verticesInLayer[toVertex_layer][toVertex_key]
        return self.getEdgeProb_vertexReferences(fromVertex, toVertex, edgeLabel)

    def getEdgeProb_vertexReferences(self, fromVertex, toVertex, edgeLabel):
        edgeKey = Edge(fromVertex, toVertex, edgeLabel).getKey()
    
        assert(edgeKey in fromVertex.outgoingEdges)
        assert(edgeKey in toVertex.incomingEdges)
        edge = fromVertex.outgoingEdges[edgeKey]
        return edge.edgeProb

    def __getVertexAndAddIfNeeded(self, stateId, verticalPosInLayer, layer):
        possibleNewVertex = Vertex(stateId, verticalPosInLayer, layer)
        vertexKey = possibleNewVertex.getKey()

        if vertexKey not in self.verticesInLayer[layer]:
            self.verticesInLayer[layer][vertexKey] = possibleNewVertex

        return self.verticesInLayer[layer][vertexKey]

    def __getEdgeAndAddIfNeeded(self, fromVertex, toVertex, edgeLabel, edgeProb = 0.0):
        possibleNewEdge = Edge(fromVertex, toVertex, edgeLabel, edgeProb)
        edgeKey = possibleNewEdge.getKey()

        if edgeKey not in fromVertex.outgoingEdges:
            assert(edgeKey not in toVertex.incomingEdges)
            fromVertex.outgoingEdges[edgeKey] = possibleNewEdge
            toVertex.incomingEdges[edgeKey] = possibleNewEdge

        assert(fromVertex.outgoingEdges[edgeKey] is toVertex.incomingEdges[edgeKey])

        return fromVertex.outgoingEdges[edgeKey]

    def __str__(self):
        s = "The input alphabet size is " + str(2) + "\n"
        s += "The number of layers is " + str(self.layers) + "\n"
        s += "The number of vertices in each layers is: \n"

        for l in range(self.layers):
            s += str(len(self.verticesInLayer[l]))
            if l < self.layers - 1:
                s += ", "
            else:
                s += "\n"
        
        for l in range(self.layers):
            s += "For layer " + str(l) + ", these vertices are:\n"

            for (vertexKey, vertex) in self.verticesInLayer[l].items():
                s += vertex.toString(printEdges=True) + "\n"

        return s

def buildTrellis_uniformInput_deletion(receivedWord, codewordLength, deletionProb, trimmedZerosAtEdges=False):
    trellis = BinaryTrellis(codewordLength)
    deletionCount = codewordLength - len(receivedWord)
    inputProb = [0.5, 0.5]

    vertex_stateId = 0
    vertex_verticalPosInLayer = 0
    vertex_layer = 0
    vertexProb = 1.0


    trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

    vertex_verticalPosInLayer = len(receivedWord) 
    vertex_layer = codewordLength

    trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

    if trimmedZerosAtEdges == True:
        assert( len(receivedWord) == 0 or (receivedWord[0] == 1 and receivedWord[-1] == 1))

    for l in range(codewordLength):
        vpos_min = max(0, l - deletionCount)
        vpos_max = min(l, len(receivedWord))

        for vpos in range(vpos_min, vpos_max+1):
            if vpos < len(receivedWord): # then we can have a non-deletion event
                fromVertex_stateId = 0
                fromVertex_verticalPosInLayer = vpos
                fromVertex_layer = l
                toVertex_stateId = 0
                toVertex_verticalPosInLayer = vpos+1 # a non-deletion
                toVertex_layer = l+1
                edgeLabel = receivedWord[vpos]
                probToAdd = inputProb[edgeLabel] * (1.0 - deletionProb)

                trellis.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, probToAdd)
            if l - vpos < deletionCount: # then we can have a deletion event
                for edgeLabel in range(2):
                    fromVertex_stateId = 0
                    fromVertex_verticalPosInLayer = vpos
                    fromVertex_layer = l
                    toVertex_stateId = 0
                    toVertex_verticalPosInLayer = vpos # a deletion
                    toVertex_layer = l+1

                    if trimmedZerosAtEdges == False or edgeLabel == 1 or (vpos > 0 and vpos < len(receivedWord)):
                        probToAdd = inputProb[edgeLabel] * deletionProb
                    else:
                        probToAdd = inputProb[edgeLabel] 
                    
                    trellis.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, probToAdd)

    return trellis

def addDeletionGuardBands(codeword, n, n0, xi):
    """Add deletion guard bands, according to equations (68) and (69) in the deletions paper"""

    if n <= n0:
        return codeword
    else:
        assert(len(codeword) % 2 == 0)

        ln = math.floor(2**((1-xi)*(n-1)))

        print( " codeword = ", codeword, " n = ", n, " ln =", ln )

        leftPartOfCodeword = codeword[0:len(codeword)//2]
        rightPartOfCodeword = codeword[len(codeword)//2:len(codeword)]

        leftAfterGuardBandsAdded = addDeletionGuardBands(leftPartOfCodeword, n-1, n0, xi)
        rightAfterGuardBandsAdded = addDeletionGuardBands(rightPartOfCodeword, n-1, n0, xi)

        guardBand = []

        for i in range(ln):
            guardBand.append(0)

        allTogether = leftAfterGuardBandsAdded + guardBand + rightAfterGuardBandsAdded

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
