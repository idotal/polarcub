import random

import numpy as np
import scipy.special

import VectorDistribution


class Vertex():
    def __init__(self, stateId=-1, verticalPosInLayer=-1, layer=-1, vertexProb=-1.0):
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
        assert (self.stateId >= 0)
        assert (self.verticalPosInLayer >= 0)
        assert (self.layer >= 0)

    def getKey(self):
        self.sanityCheck()
        return (self.stateId, self.verticalPosInLayer, self.layer)

    def toString(self, printEdges=True):
        s = "* " + str(self.getKey()) + ": stateId = " + str(self.stateId) + ", verticalPosInLayer = " + str(
            self.verticalPosInLayer) + ", layer = " + str(self.layer) + ", vertexProb = " + str(self.vertexProb) + "\n"

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
        assert (self.fromVertex is not None)
        assert (self.toVertex is not None)
        assert (self.edgeLabel != -1)
        assert (self.fromVertex.layer + 1 == self.toVertex.layer)
        assert (self.edgeProb >= 0.0 and self.edgeProb <= 1.0)

    def getKey(self):
        self.sanityCheck()
        return (
        self.fromVertex.stateId, self.fromVertex.verticalPosInLayer, self.fromVertex.layer, self.toVertex.stateId,
        self.toVertex.verticalPosInLayer, self.edgeLabel)

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
        assert (length > 0)
        self.length = length
        self.layers = length + 1
        self.verticesInLayer = []  # vericesInLayer[layer] is a dictionary accessed by vertexKey. So, vericesInLayer[layer][vertexKey].

        for l in range(self.layers):
            self.verticesInLayer.append({})

    def __len__(self):
        return self.length

    def setVertexProb(self, vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb):
        vertex = self.__getVertexAndAddIfNeeded(vertex_stateId, vertex_verticalPosInLayer, vertex_layer)
        vertex.vertexProb = vertexProb

    def addToEdgeProb(self, fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId,
                      toVertex_verticalPosInLayer, toVertex_layer, edgeLabel, probToAdd):
        fromVertex = self.__getVertexAndAddIfNeeded(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer)
        toVertex = self.__getVertexAndAddIfNeeded(toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer)
        self.addToEdgeProb_vertexReferences(fromVertex, toVertex, edgeLabel, probToAdd)

    def addToEdgeProb_vertexReferences(self, fromVertex, toVertex, edgeLabel, probToAdd):
        edge = self.__getEdgeAndAddIfNeeded(fromVertex, toVertex, edgeLabel)
        edge.edgeProb += probToAdd

    def getEdgeProb(self, fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer, toVertex_stateId,
                    toVertex_verticalPosInLayer, toVertex_layer, edgeLabel):
        fromVertex_key = Vertex(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer).getKey()
        toVertex_key = Vertex(toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer).getKey()

        fromVertex = self.verticesInLayer[fromVertex_layer][fromVertex_key]
        toVertex = self.verticesInLayer[toVertex_layer][toVertex_key]
        return self.getEdgeProb_vertexReferences(fromVertex, toVertex, edgeLabel)

    def getEdgeProb_vertexReferences(self, fromVertex, toVertex, edgeLabel):
        edgeKey = Edge(fromVertex, toVertex, edgeLabel).getKey()

        assert (edgeKey in fromVertex.outgoingEdges)
        assert (edgeKey in toVertex.incomingEdges)
        edge = fromVertex.outgoingEdges[edgeKey]
        return edge.edgeProb

    def __getVertexAndAddIfNeeded(self, stateId, verticalPosInLayer, layer):
        possibleNewVertex = Vertex(stateId, verticalPosInLayer, layer)
        vertexKey = possibleNewVertex.getKey()

        if vertexKey not in self.verticesInLayer[layer]:
            self.verticesInLayer[layer][vertexKey] = possibleNewVertex

        return self.verticesInLayer[layer][vertexKey]

    def __getEdgeAndAddIfNeeded(self, fromVertex, toVertex, edgeLabel, edgeProb=0.0):
        possibleNewEdge = Edge(fromVertex, toVertex, edgeLabel, edgeProb)
        edgeKey = possibleNewEdge.getKey()

        if edgeKey not in fromVertex.outgoingEdges:
            assert (edgeKey not in toVertex.incomingEdges)
            fromVertex.outgoingEdges[edgeKey] = possibleNewEdge
            toVertex.incomingEdges[edgeKey] = possibleNewEdge

        assert (fromVertex.outgoingEdges[edgeKey] is toVertex.incomingEdges[edgeKey])

        return fromVertex.outgoingEdges[edgeKey]

    def __str__(self):
        return self.toString()

    def toString(self):
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

    def minusTransform(self):
        return self.__miusPlusTransform()

    def plusTransform(self, decisionVector):
        return self.__miusPlusTransform(decisionVector)

    def __miusPlusTransform(self, decisionVector=None):
        """If decisionVector = None, apply a minus transform, else apply a plus transform
        """

        newTrellis = BinaryTrellis(self.length // 2)

        if decisionVector is not None:
            assert (len(decisionVector) == self.length // 2)

        # Copy vertices at start and end
        for (vkey, v) in self.verticesInLayer[0].items():
            newTrellis.setVertexProb(v.stateId, v.verticalPosInLayer, 0, v.vertexProb)

        for (vkey, v) in self.verticesInLayer[self.length].items():
            newTrellis.setVertexProb(v.stateId, v.verticalPosInLayer, self.length // 2, v.vertexProb)

        # Apply the transform
        for middleVertexLayer_parentTrellis in range(1, self.layers, 2):
            fromVertexLayer_parentTrellis = middleVertexLayer_parentTrellis - 1
            toVertexLayer_parentTrellis = middleVertexLayer_parentTrellis + 1

            for (middleVertexKey_parentTrellis, middleVertex_parentTrellis) in self.verticesInLayer[
                middleVertexLayer_parentTrellis].items():
                for (incomingEdgeKey_middleVertex_parentTrellis,
                     incomingEdge_middleVertex_parentTrellis) in middleVertex_parentTrellis.incomingEdges.items():
                    for (outgoingEdgeKey_middleVertex_parentTrellis,
                         outgoingEdge_middleVertex_parentTrellis) in middleVertex_parentTrellis.outgoingEdges.items():
                        # u --(p0,x0)--> w --(p1,x1)--> v
                        w = middleVertex_parentTrellis  # For readabilty, not really used
                        u = incomingEdge_middleVertex_parentTrellis.fromVertex
                        v = outgoingEdge_middleVertex_parentTrellis.toVertex
                        x0 = incomingEdge_middleVertex_parentTrellis.edgeLabel
                        x1 = outgoingEdge_middleVertex_parentTrellis.edgeLabel
                        p0 = incomingEdge_middleVertex_parentTrellis.edgeProb
                        p1 = outgoingEdge_middleVertex_parentTrellis.edgeProb

                        newEdgeProbToAdd = p0 * p1
                        minusEdgeLabel = 1 if x0 != x1 else 0

                        if decisionVector is None:
                            newTrellis.addToEdgeProb(u.stateId, u.verticalPosInLayer, u.layer // 2, v.stateId,
                                                     v.verticalPosInLayer, v.layer // 2, minusEdgeLabel,
                                                     newEdgeProbToAdd)
                        else:
                            if minusEdgeLabel != decisionVector[u.layer // 2]:
                                continue

                            plusEdgeLabel = x1
                            newTrellis.addToEdgeProb(u.stateId, u.verticalPosInLayer, u.layer // 2, v.stateId,
                                                     v.verticalPosInLayer, v.layer // 2, plusEdgeLabel,
                                                     newEdgeProbToAdd)

        return newTrellis

    def calcMarginalizedProbabilities(self, normalize=True):
        assert (len(self) == 1)

        marginalizedProbs = np.zeros(2)

        if normalize == True:
            s = 0.0
            for (vkey, v) in self.verticesInLayer[0].items():
                for (outgoingEdgeKey, outgoingEdge) in v.outgoingEdges.items():
                    s += v.vertexProb * outgoingEdge.edgeProb * outgoingEdge.toVertex.vertexProb
        else:
            s = 1.0

        for (vkey, v) in self.verticesInLayer[0].items():
            for (outgoingEdgeKey, outgoingEdge) in v.outgoingEdges.items():
                x = outgoingEdge.edgeLabel
                marginalizedProbs[x] += v.vertexProb * outgoingEdge.edgeProb * outgoingEdge.toVertex.vertexProb / s

        return marginalizedProbs

    def calcNormalizationVector(self):
        normalization = np.zeros(self.length)
        tempProbs = np.zeros(2)

        for i in range(self.length):
            tempProbs[0] = tempProbs[1] = 0.0

            for (vkey, v) in self.verticesInLayer[i].items():
                for (outgoingEdgeKey, outgoingEdge) in v.outgoingEdges.items():
                    x = outgoingEdge.edgeLabel
                    tempProbs[x] += outgoingEdge.edgeProb

            normalization[i] = np.maximum(tempProbs[0], tempProbs[1])

        return normalization

    # normalize according to the above described vector
    def normalize(self, normalization):
        for i in range(self.length):
            t = normalization[i]
            assert (t >= 0)
            if t == 0:
                t = 1

            for (vkey, v) in self.verticesInLayer[i].items():
                for (outgoingEdgeKey, outgoingEdge) in v.outgoingEdges.items():
                    outgoingEdge.edgeProb /= t


def buildTrellis_uniformInput_deletion(receivedWord, codewordLength, deletionProb, trimmedZerosAtEdges,
                                       numberOfOnesToAddAtBothEndsOfGuardbands):
    trellis = BinaryTrellis(codewordLength)

    deletionCount = codewordLength + 2 * numberOfOnesToAddAtBothEndsOfGuardbands - len(receivedWord)

    inputProb = [0.5, 0.5]

    vertex_stateId = 0  # only one state in the input process
    vertex_layer = 0  # start with the first layer

    if numberOfOnesToAddAtBothEndsOfGuardbands > 0:
        assert (trimmedZerosAtEdges == True)

    if numberOfOnesToAddAtBothEndsOfGuardbands > 0:
        # See longer explanation below about collapsing the leftmost layer of an auxiliary trellis into the layer to its right.

        # Commented out code for the case of adding a single 1 at the edge of each guard band
        #
        # # Add the vertex corresponding to "the 1 from the guard band was deleted"
        # vertex_verticalPosInLayer = 0
        # vertexProb = deletionProb 
        # trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)
        #
        # if len(receivedWord) > 0: # Add the vertex corresponding to "the 1 from the guard band was not deleted"
        #     vertex_verticalPosInLayer = 1 # the first y corresponds to the 1 from the guard band, so start with the second y
        #     vertexProb = 1.0 - deletionProb 
        #     trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

        for i in range(1 + min(numberOfOnesToAddAtBothEndsOfGuardbands, len(receivedWord))):
            vertex_verticalPosInLayer = i
            vertexProb = scipy.special.comb(numberOfOnesToAddAtBothEndsOfGuardbands, i, exact=True) * (
                        (1.0 - deletionProb) ** i) * (deletionProb ** (numberOfOnesToAddAtBothEndsOfGuardbands - i))
            trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

    else:
        # add a single vertex
        vertex_verticalPosInLayer = 0
        vertexProb = 1.0
        trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

    # Now add the vertices in the last layer (either a single vertex or two vertices)
    vertex_layer = codewordLength

    if numberOfOnesToAddAtBothEndsOfGuardbands > 0:
        # See longer explanation below about collapsing the rightmost layer of an auxiliary trellis into the layer to its left.

        # Commented out code for the case of adding a single 1 at the edge of each guard band
        #
        # Add the vertex corresponding to "the 1 from the guard band was deleted"
        # vertex_verticalPosInLayer = len(receivedWord) 
        # vertexProb = deletionProb 
        # trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)
        #
        # if len(receivedWord) > 0: # Add the vertex corresponding to "the 1 from the guard band was not deleted"
        #     vertex_verticalPosInLayer = len(receivedWord) - 1 # the last y corresponds to the 1 from the guard band, so end with the penultimate y
        #     vertexProb = 1.0 - deletionProb
        #     trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

        for i in range(len(receivedWord),
                       len(receivedWord) - min(numberOfOnesToAddAtBothEndsOfGuardbands, len(receivedWord)) - 1, -1):
            vertex_verticalPosInLayer = i
            j = len(receivedWord) - i
            vertexProb = scipy.special.comb(numberOfOnesToAddAtBothEndsOfGuardbands, j, exact=True) * (
                        (1.0 - deletionProb) ** j) * (deletionProb ** (numberOfOnesToAddAtBothEndsOfGuardbands - j))
            trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)
    else:
        # add a single vertex
        vertex_verticalPosInLayer = len(receivedWord)
        vertexProb = 1.0  # not really needed, but for readability
        trellis.setVertexProb(vertex_stateId, vertex_verticalPosInLayer, vertex_layer, vertexProb)

    if trimmedZerosAtEdges == True:
        assert (len(receivedWord) == 0 or (receivedWord[0] == 1 and receivedWord[-1] == 1))

    for l in range(codewordLength):  # add edges from layer l to layer l+1
        # min and max vertical positions in layer l
        if numberOfOnesToAddAtBothEndsOfGuardbands > 0:
            # Explanation for adding a single 1 at the ends of the guard bands. Need to generalize...
            #
            # The simplest way to think of what is going for this case is to write a trellis for the two forced '1' inputs from the guard bands.
            # Thus, the first and last inputs are forced to be '1'. Since the first input is known to be 1, we can collapse the first layer into the second.
            # Likewise, since the last input is known to be  '1', we can collapse the last layer into the penultimate layer.
            # Thus, layer l of the new trellis is layer l' = l+1 of the original layer
            #
            # vpos_min = max(0, l + 1 - deletionCount) # Look at the simpler case below, write l' in place of l, and then use l'=l+1
            # vpos_max = min(l + 1, len(receivedWord))

            # The generalization: l' = l + numberOfOnesToAddAtBothEndsOfGuardbands
            vpos_min = max(0,
                           l + numberOfOnesToAddAtBothEndsOfGuardbands - deletionCount)  # Look at the simpler case below, write l' in place of l, and then use l'=l+numberOfOnesToAddAtBothEndsOfGuardbands
            vpos_max = min(l + numberOfOnesToAddAtBothEndsOfGuardbands, len(receivedWord))
        else:
            vpos_min = max(0, l - deletionCount)
            vpos_max = min(l, len(receivedWord))

        for vpos in range(vpos_min, vpos_max + 1):
            if vpos < len(receivedWord):  # then we can have a non-deletion event
                fromVertex_stateId = 0
                fromVertex_verticalPosInLayer = vpos
                fromVertex_layer = l
                toVertex_stateId = 0
                toVertex_verticalPosInLayer = vpos + 1  # a non-deletion
                toVertex_layer = l + 1
                edgeLabel = receivedWord[vpos]
                probToAdd = inputProb[edgeLabel] * (1.0 - deletionProb)

                trellis.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer,
                                      toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel,
                                      probToAdd)
            # if ((numberOfOnesToAddAtBothEndsOfGuardbands > 0 and l + 1 - vpos < deletionCount) or (numberOfOnesToAddAtBothEndsOfGuardbands == 0 and l - vpos < deletionCount)): # then we can have a deletion event
            if l + 1 + numberOfOnesToAddAtBothEndsOfGuardbands - deletionCount <= vpos:  # then we can have a deletion event
                for edgeLabel in range(2):
                    fromVertex_stateId = 0
                    fromVertex_verticalPosInLayer = vpos
                    fromVertex_layer = l
                    toVertex_stateId = 0
                    toVertex_verticalPosInLayer = vpos  # a deletion
                    toVertex_layer = l + 1

                    if trimmedZerosAtEdges == False or edgeLabel == 1 or (vpos > 0 and vpos < len(receivedWord)):
                        probToAdd = inputProb[edgeLabel] * deletionProb
                    else:
                        probToAdd = inputProb[edgeLabel]

                    trellis.addToEdgeProb(fromVertex_stateId, fromVertex_verticalPosInLayer, fromVertex_layer,
                                          toVertex_stateId, toVertex_verticalPosInLayer, toVertex_layer, edgeLabel,
                                          probToAdd)

    return trellis


def deletionChannelSimulation(codeword, p, seed, randomNumberGenerator=None):
    N = len(codeword)
    receivedWord = []

    if randomNumberGenerator is not None:
        assert (seed is None)
    else:
        if seed is None:
            seed = 200  # just pick a number that is not 0
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(seed)

    for i in range(N):
        r = randomNumberGenerator.random()

        if r < p:
            pass
        else:
            receivedWord.append(codeword[i])

    return receivedWord
