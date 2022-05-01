# We will implement the heap as an array, so that the conversion to C++ is easier.
# Note that we could have just used pointers to left child, right child, and parent.

class LinkedListHeap:

    def __init__(self, keyList=None, dataList=None):
        # check if this code duplication can be avoided
        self._head = None
        self._tail = None
        self._heapArray = []

        assert (keyList == None and dataList == None) or len(dataList) == len(keyList)

        # no data given
        if dataList == None and keyList == None:
            return

        # TODO: This is O(n log n). We could use makeheap and get down to O(n)
        for key, data in zip(keyList, dataList):
            self.insertAtTail(key, data)

    def __str__(self):
        s = ""
        # s = "elementsCount = " + str(self.elementsCount)
        # s = s + " _head = " + str(self.head)
        # s = s + " _tail = " + str(self.tail) + " "

        s = s + "head -> "

        element = self._head
        while (element != None):
            # s = s + "(mem = " + str(element) + " indexInArray = " + str(element.indexInArray) + ", key = " + str(element.key) + ", data = " + str(element.data) + ")"
            s = s + "(indexInArray = " + str(element.indexInArray) + ", key = " + str(element.key) + ", data = " + str(
                element.data) + ")"

            if element != self._tail:
                s = s + " <-> "

            element = element.rightElementInList

        s = s + " <- tail"
        return s

    def numberOfElements(self):
        return len(self._heapArray)

    def getHeapMin(self):
        return self._heapArray[0]

    def extractHeapMin(self):
        element = self._heapArray[0]

        # remove from linked list
        left = element.leftElementInList
        right = element.rightElementInList

        if left != None:
            left.rightElementInList = right

        if right != None:
            right.leftElementInList = left

        # update heap
        newElement = self._heapArray.pop()
        self._heapArray[0] = newElement
        newElement.indexInArray = 0
        self._propagateDown(newElement)

        return element

    def updateKey(self, linkedListHeapElement, newKey):
        oldKey = linkedListHeapElement.key
        linkedListHeapElement.key = newKey

        if oldKey < newKey:
            self._propagateDown(linkedListHeapElement)
        elif oldKey > newKey:
            self._propagateUp(linkedListHeapElement)

    def insertAtTail(self, key, data):
        element = LinkedListHeapElement()

        element.key = key
        element.data = data
        element.indexInArray = len(self._heapArray)

        # to our left is the old tail
        element.leftElementInList = self._tail

        # now we are the new tail
        self._tail = element

        # are we also the new head?
        if self._head == None:
            self._head = element

        # if there is an element to the left of us, update its right element to point to us
        if element.leftElementInList != None:
            element.leftElementInList.rightElementInList = element

        # add ourselves to the array
        self._heapArray.append(element)

        # update the heap
        self._propagateUp(element)

    def returnData(self):
        data = []

        element = self._head

        while element != None:
            data.append(element.data)
            element = element.rightElementInList

        return data

    def _swap(self, linkedListHeapElement1, linkedListHeapElement2):
        tempIndex = linkedListHeapElement1.indexInArray
        linkedListHeapElement1.indexInArray = linkedListHeapElement2.indexInArray
        linkedListHeapElement2.indexInArray = tempIndex

        self._heapArray[linkedListHeapElement1.indexInArray] = linkedListHeapElement1
        self._heapArray[linkedListHeapElement2.indexInArray] = linkedListHeapElement2

    def _propagateUp(self, linkedListHeapElement):
        element = linkedListHeapElement

        while True:
            parentIndex = indexOfParentInArray(element.indexInArray)

            if (parentIndex == -1):  # we are the root
                break

            parent = self._heapArray[parentIndex]

            if (parent.key < element.key):
                break

            self._swap(element, parent)

    def _propagateDown(self, linkedListHeapElement):
        element = linkedListHeapElement

        while True:
            leftChildIndex = indexOfLeftChildInArray(element.indexInArray)
            rightChildIndex = indexOfRightChildInArray(element.indexInArray)

            minKey = element.key
            minChild = None

            if leftChildIndex < len(self._heapArray):
                leftChild = self._heapArray[leftChildIndex]
                if leftChild.key < minKey:
                    minKey = leftChild.key
                    minChild = leftChild

            if rightChildIndex < len(self._heapArray):
                rightChild = self._heapArray[rightChildIndex]
                if rightChild.key < minKey:
                    minKey = rightChild.key
                    minChild = rightChild

            if minChild == None:  # we are the root
                break

            self._swap(element, minChild)


class LinkedListHeapElement:
    def __init__(self):
        # set the following three fields to a value that will stand out if not initialized properly later
        self.indexInArray = None
        self.key = None
        self.data = None

        # set to a default value that makes sense
        self.leftElementInList = None  #
        self.rightElementInList = None


def indexOfLeftChildInArray(i):
    return 2 * (i + 1) - 1


def indexOfRightChildInArray(i):
    return 2 * (i + 1) + 1 - 1


def indexOfParentInArray(i):
    return (i + 1) // 2 - 1


def updateData_degrade(dataLeft, dataCenter, dataRight):  # merge dataLeft and dataCenter
    dataLeft[0] += dataCenter[0]
    dataLeft[1] += dataCenter[1]
