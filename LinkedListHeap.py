# We will implement the heap as an array, so that the conversion to c++ is faster.
# Note that we could have just used pointers to left child, right child, and parent.

class LinkedListHeap:
    def __init__(self):
        self.elementsCount = 0
        self.head = None
        self.tail = None
        self.heapArray = []

    def updateKey(self, linkedListHeapElement, newKey):
        oldKey = linkedListHeapElement.key
        linkedListHeapElement.key = newKey

        if oldKey < newKey:
            self.propagateDown(linkedListHeapElement)
        elif oldKey > newKey:
            self.propagateUp(linkedListHeapElement)

    def insertAtTail(self, key, data):
        element = LinkedListHeapElement()

        element.key = key
        element.data = data
        element.indexInArray = elementsCount
        self.elementsCount += 1

        # to our left is the old tail
        element.leftElementInList = self.tail

        # now we are the new tail
        self.tail = element

        # are we also the new head?
        if self.head == None:
            self.head = element

        # if there is an element to the left of us, update its right element to point to us
        if element.leftElementInList != None:
           element.leftElementInList.rightElementInList = element

        # add ourselves to the array
        self.heapArray.append(element)

    def swap(self, linkedListHeapElement1, linkedListHeapElement2):
        tempIndex = linkedListHeapElement1.indexInArray
        linkedListHeapElement1.indexInArray = linkedListHeapElement2.indexInArray
        linkedListHeapElement2.indexInArray = tempIndex

        self.heapArray[linkedListHeapElement1.indexInArray] = linkedListHeapElement1
        self.heapArray[linkedListHeapElement2.indexInArray] = linkedListHeapElement2

    def propagateUp(self, linkedListHeapElement):
        element = linkedListHeapElement

        while true:
            parentIndex = indexOfParentInArray(element.indexInArray)

            if ( parentIndex == -1 ): # we are the root
                break
            
            parent = self.heapArray[parent]

            if ( parent.key < element.key ):
                break

            self.swap(element,  parent)

    def propagateDown(self,linkedListHeapElement):
        element = linkedListHeapElement

        while true:
            leftChildIndex = indexOfLeftChildInArray(element.indexInArray)
            rightChildIndex = indexOfRightChildInArray(element.indexInArray)

            minKey = element.key
            minChild = None

            if leftChildIndex < elementsCount:
                leftChild = self.heapArray[leftChildIndex]
                if leftChild.key < minKey:
                    minKey = lefChild.key
                    minChild = leftChild

            if rightChildIndex < elementsCount:
                rightChild = self.heapArray[rightChildIndex]
                if rightChild.key < minKey:
                    minKey = rightChild.key
                    minChild = rightChild

            if minChild == None: # we are the root
                break

            self.swap(element,  minChild)

class LinkedListHeapElement:
    def __init__(self):
        # set the following three fields to a value that will stand out if not initialized properly later
        self.indexInArray = None 
        self.key = None
        self.data = None

        # set to a default value that makes sense
        self.leftElementInList = None # 
        self.rightElementInList = None

def indexOfLeftChildInArray(i):
    return 2*(i+1) - 1

def indexOfRightChildInArray(i):
    return 2*(i+1) + 1 - 1

def indexOfParentInArray(i):
    return (i+1) // 2 - 1
