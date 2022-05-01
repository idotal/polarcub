#!/usr/bin/python
import sys

import numpy as np


# Read Ido Frozen
def readFrozen(filename, count):
    f = open(filename, "r")

    mm = 0
    for l in f:
        if len(l) >= 3 and l[0] == "*" and l[1] == "*" and l[2] == " ":
            M = int(l[22:])
        elif len(l) >= 4 and l[0] == "*" and l[1] == "*" and l[2] == "*" and l[3] == " ":
            q = l[4:].split()
            count[int(q[0])] += float(q[1])
            if (int(q[0]) > mm):
                mm = int(q[0])

    f.close()

    return M, mm + 1


# Allocate max length vector
count = [0, ] * (2 ** 20)
total = 0

# Read and combine
for _arg in sys.argv[1:]:
    M, mm = readFrozen(_arg, count)
    total += M

# Design
biterrd = np.asarray(count[:mm])

# Sort into increasing order and compute cumulative sum
order = np.argsort(biterrd)
SE = biterrd[order] / total
CSE = np.cumsum(SE)

# Find best frozen bits
k = np.sum(CSE < 0.1)
print("N = ", mm, ", K = ", k, " Rate = ", k / mm)
frozenSet = set(order[k:])
numberOfTrials = total

f = open("out", "w")
s = "* Combined" + "\n"
f.write(s)

for i in frozenSet:
    f.write(str(i))
    f.write("\n")

s = "** number of trials = " + str(numberOfTrials) + "\n"
f.write(s)
s = "* (TotalVariation+errorProbability) * (number of trials)" + "\n"
f.write(s)

for i in range(mm):
    s = "*** " + str(i) + " " + str(biterrd[i]) + "\n"
    f.write(s)

f.close()
