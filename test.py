#! /usr/bin/env python3

import BinaryMemorylessDistribution

bmd = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

p = 0.1

bsc = BinaryMemorylessDistribution.makeBSC(p)
minus = bsc.minusTransform()
plus = bsc.plusTransform()

# print( minus.probs )
# print( plus.probs )
print( bsc.conditionalEntropy() )
print( minus.conditionalEntropy() )
print( plus.conditionalEntropy() )




