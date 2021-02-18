#! /usr/bin/env python3

import BinaryMemorylessDistribution

bmd = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

p = 0.1

bsc = BinaryMemorylessDistribution.makeBSC(p)
minus = bsc.minusTransform()
minusminus = minus.minusTransform()

print( bsc.conditionalEntropy() )
print( minus.conditionalEntropy() )

print( minus.probs )
minus.mergeEquivalentSymbols()
print( minus.probs )

# print( minusminus.conditionalEntropy() )

# print( minusminus.probs )
#
minusminus.mergeEquivalentSymbols()

print( minusminus.probs )



