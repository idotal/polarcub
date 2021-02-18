import BinaryMemorylessDistribution

bmd = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

p = 0.1

bmd.probs.append( [p, 1-p] )
bmd.probs.append( [1-p, p] )

print( bmd.errorProb() )

