import BinaryMemorylessDistribution

bmd = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

p = 0.1

# bmd.probs.append( [0.5*p, 0.5*(1-p)] )
# bmd.probs.append( [0.5*(1-p), 0.5*p] )

bmd.probs.append( [p, (1-p)] )
bmd.probs.append( [(1-p), p] )

print( bmd.errorProb() )
print( bmd.bhattacharyya() )
print( bmd.totalVariationDistance() )
bmd.normalize()
print( bmd.errorProb() )
print( bmd.bhattacharyya() )
print( bmd.totalVariationDistance() )

