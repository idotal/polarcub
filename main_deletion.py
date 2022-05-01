#! /usr/bin/env python3
import argparse
import random

import numpy as np

import BinaryPolarEncoderDecoder
import Guardbands
from ScalarDistributions import BinaryMemorylessDistribution
from VectorDistributions import BinaryTrellis
from VectorDistributions import CollectionOfBinaryTrellises


# for profiling


def make_xVectorDistribution_deletion_uniform(length):
    def make_xVectorDistribuiton():
        xDistribution = BinaryMemorylessDistribution.BinaryMemorylessDistribution()

        xDistribution.probs.append( [0.5,0.5] )
        
        xVectorDistribution = xDistribution.makeBinaryMemorylessVectorDistribution(length, None)
        return xVectorDistribution

    return make_xVectorDistribuiton


def make_codeword_addDeletionGuardBands(xi, n, n0, numberOfOnesToAddAtBothEndsOfGuardbands):
    def make_codeword(encodedVector):
        return Guardbands.addDeletionGuardBands(encodedVector, n, n0, xi, numberOfOnesToAddAtBothEndsOfGuardbands)

    return make_codeword


def make_simulateChannel_deletion(p, seed=None):
    randomNumberGenerator = random.Random()

    if seed is not None:
        randomNumberGenerator.seed(seed)

    def simulateChannel(codeword):
        return BinaryTrellis.deletionChannelSimulation(codeword, p, seed=None,
                                                       randomNumberGenerator=randomNumberGenerator)

    return simulateChannel


def make_xyVectorDistribution_deletion(deletionProbability, xi, n, n0, numberOfOnesToAddAtBothEndsOfGuardbands):
    codewordLength = 2 ** n

    def make_xyVectorDistribution(receivedWord, verbosity=0):
        return CollectionOfBinaryTrellises.buildCollectionOfBinaryTrellises_uniformInput_deletion(receivedWord,
                                                                                                  deletionProbability,
                                                                                                  xi, n, n0,
                                                                                                  numberOfOnesToAddAtBothEndsOfGuardbands,
                                                                                                  verbosity)

    return make_xyVectorDistribution


def main():
    ap = argparse.ArgumentParser(description="polar encoder/decoder for the deletion channel")

    ap.add_argument("-pd", "--deletion-probability", type=float, default=0.1,
                    help="The deletion probability of the channel. Default is 0.1.")
    ap.add_argument("-pe", "--error-probability", type=float, default=0.1,
                    help="Upper bound on error probability, when constructing the frozen set. Default is 0.1.")
    ap.add_argument("-n", "--n", type=int, required=True, help="The total number of polarization steps.")
    ap.add_argument("-n0", "--n0", type=int,
                    help="The number of slow (trellis) polarization steps. The default is n//3.")
    ap.add_argument("-g", "--genie-simulations", nargs="?", const=8000, type=int,
                    help="Perform genie encoding/decoding trials to find the frozen set. Default is 8000.")
    ap.add_argument("-e", "--encoding-decoding-simulations", type=int, nargs="?", const=8000,
                    help="Perform encoding/decoding trials to test the code. Default is 8000.")
    ap.add_argument("--xi", type=float, default=0.1,
                    help="The paramter through which the length of the guard bands is defined. Default is 0.1.")
    ap.add_argument("--ones-added-at-end-of-guard-band", type=int, default=0,
                    help="Add a sequence of ones at the start and end of a guard band (and also at the start and end of the codeword). Defult is 0 (all-zero guardbands).")
    ap.add_argument("-f", "--frozen-bits-file", help="Read/write the frozen bits from/to a file.")
    ap.add_argument("-cs", "--channel-seed", type=int, default=100,
                    help="Seed for simulating the channel. Default is 100.")
    ap.add_argument("-crs", "--common-randomness-seed", type=int, default=200,
                    help="Seed for common randomness between encoder and decoder. Default is 200.")
    ap.add_argument("-gs", "--genie-seed", type=int, default=300,
                    help="Seed for letting the genie pick a different common randomness for each encoding/decoding run. Default is 300.")
    ap.add_argument("-is", "--information-seed", type=int, default=400,
                    help="Seed for picking the random bits to encode. Default is 400.")

    args = vars(ap.parse_args())

    # print(args)

    deletionProbability = args['deletion_probability']
    numberOfGenieTrials = args['genie_simulations']
    numberOfEncodingDecodingTrials = args['encoding_decoding_simulations']
    n = args['n']
    N = 2 ** n

    n0 = args['n0'] if args['n0'] is not None else n // 3

    numberOfOnesToAddAtBothEndsOfGuardbands = args['ones_added_at_end_of_guard_band']
    filename = args['frozen_bits_file']

    # guardband parameters
    xi = args['xi']

    print("n =", n, ", n0 =", n0, ", xi =", xi, ", deletion probability =", deletionProbability)
    if filename is not None:
        print("frozen bits file = ", filename)

    upperBoundOnErrorProbability = args['error_probability']

    make_xVectorDistribuiton = make_xVectorDistribution_deletion_uniform(N)
    make_codeword = make_codeword_addDeletionGuardBands(xi, n, n0, numberOfOnesToAddAtBothEndsOfGuardbands)
    print("codeword length = ", len(make_codeword(np.zeros(2 ** n, dtype=np.int64))))

    simulateChannel = make_simulateChannel_deletion(deletionProbability, seed=args['channel_seed'])
    make_xyVectorDistribution = make_xyVectorDistribution_deletion(deletionProbability, xi, n, n0,
                                                                   numberOfOnesToAddAtBothEndsOfGuardbands)

    if numberOfGenieTrials is not None:
        print("performing", numberOfGenieTrials, "genie encoding/decoding trials to find a frozen set with WER at most",
              upperBoundOnErrorProbability)

        trustXYProbs = False if n > n0 else True

        frozenSet = BinaryPolarEncoderDecoder.genieEncodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword,
                                                                    simulateChannel, make_xyVectorDistribution,
                                                                    numberOfGenieTrials, upperBoundOnErrorProbability,
                                                                    genieSeed=args['genie_seed'],
                                                                    trustXYProbs=trustXYProbs, filename=filename)

    if numberOfEncodingDecodingTrials is not None:
        verbosity = 0

        if filename is not None:
            frozenSet = readFrozenSetFromFile(filename)

        print("performing", numberOfEncodingDecodingTrials, "encoding/decoding trials")
        BinaryPolarEncoderDecoder.encodeDecodeSimulation(N, make_xVectorDistribuiton, make_codeword, simulateChannel,
                                                   make_xyVectorDistribution, numberOfEncodingDecodingTrials, frozenSet,
                                                   commonRandomnessSeed=args['common_randomness_seed'],
                                                   randomInformationSeed=args['information_seed'], verbosity=verbosity)


def readFrozenSetFromFile(filename):
    frozenSet = set()
    f = open(filename, "r")

    for l in f:
        if l[0] == "*":
            continue
        else:
            frozenSet.add(int(l))
    f.close()
    return frozenSet


main()
# cProfile.run('main()')
