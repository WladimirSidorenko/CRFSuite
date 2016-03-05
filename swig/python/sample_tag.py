#!/usr/bin/env python
# -*- mode: python; coding: utf-8; -*-

##################################################################
# Imports
from crfsuite import Attribute, Item, ItemSequence, Tagger
import sys


##################################################################
# Constants
LINCHAIN = "1d"
SEMIM = "semim"
TREE = "tree"
MTYPE2INT = {LINCHAIN: 1, TREE: 3, SEMIM: 4}


##################################################################
# Methods
def instances(fi):
    """Iterate over instances in the provided file.

    Args:
    fi (FileInput): input file stream

    Yields:
    crfsuite instances

    """
    item_seen = False
    xseq = ItemSequence()
    for line in fi:
        line = line.strip('\n')
        if not line:
            # An empty line presents an end of a sequence.
            if item_seen:
                yield xseq
            xseq = ItemSequence()
            item_seen = False
            continue

        # Split the line on TAB characters.
        fields = line.split('\t')
        item_seen = True
        item = Item()
        for field in fields[1:]:
            p = field.rfind(':')
            if p == -1:
                # Unweighted (weight=1) attribute.
                item.append(Attribute(field))
            else:
                # Weighted attribute
                item.append(Attribute(field[:p], float(field[p+1:])))
        # Append the item to the item sequence.
        xseq.append(item)
    if item_seen:
        yield xseq

##################################################################
# Main
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=
                                     "Script for testing CRF models.")
    parser.add_argument("-m", "--model",
                        help="model in which to store the file", type=str,
                        default="")
    parser.add_argument("-t", "--type",
                        help="type of graphical model to use",
                        type=str, default=LINCHAIN,
                        choices=(LINCHAIN, TREE, SEMIM))
    parser.add_argument("files", help="input files", nargs='*',
                        type=argparse.FileType('r'),
                        default=[sys.stdin])
    args = parser.parse_args()

    # Create a tagger object.
    tagger = Tagger()

    # Load the model to the tagger.
    # the second argumend specifies the model type (1 - 1d, 3 - tree, 4 -
    # semim)
    if not tagger.open(args.model, MTYPE2INT[args.type]):
        raise RuntimeError("Could not load model file.")

    for ifile in args.files:
        for xseq in instances(ifile):
            # Tag the sequence.
            tagger.set(xseq)
            # Obtain the label sequence predicted by the tagger.
            yseq = tagger.viterbi()
            # Output the probability of the predicted label sequence.
            # print tagger.probability(yseq)
            for t, y in enumerate(yseq):
                # Output the predicted labels with their marginal
                # probabilities.
                if args.type == SEMIM:
                    print '%s' % (y)
                else:
                    print '%s:%f' % (y, tagger.marginal(y, t))
            print
