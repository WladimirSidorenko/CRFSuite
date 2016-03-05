#!/usr/bin/env python
# -*- mode: python; coding: utf-8; -*-

##################################################################
# Imports
from __future__ import print_function

import crfsuite
import codecs
import os
import sys

##################################################################
# Variables and Constants
ENCODING = "UTF-8"


##################################################################
# Class
# Inherit crfsuite.Trainer to implement message() function, which receives
# progress messages from a training process.
class Trainer(crfsuite.Trainer):
    def message(self, s):
        # Simply output the progress messages to STDOUT.
        sys.stdout.write(s)


##################################################################
# Methods
def instances(fi):
    xseq = crfsuite.ItemSequence()
    yseq = crfsuite.StringList()

    for line in fi:
        line = line.strip('\n')
        if not line:
            # An empty line presents the end of a sequence.
            yield xseq, yseq
            xseq = crfsuite.ItemSequence()
            yseq = crfsuite.StringList()
            continue

        # Split the line with TAB characters.
        fields = line.split('\t')
        # Append attributes to the item.
        item = crfsuite.Item()
        for field in fields[1:]:
            p = field.rfind(':')
            if p == -1:
                # Unweighted (weight=1) attribute.
                item.append(crfsuite.Attribute(field))
            else:
                # Weighted attribute
                item.append(crfsuite.Attribute(field[:p], float(field[p+1:])))

        # Append the item to the item sequence.
        xseq.push_back(item)
        # Append the label to the label sequence.
        yseq.push_back(fields[0])

##################################################################
# Main
if __name__ == '__main__':
    """Train CRF model on the given dataset.

    Args:
    argv (list(str)): command line arguments

    Returns:
    (void):

    """
    import argparse
    parser = argparse.ArgumentParser(description=
                                     "Script for training CRF models.")
    parser.add_argument("--help-params", help="output CRFSuite parameters",
                        action="store_true")
    parser.add_argument("-a", "--algorithm",
                        help="type of graphical model to use",
                        type=str, default="lbfgs", choices=("lbfgs", "l2sgd",
                                                            "ap", "pa",
                                                            "arow"))
    parser.add_argument("-m", "--model",
                        help="model in which to store the file", type=str,
                        default="")
    parser.add_argument("-t", "--type",
                        help="type of graphical model to use",
                        type=str, default="1d",
                        choices=("1d", "tree", "semim"))
    parser.add_argument("-v", "--version",
                        help="output CRFSuite version", action="store_true")
    parser.add_argument("files", help="input files", nargs='*',
                        type=argparse.FileType('r'),
                        default=[sys.stdin])
    args = parser.parse_args()

    # This demonstrates how to obtain the version string of CRFsuite.
    if args.version:
        print(crfsuite.version())
        sys.exit(0)
    # Create a Trainer object.
    trainer = Trainer()
    # Use L2-regularized SGD and 1st-order dyad features.
    if not trainer.select(str(args.algorithm), str(args.type)):
        raise Exception("Could not initialize trainer.")

    if args.help_params:
        for name in trainer.params():
            print(' '.join([name, trainer.get(name), trainer.help(name)]))
    else:
        if args.model:
            mdir = os.path.dirname(args.model)
            if mdir == "":
                pass
            elif os.path.exists(mdir):
                if not os.path.isdir(mdir) or not os.access(mdir, os.R_OK):
                    print("Can't write to directory '{:s}'.".format(mdir),
                          file=sys.stderr)
            else:
                os.makedirs(mdir)

        # Set the coefficient for L2 regularization to 0.1
        # trainer.set('c2', '0.1')

        # read training instances
        for ifile in args.files:
            for xseq, yseq in instances(ifile):
                trainer.append(xseq, yseq, 0)
        # print("Dataset read...", file=sys.stderr)

        # Start training; the training process will invoke trainer.message()
        # to report the progress.
        trainer.train(str(args.model), -1)
        # print("Model trained...", file=sys.stderr)
