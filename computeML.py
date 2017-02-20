#! /usr/bin/env python3
'''
Niema Moshiri 2016
ECE 286 Homework 2

Compute the likelihood of a tree given GTR parameters.
'''
import argparse, dendropy
from common import parseFASTA,gtr2matrix,normalizeGTR,L

# parse arguments
def parseArgs():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree',  required=True, type=argparse.FileType('r'), help="Tree")
    parser.add_argument('-p', '--gtrparams',  required=True, type=argparse.FileType('r'), help="GTR Parameters")
    parser.add_argument('-d', '--seqs',  required=True, type=argparse.FileType('r'), default=None, help="Sequences (FASTA)")
    parser.add_argument('-o', '--out',  required=True, type=argparse.FileType('w'), help="Output File")
    args = parser.parse_args()
    pi, gtr = [[float(i) for i in line.split()] for line in args.gtrparams]
    args.pi = {'A':pi[0], 'C':pi[1], 'G':pi[2], 'T':pi[3]}
    args.gtr = {'CT':gtr[0], 'AT':gtr[1], 'GT':gtr[2], 'AC':gtr[3], 'CG':gtr[4], 'AG':gtr[5]}
    normalizeGTR(args.gtr)
    return args

# main function
if __name__ == "__main__":
    args = parseArgs()
    tree = dendropy.Tree.get(file=args.tree, schema='newick')
    seqs = parseFASTA(args.seqs)
    args.out.write(str(L(tree, seqs, args.pi, gtr2matrix(args.gtr,args.pi))) + '\n')
