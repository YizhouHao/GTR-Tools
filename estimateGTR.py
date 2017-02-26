#! /usr/bin/env python3
'''
Niema Moshiri 2016
ECE 286 Homework 2

Estimate GTR parameters on multiple DNA sequences
'''
import argparse, dendropy
from common import parseFASTA,gtr2matrix,matrix2gtr,L
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# parse arguments
def parseArgs():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=argparse.FileType('r'), help="Tree")
    parser.add_argument('-d', '--seqs', required=True, type=argparse.FileType('r'), help="Sequences (FASTA)")
    parser.add_argument('-i', '--maxit', required=False, type=int, default=None, help="Maximum number of optimization iterations")
    parser.add_argument('-o', '--out', required=True, type=argparse.FileType('w'), help="Output File")
    args = parser.parse_args()
    assert args.maxit is None or args.maxit > 0, "Maximum number of optimization iterations must be positive"
    return args

# dummy negative likelihood function for optimization (minimization)
# x[0] = pi_A, x[1] = pi_C, x[2] = pi_G, x[3] = r_AC, x[4] = r_AG,  x[5] = r_AT, x[6] = r_CG, x[7] = r_CT, x[8] = r_GT
def f(x, tree=None, seqs=None):
    assert tree is not None, "tree cannot be None"
    assert seqs is not None, "seqs cannot be None"
    pi = {'A':x[0], 'C':x[1], 'G':x[2], 'T':(1-x[0]-x[1]-x[2])}
    if sum(pi.values()) > 1:
        return float('inf')
    for n in pi:
        if pi[n] <= 0:
            return float('inf')
    gtr = {'AC':x[3], 'AG':x[4], 'AT':x[5], 'CG':x[6], 'CT':x[7], 'GT':x[8]}
    return -1*L(tree,seqs,pi,gtr2matrix(gtr,pi))

# compute the maximum-likelihood GTR parameters
def MLGTR(tree, seqs, maxit=None):
    x0 = [0]*9 # this will be the estimate for the parameters x
    bounds = [(0,1)]*9 # all 9 parameters are probabilities
    nucFreqs = {'A':0.0,'C':0.0,'G':0.0,'T':0.0}
    for s in seqs.values():
        for c in s:
            nucFreqs[c] += 1.
    tot = sum(nucFreqs.values())
    for n in 'ACGT':
        nucFreqs[n] /= tot
    x0[0] = nucFreqs['A']
    x0[1] = nucFreqs['C']
    x0[2] = nucFreqs['G']
    # default start for R is Jukes-Cantor
    x0[3] = 1./(3*nucFreqs['C'])
    x0[4] = 1./(3*nucFreqs['G'])
    x0[5] = 1./(3*nucFreqs['T'])
    x0[6] = 1./(3*nucFreqs['G'])
    x0[7] = 1./(3*nucFreqs['T'])
    x0[8] = 1./(3*nucFreqs['T'])
    if maxit is None:
        result = minimize(f, x0, bounds=bounds, args=(tree,seqs), method='SLSQP')
    else:
        result = minimize(f, x0, bounds=bounds, args=(tree,seqs), method='SLSQP', options={'maxiter':maxit})
    x = result.x
    pi = {'A':x[0], 'C':x[1], 'G':x[2], 'T':(1-x[0]-x[1]-x[2])}
    gtr = {'AC':x[3], 'AG':x[4], 'AT':x[5], 'CG':x[6], 'CT':x[7], 'GT':x[8]}
    return pi,gtr

# main function
if __name__ == "__main__":
    args = parseArgs()
    tree = dendropy.Tree.get(file=args.tree, schema='newick')
    seqs = parseFASTA(args.seqs)
    pi,R = MLGTR(tree, seqs, maxit=args.maxit)
    args.out.write(' '.join([str(pi[c]) for c in 'ACGT']) + '\n')
    args.out.write(str(R['CT']) + ' ' + str(R['AT']) + ' ' + str(R['GT']) + ' ' + str(R['AC']) + ' ' + str(R['CG']) + ' ' + str(R['AG']) + '\n')
