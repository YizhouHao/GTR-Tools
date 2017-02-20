#! /usr/bin/env python3
'''
Niema Moshiri 2016
ECE 286 Homework 2

Estimate GTR parameters from pairs of DNA sequences
'''
import argparse
from common import parseFASTA,matrix2gtr,normalizeGTR
from math import log
from numpy import matrix
from scipy.linalg import logm
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# parse arguments
def parseArgs():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--branchlength',  required=True, type=float, help="Branch Length")
    parser.add_argument('-d', '--seqs',  required=True, type=argparse.FileType('r'), default=None, help="Pair of Sequences (FASTA)")
    parser.add_argument('-o', '--out',  required=True, type=argparse.FileType('w'), help="Output File")
    return parser.parse_args()

# negative likelihood of seeing two sequences s and t
# x[0] = pi_A, x[1] = pi_C, x[2] = pi_G, x[3] = P_AC, x[4] = P_AG, x[5] = P_AT, x[6] = P_CG, x[7] = P_CT, x[8] = P_GT
def L(x, s=None, t=None):
    assert s is not None, "Must specify a sequence for s"
    assert t is not None, "Must specify a sequence for t"
    assert len(s) == len(t), "The two strings must have equal lengths"
    out = 0
    for i in range(len(s)):
        if s[i] == 'A':
            out += log(x[0])
            if t[i] == 'A':
                out += log(1-x[3]-x[4]-x[5])
            elif t[i] == 'C':
                out += log(x[3])
            elif t[i] == 'G':
                out += log(x[4])
            else:
                out += log(x[5])
        elif s[i] == 'C':
            out += log(x[1])
            if t[i] == 'A':
                out += log(x[3])
            elif t[i] == 'C':
                out += log(1-x[3]-x[6]-x[7])
            elif t[i] == 'G':
                out += log(x[6])
            else:
                out += log(x[7])
        elif s[i] == 'G':
            out += log(x[2])
            if t[i] == 'A':
                out += log(x[4])
            elif t[i] == 'C':
                out += log(x[6])
            elif t[i] == 'G':
                out += log(1-x[4]-x[6]-x[8])
            else:
                out += log(x[8])
        else:
            out += log(1-x[0]-x[1]-x[2])
            if t[i] == 'A':
                out += log(x[5])
            elif t[i] == 'C':
                out += log(x[7])
            elif t[i] == 'G':
                out += log(x[8])
            else:
                out += log(1-x[5]-x[7]-x[8])
    return -1*out # negate because we will minimize

# constraints for optimization
def conPi(x, s=None, t=None):
    if x[0] + x[1] + x[2] < 1:
        return 0
    return 1
def conP(x, s=None, t=None):
    if x[3] + x[4] + x[5] < 1 and x[3] + x[6] + x[7] < 1 and x[4] + x[6] + [8] < 1 and x[5] + x[7] + x[8] < 1:
        return 0
    return 1

# estimate maximum likelihood GTR parameters from pair of sequences and branch length
def MLGTR(s, t, bl):
    assert len(s) == len(t), "The two strings must have equal lengths"
    x0 = [None]*9 # this will be the estimate for the parameters x
    bounds = [(0,1)]*9 # all 9 parameters are probabilities
    constraints = [{'type':'eq','fun':conPi}, {'type':'eq','fun':conP}]
    nucFreqs = {'A':0.0,'C':0.0,'G':0.0,'T':0.0}
    for i in range(len(s)):
        nucFreqs[s[i]] += 1
        nucFreqs[t[i]] += 1
    for c in 'ACGT':
        nucFreqs[c] /= (2*len(s))
    x0[0] = nucFreqs['A']
    x0[1] = nucFreqs['C']
    x0[2] = nucFreqs['G']
    transitionFreqs = {'AA':0.0, 'AC':0.0, 'AG':0.0, 'AT':0.0, 'CC':0.0, 'CG':0.0, 'CT':0.0, 'GG':0.0, 'GT':0.0, 'TT':0.0}
    for i in range(len(s)):
        transitionFreqs[''.join(sorted([s[i],t[i]]))] += 1
    transitionAsum = transitionFreqs['AA'] + transitionFreqs['AC'] + transitionFreqs['AG'] + transitionFreqs['AT']
    transitionCsum = transitionFreqs['AC'] + transitionFreqs['CC'] + transitionFreqs['CG'] + transitionFreqs['CT']
    transitionGsum = transitionFreqs['AG'] + transitionFreqs['CG'] + transitionFreqs['GG'] + transitionFreqs['GT']
    transitionTsum = transitionFreqs['AT'] + transitionFreqs['CT'] + transitionFreqs['GT'] + transitionFreqs['TT']
    x0[3] = transitionFreqs['AC'] / transitionAsum
    x0[4] = transitionFreqs['AG'] / transitionAsum
    x0[5] = transitionFreqs['AT'] / transitionAsum
    x0[6] = transitionFreqs['CG'] / transitionCsum
    x0[7] = transitionFreqs['CT'] / transitionCsum
    x0[8] = transitionFreqs['GT'] / transitionGsum
    x = minimize(L, x0, bounds=bounds, constraints=constraints, args=(s,t), method='SLSQP').x
    pi = {'A':x[0], 'C':x[1], 'G':x[2], 'T':(1-x[0]-x[1]-x[2])}
    P = {'AC':x[3], 'AG':x[4], 'AT':x[5], 'CG':x[6], 'CT':x[7], 'GT':x[8]}
    P = [[None]*4 for _ in range(4)]
    P[0][0] = (1-x[3]-x[4]-x[5])
    P[0][1] = x[3]
    P[0][2] = x[4]
    P[0][3] = x[5]
    P[1][0] = x[3]
    P[1][1] = (1-x[3]-x[6]-x[7])
    P[1][2] = x[6]
    P[1][3] = x[7]
    P[2][0] = x[4]
    P[2][1] = x[6]
    P[2][2] = (1-x[4]-x[6]-x[8])
    P[2][3] = x[8]
    P[3][0] = x[5]
    P[3][1] = x[7]
    P[3][2] = x[8]
    P[3][3] = (1-x[5]-x[7]-x[8])
    Pmat = matrix(P)
    Rmat = logm(Pmat)/bl
    R = matrix2gtr(Rmat)
    normalizeGTR(R)
    return pi,R

# main function
if __name__ == "__main__":
    args = parseArgs()
    seqs = parseFASTA(args.seqs)
    if len(seqs) != 2:
        print("ERROR: Must specify exactly 2 sequences in -d argument")
        exit(-1)
    IDs = sorted(list(seqs.keys()))
    pi,R = MLGTR(seqs[IDs[0]], seqs[IDs[1]], args.branchlength)
    args.out.write(' '.join([str(pi[c]) for c in 'ACGT']) + '\n')
    args.out.write(str(R['CT']) + ' ' + str(R['AT']) + ' ' + str(R['GT']) + ' ' + str(R['AC']) + ' ' + str(R['CG']) + ' ' + str(R['AG']) + '\n')
