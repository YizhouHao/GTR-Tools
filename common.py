#! /usr/bin/env python3
'''
Common functions to use across the alignment tools

Niema Moshiri 2017
'''
from decimal import *
from math import log
from random import uniform
from numpy import matrix
from scipy.linalg import expm

# parse FASTA file (single-line or multi-line)
def parseFASTA(f):
    seqs = {}
    currID = None
    currSeq = ''
    for line in f:
        l = line.strip()
        if len(l) == 0:
            continue
        if l[0] == '>':
            if currID is not None:
                seqs[currID] = currSeq
                currSeq = ''
            currID = l[1:]
        else:
            currSeq += l
    seqs[currID] = currSeq
    return seqs

# roll a weighted die (keys = faces, values = probabilities)
def roll(die):
    faces = sorted(die.keys())
    probs = [die[key] for key in faces]
    cdf = [probs[0]]
    while len(cdf) < len(probs):
        cdf.append(cdf[-1] + probs[len(cdf)])
    num = uniform(0, 1)
    index = 0
    while cdf[index] < num:
        index += 1
    return faces[index]

# generate a random string of length k given stationary probabilities pi
def genString(k, pi):
    return ''.join([roll(pi) for _ in range(k)])

# convert GTR parameter dictionary and stationary vector to numpy rate matrix (0 = A, 1 = C, 2 = G, 3 = T)
def gtr2matrix(gtr, pi):
    nuc = ['A','C','G','T']
    R = [[None,None,None,None] for _ in range(4)]
    # row A
    R[0][1] = gtr['AC']*pi['C']
    R[0][2] = gtr['AG']*pi['G']
    R[0][3] = gtr['AT']*pi['T']
    R[0][0] = -1*(R[0][1]+R[0][2]+R[0][3])
    # row C
    R[1][0] = gtr['AC']*pi['A']
    R[1][2] = gtr['CG']*pi['G']
    R[1][3] = gtr['CT']*pi['T']
    R[1][1] = -1*(R[1][0]+R[1][2]+R[1][3])
    # row G
    R[2][0] = gtr['AG']*pi['A']
    R[2][1] = gtr['CG']*pi['C']
    R[2][3] = gtr['GT']*pi['T']
    R[2][2] = -1*(R[2][0]+R[2][1]+R[2][3])
    # row T
    R[3][0] = gtr['AT']*pi['A']
    R[3][1] = gtr['CT']*pi['C']
    R[3][2] = gtr['GT']*pi['G']
    R[3][3] = -1*(R[3][0]+R[3][1]+R[3][2])
    # normalize such that sum(vi*pi) = 1
    norm = -1*sum([R[i][i]*pi[nuc[i]] for i in range(4)])
    for i in range(4):
        for j in range(4):
            R[i][j] /= norm
    return matrix(R)

# convert numpy matrix to GTR parameter dictionary (0 = A, 1 = C, 2 = G, 3 = T)
def matrix2gtr(R, pi):
    return {'AC':float(R[0][1])/pi['C'], 'AG':float(R[0][2])/pi['G'], 'AT':float(R[0][3])/pi['T'], 'CG':float(R[1][2])/pi['G'], 'CT':float(R[1][3])/pi['T'], 'GT':float(R[2][3])/pi['T']}

# compute maximum likelihood of tree given sequence data and GTR parameters
def L(tree, seqs, pi, R):
    nuc = ['A','C','G','T']
    for seq in seqs.values():
        k = len(seq)
        break
    # Felsenstein pruning algorithm
    for u in tree.postorder_node_iter():
        if u.is_leaf():
            s = seqs[u.taxon.label]
            u.L = [{'A':Decimal(0),'C':Decimal(0),'G':Decimal(0),'T':Decimal(0)} for i in range(k)]
            for i in range(k):
                u.L[i][s[i]] = Decimal(1)
        else:
            u.L = [{'A':Decimal(1),'C':Decimal(1),'G':Decimal(1),'T':Decimal(1)} for i in range(k)]
            for c in u.child_node_iter():
                Pmat = expm(R*c.edge_length)
                P = {'A':{}, 'C':{}, 'G':{}, 'T':{}}
                for i in range(4):
                    for j in range(4):
                        P[nuc[i]][nuc[j]] = Decimal(Pmat[i,j])
                for i in range(k):
                    for s in nuc:
                        u.L[i][s] *= sum([P[s][x] * c.L[i][x] for x in nuc])
    return float(sum([sum([Decimal(pi[s])*tree.seed_node.L[i][s] for s in nuc]).ln() for i in range(k)]))
