#! /usr/bin/env python3
'''
Common functions to use across the alignment tools

Niema Moshiri 2017
'''
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

# normalize GTR parameter dictionary such that R_AG has a rate of 1
def normalizeGTR(gtr):
    ag = gtr['AG']
    for key in gtr:
        gtr[key] /= ag

# convert GTR parameter dictionary and stationary vector to numpy rate matrix (0 = A, 1 = C, 2 = G, 3 = T)
def gtr2matrix(gtr, pi):
    normalizeGTR(gtr)
    x1 = float(gtr['AC'])
    x2 = float(gtr['AG'])
    x3 = float(gtr['AT'])
    x4 = float(gtr['CG'])
    x5 = float(gtr['CT'])
    x6 = float(gtr['GT'])
    p1 = float(pi['A'])
    p2 = float(pi['C'])
    p3 = float(pi['G'])
    p4 = float(pi['T'])
    R = [[None,None,None,None] for _ in range(4)]
    R[0][0] = -1*(x1+x2+x3)
    R[0][1] = x1
    R[0][2] = x2
    R[0][3] = x3
    R[1][0] = (p1*x1)/p2
    R[1][1] = -1*(R[1][0]+x4+x5)
    R[1][2] = x4
    R[1][3] = x5
    R[2][0] = (p1*x2)/p3
    R[2][1] = (p2*x4)/p3
    R[2][2] = -1*(R[2][0]+R[2][1]+x6)
    R[2][3] = x6
    R[3][0] = (p1*x3)/p4
    R[3][1] = (p2*x5)/p4
    R[3][2] = (p3*x6)/p4
    R[3][3] = -1*(R[3][0]+R[3][1]+R[3][2])
    return matrix(R)

# convert numpy matrix to GTR parameter dictionary (0 = A, 1 = C, 2 = G, 3 = T)
def matrix2gtr(R):
    out = {'AC':float(R[0][1]), 'AG':float(R[0][2]), 'AT':float(R[0][3]), 'CG':float(R[1][2]), 'CT':float(R[1][3]), 'GT':float(R[2][3])}
    normalizeGTR(out)
    return out

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
            u.L = [{'A':0.,'C':0.,'G':0.,'T':0.} for i in range(k)]
            for i in range(k):
                u.L[i][s[i]] = 1.
        else:
            u.L = [{'A':1.,'C':1.,'G':1.,'T':1.} for i in range(k)]
            for c in u.child_node_iter():
                Pmat = expm(R*c.edge_length)
                P = {'A':{}, 'C':{}, 'G':{}, 'T':{}}
                for i in range(4):
                    for j in range(4):
                        P[nuc[i]][nuc[j]] = Pmat[i,j]
                for i in range(k):
                    for s in nuc:
                        prob = 0.
                        for x in nuc:
                            prob += (P[s][x]*c.L[i][x])
                        u.L[i][s] *= prob
                    if sum(u.L[i].values()) == 0:
                        return float('-inf')
    out = 0.
    for i in range(k):
        pr = 0.
        for s in nuc:
            pr += (pi[s] * tree.seed_node.L[i][s])
        out += log(pr)
    return out
