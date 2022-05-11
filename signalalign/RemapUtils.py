from multiprocessing import Pool
import numpy as np
import numba
from numba.typed import List
from numba import f4
from numba.typed import Dict
import multiprocessing
from multiprocessing import Pool
from functools import partial
import time
import copy
from numba.typed import Dict
from numba import types
import gc
import ruptures as rpt
import py2bit
from Bio.Seq import Seq

TraceThres = 5  # adhoc threshold to restrict change point
TraceThresDiff = 1  # adhoc threshold to restrict change point
penalty_value = 100  # beta

SEGMENT_FROM_BASECALLER = 1
SEGMENT_FROM_CHANGE_POINT = 2
SEGMENT_ALL = 3

DIOGONAL_MOVE = 1
HORIZONTAL_MOVE = 2
SKIP_MOVE_BASE = 10

BONUS_FOR_EACH_SEGMENT = 1
IndelPenalty = 15
IndexExtentionpenalty = 2
MAXDEL_SIZE = 9
LOW_THRES_FOR_DEFULT_TRANS_PROP = 0.6
count = 0

import numpy as np

@numba.jit(nopython=True)
def toIndex(possiblemove, segment_type):
    nl = []
    nl.append(0)
    idx = 1
    for v in possiblemove:

        if segment_type == SEGMENT_ALL:
            if v != 0:
                nl.append(idx)
        else:
            if v == segment_type:
                nl.append(idx)

        idx += 1
    nl = np.array(nl, dtype='int16')
    return nl

def toMap(possiblemove):

    posmap = {}
    idx = 0
    idxbase = 0
    for v in possiblemove:

        if v == SEGMENT_FROM_BASECALLER:
            posmap[idxbase] = idx
            idxbase +=1
        if v > 0:
            idx += 1

    return posmap


def toTraceMap(possiblemove):

    posmap = {}
    idx = 0
    idxbase = 0
    for v in possiblemove:

        if v == SEGMENT_FROM_BASECALLER:
            posmap[idxbase] = idx
            idxbase +=1
        idx += 1

    return posmap

import py2bit
def getRefG(chr, s, e,ref2b):

    tb = py2bit.open(ref2b)
    ret = tb.sequence(chr, s, e)
    tb.close()
    return ret

def revcon(conpactTrace):

    compactTraceRev = np.zeros((len(conpactTrace), 4), dtype='float32')
    #
    for n in range(len(conpactTrace)):

        idx = len(conpactTrace)-1-n
        t = conpactTrace[idx]
        compactTraceRev[n][0] = t[3]
        compactTraceRev[n][1] = t[2]
        compactTraceRev[n][2] = t[1]
        compactTraceRev[n][3] = t[0]


    return compactTraceRev

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U': 'A'}
def revconSeq(seq):

    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement
####
delmargin = 50
margin = 10
def consumeref(cigarOparator,cigarlen,delmargin):
    ref_consuming_ops = (cigarOparator == 0) or (cigarOparator == 2 and cigarlen < delmargin) \
                        or (cigarOparator == 7) or (cigarOparator == 8)  # M or D or = or X
    return ref_consuming_ops

def consumeread(cigarOparator):
    ref_consuming_ops = (cigarOparator == 0) or (cigarOparator == 1) \
                        or (cigarOparator == 7) or (cigarOparator == 8)  # M or D or = or X
    return ref_consuming_ops

def adjustBoundary(exons,es,ee,distance_margin):

    for exon in exons:

        exon.st_adjusted = getColsest(exon.start,es,distance_margin)
        exon.end_adjusted = getColsest(exon.end, ee,distance_margin)

    return exons

def getColsest(pos,exonStartOrEnds,distance_margin):

    if pos in exonStartOrEnds:
        return pos
    near = searchNear(pos,exonStartOrEnds)
    if abs(near-pos) < distance_margin:
        return near
    return pos

def searchNear(pos,dlist):

    if len(dlist) == 0:
        return 0
    elif len(dlist) == 1:
        return dlist[0]

    dlist = list(map(lambda n: n-pos, dlist))
    idx = np.abs(np.asarray(dlist)).argmin()
    return dlist[idx]

def toCigarStr(cigar):

    cst =""
    for cl,co in cigar:
        cst = cst + str(cl) + toCO(co)

    return cst

def toCO(co):

    if co == 0:
        return "M"
    elif co == 1:
        return "I"
    elif co == 2:
        return "D"
    elif co == 3:
        return "N"

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def getLocalGenome(localgenome,s,e):

    st = localgenome[s:e]
    return st

def getNear(dic,pos,acc=True):

    if pos in dic:
        return dic[pos]
    keys = list(dic.keys())
    key = closest(keys, pos)
    return dic[key]

def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]

def getBase(trace,b4,bb,strand):

    nar = getTrace(trace, b4, bb,strand)
    nlist = [nar[0],nar[1],nar[2],nar[3]]
    midx = nlist.index(max(nlist))
    label = ["A","C","G","T"]
    nuc = label[midx]

    return nuc

def toCompact(trace,possiblemove_idx,frombasecaller_idx):

    intervallen = len(possiblemove_idx)
    tracelen = len(trace)
    compactTrace = np.zeros((intervallen, 5), dtype='float32')

    idx = 0
    lastbaecallidx = 0
    a, c, g, u = 0, 0, 0, 0
    for m in range(0, tracelen):

        compactTrace[idx][4] += 1  # count for length
        if (m > 0) and (m in possiblemove_idx):

            total = a + c + g + u
            if total > 0:
                a = a / total
                c = c / total
                g = g / total
                u = u / total
            compactTrace[idx][0] = a
            compactTrace[idx][1] = c
            compactTrace[idx][2] = g
            compactTrace[idx][3] = u
            idx += 1
            #print(a,c,g,u)
            # invS = (1 - H([a, c, g, u]))
            # compactTraceInvEntropy.append(invS)
            a, c, g, u = 0, 0, 0, 0

        if (m > 0) and (m in frombasecaller_idx):

            numsep = (idx - 1) - lastbaecallidx
            if numsep > 1:
                totalw = 0
                for t in range(lastbaecallidx,idx):
                    totalw += compactTrace[t][4]
                for t in range(lastbaecallidx, idx):
                    compactTrace[t][4] = compactTrace[t][4] / totalw
            else:
                compactTrace[lastbaecallidx][4] = 1.0

            lastbaecallidx = idx

        a_trace = trace[m]
        _a = max(a_trace[0], a_trace[4])
        _c = max(a_trace[1], a_trace[5])
        _g = max(a_trace[2], a_trace[6])
        _u = max(a_trace[3], a_trace[7])

        a = a + _a
        c = c + _c
        g = g + _g
        u = u + _u

    return compactTrace

@numba.jit(nopython=True)
def getStateProb(statesArray, compactTrace, n, m):


    aref = statesArray[n]
    refbase = aref

    a_trace = compactTrace[m]
    scA = a_trace[0]
    scC = a_trace[1]
    scG = a_trace[2]
    scU = a_trace[3]
    weight = a_trace[4]
    if weight >= 1:
        weight = 1

    total = scA + scC + scG + scU

    if refbase == 'A': score = scA
    if refbase == 'C': score = scC
    if refbase == 'G': score = scG
    if refbase == 'T' or refbase == 'U': score = scU
    minusscore = 0
    if total > 0:
        plusscore = score / float(total)
        minusscore = (float(total)-score) / float(total)

        score = plusscore - minusscore

        # no minus score for AG mismatch
        if refbase == 'A' and scG > scA:
            score = scA / float(total)
        # no minus score for GA mismatch
        if refbase == 'G' and scA > scG:
            score = scG / float(total)
        # no minus score for CU mismatch
        if refbase == 'T' or refbase == 'U' and scC > scU:
            score = scU / float(total)
        # no minus score for UC mismatch
        if refbase == 'C' and scU > scC:
            score = scC / float(total)
        score = weight * score
        return score

    return 0

@numba.jit(nopython=True)
def getTrasProb(compactTrace, m):
    b4_trace = compactTrace[m - 1]
    after_trace = compactTrace[m]
    tProb = (1 - np.dot(b4_trace, after_trace))  # make dot product as transition prob
    if tProb < LOW_THRES_FOR_DEFULT_TRANS_PROP:
        tProb = LOW_THRES_FOR_DEFULT_TRANS_PROP

    return tProb

def getTrace(trace,b4,bb,strand):

    acgu = np.zeros(4,dtype='int32')
    if strand == False:
        b4 = b4-1
        bb = bb-1

    for n in range(b4,bb):

        a_trace = trace[n]
        if strand == True:
            acgu[0] += a_trace[0]
            acgu[1] += a_trace[1]
            acgu[2] += a_trace[2]
            acgu[3] += a_trace[3]
        else:
            acgu[0] += a_trace[3]
            acgu[1] += a_trace[2]
            acgu[2] += a_trace[1]
            acgu[3] += a_trace[0]

    return acgu

def getSeq(trace):

    s = []
    for t in trace:
        idx = np.argmax(t)
        a = ""
        if idx == 0 :
            a = "A"
        elif idx == 1:
            a = "C"
        elif idx == 2:
            a = "G"
        else:
            a = "T"
        s.append(a)

    return "".join(s)