import signalalign.RemapUtils as ru
import numba
import ruptures as rpt


def getLoaclInterval(exons, s):
    localitvl = []
    for exon in exons:
        tp = (exon.st_adjusted - s, exon.end_adjusted - s)
        localitvl.append(tp)

    return localitvl


def applyflipplopViterbi(chr, strand, s, e, es, ee, ref2b, reads):
    localgenome = ru.getRefG(chr, s, e, ref2b)
    # print(localgenome)
    retlist = []
    for read in reads:
        read = flipplopViterbiEach(read, localgenome, s, es, ee)
        retlist.append(read)

    return retlist


delmargin = 50
margin = 10


def flipplopViterbiEach(lgenome, chrom, strand, orgcigar,r_st, r_en, q_st, q_en, trace, move):
    # print("")
    possiblemove = addPossibleChangePoint(trace, move)
    possiblemove_idx = ru.toIndex(possiblemove, SEGMENT_ALL)
    frombasecaller_idx = ru.toIndex(possiblemove, SEGMENT_FROM_BASECALLER)
    compactTrace = ru.toCompact(trace, possiblemove_idx, frombasecaller_idx)

    # if strand == -1:
    #     possiblemove = possiblemove[::-1]
    #     trace = trace[::-1]
    #     compactTrace = ru.revcon(compactTrace)
    #     possiblemove_idx = ru.toIndex(possiblemove, SEGMENT_ALL)

    compactTracePositionMap = ru.toMap(possiblemove)
    # print("lgenome",lgenome)
    seq, cigar, left, traceoffset, traceboundary = viterbi(lgenome, compactTrace, trace, compactTracePositionMap,
                                                           strand,orgcigar, q_st, q_en,
                                                           possiblemove_idx, frombasecaller_idx)

    return seq, cigar, left, traceoffset, traceboundary, frombasecaller_idx, possiblemove_idx


import signalalign.OutputUtils as ou
from numba.typed import List
import pysam

def correctCigar(targetPos, cigar):

    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 3:  # N

            refpos = refpos + cigarlen

        elif cigaroprator == 0 or cigaroprator == 4:  # match or S softclip was not correted so treat as M

            if refpos + cigarlen > targetPos:
                return relpos + (targetPos - refpos)

            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del

            refpos = refpos + cigarlen

        elif cigaroprator == 1:  # Ins

            if relpos == 0:
                if targetPos <= cigarlen:
                    return 0

            relpos = relpos + cigarlen

    return 0


def viterbi(lgenome, compactTrace, trace, compactTracePositionMap, strand,orgcigar, q_st, q_en, possiblemove_idx,
            frombasecaller_idx):

    ctstart = ru.getNear(compactTracePositionMap, q_st)
    ctend = ru.getNear(compactTracePositionMap, q_en)

    #
    ctlen = len(compactTrace)
    ctstart = 0
    compactTraceExon = compactTrace[ctstart:ctend]
    strand = (strand == 1)

    m_range = List()
    for n in range(q_en - q_st):
        l = correctCigar(n,orgcigar)
        m = ru.getNear(compactTracePositionMap, q_st+l) - ctstart
        m_range.append(m)

    tracebackPath, left = viterbiEach(compactTraceExon, lgenome, possiblemove_idx, frombasecaller_idx, ctstart, m_range)

    tracebackPathL = []
    nb4 = -1
    for boundary in tracebackPath:
        n = boundary[0]
        m = boundary[1]
        if n!=nb4:
            tracebackPathL.append(m)
        nb4 = n


    traceboundary = []

    b4 = 0
    seq = ''
    dl = []
    traceseq = []
    didx = 0
    traceboundary.append(b4)
    for bb in tracebackPathL:

        # print(boundary)
        if bb != b4 and didx > 0:
            dl.append(didx)

            seq = seq + ru.getBase(trace, b4, bb, strand)
            traceseq.append(ru.getTrace(trace, b4, bb, strand))

            traceboundary.append(bb)

            # if strand == True:
            #
            # else:
            #     traceboundary.append(bb-1)

        b4 = bb
        didx += 1

    seqlen = len(seq)
    cigar = ou.toCigar(tracebackPath, seqlen)
    traceoffset = ctstart
    return seq, cigar, left, traceoffset, traceboundary


def addPossibleChangePoint(trace, move):
    movec = move.T
    tt = trace.T
    for flipflop in tt:
        algo_c = rpt.KernelCPD(kernel="linear", min_size=1).fit(flipflop)
        result = algo_c.predict(pen=penalty_value)
        for idx in result:
            if idx >= len(flipflop): break
            # b4, after = int(flipflop[idx - 1]), int(flipflop[idx])
            # diff = abs(b4 - after)
            # if max(b4, after) > TraceThres and diff > TraceThresDiff:
            #     # check if order change
            if movec[idx - 1] == 0:
                movec[idx - 1] = SEGMENT_FROM_CHANGE_POINT

    cnt = 1
    density0 = len(movec) // cnt
    lastidx = 0
    for idx in range(len(movec) - 1):

        if movec[idx] == SEGMENT_FROM_CHANGE_POINT:
            # clear change point if trace order unchange in small interval
            if not changeOrder(trace, idx) and (idx - lastidx) < density0:
                movec[idx] = 0

        if movec[idx] > 0:
            cnt += 1
        lastidx = idx

    density = len(movec) // cnt
    # devide too long segment
    idx = 0
    lastidx = 0
    for m in movec:
        if m > 0:

            if (idx - lastidx) > density * 3:
                unit = (idx - lastidx) // 3
                movec[lastidx + unit] = SEGMENT_FROM_CHANGE_POINT
                movec[lastidx + (unit * 2)] = SEGMENT_FROM_CHANGE_POINT

            elif (idx - lastidx) > density * 1.5:
                middleidx = (idx + lastidx) // 2
                movec[middleidx] = SEGMENT_FROM_CHANGE_POINT

            lastidx = idx

        idx += 1

    return movec


@numba.jit(nopython=True)
def changeOrder(trace, idx):
    firstidx0, secondidx0 = getFirstSecondIdx(trace[idx])
    firstidx1, secondidx1 = getFirstSecondIdx(trace[idx + 1])
    if (firstidx0 == firstidx1) and (secondidx0 == secondidx1):
        return False
    else:
        return True


@numba.jit(nopython=True)
def getFirstSecondIdx(trace):
    max = 0
    maxidx = 0
    second = 0
    secondidx = 0
    idx = -1
    for n in trace:
        idx += 1
        if n > max:

            if second < max:
                second = max
                secondidx = maxidx
            max = n
            maxidx = idx
        elif n > second:
            second = n
            secondidx = idx

    return maxidx, secondidx


import numpy as np
import numba

TraceThres = 5  # adhoc threshold to restrict change point
TraceThresDiff = 1  # adhoc threshold to restrict change point
penalty_value = 80  # beta

SEGMENT_FROM_BASECALLER = 1
SEGMENT_FROM_CHANGE_POINT = 2
SEGMENT_ALL = 3

DIOGONAL_MOVE = 1
HORIZONTAL_MOVE = 2
SKIP_MOVE_BASE = 10

IndelPenalty = 0.8
IndexExtentionpenalty = 0.1
MAXDEL_SIZE = 50
LOW_THRES_FOR_DEFULT_TRANS_PROP = 0.8
BONUS_FOR_EACH_SEGMENT = 0.1
count = 0
bannedinterval = 80


@numba.jit(nopython=True)
def rangeCheck(n, m, m_range):
    if n > len(m_range) - 1:
        m_in_alinment = m_range[len(m_range) - 1]
    else:
        m_in_alinment = m_range[n]

    inrange = ((m >= (m_in_alinment - bannedinterval)) and (m <= (m_in_alinment + bannedinterval)))
    diff = abs(m_in_alinment - m)
    wscore = 0
    if diff < 15:
        wscore = 1 - 0.05 * diff

    border = False
    if inrange == False:
        border = diff < bannedinterval
    return inrange, wscore, border

@numba.jit(nopython=True)
def getMRange(n, m_range,intervallen):

    if n > len(m_range) - 1:
        m_in_alinment = m_range[len(m_range) - 1]
    else:
        m_in_alinment = m_range[n]
    start = m_in_alinment - bannedinterval
    end = m_in_alinment + bannedinterval
    if start < 0:
        start = 0
    if end > intervallen:
        end = intervallen
    return start,end

from numba import objmode
import time
@numba.jit(nopython=True)
def viterbiEach(compactTrace, localgenome, possiblemove_idx, frombasecaller_idx, ctstart, m_range):

    intervallen = len(compactTrace)
    genomelen = len(localgenome)

    scorematrix = np.zeros((genomelen, intervallen), dtype='float32')
    movematrix = np.ones((genomelen, intervallen), dtype='uint8') # defult DIOGONAL move
    # score matrix
    maxscore = 0
    penalty = 0.05

    with objmode(t1='f8'):
        t1 = time.perf_counter()

    (maxn, maxm) = (0, 0)

    for n in range(genomelen):

        start,end = getMRange(n, m_range,intervallen)
        mrange = range(start,end)
        # print("in range",n,start,end)
        for m in mrange:

            inrange, wscore, border = rangeCheck(n, m, m_range)
            if inrange:
                if m == 0:

                    scorematrix[n][m] = round(ru.getStateProb(localgenome, compactTrace, n, m),2)

                else:

                    stateP = ru.getStateProb(localgenome, compactTrace, n, m)
                    if n == 0:
                        unTransScore = scorematrix[n][m - 1] + stateP
                        scorematrix[n][m] = round(unTransScore,2)
                        movematrix[n][m] = HORIZONTAL_MOVE

                    unTransScore = scorematrix[n][m - 1] + stateP - penalty
                    transScore = scorematrix[n - 1][m - 1] + stateP

                    maxdscore = 0
                    maxdin = 0
                    if n > MAXDEL_SIZE:

                        maxdscore = 0
                        maxdin = 0
                        for dlen in range(2, MAXDEL_SIZE + 1):

                            deltransScore = scorematrix[n - dlen][m - 1] + stateP - IndelPenalty - (
                                        dlen * IndexExtentionpenalty)
                            if deltransScore > maxdscore:
                                maxdscore = deltransScore
                                maxdin = dlen

                    if transScore > unTransScore and transScore > maxdscore:
                        scorematrix[n][m] = round(transScore,2)
                        movematrix[n][m] = DIOGONAL_MOVE

                        # if debug:
                        #     print(n,m,"diogonal",transScore,wscore)

                    elif unTransScore > maxdscore:
                        scorematrix[n][m] = round(unTransScore,2)
                        movematrix[n][m] = HORIZONTAL_MOVE
                        # if debug:
                        #     print(n,m,"untrans",unTransScore,wscore)

                    else:
                        scorematrix[n][m] = round(maxdscore,2)
                        movematrix[n][m] = SKIP_MOVE_BASE + maxdin
                        # if debug:
                        #     print(n,m,"SKIPBASE",maxdscore,wscore)

                    if maxscore < scorematrix[n][m]:
                        maxscore = scorematrix[n][m]
                        maxm = m
                        maxn = n

            else:

                movematrix[n][m] = DIOGONAL_MOVE

    left = genomelen - 1 - maxn
    tracebackPath = []
    n = maxn
    m = maxm

    countResolvedSegment = 1
    tracebackPath.append((n, possiblemove_idx[ctstart + m]))
    while n >= 0:

        if movematrix[n][m] == DIOGONAL_MOVE:
            n = n - 1
            m = m - 1
            if n >= 0 and m >= 0:
                tracebackPath.append((n, possiblemove_idx[ctstart + m]))

            # for final score to add bonus for sereoved segment
            countResolvedSegment += 1
            # print("diogonal",n, m)

        elif movematrix[n][m] >= SKIP_MOVE_BASE:

            dellen = int(movematrix[n][m] - SKIP_MOVE_BASE)
            m = m - 1
            for l in range(1, 1 + dellen):
                if n >= 0 and m >= 0:
                    n -= 1
                    tracebackPath.append((n, possiblemove_idx[ctstart + m]))
                    # print("skipbase", n, m)


        else:
            # no state change
            m = m - 1
            # overwrite
            if m >= 0:
                tracebackPath.pop(-1)
                tracebackPath.append((n, possiblemove_idx[ctstart + m]))
            # print("hrizontal", n, m)

        if m <= 0:
            break

    tracebackPath.reverse()
    trace = None
    compactTrace = None
    scorematrix = None
    scorematrix_os = None
    movematrix = None
    # print("tracebackPath,left",left)

    with objmode(t2='f8'):
        t2 = time.perf_counter()
    print("maxscore",maxscore,genomelen,t2-t1)

    return tracebackPath, left



