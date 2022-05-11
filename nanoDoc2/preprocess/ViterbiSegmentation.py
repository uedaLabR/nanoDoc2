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


def flipplopViterbiEach(lgenome, chrom, strand, r_st, r_en, q_st, q_en, trace, move):
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
                                                           strand, q_st, q_en,
                                                           possiblemove_idx, frombasecaller_idx)
    return seq, cigar, left, traceoffset, traceboundary, frombasecaller_idx, possiblemove_idx


import signalalign.OutputUtils as ou
from numba.typed import List


def viterbi(lgenome, compactTrace, trace, compactTracePositionMap, strand, q_st, q_en, possiblemove_idx,
            frombasecaller_idx):
    ctstart = ru.getNear(compactTracePositionMap, q_st)
    # print(compactTracePositionMap)
    ctend = ru.getNear(compactTracePositionMap, q_en)
    # print("ct_starte_end",q_st,q_en,ctstart,ctend)
    #
    ctlen = len(compactTrace)
    ctstart = 0
    compactTraceExon = compactTrace[ctstart:ctend]
    strand = (strand == 1)

    m_range = List()
    for n in range(q_en - q_st):
        m = ru.getNear(compactTracePositionMap, q_st + n) - ctstart
        m_range.append(m)

    tracebackPath, left = viterbiEach(compactTraceExon, lgenome, possiblemove_idx, frombasecaller_idx, ctstart, m_range)
    traceboundary = []

    b4 = 0
    seq = ''
    dl = []
    traceseq = []
    didx = 0
    traceboundary.append(b4)
    for boundary in tracebackPath:

        bb = boundary[1]
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
            if movec[idx - 1] == 0:
                movec[idx - 1] = SEGMENT_FROM_CHANGE_POINT

    cnt = 1
    density0 = len(movec) // cnt
    lastidx = 0
    # for idx in range(len(movec) - 1):
    #
    #     if movec[idx] == SEGMENT_FROM_CHANGE_POINT:
    #         # clear change point if trace order unchange in small interval
    #         if not changeOrder(trace, idx) and (idx - lastidx) < density0:
    #             movec[idx] = 0
    #
    #     if movec[idx] > 0:
    #         cnt += 1
    #     lastidx = idx

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

IndelPenalty = 1.5
IndexExtentionpenalty = 0.4
MAXDEL_SIZE = 9
LOW_THRES_FOR_DEFULT_TRANS_PROP = 0.8
BONUS_FOR_EACH_SEGMENT = 0.1
count = 0

bannedinterval = 150


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

    return inrange, wscore


@numba.jit(nopython=True)
def viterbiEach(compactTrace, localgenome, possiblemove_idx, frombasecaller_idx, ctstart, m_range):
    intervallen = len(compactTrace)
    genomelen = len(localgenome)
    # print("genome len",genomelen)
    # print("iv len",intervallen)

    scorematrix = np.zeros((genomelen, intervallen), dtype='float32')
    movematrix = np.zeros((genomelen, intervallen), dtype='float32')
    # score matrix
    maxscore = 0
    penalty = 0.05
    untransPenaltyFordefultsegment = 0.95
    traceRatios = None

    (maxn, maxm) = (0, 0)

    for n in range(genomelen):
        for m in range(intervallen):

            inrange, wscore = rangeCheck(n, m, m_range)
            if inrange:
                if m == 0:

                    scorematrix[n][m] = ru.getStateProb(localgenome, compactTrace, n, m)

                else:

                    stateP = ru.getStateProb(localgenome, compactTrace, n, m)
                    if n == 0:
                        unTransScore = scorematrix[n][m - 1] + stateP
                        scorematrix[n][m] = unTransScore
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
                        scorematrix[n][m] = transScore
                        movematrix[n][m] = DIOGONAL_MOVE

                        # if debug:
                        #     print(n,m,"diogonal",transScore,wscore)

                    elif unTransScore > maxdscore:
                        scorematrix[n][m] = unTransScore
                        movematrix[n][m] = HORIZONTAL_MOVE
                        # if debug:
                        #     print(n,m,"untrans",unTransScore,wscore)

                    else:
                        scorematrix[n][m] = maxdscore
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
    return tracebackPath, left



