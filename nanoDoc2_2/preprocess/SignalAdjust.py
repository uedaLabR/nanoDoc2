import numpy as np

unit = 10

from numba import jit,u1,i8,f8

@jit
def decode16bit(a_trace):

    a = (a_trace & 0b1111000000000000) >> 12
    c = (a_trace & 0b0000111100000000) >> 8
    g = (a_trace & 0b0000000011110000) >> 4
    t = (a_trace & 0b0000000000001111)
    #
    return (a,c,g,t)


import numpy as np
def decode(trace):

    tp = list(map(decode16bit, trace))
    ar = np.array(tp)
    return ar

def _getTraceSeq(traceintervals,trace):

    print(trace.shape)
    trace = trace.T
    print("decode trace")
    print(trace.shape)
    sl = []
    seq = ["A","C","G","T"]
    for m in range(1,len(traceintervals)):

        b4 = traceintervals[m-1]
        idx = traceintervals[m]
        partialtrace = trace[:,b4:idx]
        su = np.sum(partialtrace, axis=1)
        maxidx = su.argmax()
        s = seq[maxidx]
        sl.append(s)

    return "".join(sl)


def reconst(traceintervals):

    r = []
    t = 0
    for m in traceintervals:
        t = t+m
        r.append(t)
    r = np.array(r)
    return r



@jit
def getTraceSeq_J(traceintervals,trace):

    seq = []
    tracelen = len(trace)
    idx = 0
    lastbaecallidx = 0
    a, c, g, u = 0, 0, 0, 0
    for m in range(1, tracelen):

       if m in traceintervals:

           maxn = max(a,c,g,u)
           if a == maxn:
               seq.append("A")
           if c == maxn:
               seq.append("C")
           if g == maxn:
               seq.append("G")
           if u == maxn:
               seq.append("T")

           a, c, g, u = 0, 0, 0, 0

       a_trace = trace[m]
       _a = max(a_trace[0], a_trace[4])
       _c = max(a_trace[1], a_trace[5])
       _g = max(a_trace[2], a_trace[6])
       _u = max(a_trace[3], a_trace[7])

       a = a + _a
       c = c + _c
       g = g + _g
       u = u + _u

    return seq



import pysam
def addsla(readseq,cigar):

    str_list = list(readseq)
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    relpos = 0
    refpos = 0
    for cigaroprator, cigarlen in a.cigar:


        if cigaroprator == 3:  # N

            # refpos = refpos + cigarlen
            pass

        elif cigaroprator == 0:  # match or S softclip was not correted so treat as M

            refpos = refpos + cigarlen
            relpos = relpos + cigarlen

        elif cigaroprator == 1:  # Ins

            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del

            for n in range(cigarlen):

                str_list.insert(relpos, '-')

            relpos = relpos + cigarlen

    return ''.join(str_list)



def toNum(move):

    cnt = 0
    l = []
    for m in move:
        if m == 1:
            l.append(cnt)
        cnt+=1
    return l

def getNadd(cigar):

    seq = []
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    for cigaroprator, cigarlen in a.cigar:
        break

    if cigaroprator == 3:
        for m in range(cigarlen):
            seq.append("N")

    return seq

def getTraceSeq(signalintervals,trace):


    trace = trace.T
    #print(trace)
    sl = []
    seq = ["A","C","G","T","A","C","G","T"]

    for m in range(1,len(signalintervals)):

        b4 = int(signalintervals[m-1])
        idx = int(signalintervals[m])

        if b4 == idx:
            sl.append("N")
            # print("N",m)
            continue

        partialtrace = trace[:,b4:idx]
        su = np.sum(partialtrace, axis=1)
        #print("su",su)
        maxidx = su.argmax()
        s = seq[maxidx]
        sl.append(s)

    return "".join(sl)


def cm(r1,r2):

    cnt = 0
    for n in range(min(len(r1),len(r2))):

        rr1 = r1[n]
        rr2 = r2[n]
        if rr1 == rr2:
           cnt = cnt+1

    return cnt

def adjustMismatchindel(read,fmerDict):

    gseq = read.refgenome
    # print("pass1")
    traceintervals = np.array(read.traceboundary)
    # print(read.traceboundary)
    # print("pass2")
    #readseq = getNadd(read.cigar_str).extend(getTraceSeq(traceintervals,read.trace))
    ns = getNadd(read.cigar_str)
    numberofFisrstN = len(ns)
    print(ns)
    readseq = getTraceSeq(traceintervals,read.trace)
    ns.extend(readseq)

    # print("pass3")
    readseqwithsla = addsla(ns,read.cigar_str)
    # print("pass4")

    # print(read.cigar_org)
    # print(read.strand)
    r = read.fastq.split('\n')[1]
    # print("traceboundary", len(traceintervals))
    # print("gseq     ",len(gseq),gseq)
    # print("read",len(r),"".join(r))
    # print(read.cigar_str,ns)
    # print("readseq", len(readseq), "".join(readseq))
    # print("readwiths",len(readseqwithsla),readseqwithsla)
    intarvals = au.findMismatchInterval(gseq,readseqwithsla)
    # print("interval",intarvals)
    # print("pass5")

    n = -2
    signalInterval = adjustWithDTW(read,intarvals,len(gseq),fmerDict,traceintervals,n)
    traceb = np.array(signalInterval) // 10
    #readseq2 = getTraceSeq(traceb, read.trace)

    print("traceinterval",len(traceintervals),traceintervals)
    print("signalinterval",len(traceb),traceb)
    print("interval", intarvals)

    # print(read.cigar_org)
    # print("traceintervals", len(traceintervals))
    # print("signalInteral",len(signalInteral))

    return signalInterval

def relPos(pos,cgl):


    lpos = pos
    cpos = 0
    isFirst = True
    for cigaroprator, cigarlen in cgl:

        if cigaroprator == 3:  # N

            if isFirst:
                cpos = cpos - cigarlen
            # refpos = refpos + cigarlen
            pass

        elif cigaroprator == 0:  # match or S softclip was not correted so treat as M

            if lpos <= cigarlen:
                return cpos + lpos
            else:
                cpos += cigarlen
                lpos = lpos - cigarlen

        elif cigaroprator == 2:  # Del

            lpos = lpos - cigarlen

        isFirst = False
    #should not come here
    return cpos + lpos


def pathToIndex(path,subsignal_len):

    pl = list(map(lambda x:x[0],path))
    arr = np.array(pl)
    diff = np.diff(arr)
    idx = np.where(diff == 1)[0]
    idx = list(map(lambda x:x+1,idx))
    l = []
    l.append(0)
    for m in idx:
        refidx = path[m][1]
        l.append(refidx*5)
    l.append(subsignal_len)
    return l


def replaceList(signalInterval ,indexes,size):


    try:

        # start = signalInterval.index(indexes[0])
        start = np.where(signalInterval == indexes[0])[0]
        if len(start) == 0:
            print("Error in replace key " ,indexes[0] )
            print(signalInterval)
            return signalInterval,False

        start = int(start[0])
        end = start + size
        # print("index",len(indexes),indexes)
        # print("size",size)
        # print("siginterval",len(signalInterval))
        # print("b4",len(signalInterval[start:end+1]),signalInterval[start:end+1])
        # print("after",len(indexes),indexes)
        signalInterval[start:end+1] = indexes
        # print("siginterval", len(signalInterval))


        return signalInterval,True

    except ValueError:
        import traceback
        traceback.print_exc()
        print("Value Error at replace list")
        pass

from scipy import signal
@jit
def binnedSignalAve(subsignal):

    idx = 0
    l = []
    while idx+5 < len(subsignal):

        m = np.mean(subsignal[idx:idx+5])
        l.append(m)
        idx+=5
    if idx < len(subsignal):
        l.append(np.mean(subsignal[idx:]))
    else:
        l.append(0.0)

    return l

import nanoDoc2_2.preprocess.AdjustUtils as au
import fastdtw as fastdtw
from scipy.spatial.distance import euclidean
import ruptures as rpt
penalty_value = 60
def adjustWithDTW(read,intarvals,leng,fmerDict,traceintervals,n):

    a = pysam.AlignedSegment()
    a.cigarstring = read.cigar_org
    cgl = []

    for cg in a.cigar:
        cgl.append(cg)

    l =[]
    for iv in intarvals:
        l.append(iv)

    # print("intarvals", l)
    # print("cgs", cgl)
    lfunc = lambda x: x *10
    signalInterval = lfunc(traceintervals)
    signalInterval = signalInterval.tolist()

    addindex = 0

    for iv in l:

        gstart = iv[0]
        gend = iv[1]

        if gstart < 3:
            continue
        if gend > len(read.refgenome)-6:
            continue

        start = relPos(gstart,cgl)
        end = relPos(gend,cgl)

        # if read.strand == True:
        #
        #     start0 = leng - end
        #     end = leng - start
        #     start = start0


        if start < 3:
            continue
        if end > len(traceintervals)-1:
            continue

        signal_ivstart = traceintervals[start] * unit
        signal_ivend = traceintervals[end] * unit

        #
        strand = True
        lgenome = read.refgenome[gstart+n:gend+n+6]

        theorymean = au.theoryMean(fmerDict, lgenome, strand)
        subsignal = read.signal[signal_ivstart:signal_ivend]

        binnedSignal = binnedSignalAve(subsignal)

        distance, path = fastdtw.fastdtw(theorymean, binnedSignal, dist=euclidean)
        print("theorymean",theorymean)

        indexes = pathToIndex(path,len(subsignal))
        # print(indexes)

        diffidx = (gend-gstart) - (end-start)
        indexes = indexes + (signal_ivstart)
        #
        size = (end-start)
        signalInterval,success = replaceList(signalInterval,indexes,size)

        if success:
            addindex = addindex + diffidx
            # print("diffidx",diffidx)

    return signalInterval
