import h5py
from statsmodels import robust
import numpy as np
import random
from scipy import interpolate
import pyarrow as pa
import pyarrow.parquet
import pandas as pd

from signalalign import ViterbiSegmantation


def equalize(bindata,binsize):


    l = len(bindata)

    if l == 0:
        return addnoise([0] * binsize)
    elif l == 1:
        return addnoise([int(bindata[0])] * binsize)
    elif l == 2:
        datal1 = [int(bindata[0])] * int(binsize / 2)
        datal2 = [int(bindata[1])] * int(binsize / 2)
        datal1.extend(datal2)
        return addnoise(datal1)
    elif l == binsize:
        return bindata
    elif l < binsize:
        rd = np.zeros(binsize)
        for m in range(len(rd)):
            idx = (int)((m/len(rd)) *len(bindata))
            if idx >= len(bindata):
                idx = len(bindata)-1
            rd[m] = bindata[idx]
        return rd
    else:
        return downsampling(bindata, binsize)

def addnoise(npa):
    npa = np.array(npa)
    noise_unit = 0.005
    noise = np.random.normal(-noise_unit, noise_unit, len(npa))
    return npa + noise

def downsampling(bindata, binsize):
    lsize = len(bindata)
    indexs = random.sample(range(0, lsize), binsize)
    ls = []
    for i in range(0, lsize):
        if i in indexs:
            ls.append(bindata[i])
    return ls


def toByteRange(x):
    y = int(x) + 128
    if y > 255:
        #        print(y)
        y = 255
    if y < 0:
        #        print(y)
        y = 0
    return y


def indexFiles(path):

    unit = 4000
    cnt = 0
    idxcnt = 0
    ref = ""
    lastref = ""
    blist = []
    indexlist = []
    with open(path) as f:
        for s_line in f:
            ref = s_line.split(",")[1]
            if ref != lastref or cnt % unit == 0:
                if len(indexlist) > 0:
                    blist.append(indexlist)
                indexlist = []
                cnt = 0
            indexlist.append(s_line)
            cnt = cnt + 1
            lastref = ref
        if len(indexlist) > 0:
            blist.append(indexlist)

    return blist

import multiprocessing
from multiprocessing import Pool
import functools


def _makeParquet(tp,ref):

    aligner = mp.Aligner(ref, preset="map-ont")
    (indexcontent, pathout, ref) = tp

    wholedata = []
    chrom = "NC_000913.3"
    lenlimit = 6000
    cnt_1 = 0
    for s_line in indexcontent:
        fdata = s_line.split(",")
        onerecord = getRecord(fdata)
        if onerecord != None or withinRange(onerecord,lenlimit):

            onerecord = alginandremap(onerecord, aligner)
            if onerecord != None:
                wholedata.append(onerecord)
                cnt_1 +=1

        if cnt_1 % 10 ==0:
            print(cnt_1)



    df = pd.DataFrame(wholedata,
                      columns=['mapped_chrom', 'mapped_strand', 'mapped_start', 'mapped_end', 'clipped_bases_start',
                               'clipped_bases_end', 'fastq', 'chrom', 'st', 'ed', 'cigar','md', 'move','trace','signal',
                               'sgmean', 'sgst', 'sglen', 'bases'
                               ])
    df.to_parquet(pathout)

def withinRange(onerecord,lenlimit):

    mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, ctg, st, ed, \
    cigar, md, move, trace, signal, sgmean, sgst, sglen, bases, viterbistart, raw_start = onerecord
    fastq = fastq.split("\n")
    seq = fastq[1]
    return len(seq) <= lenlimit

def ciglen(cigar):

    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 0:  # match


            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del
            refpos = refpos + cigarlen
        elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
            relpos = relpos + cigarlen

    return refpos

import mappy
def alginandremap(onerecord, aligner):

    mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, ctg, st, ed, \
    cigar, md, move, trace, signal, sgmean, sgst, sglen, bases, viterbistart, raw_start = onerecord

    # tRNA_species = hit.ctg, r_st = hit.r_st, r_en = hit.r_en,
    # q_st = hit.q_st, q_en = hit.q_en, cigar_str = hit.cigar_str, MD = hit.MD, score = score

    chrom,strand, r_st, r_en, q_st, q_en, cigar_str = "", 0, 0, 0, 0, 0, ""
    fastq = fastq.split("\n")
    seq = fastq[1]

    for hit in aligner.map(seq):

        chrom = hit.ctg
        strand = hit.strand
        r_st = hit.r_st
        r_en = hit.r_en
        q_st = hit.q_st
        q_en = hit.q_en
        cigar_str = hit.cigar_str


    lgenome = aligner.seq(chrom, start=r_st, end=r_en)
    if lgenome is None:
        return None
    lgenome = mappy.revcomp(lgenome)
    seq, cigar, left, traceboundary,frombasecaller_idx,possiblemove_idx = ViterbiSegmantation.flipplopViterbiEach(lgenome,chrom,strand,r_st,r_en,q_st,q_en,trace,move)

    print("information")
    print(r_st,r_en)
    print(r_st, r_st + ciglen(cigar))
    mapped_start = r_st
    mapped_end = r_en
    clipped_bases_start = 0
    clipped_bases_end = q_st
    traceboundary = [n * 10 for n in traceboundary]
    sgst = traceboundary
    #print("sgst",sgst)
    signal = signal[viterbistart:]
    trace = trace.flatten()
    return mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, ctg, st, ed, \
    cigar, md, move, trace, signal, sgmean, sgst, sglen, bases


from matplotlib import pyplot as plt
from matplotlib import gridspec
import numba

import pysam
def getrelpos(cigar,pos):

    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 0:  # match

            if refpos + cigarlen > pos:
                return relpos + (pos - refpos)

            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del
            refpos = refpos + cigarlen
        elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
            relpos = relpos + cigarlen

    return 0

import itertools
base_corresponding_table = {0:'A',2:'G',1:'C',3:'T',4:'A-',6:'G-',5:'C-',7:'T-'}
base_color = {'A':'#228b22','T':'#db7093','G':'#ff8c00','C':'#4169e1','A-':'#228b22','T-':'#db7093','G-':'#ff8c00','C-':'#4169e1'}
def plotboth(signal_t, signal_v,traceboundary,trace,move,lgenome,cigar):

    limit = 25000
    fig = plt.figure(figsize=(360, 20))

    gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 0.5])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    # ax3 = fig.add_subplot(gs[2])
    # ax4 = fig.add_subplot(gs[3])
    # ax5 = fig.add_subplot(gs[4])


    #
    # plot signal
    ax2.plot(signal_t, linewidth=2)
    # ax5.plot(signal_v, linewidth=2)

    # plot trace
    traces = tuple(zip(*trace))
    for i in range(len(traces)):
        temp = [[trace] * 10 for trace in traces[i]]  # extend the index to match the index of signal
        trace = list(itertools.chain.from_iterable(temp))
        lstyle = "solid" if i <= 3 else "dashed"
        ax1.plot(trace, color=base_color[base_corresponding_table[i]], linestyle=lstyle, linewidth=1.0)
        #ax2.plot(trace, color=base_color[base_corresponding_table[i]], linestyle=lstyle, linewidth=1.0)
        #ax3.plot(trace, color=base_color[base_corresponding_table[i]], linestyle=lstyle, linewidth=1.0)

    # plot segmentation by guppy
    move_positions = move
    print(move)

    # for n in range(len(move_positions)):
    #     if move[n] == 1:
    #         ax1.axvline(x=n*10, ymin=0, ymax=1, color='black', alpha=0.2)

    print(traceboundary)
    # plot segmentation by viterbi algorithm
    segmentation_positions = traceboundary
    cnt = 0
    for index in segmentation_positions:
        color = 'black'
        ax1.axvline(x=index*10, ymin=0, ymax=1, color=color, alpha=0.2)
    for n in range(min(len(lgenome),len(segmentation_positions))):
        nuc = lgenome[n]
        relpos = getrelpos(cigar,n)
        index = segmentation_positions[relpos]*10
        ax1.text(index, 260, nuc, size=7, color=base_color[nuc], ha='center')


    return fig


def normalise(signal):

    med = np.median(signal)  # take median
    mad = robust.mad(signal)  # take mad
    signal = ((signal - med) / mad)  # normalize
    signal = (((signal - 1) / 4) * 128)  # fit into byte
    signal = signal.astype(np.int8)  # type conversion
    return signal


def jointByteArray(binnedlist):
    elist = []
    for ll in binnedlist:
        ll = list(map(toByteRange, ll))
        elist.extend(ll)
    return bytes(elist)


def getFastQ(f, key):
    if key in f:
        group = f[key]
        datalist = list(group.values())
        fastq = None
        trace = None
        move = None
        for d in datalist:
            if "Fastq" in str(d):
                fastq = d
            if "Trace" in str(d):
                trace = d
            if "Move" in str(d):
                move = d

        fastq = fastq[()]
        trace = trace[()]
        move = move[()]
        try:
            fastq = fastq.decode('utf-8')
        except AttributeError:
            pass
        return fastq,trace,move
    return None

import mappy as mp
def getRecord(fdata):
    binsize = 60

    with h5py.File(fdata, 'r') as f:

        group = f['/Analyses/RawGenomeCorrected_000/BaseCalled_template']
        datalist = list(group.values())
        alignment = datalist[0]

        event = datalist[1]
        allevent = event[()]
        raw_start = event.attrs.get('read_start_rel_to_raw')
        dlen = len(allevent)-1
        end = allevent[dlen][2]

        sgmean = []
        sgst = []
        sglen = []
        bases = []
        seq = ""
        for n in range(len(allevent)):

            mean = allevent[n][0]
            start = allevent[n][2]
            size = allevent[n][3]
            b = allevent[n][4]

            sgmean.append(mean)
            sgst.append(start)
            sglen.append(size)
            bases.append(b)
            seq = seq + str(b).replace("b","").replace("'","")


        group = f['/Raw/Reads/']
        readid = list(group.values())[0].name
        signalkey = readid + "/Signal"
        signal = f[signalkey][()]

        scaling = 1163.0 / 8192.0
        offset = 10
        signal = np.array(scaling * (signal + offset), dtype=np.float32)

        #signal = normalise(signal)
        #signal = signal[raw_start:end] #upto tombo segmantation

        mapped_chrom = alignment.attrs.get('mapped_chrom')
        mapped_strand = alignment.attrs.get('mapped_strand')
        mapped_start = alignment.attrs.get('mapped_start')
        mapped_end = alignment.attrs.get('mapped_end')

        clipped_bases_start = alignment.attrs.get('clipped_bases_start')
        clipped_bases_end = alignment.attrs.get('clipped_bases_end')

        fastqtp =  getFastQ(f, '/Analyses/Basecall_1D_001/BaseCalled_template')
        if fastqtp == None:
            fastq,trace,move = getFastQ(f, '/Analyses/Basecall_1D_000/BaseCalled_template')
        else:
            fastq, trace, move = fastqtp

        ctg,st,ed,cigar,md = "",0,0,"",""
        trace = trace[::-1].astype(np.int16)

        viterbistart = len(signal) - 10 * len(trace)

        move = move[::-1].astype(np.int16)
        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,move,trace,signal,sgmean,sgst,sglen,bases,viterbistart,raw_start)
        return tp

    return None


def toSeq(fastq):

    fa = fastq.split("\n")
    return fa[1]

from statistics import mean
def getMeans(signal,traceboundary):

    means = []
    tbb4 = 0
    for tb in traceboundary:

        if tb == 0:
            continue
        if tb*10 > len(signal):
            break
        means.append(mean(signal[tbb4*10:tb*10]))
        tbb4 = tb

    return means

def getMeansT(signal,traceboundary):

    means = []
    tbb4 = 0
    for tb in traceboundary:

        if tb == 0:
            continue
        if tb > len(signal):
            break
        means.append(mean(signal[tbb4:tb]))
        tbb4 = tb

    return means


def theoryMean(fmercurrent,lgenome):

    a = {}
    with open(fmercurrent) as f:
        cnt = 0
        for line in f:
            if cnt > 0:
                data = line.split()
                a[data[0]] = float(data[1])
            cnt = cnt+1

    means = []
    rg = lgenome[::-1]
    for n in range(0,len(rg)-5):

       fmer = rg[n:n+5]
       cv = a[fmer]
       means.append(cv)


    means.reverse()
    return means

def moda(a,shift):

    if shift ==0:
        return a
    elif shift > 0:
        return a[shift:]
    else:
        a = list(a)
        for n in range(abs(shift)):
            a.insert(0, 0.5)
        return np.array(a)

def predictShift(a,b):

    maxn = 0
    max = 0
    ar = None
    br = None
    for n in range(-3,3):

        bmod = moda(b,n)
        if len(a) > len(bmod):
            a = a[0:len(bmod)]
        if len(bmod) > len(a):
            bmod = bmod[0:len(a)]

        v = np.dot(a,bmod)
        if v > max:
            max = v
            maxn = n
            ar = a
            br = bmod

    return maxn,ar,br

def calcNormalizeScaleLMS(signalmeans,theorymean):

    y1 = 0
    y2 = 0
    y2_pow2 = 0
    y1y2 = 0
    n_len = len(signalmeans)
    for n in range(len(signalmeans)):

        y1 = y1 + theorymean[n]
        y2 = y2 + signalmeans[n]
        y1y2 = y1y2 + y1*y2
        y2_pow2 = y2_pow2 + y2*y2

    a = [[y2_pow2, y2], [y2, n_len]]
    ainv = np.linalg.inv(a)
    b = np.array([[y1y2], [y1]])
    return np.dot(ainv, b)

def toMatrix(frombasecaller_idx,traceboundary):

    maxn = max(len(frombasecaller_idx),len(traceboundary))+1
    maxm = max(max(frombasecaller_idx),max(traceboundary))+1
    scorematrix = np.zeros((maxn, maxm), dtype='float32')

    c = 0
    for idx in frombasecaller_idx:
        scorematrix[c][idx] = 1.0
        c+=1
    c = 0
    for idx in traceboundary:
        scorematrix[c][idx] = 2.0
        c+=1

    return scorematrix


def traceseq(compactTrace):

    ret = []
    for idx in compactTrace:
        idx = idx[0:4]
        midx = np.argmax(idx)
        print(idx,midx)
        label = ["A", "C", "G", "T"]
        nuc = label[midx]
        ret.append(nuc)
    return ret

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import signalalign.RemapUtils as ru

if __name__ == "__main__":

   #path = '/data/nanopore/rRNA/1623_ivt-multi/singleFast5/137/462bc3b0-c01e-4404-b289-aa22066df3c9.fast5'
   #path = '/data/nanopore/rRNA/1623_ivt-multi/singleFast5/273/1ea4c5e1-9ded-4873-9e18-0fb706e20e16.fast5'
   path = '/data/nanopore/rRNA/1623_ivt-multi/singleFast5/118/c2610804-dcb2-4305-9fb6-ffc0a64a74f4.fast5'
   pathout = '/data/nanopore/rRNA/test'
   ref ="/data/nanopore/reference/NC000913.fa"
   indexf = "/data/nanopore/rRNA/test/index.txt"
   path_w = pathout + "/sampleplingplan.txt"

   fmercurrent ="/data/nanopore/signalStatRNA180.txt"

   samplesize = 400

   p = getRecord(path)
   print(p)
   aligner = mp.Aligner(ref, preset="map-ont")

   mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, ctg, st, ed, \
    cigar, md, move, trace, signal, sgmean, sgst, sglen, bases, viterbistart, raw_start = p

    # tRNA_species = hit.ctg, r_st = hit.r_st, r_en = hit.r_en,
    # q_st = hit.q_st, q_en = hit.q_en, cigar_str = hit.cigar_str, MD = hit.MD, score = score

   chrom,strand, r_st, r_en, q_st, q_en, cigar_str = "", 0, 0, 0, 0, 0, ""
   fastq = fastq.split("\n")
   seq = fastq[1]

   for hit in aligner.map(seq):

       chrom = hit.ctg
       strand = hit.strand
       r_st = hit.r_st
       r_en = hit.r_en
       q_st = hit.q_st
       q_en = hit.q_en
       cigar_str = hit.cigar_str

   print(chrom,strand,r_st,r_en,q_st,q_en,cigar_str)
   lgenome = aligner.seq(chrom, start=r_st, end=r_en)
   print(len(lgenome),lgenome)
   print(len(seq),seq.replace('U','T'))
   # lgenome = lgenome[0:50]
   # trace = trace[0:1500]
   # move = move[0:1500]
   # r_en = r_st +50
   # q_en = q_st +50
   seq, cigar, left,traceoffset, traceboundary, frombasecaller_idx,possiblemove_idx = ViterbiSegmantation.flipplopViterbiEach(lgenome,
                                                                                                                chrom,
                                                                                                                strand,
                                                                                                                r_st,
                                                                                                                r_en,
                                                                                                                q_st,
                                                                                                                q_en,
                                                                                                                trace,
                                                                                                                move)

   print(len(seq), seq.replace('U', 'T'))
   print(len(signal),raw_start, viterbistart)
   print("information")
   print(r_st, r_en)
   print(r_st, r_st + ciglen(cigar))
   signal_t = signal
   signal_v = signal[viterbistart:]
   signal_t = signal_t[::-1]  # reverse for rna
   signal_t = signal_t[raw_start:]
   signal_v = signal_v[::-1]  # reverse for rna
   # frombasecaller_idx = frombasecaller_idx[0:50]
   # signal_t = signal_t[0:15000]
   # signal_v = signal_v[0:15000]
   print(len(signal_t),len(signal_v))
   signalmeans = getMeans(signal_v,traceboundary)
   #signalmeans = getMeans(signal_v, frombasecaller_idx)

   tombosignalmean = getMeansT(signal_t,sgst)
   theorymean = theoryMean(fmercurrent,lgenome)

   print(cigar, left,traceboundary)
   print(signalmeans)
   print(theorymean)
   print("sgst",len(sgst),sgst)
   print("traceboundary", len(traceboundary),traceboundary)
   print("frombasecaller_idx",len(frombasecaller_idx),frombasecaller_idx)

   print(chrom,r_st,cigar,seq.replace('U', 'T'))

   shift,signalmeans,theorymean = predictShift(signalmeans,theorymean)
   scaleshift = calcNormalizeScaleLMS(signalmeans,theorymean)
   print(scaleshift)
   a =  scaleshift[0][0]
   b =  scaleshift[1][0]
   signal_t = signal_t * a + b
   signal_v = signal_v * a + b

   signal_t = np.clip(signal_t, 60, 160)
   signal_t = np.clip(signal_v, 60, 160)
   signal_t = signal_t -60  / 128
   signal_v = signal_v - 60 / 128
   fig = plotboth(signal_t, signal_v,traceboundary,trace,move,lgenome,cigar)
   fig.savefig("/data/nanopore/testnormalize_trace.png")

   # print(type(scorematrix))
   # print(scorematrix)
   # plt.figure()
   # sns.heatmap(scorematrix)
   # plt.savefig('/data/nanopore/scoreMatrix.png')
   #
   # sc = toMatrix(frombasecaller_idx,traceboundary)
   # print(sc)
   # plt.figure()
   # sns.heatmap(sc)
   # plt.savefig('/data/nanopore/traceMatrix.png')
   # plt.close('all')
   #
   # compactTrace = ru.toCompact(trace, possiblemove_idx,frombasecaller_idx)
   # compactTrace = compactTrace[0:90]
   # traceseq = traceseq(compactTrace)
   # df = pd.DataFrame(data=scorematrix,columns=traceseq, index=list(lgenome))
   # df.to_csv('/data/nanopore/score.csv')
