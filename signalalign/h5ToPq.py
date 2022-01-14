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
def makeParquet(indexf,ref, pathout,  thread):




    indexlist = indexFiles(indexf)
    blist = []
    idx = 1
    for d in indexlist:

        pathoutpq = pathout + "/" + str(idx) + ".pq"
        tp = (d, pathoutpq, ref)
        print(tp)
        blist.append(tp)
        idx = idx + 1

    size = len(blist)
    #    print(size)
    # graph_manager = GraphManager('/share/data/IVT/wt.pdf')
    ncore = 8
    _makeParquet_2 = functools.partial(_makeParquet, ref=ref)
    with Pool(ncore) as p:
        p.map(_makeParquet_2, blist)


    # graph_manager.save()

    # with Pool(thread) as p:
    #      p.map(_makeParquet, blist)


#        tqdm(imap,total=size)


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


    # print("tracelen",len(trace))
    # print("siglen", len(signal))

    #
    lgenome = aligner.seq(chrom, start=r_st, end=r_en)
    if lgenome is None:
        return None
    lgenome = mappy.revcomp(lgenome)
    seq, cigar, left, traceboundary,frombasecaller_idx,possiblemove_idx = ViterbiSegmantation.flipplopViterbiEach(lgenome,chrom,strand,r_st,r_en,q_st,q_en,trace,move)
    # print(seq)
    # print(cigar)
    # print(chrom, strand,q_st,q_en, r_st, r_en, mapped_start, mapped_end)

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

    # viterbimean = []
    # tombomean = []
    # b4 = 0
    # for b in traceboundary:
    #     if b4!=0:
    #         me = np.mean(signal[b4*10:b*10])
    #         viterbimean.append(me)
    #     b4 = b
    #
    # idx = 0
    #
    # for s in sgst:
    #
    #     sl = sglen[idx]
    #     me = np.mean(signaltommbo[s:s+sl])
    #     tombomean.append(me)
    #     idx +=1

    # print("tb",len(traceboundary),traceboundary)
    # print("sgst",len(sgst),sgst)
    #
    # fig = plot100(traceboundary,sgst,signal,trace,viterbistart, raw_start,frombasecaller_idx,possiblemove_idx)
    # fig.savefig("/data/nanopore/img.png")

    # print("vm",viterbimean)
    # print("tm", tombomean)

from matplotlib import pyplot as plt
from matplotlib import gridspec
import itertools
base_corresponding_table = {0:'A',2:'G',1:'C',3:'U',4:'A-',6:'G-',5:'C-',7:'U-'}
base_color = {'A':'#228b22','U':'#db7093','G':'#ff8c00','C':'#4169e1','A-':'#228b22','U-':'#db7093','G-':'#ff8c00','C-':'#4169e1'}
def plot100(traceboundary,sgst,signal,trace,viterbistart, raw_start, frombasecaller_idx ,possiblemove_idx):

    limit = 25000

    fig = plt.figure(figsize=(180, 15))
    signal = signal[viterbistart:viterbistart+limit]
    trace = trace[0:int(limit/10)]
    gs = gridspec.GridSpec(5, 1, height_ratios=[0.2, 0.2, 0.2,0.2,0.2])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax2)
    ax4 = fig.add_subplot(gs[3], sharex=ax3)
    ax5 = fig.add_subplot(gs[4], sharex=ax4)
    #
    # plot signal
    ax1.plot(signal, linewidth=1)
    ax2.plot(signal, linewidth=1)

    # plot trace
    traces = tuple(zip(*trace))
    for i in range(len(traces)):

        temp = [[trace] * 10 for trace in traces[i]]  # extend the index to match the index of signal
        trace = list(itertools.chain.from_iterable(temp))
        lstyle = "solid" if i <= 3 else "dashed"
        ax3.plot(trace, color=base_color[base_corresponding_table[i]], linestyle = lstyle, linewidth=1.0)
        ax4.plot(trace, color=base_color[base_corresponding_table[i]], linestyle=lstyle, linewidth=1.0)
        ax5.plot(trace, color=base_color[base_corresponding_table[i]], linestyle=lstyle, linewidth=1.0)


    for index in sgst:

        pos = index  - (viterbistart- raw_start)
        if pos > 0 and pos < limit:
            ax1.axvline(x=pos, ymin=0, ymax=1, color="black", alpha=0.2)

    for index in traceboundary:

        if index*10 < limit:
            ax2.axvline(x=index*10, ymin=0, ymax=1, color="black", alpha=0.2)

    for index in traceboundary:

        if index*10 < limit:
            ax3.axvline(x=index*10, ymin=0, ymax=1, color="black", alpha=0.2)

    for index in frombasecaller_idx:

        if index*10 < limit:
            ax4.axvline(x=index*10, ymin=0, ymax=1, color="black", alpha=0.2)

    for index in possiblemove_idx:

        if index*10 < limit:
            ax5.axvline(x=index*10, ymin=0, ymax=1, color="black", alpha=0.2)

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

    with h5py.File(fdata[0], 'r') as f:

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

        signal = signal[::-1] #reverse for rna
        signal = normalise(signal)
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

        #trace = trace.T
        # print("viterbi start",viterbistart,"rawstrt",raw_start)
        # print("end",end)
        # print("bases",seq)
        # print("fasta",fastq.split("\n")[1].replace("U","T"))
        #trace = trace.flatten()

        move = move[::-1].astype(np.int16)
        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,move,trace,signal,sgmean,sgst,sglen,bases,viterbistart,raw_start)
        return tp

    return None


def toSeq(fastq):

    fa = fastq.split("\n")
    return fa[1]


if __name__ == "__main__":

   path = '/data/nanopore/rRNA/1623_ivt-multi/singleFast5'
   pathout = '/data/nanopore/rRNA/test'
   ref ="/data/nanopore/reference/NC000913.fa"
   indexf = "/data/nanopore/rRNA/test/index.txt"
   path_w = pathout + "/sampleplingplan.txt"

   samplesize = 400

   makeParquet(indexf,ref, pathout, 8)

