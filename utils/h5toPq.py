import h5py
from statsmodels import robust
import numpy as np
import random
from scipy import interpolate
import pyarrow as pa
import pyarrow.parquet
import pandas as pd
# import mappy as mp
import tqdm
from functools import cmp_to_key
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial

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
    for tp in blist:

       _makeParquet(tp,aligner)

    # graph_manager.save()

    # with Pool(thread) as p:
    #      p.map(_makeParquet, blist)


#        tqdm(imap,total=size)


def _makeParquet(tp):

    (indexcontent, pathout, ref) = tp


    wholedata = []

    cnt_1 = 0
    for s_line in indexcontent:
        fdata = s_line.split(",")
        onerecord = getRecord(fdata)


        if onerecord != None:


            wholedata.append(onerecord)
            cnt_1 = cnt_1 + 1
            if cnt_1 % 20 == 0:
                print(cnt_1)
                #break
        # (m
        #
        # apped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, bytedata,binsizelist)

    # tp = (
    #     mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, chrom,
    #     r_st, r_en,
    #     cigar, bytedata, tracedata, binsizelist)

    # tp = (
    # mapped_chrom, mapped_strand, mapped_start, mapped_end, clipped_bases_start, clipped_bases_end, fastq, ctg, st, ed,
    # cigar, md, trace, move, signal, allevent)
    df = pd.DataFrame(wholedata,
                      columns=['mapped_chrom', 'mapped_strand', 'mapped_start', 'mapped_end', 'clipped_bases_start',
                               'clipped_bases_end', 'fastq', 'chrom', 'st', 'ed', 'cigar','md', 'move','trace','signal',
                               'sgmean', 'sgst', 'sglen', 'bases'
                               ])
    df.to_parquet(pathout)






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

        #print(allevent.shape)
        sgmean = []
        sgst = []
        sglen = []
        bases = []
        for n in range(len(allevent)):

            mean = allevent[n][0]
            start = allevent[n][2]
            size = allevent[n][3]
            b = allevent[n][4]

            sgmean.append(mean)
            sgst.append(start)
            sglen.append(size)
            bases.append(b)


        group = f['/Raw/Reads/']
        readid = list(group.values())[0].name
        signalkey = readid + "/Signal"
        signal = f[signalkey][()]

        signal = signal[::-1] #reverse for rna
        signal = signal[raw_start:end] #upto tombo segmantation
        signal = normalise(signal)

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
        trace = trace.T
        trace = trace.flatten()

        move = move[::-1].astype(np.int16)
        tp = (mapped_chrom,mapped_strand,mapped_start,mapped_end,clipped_bases_start,clipped_bases_end,fastq,ctg,st,ed,cigar,md,move,trace,signal,sgmean,sgst,sglen,bases)
        return tp

    return None


def toSeq(fastq):

    fa = fastq.split("\n")
    return fa[1]




if __name__ == "__main__":

   path = '/share/data/IVT/m6AIvt/workspace'
   pathout = '/data/nanopore/IVT/pqtest'
   ref = "/share/reference/Curlcake.fa"
   indexf = "/data/nanopore/IVT/pqtest/index.txt"
   path_w = pathout + "/sampleplingplan.txt"

   samplesize = 400

   makeParquet(indexf,ref, pathout, 8)

