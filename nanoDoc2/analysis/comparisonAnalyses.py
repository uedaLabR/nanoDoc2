import glob
from sklearn.model_selection import train_test_split

import tensorflow as tf  # add
import numpy as np
from tensorflow.keras.layers import GlobalAveragePooling1D
import numpy as np
from tensorflow.keras import Model

from nanoDoc2.network import cnnwavenet_decfilter

DATA_LENGTH_UNIT = 60
DATA_LENGTH = 40
from numba import jit,u1,i8,f8
from nanoDoc2.utils.PqTraceFileReader import PqReader

def getModel():

    shape1 = (None, DATA_LENGTH, 4)
    num_classes_org = 1024
    # with tf.device("/cpu:0"):
    #model = cnnwavenettrace.build_network(shape=shape1, num_classes=num_classes_org)
    model = cnnwavenet_decfilter.build_network(shape=shape1, num_classes=num_classes_org)
    flat = GlobalAveragePooling1D()(model.layers[-3].output)
    model_t = Model(inputs=model.input, outputs=flat)
    model_t.summary()
    return model_t


def getSD(cnt, coeffA, coeffB):
    # y=b*(x**a) to approximate std value for the score
    sd = coeffB * (cnt ** coeffA)
    minv = 0.0001
    if sd < minv:
        sd = minv
    return sd


def callMinusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt,
                    start,
                    end, refpq, targetpq, minreadlen):
    n = end
    idx = 0
    res = None #faiss.StandardGpuResources()
    strand = "-"
    while n > start:
        subs = seq[idx:idx + 5]
        cnt, cntref = eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr,
                                  model_t, fw, chrom,
                                  chromtgt,res)

        n = n - 1
        idx = idx + 1


def callPlusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt,
                   start,
                   end, refpq, targetpq, minreadlen):
    strand = "+"
    res = None #faiss.StandardGpuResources()
    for n in range(start, end-5):
        subs = seq[(n - start):(n - start) + 5]
        cnt, cntref = eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr,
                                  model_t, fw, chrom,
                                  chromtgt,res)


def getFormat(traces):

    traces = np.array(traces)
    traces = traces/255
    DATA_LENGTH_UNIT = 8 * 5
    traces = np.reshape(traces, (-1, DATA_LENGTH_UNIT, 4))

    return traces


import faiss

def getDist(a, b,res):

    d = a.shape[1]
    #flat_config = faiss.GpuIndexFlatConfig()
    #index = faiss.GpuIndexFlatL2(res, d, flat_config)   # build the index
    index = faiss.IndexFlatL2(a.shape[1])
    index.add(a)  # add vectors to the index
    dists, result = index.search(b, k=10)  # actual search

    #index = faiss.GpuIndexFlatL2(res, d, flat_config)   # build the index
    index = faiss.IndexFlatL2(a.shape[1])
    index.add(b)  # add vectors to the index
    dists2, result2 = index.search(a, k=10)  # actual search

    dists = np.mean(dists, axis=1)
    dists2 = np.mean(dists2, axis=1)

    dists = np.clip(dists, 0, 3000)
    dists2 = np.clip(dists2, 0, 3000)

    return dists, dists2

def getDistOneside(a, b,res):

    d = a.shape[1]
    flat_config = faiss.GpuIndexFlatConfig()
    #index = faiss.GpuIndexFlatL2(res, d, flat_config)
    index = faiss.IndexFlatL2(a.shape[1])  # build the index
    index.add(a)  # add vectors to the index
    dists, result = index.search(b, k=10)  # actual search
    dists = np.mean(dists, axis=1)
    dists = np.clip(dists, 0, 3000)

    return dists

import os.path
import matplotlib.pyplot as plt
import nanoDoc2.utils.nanoDocUtils as nanoDocUtils
def eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw,
                chrom,
                chromtgt,res):

    weight_path = wfile + "/" +str(subs) + "/model_t_ep_1.h5"
    if not os.path.isfile(weight_path):
        #     no 6mer found
        print(weight_path)
        infos = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(n, str(subs), 0, 0, 0, 0)
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (0, 0)

    # target signal
    rawdatas, cnt = targetpr.getRowData(chromtgt, strand, n, uplimit)
    # reference signal
    refdatas, cntref = refpr.getRowData(chrom, strand, n, cnt * 3)

    #    reference start or end, or nodepth
    if (cnt < 10 or cntref < 10 or (rawdatas is None) or (refdatas is None)):
        infos = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(n, str(subs), cnt, cnt, 0, 0)
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (cnt, cntref)

    model_t.load_weights(weight_path)
    rawdatas, cnt = nanoDocUtils.reducesize(rawdatas, cntref // 2)  # raw data have at least half size as reference data

    refdata1, refdata2 = train_test_split(refdatas, test_size=(cntref // 3))
    refdata1, refdata3 = train_test_split(refdata1, test_size=(cntref // 3))
    #ref for stat comparison
    refdata1 = getFormat(refdata1)
    refdata2 = getFormat(refdata2)
    refdata3 = getFormat(refdata3)

    rowdata = getFormat(rawdatas)

    xref = model_t.predict(refdata1)
    xref2 = model_t.predict(refdata2)
    xref3 = model_t.predict(refdata3)
    xrow = model_t.predict(rowdata)

    dist = getDistOneside(xrow, xref, res)
    distref = getDistOneside(xref2, xref, res)
    distref2 = getDistOneside(xref3, xref, res)
    ratioraf = distref2/distref

    querythres = np.percentile(ratioraf, 95)
    hiQueryValueIndex = np.where(ratioraf > querythres)

    ratio = dist/distref
    thres = np.percentile(ratioraf, 99)
    thres = 1.38


    hiValueIndex = np.where(ratio > thres)
    diffset = np.setdiff1d(hiValueIndex,hiQueryValueIndex)
    # countmod = np.count_nonzero(ratio > thres)
    score = len(diffset) / (cnt * 0.95)

    if cnt < 30:
        score = 0
    if score < 0:
        score = 0

    scoreDisplay = '{:.7f}'.format(score)
    infos = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(n, str(subs), cnt, cntref,thres,scoreDisplay)

    print(infos)
    fw.writelines(infos + "\n")
    fw.flush()
    return (cnt, cntref)


def modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 2000):


    if chrom == "":
        chrom = nanoDocUtils.getFirstChrom(ref)
        chromtgt = chrom
        print("modcallinit", chrom)
    seq = nanoDocUtils.getSeq(ref, chrom, start, end, strand)
    if end < 0:
        end = len(seq)

    coeffA, coeffB, uplimitb, takeparcentile = nanoDocUtils.readParam(paramf)

    refpr = PqReader(refpq, ref)
    targetpr = PqReader(targetpq,ref)
    model_t = getModel()

    fw = open(out, mode='w')

    if strand == "-":
        callMinusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                        chromtgt,
                        start, end, refpq, targetpq, minreadlen)
    else:
        callPlusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                       chromtgt,
                       start, end, refpq, targetpq, minreadlen)

    fw.close()


import sys

if __name__ == '__main__':

    #    wfile = "/groups2/gac50430/nanopore/dataset4DL/weight5merm6A/"
    #    paramf = "/groups2/gac50430/nanopore/shell/modcall/param.txt"
    #    ref ="/groups2/gac50430/nanopore/reference/NC000913.fa"
    #    refpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaIvt"
    #    targetpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaNative"
    #    out = "/groups2/gac50430/nanopore/detection/ecoli/23S\m6aplus.txt"
    #    chrom = "NC_000913.3"
    #    start = 4037519+1600
    #    end =  4037519+1730

    # wfile = "/data/nanopore/IVT/doc_weight_bk2"
    # paramf = "/data/param20.txt"
    # ref = "/data/nanopore/reference/NC000913.fa"
    # refpq = "/data/nanopore/rRNA/1623_ivtpq"
    # targetpq = "/data/nanopore/rRNA/1623_nativepq"
    # # out = "/data/nanopore/rRNA/16S_test.txt"
    # out = "/data/nanopore/rRNA/23S_test.txt"
    # chrom = "NC_000913.3"
    # chromtgt = "NC_000913.3"
    # # start = 4035531
    # # end = start+1541
    # start = 4037519
    # end = 4040423

    #    modCall(wfile,paramf,ref,refpq,targetpq,out,chrom,start,end)
    wfile = sys.argv[1]
    paramf = sys.argv[2]
    ref = sys.argv[3]
    refpq = sys.argv[4]
    targetpq = sys.argv[5]
    out = sys.argv[6]
    chrom = sys.argv[7]
    start = int(sys.argv[8])
    end = int(sys.argv[9])
    strand = sys.argv[10]
    #    minreadlen = 700
    minreadlen = 200
    #    if len(sys.argv) > 11 :
    #        minreadlen = int(sys.argv[11])
    chromtgt = chrom
    # # for covid analysis
    # if "england" in out:
    #     chromtgt = "hCoV-19/England/02/2020|EPI_ISL_407073"
    # if "austraria" in out:
    #     chromtgt = "MT007544.1"
    modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)

