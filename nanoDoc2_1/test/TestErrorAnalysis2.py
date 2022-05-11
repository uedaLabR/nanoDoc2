import glob
from sklearn.model_selection import train_test_split

import tensorflow as tf  # add
import numpy as np
from tensorflow.keras.layers import GlobalAveragePooling1D
import numpy as np
from tensorflow.keras import Model

from nanoDoc2.network import cnnwavenet_decfilter

DATA_LENGTH_UNIT = 60
DATA_LENGTH = 1024
from numba import jit,u1,i8,f8
from nanoDoc2_1.utils.PqFile6merReader import PqReader
from nanoDoc2_1.network import CnnWavenetDecDimention
from nanoDoc2.graph.GraphManager import GraphManager

def getModel():

    num_classes_org = 4078
    shape1 = (None, DATA_LENGTH, 1)
    model = CnnWavenetDecDimention.build_network(shape=shape1, num_classes=num_classes_org)
    flat = model.layers[-11].output
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
                    end):
    n = end
    idx = 0
    res = None #faiss.StandardGpuResources()
    strand = "-"
    while n > start+5:
        subs = seq[idx:idx + 6]
        cnt, cntref = eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr,
                                  model_t, fw, chrom,
                                  chromtgt,res)

        n = n - 1
        idx = idx + 1

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
    #print(tp)
    ar = np.array(tp)
    # a = ar[:, 0]
    # c = ar[:, 1]
    # g = ar[:, 2]
    # t = ar[:, 3]
    return ar

base_corresponding_table = {0:'A',2:'G',1:'C',3:'T',4:'A-',6:'G-',5:'C-',7:'U-'}
base_color = {'A':'#228b22','T':'#db7093','G':'#ff8c00','C':'#4169e1','A-':'#228b22','U-':'#db7093','G-':'#ff8c00','C-':'#4169e1'}
from matplotlib import gridspec
def plotGraph(traces,signals,infos,path):

    gm = GraphManager(path)

    for n in range(len(signals)):

        info = infos[n]
        if "1N" not in info:
            continue

        fig = plt.figure(figsize=(40, 20))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.5])

        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.plot(signals[n])

        trace = traces[n]
        trace = decode(trace)
        trace = np.array(trace).T
        #print(trace)

        print(info)
        plt.title(info)

        i = 0
        for atrace in trace:

            ax2.plot(atrace, color=base_color[base_corresponding_table[i]],linewidth=1)
            i +=1

        gm.add_figure(fig)

    gm.save()

def callPlusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom, chromtgt,
                   start,
                   end):
    strand = "+"
    poss = [4035552,4035546]
    res = None #faiss.StandardGpuResources()
    for n in range(start, end-6):

        if n in poss:

            subs = seq[(n - start):(n - start) + 6]
            rtraces, rawdatas, cnt, rinfos, retraces,refdatas, cntref,reinfos\
                = eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr,
                                      model_t, fw, chrom,
                                      chromtgt,res)
            if n == 4035546:
                print(len(rtraces),cnt)
                plotGraph(rtraces, rawdatas,rinfos,"/data/nanopore/nanoDoc2_1/test/wt.pdf")
                print(len(retraces), cntref)
                plotGraph(retraces,refdatas,reinfos,"/data/nanopore/nanoDoc2_1/test/ref.pdf")




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
    minv = 0.000001 # to prevent zero dev
    dists = np.clip(dists, minv, 3000)
    dists2 = np.clip(dists2, minv, 3000)

    return dists, dists2

def getDistOneside(a, b,res):

    d = a.shape[1]
    flat_config = faiss.GpuIndexFlatConfig()
    #index = faiss.GpuIndexFlatL2(res, d, flat_config)
    index = faiss.IndexFlatL2(a.shape[1])  # build the index
    index.add(a)  # add vectors to the index
    dists, result = index.search(b, k=10)  # actual search
    dists = np.mean(dists, axis=1)
    minv = 0.000001  # to prevent zero dev
    dists = np.clip(dists, minv, 3000)

    return dists

def lof(data):

    if len(data) < 5:
        return 0

    # LOF
    # print("lof",data)
    # print(data.shape)
    index = faiss.IndexFlatL2(data.shape[1])  # build the index
    index.add(data)  # add vectors to the index
    dists, result = index.search(data, k=5)  # actual search

    mincols = dists.argmin(axis=1)
    dists = np.delete(dists, mincols, axis=1)
    dists = np.mean(dists, axis=1)
    dists = np.clip(dists, 0, 5000)
    #
    score = np.mean(dists)

    return score


import os.path
import matplotlib.pyplot as plt
import nanoDoc2.utils.nanoDocUtils as nanoDocUtils
from scipy.stats import chi2_contingency
import pandas as pd
import math


def getScoreFromDistanceComparison(xref, xrow, res):

    xref1, xref2 = train_test_split(xref, test_size=(len(xref) // 2))
    xrow1, xrow2 = train_test_split(xrow, test_size=(len(xrow) // 2))

    dist = getDistOneside(xref1, xrow1, res)
    dist2 = getDistOneside(xref2, xrow1, res)
    dist_self = getDistOneside(xrow2, xrow1, res)

    ratio = dist/ dist_self
    ratio2 = dist2 / dist_self

    thres = 1.8
    #threslow = 1 / thres
    cnt = len(xref1)

    hiValueIndex = np.where(ratio > thres)[0]
    hiValueIndex2 = np.where(ratio2 > thres)[0]

    #lowValueIndex = np.where(ratio < threslow)[0]
    #lowValueIndex2 = np.where(ratio2 < threslow)[0]

    hivIdx = np.intersect1d(hiValueIndex, hiValueIndex2)
    #lowIdx = np.intersect1d(lowValueIndex, lowValueIndex2)

    score = len(hivIdx) / cnt

    #print(hiValueIndex, cnt, score)

    if cnt < 30:
        score = 0
    if score < 0:
        score = 0

    return score


def getScoreFromCluster(xref, xrow, niter):

    x = np.concatenate([xref, xrow])
    print(x.shape)
    nn, d = x.shape
    mindatasizeFor3centroid = 120
    scorenormalizefactor = 75

    k = 3 # k=3 for clustering one for unmod, one for mod, other
    if nn < mindatasizeFor3centroid: #if data size < 120 to small to 3 centroids
        k = 2

    kmeans = faiss.Kmeans(d=d, k=k, niter=niter)
    kmeans.train(x)
    # index = faiss.IndexFlatL2(16)
    D, I = kmeans.index.search(x, 1)
    Iref, Irow = np.split(I, 2)


    maxscore = 0
    for m in range(k):
        countref = np.count_nonzero(Iref == m)
        countrow = np.count_nonzero(Irow == m)

        if countref > 0 and countrow > 0:

            tc = (countref + countrow) // 2
            df = pd.DataFrame([[tc, tc], [countref, countrow]])
            chi2, p, dof, expected = chi2_contingency(df, correction=False)

            score = 0
            if abs(1-p) < 0.00001:
                score = 0
            elif p > 0:
                score = -1 * math.log10(p) / scorenormalizefactor
            print("ountref,countrow", countref, countrow,score)
            # print(p, score)
            if score > maxscore:
                maxscore = score

    if maxscore > 1:
        maxscore = 1

    maxscore = maxscore
    return maxscore

def getFormat(dlist):

    signal = np.array(dlist)
    signal = signal.astype('float16') / 255
    train_x = np.reshape(signal, (-1, DATA_LENGTH, 1))

    return train_x

def eachProcess(wfile, n, subs, strand, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw,
                chrom,
                chromtgt,res):

    weight_path = wfile + "/" +str(subs) + "/model_t_ep_2.h5"
    if not os.path.isfile(weight_path):
        #     no 6mer found
        print(weight_path)
        infos = "{0}\t{1}\t{2}\t{3}\t{4}".format(n, str(subs), 0, 0, 0, 0)
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (0, 0)

    try:
        # target signal
        rtraces, rawdatas, cnt, rinfos = targetpr.getRowData(chromtgt, strand, n, uplimit)
        # reference signal
        retraces,refdatas, cntref, reinfos = refpr.getRowData(chrom, strand, n, cnt)
    except:
        infos = "{0}\t{1}\t{2}\t{3}\t{4}".format(n, str(subs), 0, 0,  0)
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (0,0)

    #    reference start or end, or nodepth
    if (cnt < 5 or cntref < 5 or (rawdatas is None) or (refdatas is None)):
        infos = "{0}\t{1}\t{2}\t{3}\t{4}".format(n, str(subs), cnt, cnt,  0)
        print(infos)
        fw.writelines(infos + "\n")
        fw.flush()
        return (cnt, cntref)

    model_t.load_weights(weight_path)
    if cntref < cnt:
        rawdatas, cnt = nanoDocUtils.reducesize(rawdatas, cntref)  # raw data have the same size as reference data


    refdata = getFormat(refdatas)
    rowdata = getFormat(rawdatas)
    #
    #
    #k = 3
    niter = 20
    xref = model_t.predict(refdata)
    xrow = model_t.predict(rowdata)

    score = getScoreFromCluster(xref, xrow, niter)
    scoreDisplay = '{:.7f}'.format(score)

    # score2  = getScoreFromDistanceComparison(xref, xrow, res)
    # scoreDisplay2 = '{:.7f}'.format(score2)

    infos = "{0}\t{1}\t{2}\t{3}\t{4}".format(n, str(subs), cnt, cntref,scoreDisplay)
    # print(n-4037519)
    print(infos)
    fw.writelines(infos + "\n")
    fw.flush()
    return (rtraces, rawdatas, cnt, rinfos, retraces,refdatas, cntref,reinfos)


def modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 200):


    if chrom == "":
        chrom = nanoDocUtils.getFirstChrom(ref)
        chromtgt = chrom
        print("modcallinit", chrom)
    seq = nanoDocUtils.getSeq(ref, chrom, start, end, strand)
    if end < 0:
        end = len(seq)

    coeffA, coeffB, uplimitb, takeparcentile = nanoDocUtils.readParam(paramf)

    refpr = PqReader(refpq, ref,minreadlen,strand,IndelStrict=True) # Assume
    targetpr = PqReader(targetpq,ref,minreadlen,strand,IndelStrict=True)
    model_t = getModel()

    fw = open(out, mode='w')

    if strand == "-":
        callMinusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                        chromtgt,
                        start, end)
    else:
        callPlusStrand(wfile, coeffA, coeffB, uplimit, takeparcentile, seq, refpr, targetpr, model_t, fw, chrom,
                       chromtgt,
                       start, end)

    fw.close()


import sys

if __name__ == '__main__':

   # wfile = "/groups2/gac50430/nanopore/dataset4DL/weight5merm6A/"
   # paramf = "/groups2/gac50430/nanopore/shell/modcall/param.txt"
   # ref ="/groups2/gac50430/nanopore/reference/NC000913.fa"
   # refpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaIvt"
   # targetpq = "/groups2/gac50430/nanopore/equalbinnedpq/ecrRnaNative"
   # out = "/groups2/gac50430/nanopore/detection/ecoli/23S\m6aplus.txt"
   # chrom = "NC_000913.3"
   # start = 4037519+1600
   # end =  4037519+1730

    wfile = "/data/nanopore/nanoDoc2_1/weight/docweight"
    paramf = "/data/param20.txt"
    ref = "/data/nanopore/reference/NC000913.fa"
    refpq = "/data/nanopore/nanoDoc2_1/1623_ivt"
    targetpq = "/data/nanopore/nanoDoc2_1/1623_wt"
    # out = "/data/nanopore/rRNA/16S_test.txt"
    out = "/data/nanopore/nanoDoc2_1/error2.txt"
    chrom = "NC_000913.3"
    chromtgt = "NC_000913.3"

    start = 4035545
    end = 4035561
    strand = "+"


    minreadlen = 400
    #    if len(sys.argv) > 11 :
    #        minreadlen = int(sys.argv[11])
    chromtgt = chrom
    # # for covid analysis
    # if "england" in out:
    #     chromtgt = "hCoV-19/England/02/2020|EPI_ISL_407073"
    # if "austraria" in out:
    #     chromtgt = "MT007544.1"
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)

