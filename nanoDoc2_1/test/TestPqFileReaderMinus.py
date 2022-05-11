from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2_1.utils.PqFile6merReader import PqReader
from matplotlib import pyplot as plt
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import pandas as pd
import mappy as mp
from numba import jit

def getPlot(data,pos):


    plt.suptitle(str(pos))
    fig = plt.plot(data, linewidth=1)
    return fig

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
def plotGraph(traces,signals,rseq,infos,fp):

    gm = GraphManager(fp)

    for n in range(len(signals)):

        fig = plt.figure(figsize=(40, 20))
        info = infos[n]
        plt.title(info + " " + str(rseq))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.5])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.plot(signals[n])

        trace = traces[n]
        trace = decode(trace)
        trace = np.array(trace).T
        print(trace)

        i = 0
        for atrace in trace:

            ax2.plot(atrace, color=base_color[base_corresponding_table[i]],linewidth=1)
            i +=1

        gm.add_figure(fig)


    gm.save()

import mappy as mp
if __name__ == "__main__":

    # ref = "/data/nanopore/reference/Yeast_sk1.fa"
    # refpq = "/data/nanopore/nanoDoc2_1/1825_native"
    # out = "/data/nanopore/nanoDoc2_1/18S_test.txt"
    # #path = '/data/nanopore/nanoDoc2/testSARSCOV2'
    # #margin = 1000
    # minreadlen = 10
    # strand = False
    # chrom = "chr12"
    # chromtgt = "chr12"
    # start = 455938
    # end = 457732

    start = 451786
    end = 455181

    wfile = "/data/nanopore/nanoDoc2_1/weight/docweight"
    paramf = "/data/param20.txt"
    ref = "/data/nanopore/reference/NC000913.fa"
    refpq = "/data/nanopore/nanoDoc2_1/1623_ivt"
    targetpq = "/data/nanopore/nanoDoc2_1/1623_wt"
    out = "/data/nanopore/nanoDoc2_1/error.txt"

    chrom = "NC_000913.3"
    chromtgt = "NC_000913.3"
    minreadlen = 200
    strand = True

    start = 4035570
    end = 4035580

    refpr = PqReader(refpq, ref,minreadlen,strand,start,end,IndelStrict=True) # Assume

    a = mp.Aligner(ref)
    poss = [4035575+1249,4035575+1270,4035575+1290]
    # = range(start,end)
    for pos in poss:
        rseq = a.seq(chrom, start=pos-6, end=pos)
        rseq = mp.revcomp(rseq)
        print(rseq)
        # depth = refpr.getDepth(chrom,pos,strand)
        # print(pos,depth,rseq)

        # load or build index
        traces, signals, sampledlen,infos = refpr.getRowData(chrom, True, pos,takecnt=100)
        # print(traces)
        # print(signals)

        print(sampledlen)

        fp = "/data/nanopore/nanoDoc2_1/figs/testSignalTraceminus"+str(pos)+".pdf"
        plotGraph(traces, signals,rseq,infos,fp)

