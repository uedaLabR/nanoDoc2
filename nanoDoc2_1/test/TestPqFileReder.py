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
def plotGraph(traces,signals):

    gm = GraphManager("/data/nanopore/nanoDoc2/testSignalTrace.pdf")

    for n in range(len(signals)):

        fig = plt.figure(figsize=(40, 20))
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
        # id = ids[n]
        # print(n)
        # fig = plt.figure(figsize=(40, 20))
        # fig.suptitle(id)
        # gs = gridspec.GridSpec(4, 1, height_ratios=[0.25, 0.25,0.25, 0.25])
        # ax1 = fig.add_subplot(gs[0])
        # ax2 = fig.add_subplot(gs[1])
        # ax3 = fig.add_subplot(gs[2])
        # ax4 = fig.add_subplot(gs[3])
        #
        # trace = traces[n]
        # if len(trace) > 0:
        #     i = 0
        #     for atrace in trace:
        #
        #         ax1.plot(atrace, color=base_color[base_corresponding_table[i]],linewidth=1)
        #         i +=1
        #
        #     ax2.plot(signals[n])
        #
        #     i = 0
        #     unt_trace = untrimtraces[n]
        #     print(unt_trace)
        #     for atrace in unt_trace:
        #         ax3.plot(atrace, color=base_color[base_corresponding_table[i]], linewidth=1)
        #         i += 1
        #
        #     ax4.plot(untrimsignals[n])
        #     gm.add_figure(fig)

    gm.save()

import mappy as mp
if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2_1/CurlcakeIVT'
    ref = "/data/nanopore/reference/Curlcake.fa"


    wfile = "/data/nanopore/nanoDoc2_1/weight/docweight"
    paramf = "/data/param20.txt"
    ref = "/data/nanopore/reference/NC000913.fa"
    refpq = "/data/nanopore/nanoDoc2_1/1623_ivt"
    targetpq = "/data/nanopore/nanoDoc2_1/1623_wt"
    out = "/data/nanopore/nanoDoc2_1/error.txt"

    chrom = "NC_000913.3"
    chromtgt = "NC_000913.3"
    start = 4035531
    end = 4035561

    #path = '/data/nanopore/nanoDoc2/testSARSCOV2'
    fr = PqReader(targetpq,ref,True,50)
    # load or build index


    # pathout = '/data/nanopore/nanoDoc2/testfig/test.pdf'
    # output_file = '/data/nanopore/nanoDoc2/nanodoctest.pq'

    #gm = GraphManager(pathout)
    datab =[]
    #depth = fr.getRowData("cc6m_2244_t7_ecorv", True, 2000)
    pos = start
    a = mp.Aligner(ref)
    rseq = a.seq(chrom,start=pos,end=pos+6)
    print(rseq)


    traces,signals,sampledlen = fr.getRowData(chrom, True, pos, takecnt=50)
    # print(len(untrimtraces))
    # print(len(untrimsignals))
    plotGraph(traces,signals)

