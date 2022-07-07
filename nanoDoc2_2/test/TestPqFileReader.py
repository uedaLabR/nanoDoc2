from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2_2.utils.Pq5merReader import PqReader
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
def plotGraph(signals,traces,path):

    gm = GraphManager(path)

    for n in range(len(signals)):

        fig = plt.figure(figsize=(40, 20))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 0.5])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.plot(signals[n])

        trace = traces[n]
        trace = decode(trace)
        trace = np.array(trace).T
        print(trace)

        i = 0
        for atrace in trace:
            ax2.plot(atrace, color=base_color[base_corresponding_table[i]], linewidth=1)
            i += 1

        gm.add_figure(fig)

    gm.save()


import mappy as mp
if __name__ == "__main__":

    path = '/data/nanoDoc2_2/varidate/CurlcakeIVT'
    ref = "/data/nanopore/reference/Curlcake.fa"
    fp = "/data/nanoDoc2_2/test.pdf"

    minreadlen = 200
    strand = True
    start = 200
    end = 500
    fr = PqReader(path, ref, minreadlen, strand, start, end, IndelStrict=True)
    chrom = "cc6m_2244_t7_ecorv"
    pos = 310
    signals,traces,sampledlen,infos = fr.getRowData(chrom, True, pos, takecnt=50)
    # print(signals)
    plotGraph(signals,traces,fp)


