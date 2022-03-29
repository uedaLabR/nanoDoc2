import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO
from operator import itemgetter, attrgetter
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import numpy as np
import mappy as mp
from nanoDoc2_1.utils.PqFileReader import PqReader
from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2_1.utils.PqFileReader import PqReader
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
def plotGraph(fmerlist,traces,signals):

    gm = GraphManager("/data/nanopore/nanoDoc2/traceview.pdf")

    for n in range(len(signals)):

        fig = plt.figure(figsize=(40, 20))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.5])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax1.plot(signals[n])
        fmer = fmerlist[n]
        trace = traces[n]
        trace = decode(trace)
        trace = np.array(trace).T
        plt.title(fmer)
        print(trace)

        i = 0
        for atrace in trace:

            ax2.plot(atrace, color=base_color[base_corresponding_table[i]],linewidth=1)
            i +=1

        gm.add_figure(fig)

    gm.save()



if __name__ == "__main__":

    pqf = "/data/nanopore/nanoDoc2_1/1200signal.pq"
    df = pq.read_table(pqf).to_pandas()
    labelidx = df["fmer"].unique().tolist()
    fmerlist = []
    traceslist = []
    signallist = []
    cnt = 0
    regcnt = 0
    for idx, row in df.iterrows():

        print(row)
        flg = labelidx.index(row[1])
        fmer = row[1]
        trace = np.array(row[2])
        signal = np.array(row[3])
        if cnt%240 == 0:
            fmerlist.append(fmer)
            traceslist.append(trace)
            signallist.append(signal)
            regcnt +=1
        cnt +=1
        if regcnt > 200:
            break

    plotGraph(fmerlist,traceslist, signallist)