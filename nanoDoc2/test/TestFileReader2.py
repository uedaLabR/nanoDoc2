from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2.utils.PqFileReader import PqReader
from matplotlib import pyplot as plt
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import pandas as pd


def getPlot(data,pos):


    plt.suptitle(str(pos))
    fig = plt.plot(data, linewidth=1)
    return fig

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

import pysam
def correctCigar(targetPos,cigar):

    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 0:  # match

            if refpos + cigarlen > targetPos:
                return relpos + (targetPos - refpos)

            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del
            refpos = refpos + cigarlen
        elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
            relpos = relpos + cigarlen

    return 0

from matplotlib import pyplot as plt
from matplotlib import gridspec
import itertools
base_corresponding_table = {0:'A',2:'G',1:'C',3:'T',4:'A-',6:'G-',5:'C-',7:'U-'}
base_color = {'A':'#228b22','T':'#db7093','G':'#ff8c00','C':'#4169e1','A-':'#228b22','U-':'#db7093','G-':'#ff8c00','C-':'#4169e1'}
from matplotlib import gridspec
def plotGraph(data):

    gm = GraphManager("/data/nanopore/nanoDoc2/test2")
    cnt = 0
    for index, row in data:

        chr = row['chr']
        start = row['start']
        end = row['end']
        cigar = row['cigar']
        print(start,end,cigar)
        genomeref = row['genome']
        offset = row['offset']
        traceintervals = row['traceintervals']
        trace = row['trace']
        signal = row['signal']
        trace = decode(trace)
        fig = plt.figure(figsize=(250, 10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5, 0.5])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        tracestart = traceintervals[0]
        idx = 0
        for ti in traceintervals:
            pos = ti-tracestart
            ax1.axvline(x=pos, ymin=0, ymax=15, color='black', alpha=0.2)

        for n in range(len(genomeref)):
            idx = correctCigar(n,cigar)
            base = genomeref[n]
            pos = traceintervals[idx] - tracestart
            ax1.text(pos, 16, base, size=7, color=base_color[base], ha='center')


        # cmap =['#228b22','#4169e1','#ff8c00','#db7093']
        # ax1.plot(trace,linewidth=1)
        for i in range(4):

            atrace = trace[:,i]
            ax1.plot(atrace, color=base_color[base_corresponding_table[i]],linewidth=1)


        ax2.plot(signal, linewidth=1)

        gm.add_figure(fig)
        cnt = cnt+1
        if cnt > 10:
            break

    gm.save()



if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    #path = '/data/nanopore/nanoDoc2/testSARSCOV2'
    fr = PqReader(path,4000)

    # pathout = '/data/nanopore/nanoDoc2/testfig/test.pdf'
    # output_file = '/data/nanopore/nanoDoc2/nanodoctest.pq'

    #gm = GraphManager(pathout)
    datab =[]
    #depth = fr.getRowData("cc6m_2244_t7_ecorv", True, 2000)
    pos = 300
    #traces,traceItvs,signals,takecnt = fr.getRowData("cc6m_2244_t7_ecorv", True, pos, takecnt=10)
    data = fr.getRowSequence("cc6m_2244_t7_ecorv", True, pos, takecnt=-1)

    plotGraph(data)

    # print(datab)
    # pschema = schema(
    #         [
    #             ('key', string()),
    #             ('signal', list_(float32()))
    #         ]
    # )
    # r = []
    # tp = ("ATGCA", datab)
    # r.append(tp)
    # print(r)
    # df = pd.DataFrame(r,
    #                   columns=['key', 'signal'])
    # pyarrow_table = Table.from_pandas(df, pschema)
    # pq.write_table(
    #     pyarrow_table,
    #     output_file,
    #     row_group_size=4000,
    #     compression='snappy',
    #     flavor=['spark'],
    # )



    #gm.save()