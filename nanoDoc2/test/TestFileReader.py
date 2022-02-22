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

from matplotlib import gridspec
def plotGraph(traces,traceItvs,signals):

    gm = GraphManager("/data/nanopore/nanoDoc2/test.pdf")
    l = min(len(traces),len(signals))

    print("traceItv",traceItvs)

    for n in range(l):
        print(n)
        fig = plt.figure(figsize=(40, 10))
        gs = gridspec.GridSpec(2, 1, height_ratios=[0.5,0.5])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        trace = traces[n]
        trace = decode(trace)
        # print("trace",trace)
        signal = signals[n]
        # print("signal", signal)
        # print(signal)
        traceitv = traceItvs[n]
        ax1.plot(trace, linewidth=1)
        tracestart = traceitv[0]
        for ti in traceitv:
            pos = ti-tracestart
            ax1.axvline(x=pos, ymin=0, ymax=15, color='black', alpha=0.2)

        ax2.plot(signal, linewidth=1)

        gm.add_figure(fig)

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
    pos = 1050
    traces,traceItvs,signals,takecnt = fr.getRowData("cc6m_2244_t7_ecorv", True, pos, takecnt=50)
    plotGraph(traces,traceItvs,signals)

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