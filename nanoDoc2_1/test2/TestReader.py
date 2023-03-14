from nanoDoc2_1.test2.PqTraceSeqReader import TraceSeqReader
from matplotlib import pyplot as plt
from nanoDoc2.graph.GraphManager import GraphManager
from numba import jit

def plotGraph(traces,infos,pathout,seq):

    gm = GraphManager(pathout)

    for n in range(len(traces)):

        fig = plt.figure(figsize=(40, 20))
        value = traces[n]
        info = str(infos[n])
        print(value)
        plt.title(seq+" "+info)
        plt.plot(value)
        #
        gm.add_figure(fig)

        if n == 100:
            break


    gm.save()

import mappy as mp
import pysam
if __name__ == "__main__":


    print("start")
    path = "/data/nanopore/nanoTune/test/rna_pq_ck"
    start = 1
    end = 1500
    chrom = "cc6m_2709_t7_ecorv"
    ref = "/data/nanopore/reference/Curlcake.fa"

    #path = '/data/nanopore/nanoDoc2/testSARSCOV2'
    #fr = PqReader(targetpq,ref,True,50)
    minreadlen = 200
    strand = True
    findexs = [0,1,2,3,4,5,6,7,8,9,10]
    maxreads = 500

    fasta = pysam.Fastafile(ref)
    pos = 100
    seq = fasta.fetch(chrom, pos, pos+6).upper()
    print(seq)
    seqall = fasta.fetch(chrom, 1, 2000).upper()
    # print(seqall)
    # sets = set()
    # for n in range(1,4500):
    #     smer = seqall[n-1:n+5]
    #     if smer == "CGAGAC":
    #         print(smer,n)
    #     # sets.add(smer)
    #     #
    #     # if smer == seq:
    #     #     print(n)
    # CGAGAC
    # 1
    # CGAGAC
    # 1119
    # CGAGAC
    # 1252

    fr = TraceSeqReader(chrom,strand,start, end, path,findexs,maxreads,minreadlen)
    # load or build index
    traces,infos = fr.getRowData(pos,takecnt=maxreads)
    # print(len(traces))
    # print(traces)
    pathout = "/data/nanopore/nanoTune/test/test_1.pdf"
    plotGraph(traces,infos,pathout,seq)

    # traces,infos = fr.getRowData(1252,takecnt=maxreads)
    # # print(len(traces))
    # # print(traces)
    # pathout = "/data/nanopore/nanoTune/test/test2.pdf"
    # plotGraph(traces,infos,pathout,seq)




