import nanoDoc2_1.preprocess.ViterbiSegmentation  as vs
import nanoDoc2_1.preprocess.SignalNormalization as ss

import nanoDoc2_2.preprocess.SignalAdjust as sn
def preprocess(read,fmerDict):

    #Viterbi segmentation
    lgenome = read.refgenome
    chrom = read.chrom
    strand = read.strand
    orgcigar = read.cigar_str
    r_st = read.r_st
    r_en = read.r_en
    q_st = read.q_st
    q_en = read.q_en
    trace = read.trace
    move = read.move
    # print("cigar1",orgcigar )
    print("flipplopViterbi")
    seq, cigar, left,traceoffset, traceboundary, frombasecaller_idx,possiblemove_idx = vs.flipplopViterbiEach(lgenome,
                                                                                                  chrom,
                                                                                                  strand,
                                                                                                  orgcigar,
                                                                                                  r_st,
                                                                                                  r_en,
                                                                                                  q_st,
                                                                                                  q_en,
                                                                                                  trace,
                                                                                                  move)


    # print(len(read.trace))

    read.cigar_str = cigar
    read.cigar_org = cigar
    read.settraceboundary(traceboundary)
    print("normalize")
    read.normSignal = ss.normalizeSignal(read, traceboundary, fmerDict)
    #print("adjust start")
    print("adjust start")
    read.signalboundary = sn.adjustMismatchindel(read,fmerDict)
    # print(len(read.traceboundary))
    # print(len(read.signalboundary))

    print("adjust end")
    return read

