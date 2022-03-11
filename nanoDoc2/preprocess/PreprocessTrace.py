import nanoDoc2.preprocess.ViterbiSegmentation  as vs
import nanoDoc2.preprocess.SignalNormalization as ss


def preprocess(read,fmerDict):

    #Viterbi segmentation
    lgenome = read.refgenome
    chrom = read.chrom
    strand = read.strand
    r_st = read.r_st
    r_en = read.r_en
    q_st = read.q_st
    q_en = read.q_en
    trace = read.trace
    move = read.move

    seq, cigar, left,traceoffset, traceboundary, frombasecaller_idx,possiblemove_idx = vs.flipplopViterbiEach(lgenome,
                                                                                                  chrom,
                                                                                                  strand,
                                                                                                  r_st,
                                                                                                  r_en,
                                                                                                  q_st,
                                                                                                  q_en,
                                                                                                  trace,
                                                                                                  move)

    # print(traceboundary)
    # print(len(read.trace))

    read.cigar_str = cigar
    read.settraceboundary(traceboundary)
    read.normSignal = ss.normalizeSignal(read, traceboundary, fmerDict)
    return read

