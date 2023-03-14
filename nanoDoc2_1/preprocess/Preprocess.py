import nanoDoc2_1.preprocess.ViterbiSegmentation  as vs
import nanoDoc2_1.preprocess.SignalNormalization as ss

# import pysam
# def devideRange(cigar, r_st, r_en, q_st, q_en, strand,devidesize = 400,min=200):
#
#     a = pysam.AlignedSegment()
#     a.cigarstring = cigar
#     refpos = 0
#     relpos = 0
#     cidx = 0
#     mintilelen = 25
#     cl = []
#     for cigaroprator, cigarlen in a.cigar:
#
#         if cigaroprator == 0 and cigarlen > mintilelen:
#
#            cl.append()
#
#
#         if cigaroprator == 3 and cidx > 0:  # N
#
#             refpos = refpos + cigarlen
#
#         elif cigaroprator == 0 or cigaroprator == 4:  # match or S softclip was not correted so treat as M
#
#             if refpos + cigarlen > targetPos:
#                 return relpos + (targetPos - refpos)
#
#             relpos = relpos + cigarlen
#             refpos = refpos + cigarlen
#
#         elif cigaroprator == 2:  # Del
#
#             refpos = refpos + cigarlen
#
#         elif cigaroprator == 1:  # Ins
#
#             if relpos == 0:
#                 if targetPos <= cigarlen:
#                     return 0
#
#             relpos = relpos + cigarlen
#
#
#





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
    print("cigar1",orgcigar )

    # ret = devideRange(orgcigar, r_st, r_en, q_st, q_en, devidesize = 400,min=200)
    # if ret == None:

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

    # else:
    #
    #     for r_st,r_en,q_st,q_en in ret:
    #
    #         seq, cigar, left, traceoffset, traceboundary, frombasecaller_idx, possiblemove_idx = vs.flipplopViterbiEach(
    #             lgenome,
    #             chrom,
    #             strand,
    #             orgcigar,
    #             r_st,
    #             r_en,
    #             q_st,
    #             q_en,
    #             trace,
    #             move)





    print(traceboundary)
    # print(len(read.trace))

    read.cigar_str = cigar
    read.settraceboundary(traceboundary)
    read.normSignal = ss.normalizeSignal(read, traceboundary, fmerDict)
    return read

