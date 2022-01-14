
def toPqRead(read):

    mr = read.getMapResult()

    reference_id = mr.ctg
    strand = mr.strand
    r_st = mr.r_st
    r_en = mr.r_en
    q_st = mr.q_st
    q_en = mr.q_en
    cigar_old = mr.cigar

    rst = read.r_st_re
    cigar = read.cigar
    traceboundary = read.traceboundary
    traceseq = read.traceseq
    seq = read.reseq
    trace = read.trace
    signal = read.signal


    tp = (read.read_id,reference_id,strand,rst,cigar,seq,traceboundary,traceseq,r_st,r_en,q_st,q_en,cigar_old,read.fastq,trace,signal)
    return tp

def toCigar(traceback_path,seqlen):

    nprev = -1
    mprev = -1
    cnttotal = 0
    cigigarlist = []
    countdel = 0
    countmatchmismach = 0
    init = True
    for n, m in traceback_path:

        if init:
            if n>0:
                cigigarlist.append(str(n) + "N")
            init = False

        if m == mprev:  # delation
            countdel += 1
            if (countmatchmismach > 0):

                if (cnttotal + countmatchmismach) > seqlen:
                    mm = seqlen - cnttotal
                else:
                    mm = countmatchmismach
                cnttotal += mm
                cigigarlist.append(str(mm) + "M")

            countmatchmismach = 0

        else:  # match or mismatch
            countmatchmismach += 1
            if (countdel > 0):
                cigigarlist.append(str(countdel) + "D")
            countdel = 0
        mprev = m

    if (countmatchmismach > 0):
        if (cnttotal + countmatchmismach) > seqlen:
            mm = seqlen - cnttotal
        else:
            mm = countmatchmismach
        cnttotal += mm
        cigigarlist.append(str(mm) + "M")
    if (countdel > 0):
        cigigarlist.append(str(countdel) + "D")
    cigarstring = "".join(cigigarlist)

    return cigarstring