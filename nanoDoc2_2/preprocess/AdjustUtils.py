
unit = 10
import pysam
def getrelpos(cigar,pos):
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 0:  # match

            if refpos + cigarlen > pos:
                return relpos + (pos - refpos)

            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del
            refpos = refpos + cigarlen
        elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
            relpos = relpos + cigarlen

    return 0

from statistics import mean
def getMeans(signal,traceboundary,cigar,seqlen):

    means = []
    for n in range(0, seqlen - 5):

        relpos = getrelpos(cigar,n)
        if relpos+1 < len(traceboundary):
            start = traceboundary[relpos] * unit
            end = traceboundary[relpos+1] * unit
            if end < len(signal):
                subsignal = signal[start:end]
                if len(subsignal) > 0:
                    means.append(mean(subsignal))


    return means


def theoryMean(fmerDict,lgenome,strand):

    means = []
    # rg = lgenome
    # for n in range(0, len(rg) - 5):
    #     fmer = rg[n:n + 5]
    #     if "N" in fmer:
    #         fmer = fmer.replace('N', 'A')
    #     cv = fmerDict[fmer]
    #     means.append(cv)

    if strand == "-":

        rg = lgenome
        for n in range(0, len(rg) - 5):
            fmer = rg[n:n + 5]
            if "N" in fmer:
                fmer = fmer.replace('N', 'A')
            cv = fmerDict[fmer]
            means.append(cv)

    else:

        #plus strand
        rg = lgenome[::-1]
        for n in range(0,len(rg)-5):

           fmer = rg[n:n+5]
           if "N" in fmer:
               fmer = fmer.replace('N', 'A')
           cv = fmerDict[fmer]
           means.append(cv)

        means.reverse()

    return means

def getCurrentDict(fmercurrent):
    a = {}
    with open(fmercurrent) as f:
        cnt = 0
        for line in f:
            if cnt > 0:
                data = line.split()
                a[data[0]] = float(data[1])
            cnt = cnt + 1
    return a
def countmatch(g,r):

    cnt = 0
    for m in range(5):
        if g[m] == r[m]:
            cnt += 1
    return cnt


def findMismatchInterval(gseq,readseq):

    ml = []
    l = min(len(gseq),len(readseq))
    for n in range(l-5):
        m = countmatch(gseq[n:n+5],readseq[n:n+5])
        if m <= 2:
            ml.append(n)
        elif readseq[n+2] =='-':
            ml.append(n)

    return interval_extract(ml)


def interval_extract(list):
    if len(list) ==0:
        return []

    list = sorted(set(list))
    range_start = previous_number = list[0]

    for number in list[1:]:
        if number == previous_number + 1:
            previous_number = number
        else:
            yield [range_start+1, previous_number+4]
            range_start = previous_number = number
    yield [range_start+1, previous_number+4]

# gseq =     "CAAGAAGAAGAAGGACGCTGGAAAGTCGGCCAAGAAAGACAAAGACCCAGTGAACAAATCCGGGGGCAAGGCCAAAAAGAAGAAGTGGTCCAAAGGCAAAGTTCGGGACAAGCTCAATAACTTAGTCTTGTTTGACAAAGCTACCTATGATAAACTCTGTAAGGAAGTTCCCAACTATAAACTTATAACCCCAGCTGTGGTCTCTGAGAGACTGAAGATTCGAGGCTCCCTGGCCAGGGCAGCCCTTCAGGAGCTCCTTAGTAAAGGACTTATCAAACTGGTTTCAAAGCACAGAGCTCAAGTAATTTACACCAGAAATACCAAGGGTGGAGATGCTCCAGCTGCTGGTGAAGATGCATGAATAGGTCCAACCAGCTGTACATTTGGAAAAATAAAACT"
# traceseq = "TAAGAAGAAGAAGGACGCTGGAAAGTCGGCCAAGAAAGACAAAGACCCAGTGAACAAATCCGGGGGCAAGGCCAAAAAGAAGAAGTGGTCCAAAGGCAAAGTTCGGGACAAGCTCAATAACTTAGTTTCGTTTGACAAAGCTACCTATGATAAACTCTGTTGGGAAGTTCCCAACTATAAACTTATAACCCCAGCTGTGGTCTTTGAGAGACTGAAGATTCGAA---CCCTGGCCAGGGCAGCCCTTCAGGAGCTCCTTAGTAAAGGACCTCCCAAACTGGTTTCAAAGCACAGAGCTCAAGTAATTCATACCAGAAATTCTAAGGGTGGAGATGCTCCAGCTGCTGGTGAAGATGCATGA-TAGGTCCAACCAGCTGTACATTTGGAAAAATAAAAC"

gseq =     "TGCATGAATAGGTCCAACCAGCTGTACATTTGGAAAAATAAAAC"
traceseq = "TGCATGA-TAGGTCCAACCAGCTGTACATTTGGAAAAATAAAAC"


itvl = findMismatchInterval(gseq, traceseq)

for i in itvl:
    print(i)

print(itvl)

