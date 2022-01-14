import csv
import py2bit
from Bio.Seq import Seq
import pysam

binsize = 10**7

def toIntlist(s):
    l = s.rstrip(',').split(",")
    return list(map(int, l))


def dB(dbref):
    dc = {}
    with open(dbref, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        cnt = 0
        for row in reader:
            id = row[0]
            chr = row[1]
            strand = row[2]
            start = int(row[3])
            end = int(row[4])
            numexon = row[7]
            es = toIntlist(row[8])
            ee = toIntlist(row[9])
            # print((id, chr, strand, numexon, es, ee))
            dc[id] = (chr, strand, start, end, numexon, es, ee)
    return dc


def dBbyChr(dbref):
    dc = {}
    with open(dbref, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        cnt = 0
        for row in reader:
            id = row[0]
            chr = row[1]
            strand = row[2]
            start = int(row[3])
            end = int(row[4])
            numexon = row[7]
            es = toIntlist(row[8])
            ee = toIntlist(row[9])

            sl = dc.get(chr)
            if sl is None:
                sl = []
                dc[chr] = sl
            tp = (chr, strand, start, end, numexon, es, ee)
            sl.append(tp)

    return dc

def cluster(dbref):
    with open(dbref, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        cnt = 0
        chrprev = ""
        plusTs = []
        minusTs = []
        for row in reader:
            id = row[0]
            chr = row[1]
            strand = row[2]
            start = int(row[3])

            end = int(row[4])
            tp = (chr, start, end)
            if strand == "+":
                plusTs.append(tp)
            else:
                minusTs.append(tp)
        plusCluster = toCluster(plusTs)
        minusCluster = toCluster(minusTs)

    return (plusCluster, minusCluster)


def toCluster(tslist):

    clusterDict = {}
    cluster = []

    tslist = sorted(tslist, key=lambda k: (k[0], k[1]))  # sort by chr and start

    for chr, start, end in tslist:

        if len(cluster) == 0:
            cluster.append((chr, start, end))
            continue

        #
        lastidx = len(cluster) - 1
        prev_chr, prev_start, prev_end = cluster[lastidx]

        if prev_chr != chr:
            clusterDict[prev_chr] = convertF(cluster)
            cluster = []
            cluster.append((chr, start, end))
            continue

        if intercect(prev_chr, prev_start, prev_end, chr, start, end):
            # expand
            cluster[lastidx] = (prev_chr, min(prev_start, start), max(prev_end, end))
        else:
            cluster.append((chr, start, end))

    clusterDict[chr] = convertF(cluster)
    return clusterDict

def getBinKey(strand, chr,r_st, r_en,clusterAnnotation):


    if clusterAnnotation is None:
        return keyByBin(strand, chr, r_st, r_en)

    if strand == 1:
        cl = clusterAnnotation[0]
    else:
        cl = clusterAnnotation[1]

    if chr not in cl:
        return None
    (startlist,endlist) = cl[chr]


    cpos = searchCluster(r_st, r_en,startlist,endlist)
    if cpos is not None:

        key = "ts" + ":" + chr + ":" + str(strand) + ":" + str(cpos[0]) + ":" + str(cpos[1])

    else:
        key = keyByBin(strand, chr, r_st, r_en)

    return key

import bisect
def searchCluster(r_st, r_en,startlist,endlist):

    indexs = bisect.bisect_left(startlist, r_st)-1
    indexe = bisect.bisect_left(startlist, r_en) - 1

    if indexs < 0 and indexe < 0 :
        return None ##no hit


    # take most overlapped transcript definition
    lidx = 0
    maxoverlap = 0
    maxstart = 0
    maxend = 0

    if indexe+1 < len(startlist):
        indexe+=1


    st = 0
    for i in range(indexs,indexe):
        st =  startlist[i]
        ed = endlist[i]
        overl = overlap(r_st, r_en,st,ed)
        if overl > maxoverlap:
            maxoverlap = overl
            maxstart = st
            maxend = ed

        lidx = lidx +1
    if maxoverlap <= 0:
        return None
    return (maxstart,maxend)

def overlap(a,b,c,d):
    r = 0 if a==c and b==d else min(b,d)-max(a,c)
    if r>=0:
        return r
    else:
        return -1


def keyByBin(strand, chr,  r_st, r_en):

    genomebin = r_st // binsize
    return "g"+":"+chr + ":" + str(strand)+":" +str(genomebin)

def binToInterval(bin):

    return (binsize*bin,binsize*(bin+1))

def convertF(cluster):
    stlist = []
    edlist = []

    for chr, start, end in cluster:
        stlist.append(start)
        edlist.append(end)

    return (stlist, edlist)


def intercect(chr, st, en, chrts, st_ts, en_ts):
    if chr == chrts:
        if st <= en_ts and st_ts <= en:
            return True

    return False


def convertAndAdd(tsDb, hit, hitgenome, convertCigar):
    appended = False
    transcriptID = hit[3]
    if transcriptID in tsDb:

        tsCoord = tsDb[transcriptID]
        (chrts, strand, start, end, numexon, es, ee) = tsCoord

        for (hitG, tslist, convertlist) in hitgenome:

            (q_st, q_en, strand, chr, ctg_len, r_st, r_en, mlen, blen, mapq, primary, NM, cigar) = hitG
            if intercect(chr, r_st, r_en, chrts, start, end) and len(tslist) < 2:  # max 2 optimal hit
                tslist.append(hit)
                if convertCigar:
                    convertHit = convertCoordinate(hit, tsCoord)
                    convertlist.append(convertHit)
                appended = True

    return appended


def toTp(hit):
    return (hit.q_st, hit.q_en, hit.strand, hit.ctg, hit.ctg_len, hit.r_st, hit.r_en, hit.mlen, hit.blen, hit.mapq,
            hit.is_primary, hit.NM, hit.cigar)


def convertCoordinate(hit, tsCoord):
    (q_st, q_en, strand, ctg, ctg_len, r_st, r_en, mlen, blen, mapq, primary, nm, cigar) = hit
    (chrts, strand, start, end, numexon, es, ee) = tsCoord
    #
    st_ts = convertCoord(r_st, tsCoord)
    en_ts = convertCoord(r_en, tsCoord)
    cigar_ts = convertCigar(r_st, cigar, tsCoord)

    hitConvert = (q_st, q_en, strand, chrts, ctg_len, st_ts, en_ts, mlen, blen, mapq, primary, nm, cigar_ts)
    return hitConvert


def convertCigar(r_st, cigar, tsCoord):
    (chrts, strand, start, end, numexon, es, ee) = tsCoord

    if len(es) <= 1:
        return cigar

    idx = 0
    localpos = 0
    introns = []

    for n in range(len(es) - 1):
        s0 = es[n]
        e = ee[n]
        s = es[n + 1]
        exonlen = e - s0
        intronlen = s - e
        localpos = localpos + exonlen
        introns.append((localpos, intronlen))

    cigarOparatorN = 3
    cigarNew = []

    while len(introns) > 0:

        cigarNew = []
        (localpos, intronlen) = introns.pop(0)

        lpos = r_st
        if localpos < lpos:
            continue

        addintron = False
        for cginfo in cigar:

            cigarLen, cigarOparator = cginfo[0], cginfo[1]

            if cigarOparator > 3:
                print("CIGAR=" + cigar)

            read_consuming_ops = cigarOparator <= 1  # M or I
            if read_consuming_ops and (addintron == False):

                if lpos + cigarLen < localpos:
                    if cigarOparator == 0:  # match or mismatch
                        lpos = lpos + cigarLen
                    cigarNew.append([cigarLen, cigarOparator])
                else:

                    first = localpos - lpos + 1
                    last = cigarLen - first
                    if (first > 0):
                        cigarNew.append([first, cigarOparator])
                        cigarNew.append([intronlen, cigarOparatorN])  # N

                    if (last > 0):
                        cigarNew.append([last, cigarOparator])

                    lpos = lpos + cigarLen
                    addintron = True
                    # print("add intron",lpos,intronlen, cigarOparatorN)

            else:

                if cigarOparator == 2 and (addintron == False):  # Deletion
                    if lpos + cigarLen < localpos:

                        cigarNew.append([cigarLen, cigarOparator])
                        localpos = localpos - cigarLen

                    else:

                        first = localpos - lpos + 1
                        cigarNew.append([intronlen + first, cigarOparatorN])
                        addintron = True

                else:
                    cigarNew.append([cigarLen, cigarOparator])

        cigar = cigarNew

    return cigarNew


def convertCoord(pos, tsCoord):
    (chrts, strand, start, end, numexon, es, ee) = tsCoord
    genomic_distance = 0
    tsdistance = 0
    idx = 0
    for s in es:
        e = ee[idx]

        exonlen = e - s
        if (tsdistance + exonlen) < pos:

            genomic_distance = genomic_distance + exonlen
            if len(es) > idx + 1:
                genomic_distance = genomic_distance + (es[idx + 1] - e)
            tsdistance = tsdistance + exonlen
            idx = idx + 1
            continue
        else:
            left = pos - tsdistance
            genomic_distance = genomic_distance + left
            break

    return start + genomic_distance - 1


def filterCand(hitgenome, totalAppended):
    if totalAppended:

        ret = []
        for hit, tslist, convertlist in hitgenome:
            if len(tslist) > 0:
                ret.append((hit, tslist, convertlist))
        return ret

    else:

        if len(hitgenome) > 0:
            return hitgenome[0]
        else:
            return []


# str(hit) =
# q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str
def getMap(aligner, alignerTS, tsDb, seq, convertCigar=False):
    hitgenome = []
    hits = aligner.map(seq)
    for hit in hits:
        hit0 = toTp(hit)
        tslist = []
        convertlist = []
        hitgenome.append((hit0, tslist, convertlist))

    totalAppended = False
    # add transcript hit
    hits = alignerTS.map(seq)
    for hit in hits:
        hit0 = toTp(hit)
        appended = convertAndAdd(tsDb, hit0, hitgenome, convertCigar)
        if appended:
            totalAppended = True

    # filter with candidate with corresponding transcript
    ret = filterCand(hitgenome, totalAppended)
    return ret


# header = {'HD': {'VN': '1.0'},
#           'SQ': [{'LN': 1575, 'SN': 'chr1'},
#                  {'LN': 1584, 'SN': 'chr2'}]}
def headerFrom2bit(genomeref):
    chrtoidx = {}
    tb = py2bit.open(genomeref)
    d = tb.chroms()
    dlist = []
    idx = 0
    for chr in d:
        dlist.append({'LN': d[chr], 'SN': chr})
        chrtoidx[chr] = idx
        idx = idx + 1

    tb.close()
    header = {'HD': {'VN': '1.0'},
              'SQ': dlist}
    return header, chrtoidx


def getFlg(strand, primary):
    if strand:
        if primary:
            return 0
        else:
            return 256
    else:
        if primary:
            return 16
        else:
            return 272


def getHitToAL(readidx, id, seq, qual, chrtoidx, hit, transcriptid):
    (q_st, q_en, strand, chr, ctg_len, st, en, mlen, blen, mapq, primary, nm, cigar) = hit
    a = pysam.AlignedSegment()
    a.query_name = id + "_" + str(readidx)

    sq = seq[q_st:q_en]
    bstrand = True
    if isinstance(strand, str):
        if strand == "-":
            bstrand = False
    else:
        if strand == -1:
            bstrand = False

    if not bstrand:
        sq = str(Seq(sq).reverse_complement())
    a.query_sequence = sq
    a.flag = getFlg(bstrand, primary)
    a.reference_id = chrtoidx[chr]
    a.reference_start = st
    a.mapping_quality = mapq
    a.cigar = toTuple(cigar)

    qual = qual[q_st:q_en]
    # a.query_qualities = pysam.qualitystring_to_array(qual)
    a.query_qualities = qual
    if transcriptid is not None:
        a.tags = [("NM", nm), ("XS", transcriptid)]
    else:
        a.tags = [("NM", nm)]
    return a


def toTuple(cigar):
    ret = []
    for cg in cigar:
        ret.append((cg[1], cg[0]))
    return ret


def getAl(id, seq, qual, chrtoidx, result):
    retArray = []
    readidx = 0
    convertlist = None
    for ret in result:

        if len(ret) == 3:
            (hit, tslist, convertlist) = ret
        elif len(ret) == 0:
            continue
        else:
            hit = ret

        hital = getHitToAL(readidx, id, seq, qual, chrtoidx, hit, None)
        retArray.append(hital)
        readidx = readidx + 1

        if convertlist is not None:
            lidx = 0
            for hit in convertlist:
                transcriptid = tslist[lidx][3]
                hital = getHitToAL(readidx, id, seq, qual, chrtoidx, hit, transcriptid)
                retArray.append(hital)
                readidx = readidx + 1
                lidx = lidx + 1

    return retArray


