from sklearn.model_selection import train_test_split

def getCurrentDict(fmercurrent):
    a = {}
    with open(fmercurrent) as f:
        for line in f:
            if not line.startswith("#"):
                data = line.split()
                a[data[0]] = float(data[1])
    return a


def reducesize(data, size):
    if len(data) <= size:
        return data, len(data)

    a_train, a_test = train_test_split(data, test_size=size)
    return a_test, len(a_test)


def readParam(paramf):
    f = open(paramf)
    for line in f:
        if line.startswith("a="):
            a = float(line.replace("a=", ""))
        if line.startswith("b="):
            b = float(line.replace("b=", ""))
        if line.startswith("uplimit="):
            uplimit = int(line.replace("uplimit=", ""))
        if line.startswith("takeparcentile="):
            takeparcentile = float(line.replace("takeparcentile=", ""))

    f.close()
    return a, b, uplimit, takeparcentile

from Bio.Seq import Seq
from Bio import SeqIO

def getSeq(ref, chrtgt, start, end, strand):
    records = SeqIO.parse(ref, 'fasta')
    record = None
    for r in records:
        if r.id == chrtgt:
            record = r
    if record is None:
        return None
    seq = record.seq[start:end]
    if strand == "-":
        seq2 = Seq(str(seq))
        seq = seq2.reverse_complement()
    seq = seq.upper()
    return seq


def getFirstChrom(ref):
    records = SeqIO.parse(ref, 'fasta')
    for r in records:
        return r.id
    return ""
