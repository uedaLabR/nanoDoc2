import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO
from sklearn.model_selection import train_test_split
import os
import shutil


from src.nanoDoc2 import PqReader


def loadParquet(p, chr, idxs):
    cnt = 0
    for idx in idxs:
        idx = idx.replace('}', '')
        idx = idx.replace('{', '')
        idx = idx.replace(' ', '')
        pp = p + "/algined" + idx + ".pq"
        if not os.path.exists(pp):
            pp = p + "/" + idx + ".pq"
        table = pq.read_table(pp)
        df = table.to_pandas()
        if cnt == 0:
            totaldf = df
        else:
            totaldf = pd.concat([totaldf, df], axis=0)
        cnt = cnt + 1
    #    print(totaldf)
    return totaldf


#     table2 = pq.read_table('example.parquet')

def loadRef(ref, chr):
    record = None
    records = SeqIO.parse(ref, 'fasta')
    for r in records:
        record = r
        if r.id == chr:
            break
    return record


def getData(df1, reference, position, samplenum):
    df1 = df1[(df1.mapped_start < position) & (df1.mapped_end > position + 5)]
    train, test = train_test_split(df1, test_size=samplenum)
    print(test)
    cnt = 0
    unitwidth = 60

    #         df = pd.DataFrame(data, columns=['nucb4After','mapped_chrom','position','mapped_start','mapped_end','signal','originalsize'])
    data = []
    for index, row in test.iterrows():
        mapped_chrom = row['mapped_chrom']
        mapped_start = row['mapped_start']
        mapped_end = row['mapped_end']
        nucb4After = reference.seq[position - 1] + reference.seq[position + 5]

        relpos = position - mapped_start
        signal = row['signal'][relpos * unitwidth:(relpos + 5) * unitwidth]

        #         if cnt <= 1:
        #             isignal = list(signal)
        #             plt.figure(figsize=(20, 5))
        #             plt.plot(isignal)

        originalsize = row['originalsize'][relpos:relpos + 5]
        cnt = cnt + 1
        data.append((nucb4After, mapped_chrom, position, mapped_start, mapped_end, signal, originalsize))

    return data


def mkpq(planf, ref, p, outdir):

    f = open(planf)
    line = True
    parquetLoaded = False
    matrix = None
    reference = None
    lastreference = None
    cnt = 0
    minreadlen = 100
    refpr = PqReader.PqReader(p, minreadlen)
    uplimit = 10000

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    while line:

        line = f.readline().rstrip('\n')
        data = line.split(",")
        if len(data) < 2:
            break

        chr = data[0]
        position = int(data[1])
        samplenum = int(data[3])
        fivemer = data[4]
        print(chr)
        if chr != lastreference:
            reference = loadRef(ref, chr)
            if reference is None:
                continue

        outpath = outdir + "/" + fivemer + "/" + chr.replace('/','_') + "_" + data[1] + ".pq"
        if os.path.exists(outpath):
            continue

        strand = '+'
        rawdatas, cnt = refpr.getData(chr, strand, position, samplenum)

        if not os.path.exists(outdir + "/" + fivemer):
            os.makedirs(outdir + "/" + fivemer)

        df = pd.DataFrame(rawdatas,
                          columns=['signal', 'originalsize'])
        df.to_parquet(outpath)
        cnt = cnt + 1
        lastreference = chr

    f.close

    # merge files
    # files = []
    # for x in os.listdir(outdir):
    #     if os.path.isdir(outdir + "/" + x):
    #         files.append(x)
    #
    # totaldf = None
    # for dir in files:
    #     cnt = 0
    #     for each in os.listdir(outdir + "/" + dir):
    #         s = outdir + "/" + dir + "/" + each
    #         try:
    #             table = pq.read_table(s, columns=['signal', 'originalsize'])
    #             df = table.to_pandas()
    #             if cnt == 0:
    #                 totaldf = df
    #             else:
    #                 totaldf = pd.concat([totaldf, df], axis=0)
    #             cnt = cnt + 1
    #         except:
    #             pass
    #     outpath = outdir + "/" + dir + ".pq"
    #     print(outpath)
    #     totaldf.to_parquet(outpath)
    #     shutil.rmtree(outdir + "/" + dir)

def margeFile(indir):

    files = []
    for x in os.listdir(indir):
         if os.path.isdir(indir + "/" + x):
             files.append(x)

    totaldf = []
    files = sorted(files)
    for dir in files:
        cnt = 0
        for each in os.listdir(indir + "/" + dir):
            s = indir + "/" + dir + "/" + each
            print(s)

            table = pq.read_table(s, columns=['signal', 'originalsize'])# mistake actually originalsize is signal
            df = table.to_pandas()
            for index, row in df.iterrows():
                sixmer = dir
                tp = (sixmer,row[1])
                print(sixmer,len(row[1]),len(totaldf))
                totaldf.append(tp)


    outpath = indir + "/marged6mer.pq"
    print(outpath)

    df = pd.DataFrame(totaldf,
                      columns=['sixmer','signal'])
    df.to_parquet(outpath)

import random
def margeFile2(indir):

    files = []
    for x in os.listdir(indir):
         if os.path.isdir(indir + "/" + x):
             files.append(x)

    totaldf = []
    files = sorted(files)
    for dir in files:
        cnt = 0
        sixmerdata = []
        for each in os.listdir(indir + "/" + dir):
            s = indir + "/" + dir + "/" + each
            print(s)

            table = pq.read_table(s, columns=['signal', 'originalsize'])# mistake actually originalsize is signal
            df = table.to_pandas()
            for index, row in df.iterrows():
                sixmer = dir
                sixmeridx = files.index(sixmer)
                tp = (sixmer,sixmeridx,row[1])
                sixmerdata.append(tp)
                #print(sixmer,len(row[1]),len(totaldf))

        random.shuffle(sixmerdata)
        sixmerdata = sixmerdata[0:400]
        totaldf.extend(sixmerdata)


    outpath = indir + "/marged6mer_500.pq"
    print(outpath)

    df = pd.DataFrame(totaldf,
                      columns=['sixmer','sixmeridx','signal'])
    df.to_parquet(outpath)

import sys

if __name__ == "__main__":

    outdir = "/groups2/gac50430/nanopore/dataset4DL/fivemer"
    args = sys.argv
    path = args[1]
    main(path, outdir)
    print(path)

