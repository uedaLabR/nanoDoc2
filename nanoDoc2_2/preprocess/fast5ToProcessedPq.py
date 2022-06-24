#load mutli fast5 file require,  basecalled multu fast5
#map by mappy
#viterbi segmantaiton
#bandle to parquet file
import mappy as mp
from ont_fast5_api.fast5_interface import get_fast5_file
import glob
import multiprocessing
from multiprocessing import Pool
from functools import partial

from nanoDoc2_2.preprocess import Preprocess
from nanoDoc2_1.utils import FileIO
from nanoDoc2.utils.nanoDocRead import nanoDocRead
import nanoDoc2_1.preprocess.WriteToFile as wf
import pyarrow.parquet as pq
from nanoDoc2.utils import nanoDocUtils as nutils

def get_number_of_core(MAX_CORE:int):

    ncore = int(multiprocessing.cpu_count())
    print(ncore,MAX_CORE)
    if ncore > MAX_CORE:
        ncore = MAX_CORE
    return ncore

def get_fast5_files_in_dir(directory:str):
    return list(sorted(glob.glob(directory + '/*.fast5',recursive=True)))

import os
def preprocess(f5file,pathout,ref,ncore,qvaluethres,fmercurrent,mappyoption):

    reads = []
    ret = []
    k,w,min_chain_score,min_dp_score = mappyoption
    aligner = mp.Aligner(ref,k=k,w=w,min_chain_score=min_chain_score,min_dp_score=min_dp_score,best_n=1)
    with get_fast5_file(f5file, mode="r") as f5:
        for read in f5.get_reads():
            try:
                basecall_run = read.get_latest_analysis("Basecall_1D")
                fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
                trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
                move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")
                row_data = read.get_raw_data(scale=True)

                #
                #print(fastq)
                seq = fastq.split("\n")[1]
                chrom, strand, r_st, r_en, q_st, q_en, cigar_str = "", 0, 0, 0, 0, 0, ""
                for hit in aligner.map(seq):
                    #take tophit only
                    chrom = hit.ctg
                    strand = hit.strand
                    r_st = hit.r_st
                    r_en = hit.r_en
                    q_st = hit.q_st
                    q_en = hit.q_en
                    cigar_str = hit.cigar_str
                    print(chrom,strand,r_st,r_en)
                    read = nanoDocRead(read.read_id, chrom, strand, r_st, r_en, q_st, q_en, cigar_str, fastq, trace, move,
                                row_data)

                    #filter by qvalue score
                    if(read.mean_qscore > qvaluethres):

                        #set reference
                        lgenome = aligner.seq(chrom, start=r_st, end=r_en)
                        if strand == -1:
                            lgenome = mp.revcomp(lgenome)

                        read.setrefgenome(lgenome)
                        reads.append(read)

                    #take top hit only
                    break
            except KeyError:
                print('Key Error')


    fmerdict = nutils.getCurrentDict(fmercurrent)
    preprocess = partial(Preprocess.preprocess,fmerDict=fmerdict)
    print("finish mapping to the reference")
    with Pool(ncore) as p:
        #Viterbi Segmentation/Normalize/
        ret = p.map(preprocess , reads)

    #separate to key
    filename = os.path.basename(f5file)
    filename = os.path.splitext(filename)[0]
    wf.writeToFile2(pathout,ncore,ret,filename)
    #

    print("finish segmentation and output file")



def records(df): return df.to_records(index=False).tolist()

import pandas as pd
import shutil
def _mergeParquet(dirinfo):

    pqidx,dir = dirinfo

    files = glob.glob(dir+"/*_pre.pq")
    data = []
    for file in files:
        asfile = os.path.abspath(file)
        #df = pd.read_pickle(asfile)
        df = pd.read_parquet(asfile)
        rec = records(df)
        data.extend(rec)


    df = pd.DataFrame(data,
                      columns=['read_id','chr', 'strand', 'start', 'end', 'cigar','genome','fastq','offset','traceintervals', 'trace','signal','signalboundary'])

    df_s = df.sort_values('start')
    df_s = df_s.reset_index()
    df_s = df_s.rename(columns={"index": "read_no"})
    df_s['read_no'] = df_s.index

    file_out = dir+".pq"
    FileIO.writeToPqWithIdx(df_s, file_out)
    shutil.rmtree(dir)

def mergeParquet(pathout,ncore):

    dirlist = os.listdir(pathout)
    dirlistTohandle = []

    pqidx = 0
    for dir in dirlist:

        p = pathout+"/"+dir
        pqidx = pqidx+1
        dirlistTohandle.append((pqidx,p))


    with Pool(ncore) as p:

        p.map(_mergeParquet,dirlistTohandle)


def h5tosegmantedPq(path,pathout,ref,MAX_CORE,qvaluethres,fmercurrent,mappyoption):

    f5list = get_fast5_files_in_dir(path)
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    cnt = 0
    for f5file in f5list:

        print(cnt,f5file)
        preprocess(f5file,pathout,ref,ncore,qvaluethres,fmercurrent,mappyoption)
        cnt += 1

    #marge Parquet
    print("merge files")
    mergeParquet(pathout, ncore)




