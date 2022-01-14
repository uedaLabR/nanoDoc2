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

from nanoDoc2.preprocess import Preprocess
from nanoDoc2.utils.nanoDocRead import nanoDocRead


def get_number_of_core(MAX_CORE:int):
    ncore = multiprocessing.cpu_count()
    if ncore > MAX_CORE:
        ncore = MAX_CORE
    return ncore

def get_fast5_files_in_dir(directory:str):
    return list(sorted(glob.glob(directory + '/**/*.fast5',recursive=True)))


def preprocess(f5file,pathout,ref,fmercurrent,ncore,qvaluethres):

    reads = []
    aligner = mp.Aligner(ref,preset="map-ont")
    with get_fast5_file(f5file, mode="r") as f5:
        for read in f5.get_reads():

            basecall_run = read.get_latest_analysis("Basecall_1D")
            fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
            trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
            move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")
            row_data = read.get_raw_data()
            #
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

                read = nanoDocRead(read.read_id, chrom, strand, r_st, r_en, q_st, q_en, cigar_str, fastq, trace, move,
                            row_data)

                #filter by qvalue score
                if(read.mean_qscore > qvaluethres):

                    #set reference
                    lgenome = aligner.seq(chrom, start=r_st, end=r_en)
                    read.setrefgenome(lgenome)
                    reads.append(read)

                #take top hit only
                break


        preprocess = partial(Preprocess.preprocess,pathout=pathout,ref=ref,fmercurrent=fmercurrent)
        print("finish mapping to reference")
        with Pool(ncore) as p:
            #Viterbi Segmentation/Normalize/Write To Partial File
            p.map(preprocess , reads)
        print("finish segmentation and output file")

def h5tosegmantedPq(path,pathout,ref,fmercurrent,MAX_CORE):

    f5list = get_fast5_files_in_dir(path)
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    for f5file in f5list:

        preprocess(f5file,pathout,ref,fmercurrent,ncore)




if __name__ == "__main__":

    path = '/data/nanopore/IVT/m6aIVT/multifast5/fast5'
    pathout = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    ref = "/data/nanopore/reference/Curlcake.fa"
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"

    MAX_CORE = 8
    qvaluethres = 5
    h5tosegmantedPq(path,pathout,ref,fmercurrent,MAX_CORE,qvaluethres)