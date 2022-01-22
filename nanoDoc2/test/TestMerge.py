import os
import mappy as mp
from ont_fast5_api.fast5_interface import get_fast5_file
import glob
import multiprocessing
from multiprocessing import Pool
from functools import partial

from nanoDoc2.preprocess import Preprocess
from nanoDoc2.utils import FileIO
from nanoDoc2.utils.nanoDocRead import nanoDocRead
import nanoDoc2.preprocess.WriteToFile as wf
import pyarrow.parquet as pq
import pandas as pd
import shutil
from pyarrow import int8, int16, int32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq

def records(df): return df.to_records(index=False).tolist()

def _mergeParquet(dir):

    files = glob.glob(dir+"/*.pickle")
    data = []
    cnt = 0
    for file in files:
        asfile = os.path.abspath(file)
        df = pd.read_pickle(asfile)
        rec = records(df)
        data.extend(rec)
        cnt = cnt+1
        if cnt == 10:
            break

    df = pd.DataFrame(data,
                      columns=['chr', 'strand', 'start', 'end', 'cigar', 'fastq', 'traceboundary', 'signal'])

    df_s = df.sort_values('start')
    output_file = dir + ".pq"
    FileIO.writeToPq(df_s,output_file)
    #df_s.to_parquet(output_file, compression='snappy')
    print(pq.read_metadata(output_file).schema)







if __name__ == "__main__":

    path = '/data/nanopore/IVT/m6aIVT/multifast5/'
    pathout = '/data/nanopore/nanoDoc2/testCurlcakeIVT/cc6m_2709_t7_ecorv_1_0'
    ref = "/data/nanopore/reference/Curlcake.fa"
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"

    MAX_CORE = 20
    qvaluethres = 5
    ncore = 1
    _mergeParquet(pathout)