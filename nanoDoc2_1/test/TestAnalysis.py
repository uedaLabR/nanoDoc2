import glob
from sklearn.model_selection import train_test_split

import tensorflow as tf  # add
import numpy as np
from tensorflow.keras.layers import GlobalAveragePooling1D
import numpy as np
from tensorflow.keras import Model
from nanoDoc2.network import cnnwavenet_decfilter
import nanoDoc2_1.analysis.comparisonAnalysisKmean as  comparisonAnalysisKmean
DATA_LENGTH_UNIT = 60
DATA_LENGTH = 1024
from numba import jit,u1,i8,f8
from nanoDoc2_1.utils.PqFileReader import PqReader
from nanoDoc2_1.network import CnnWavenetDecDimention
import os

if __name__ == '__main__':


    wfile = "/data/nanopore/nanoDoc2_1/weight/docweight"
    paramf = "/data/param20.txt"
    ref = "/data/nanopore/reference/NC000913.fa"
    refpq = "/data/nanopore/nanoDoc2_1/1623_ivt"
    targetpq = "/data/nanopore/nanoDoc2_1/1623_wt"
    out = "/data/nanopore/nanoDoc2_1/16S_test.txt"
    #out = "/data/nanopore/nanoDoc2_1/23S_test.txt"
    chrom = "NC_000913.3"
    chromtgt = "NC_000913.3"
    # start = 4035531
    # end = start+1541
    start = 4035531
    end = 4037072
    strand = "+"

    #    if len(sys.argv) > 11 :
    #        minreadlen = int(sys.argv[11])
    chromtgt = chrom
    minreadlen = 300
    # # for covid analysis
    # if "england" in out:
    #     chromtgt = "hCoV-19/England/02/2020|EPI_ISL_407073"
    # if "austraria" in out:
    #     chromtgt = "MT007544.1"
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)