import glob
from sklearn.model_selection import train_test_split

import tensorflow as tf  # add
import numpy as np
from tensorflow.keras.layers import GlobalAveragePooling1D
import numpy as np
from tensorflow.keras import Model
from nanoDoc2.network import cnnwavenet_decfilter
import nanoDoc2_2.analysis.ComparisonAnalysisKmean2 as  comparisonAnalysisKmean
DATA_LENGTH_UNIT = 60
DATA_LENGTH = 1024
from numba import jit,u1,i8,f8
from nanoDoc2_1.utils.PqFileReader import PqReader
from nanoDoc2_1.network import CnnWavenetDecDimention
import os

if __name__ == '__main__':


    wfile = "/data/nanoDoc2_2/docweight"
    paramf = "/data/param20.txt"
    # ref = "/data/nanopore/reference/NC000913.fa"
    # refpq = "/data/nanopore/nanoDoc2_1/1623_ivt"
    # targetpq = "/data/nanopore/nanoDoc2_1/1623_wt"
    # out = "/data/nanopore/nanoDoc2_1/16S_score.txt"
    #
    # chrom = "NC_000913.3"
    # chromtgt = "NC_000913.3"
    # start = 4035531
    # end = 4037072

    ref = "/data/nanopore/nanoDoc2_1/testrun/ecolirRNA.fa"
    refpq = "/data/nanoDoc2_2/varidate/1623_ivt"
    targetpq = "/data/nanoDoc2_2/varidate/1623_wt"
    out = "/data/nanoDoc2_2/varidate/23S_score.txt"

    chrom = "ecoli23S"
    chromtgt = "ecoli23S"
    start = 2420
    end = 2904

    # start = 4035570
    # end = 4035580
    strand = "+"

    chromtgt = chrom
    minreadlen = 1000
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        #modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen, uplimit=500):
        comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen)