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

    wfile = "/data/nanopore/nanoDoc2_1/varidate/docweight"
    paramf = "/data/param20.txt"

    ref = "/data/project/dataFromNanoCompare/fa/oligo.fa"

    refpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo_control'
    targetpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo1'
    out = "/data/nanopore/nanoDoc2_1/varidate/Oligo1.txt"

    chrom = "control"
    chromtgt = "control"
    start = 1
    end = 95


    strand = "+"


    chromtgt = chrom
    minreadlen = 10
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
                                        uplimit=1000)


    refpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo_control'
    targetpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo2'
    out = "/data/nanopore/nanoDoc2_1/varidate/Oligo2.txt"

    chrom = "control"
    chromtgt = "control"
    start = 1
    end = 95


    strand = "+"


    chromtgt = chrom
    minreadlen = 10
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
                                        uplimit=1000)

    refpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo_control'
    targetpq = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo3'
    out = "/data/nanopore/nanoDoc2_1/varidate/Oligo3.txt"

    chrom = "control"
    chromtgt = "control"
    start = 1
    end = 95


    strand = "+"


    chromtgt = chrom
    minreadlen = 10
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
                                        uplimit=1000)