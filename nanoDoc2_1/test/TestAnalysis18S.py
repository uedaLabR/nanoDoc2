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
    #ref = "/data/nanopore/reference/S288C_reference_sequence_R1-1-1_19960731.fa"
    #ref = "/data/nanopore/reference/Yeast_sk1.fa"
    ref = "/data/nanopore/nanoDoc2_1/testrun/yeastrRNA.fa"
    refpq = "/data/nanopore/nanoDoc2_1/testrun/1825_ivt"
    targetpq = "/data/nanopore/nanoDoc2_1/testrun/1825_native"
    out = "/data/nanopore/nanoDoc2_1/testrun/18S_score.txt"

    chrom = "yeast18S"
    chromtgt = "yeast18S"
    # start = 455938
    # end = 457732
    start = 1
    end = 1794
    # start = 465070
    # end = 466869

    strand = "+"

    chromtgt = chrom
    minreadlen = 200
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)