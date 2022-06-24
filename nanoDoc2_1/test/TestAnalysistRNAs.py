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

    ref = "/data/project/dataFromNanoCompare/fa/tRNA.fa"

    refpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_phe_tRNA_synthetic'
    targetpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_phe_tRNA_biological'
    out = "/data/nanopore/nanoDoc2_1/varidate/E.coli_phe_tRNA_biological"

    chrom = "Phe"
    chromtgt = "Phe"
    # start = 455938
    # end = 457732
    start = 1
    end = 76

    # start = 460923
    # end = 464318

    strand = "+"


    chromtgt = chrom
    minreadlen = 10
    with tf.device('/CPU:0'):
        # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
                                        uplimit=500)
        #comparisonAnalysisKmean.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen)

    # refpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_fmet_tRNA_synthetic'
    # targetpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_fmet_tRNA_biological'
    # out = "/data/nanopore/nanoDoc2_1/varidate/E.coli_fmet_tRNA_biological"
    #
    # chrom = "fMet1"
    # chromtgt = "fMet1"
    # # start = 455938
    # # end = 457732
    # start = 1
    # end = 76
    #
    # # start = 460923
    # # end = 464318
    #
    # strand = "+"
    #
    # chromtgt = chrom
    # minreadlen = 30
    # with tf.device('/CPU:0'):
    #     # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
    #     # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
    #     os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    #     comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
    #                                     uplimit=500)
    #
    # refpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_lys_tRNA_synthetic'
    # targetpq = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_lys_tRNA_biological'
    # out = "/data/nanopore/nanoDoc2_1/varidate/E.coli_lys_tRNA_biological"
    #
    # chrom = "Lys"
    # chromtgt = "Lys"
    # # start = 455938
    # # end = 457732
    # start = 1
    # end = 76
    #
    # # start = 460923
    # # end = 464318
    #
    # strand = "+"
    #
    # chromtgt = chrom
    # minreadlen = 30
    # with tf.device('/CPU:0'):
    #     # os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
    #     # os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
    #     os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    #     comparisonAnalysisKmean.modCall(wfile, ref, refpq, targetpq, out, chrom, chromtgt, start, end, minreadlen,
    #                                     uplimit=500)