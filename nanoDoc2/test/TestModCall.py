from nanoDoc2.analysis import comparisonAnalyses

wfile = "/data/nanopore/IVT_Trace/weight/docweight"
paramf = "/data/param20.txt"
ref ="/data/nanopore/reference/NC000913.fa"
refpq = '/data/nanopore/nanoDoc2/1623_ivt'
targetpq = '/data/nanopore/nanoDoc2/1623_native'
#out = "/data/nanopore/rRNA/16S_test.txt"
out = "/data/nanopore/nanoDoc2/23S_test.txt"
chrom = "NC_000913.3"
chromtgt = "NC_000913.3"
# start = 4035531
# end = start+1541
start = 4037519
end = 4040423

# start = 4037519 + 1600
# end = 4037519 + 1650

# start = 4037519+1800
# end = 4037519+2200

strand = True
minreadlen = 100
#out = "/data/nanopore/rRNA/23S_1000_test7.txt
import os
import tensorflow as tf
with tf.device('/GPU:1'):

    os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
    os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"

    comparisonAnalyses.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 1000)
# out = "/data/nanopore/rRNA/23S_500.txt"
# nanoDoc2AnalysisOneSide.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 500)
# out = "/data/nanopore/rRNA/23S_250.txt"
# nanoDoc2AnalysisOneSide.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 250)
# out = "/data/nanopore/rRNA/23S_100.txt"
# nanoDoc2AnalysisOneSide.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 100)
# out = "/data/nanopore/rRNA/23S_50.txt"
# nanoDoc2AnalysisOneSide.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 50)
# out = "/data/nanopore/rRNA/23S_30.txt"
# nanoDoc2AnalysisOneSide.modCall(wfile, paramf, ref, refpq, targetpq, out, chrom, chromtgt, start, end, strand, minreadlen,uplimit = 30)