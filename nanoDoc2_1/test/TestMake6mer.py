import nanoDoc2_1.trainprep.make6merParquet as make6merParquet

if __name__ == "__main__":

    ref1 = "/data/nanopore/reference/Curlcake.fa"
    ref2 = "/data/nanopore/reference/Cov2_Korea.fa"

    path_w = "/data/nanopore/nanoDoc2/5000signal.pq"
    refs= [ref1,ref2]
    pq1 = "/data/nanopore/nanoDoc2_1/CurlcakeIVT"
    pq2 = "/data/nanopore/nanoDoc2_1/SARSCOV2"
    pqs = [pq1,pq2]
    #

    # fr = PqReader(pq1, 4000)
    # chr = "cc6m_2244_t7_ecorv"
    # data = fr.getRowData(chr, True, 30, 50)

    takeCnt = 3000
    make6merParquet.makeSamplePlan(refs,pqs, path_w,takeCnt)