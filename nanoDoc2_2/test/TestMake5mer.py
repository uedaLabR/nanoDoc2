import nanoDoc2_2.trainprep.make5merParquet as make5merParquet

if __name__ == "__main__":

    ref1 = "/data/nanopore/reference/Curlcake.fa"
    ref2 = "/data/nanopore/reference/Cov2_Korea.fa"

    path_w = "/data/nanoDoc2_2/varidate/5000signal.pq"
    refs= [ref1,ref2]
    pq1 = "/data/nanoDoc2_2/varidate/CurlcakeIVT"
    pq2 = "/data/nanoDoc2_2/varidate/SARSCOV2"
    pqs = [pq1,pq2]
    #

    # fr = PqReader(pq1, 4000)
    # chr = "cc6m_2244_t7_ecorv"
    # data = fr.getRowData(chr, True, 30, 50)

    takeCnt = 5000
    make5merParquet.makeSamplePlan(refs,pqs, path_w,takeCnt)