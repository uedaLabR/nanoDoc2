import nanoDoc2_1.trainprep.make6merParquetEach as make6merParquetEach

if __name__ == "__main__":

    ref1 = "/data/nanopore/reference/Curlcake.fa"
    ref2 = "/data/nanopore/reference/Cov2_Korea.fa"

    path_w = "/data/nanopore/nanoDoc2_1/varidate/50signal/"
    refs= [ref1,ref2]
    pq1 = "/data/nanopore/nanoDoc2_1/varidate/CurlcakeIVT"
    pq2 = "/data/nanopore/nanoDoc2_1/varidate/SARSCOV2"
    pqs = [pq1,pq2]

    # fr = PqReader(pq1, 4000)
    # chr = "cc6m_2244_t7_ecorv"
    # data = fr.getRowData(chr, True, 30, 50)

    takeCnt = 50
    make6merParquetEach.makeSamplePlan(refs,pqs, path_w,takeCnt)
    #make6merParquetEach.mergePq(path_w)