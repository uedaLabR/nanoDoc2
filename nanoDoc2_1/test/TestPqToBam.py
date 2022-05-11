import nanoDoc2_1.utils.PqToBam as PqToBam

if __name__ == "__main__":

    #ref = "/data/nanopore/reference/NC000913.fa"
    ref = "/data/nanopore/nanoDoc2_1/testrun/ecolirRNA.fa"
    #path = '/data/nanopore/nanoDoc2/1623_ivt'
    path = '/data/nanopore/nanoDoc2_1/testrun2/ecoli'
    pathout = "/data/nanopore/nanoDoc2_1/testrun2/1623_ivt_test.bam"
    # ref = "/data/nanopore/reference/Yeast_sk1.fa"
    # path = '/data/nanopore/nanoDoc2_1/bk/1825_native'
    # pathout = "/data/nanopore/nanoDoc2_1/1825_native.bam"
    PqToBam.toBam(ref,path,pathout)
