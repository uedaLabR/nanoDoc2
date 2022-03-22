
import nanoDoc2_1.preprocess.fast5ToProcessedPq as fast5ToProcessedPq
if __name__ == "__main__":

    path = '/data/nanopore/IVT/m6aIVT/multifats5_2/workspace'
    pathout = '/data/nanopore/nanoDoc2_1/CurlcakeIVT'
    ref = "/data/nanopore/reference/Curlcake.fa"

    # path = '/data/nanopore/rRNA/1623_ivt-multi/multifast5_2/workspace'
    # pathout = '/data/nanopore/nanoDoc2/1623_ivt'
    # ref = "/data/nanopore/reference/NC000913.fa"
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"


    MAX_CORE = 24
    qvaluethres = 5
    fast5ToProcessedPq.h5tosegmantedPq(path,pathout,ref,MAX_CORE,qvaluethres,fmercurrent)