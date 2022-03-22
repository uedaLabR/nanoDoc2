
import nanoDoc2_1.preprocess.fast5ToProcessedPq as fast5ToProcessedPq
if __name__ == "__main__":

    path = '/data/nanopore/IVT/koreaIVT/multifast5_2/workspace'
    pathout = '/data/nanopore/nanoDoc2_1/SARSCOV2'
    ref = "/data/nanopore/reference/Cov2_Korea.fa"
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"


    MAX_CORE = 24
    qvaluethres = 5
    fast5ToProcessedPq.h5tosegmantedPq(path,pathout,ref,MAX_CORE,qvaluethres,fmercurrent)