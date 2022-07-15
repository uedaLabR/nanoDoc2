
import nanoDoc2_2.preprocess.fast5ToProcessedPq as fast5ToProcessedPq
if __name__ == "__main__":

    path = '/data/nanopore/IVT/m6aIVT/multifats5_2/workspace'
    pathout = '/data/nanoDoc2_2/varidate/CurlcakeIVT'
    ref = "/data/nanopore/reference/Curlcake.fa"
    MAX_CORE = 24
    qvaluethres = 5
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"


    fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))


