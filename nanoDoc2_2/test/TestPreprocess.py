
import nanoDoc2_2.preprocess.fast5ToProcessedPq as fast5ToProcessedPq
if __name__ == "__main__":

    # path = '/data/nanopore/IVT/m6aIVT/multifats5_2/workspace'
    # pathout = '/data/nanoDoc2_2/varidate/CurlcakeIVT'
    # ref = "/data/nanopore/reference/Curlcake.fa"
    MAX_CORE = 24
    qvaluethres = 5
    fmercurrent = "/data/nanopore/signalStatRNA180.txt"

    path = '/data/nanopore/IVT/koreaIVT/multifast5_2/workspace'
    pathout = '/data/nanoDoc2_2/varidate/SARSCOV2'
    ref = "/data/nanopore/reference/Cov2_Korea.fa"
    # fmercurrent = "/data/nanopore/signalStatRNA180.txt"


    fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))


    # path = '/data/nanopore/rRNA/1623_ivt-multi/multifast5_2/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/1623_ivt'
    # #ref = "/data/nanopore/reference/NC000913.fa"
    # ref = "/data/nanopore/nanoDoc2_1/testrun/ecolirRNA.fa"
    # fmercurrent = "/data/nanopore/signalStatRNA180.txt"

    # path = '/data/nanopore/rRNA/1825_ivt-multi/multifast5_2/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/1825_ivt'
    # #ref = "/data/nanopore/reference/Yeast_sk1.fa"
    # ref = "/data/nanopore/nanoDoc2_1/testrun/yeastrRNA.fa"
    # fmercurrent = "/data/nanopore/signalStatRNA180.txt"

    # path = '/data/project/dataFromNanoCompare/basecalled/ncRNAMETTL3KD/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/ncRNAMETTL3KD'
    # ref = "/data/project/dataFromNanoCompare/fa/7SK.fa"
    # fmercurrent = "/data/nanopore/signalStatRNA180.txt"
    #
    # MAX_CORE = 24
    # qvaluethres = 5
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/7SKIVT/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/7SKIVT'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/ncRNAWT/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/ncRNAWT'
    # #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    #
    # # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_fmet_tRNA_biological/workspace'
    # # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_fmet_tRNA_biological'
    # # #ref = "/data/nanopore/reference/Yeast_sk1.fa"
    ref = "/data/project/dataFromNanoCompare/fa/tRNA.fa"
    # # fmercurrent = "/data/nanopore/signalStatRNA180.txt"
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))

    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_fmet_tRNA_synthetic/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_fmet_tRNA_synthetic'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_lys_tRNA_biological/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_lys_tRNA_biological'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_lys_tRNA_synthetic/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_lys_tRNA_synthetic'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))
    #

    #
    # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_phe_tRNA_synthetic/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_phe_tRNA_synthetic'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/E.coli_phe_tRNA_biological/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/E.coli_phe_tRNA_biological'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))

    # path = '/data/project/dataFromNanoCompare/basecalled/Oligo_control/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo_control'
    # #ref = "/data/nanopore/reference/Yeast_sk1.fa"
    # ref = "/data/project/dataFromNanoCompare/fa/oligo.fa"
    # fmercurrent = "/data/nanopore/signalStatRNA180.txt"
    #
    # MAX_CORE = 24
    # qvaluethres = 5
    # #
    # #
    # # # @cmd.command()
    # # # @click.option('-i', '--input', required='True')
    # # # @click.option('-o', '--output', required='True')
    # # # @click.option('-r', '--ref', required='True')
    # # # @click.option('-fm', '--fmercurrent', required='True')
    # # # @click.option('-t', '--thread', default=12)
    # # # @click.option('-qv', '--qvaluethres', default=5)
    # # # @click.option('-mo', '--mappyoption', nargs=4, type=click.Tuple([int, int, int, int]), default=(12, 10, 30, 20))
    # # #def fast5ToReSegmentedPq(input, output, ref, fmercurrent, thread, qvaluethres, mappyoption):
    # # #path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, mappyoption
    # #
    # fast5ToProcessedPq.h5tosegmantedPq(path,pathout,ref,MAX_CORE,qvaluethres,fmercurrent,(12, 10, 30, 20))
    #
    # ref = "/data/project/dataFromNanoCompare/fa/oligo.fa"
    # path = '/data/project/dataFromNanoCompare/basecalled/Oligo1/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo1'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/Oligo2/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo2'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    #
    # path = '/data/project/dataFromNanoCompare/basecalled/Oligo3/workspace'
    # pathout = '/data/nanopore/nanoDoc2_1/varidate/pq/Oligo3'
    #
    # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (12, 10, 30, 20))
    # #
    # # path = '/data/project/dataFromNanoCompare/basecalled/Oligo_control/workspace'
    # # pathout = '/data/nanopore/nanoDoc2_1/varidate/additional/Oligo_control'
    # #
    # # fast5ToProcessedPq.h5tosegmantedPq(path, pathout, ref, MAX_CORE, qvaluethres, fmercurrent, (4, 1, 15, 1))