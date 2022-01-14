import numpy as np

ascii_order = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
ascii_score_dict = {ascii_order[k]: k for k in range(len(ascii_order))}


class nanoDocRead():


    def __init__(self,read_id,chrom, strand, r_st, r_en, q_st, q_en, cigar_str,fastq,trace,move,signal):

        self.read_id = read_id
        tracelen = len(trace)
        self.adapter_signal = signal[:len(signal) - 10 * tracelen].astype(np.float64)
        self.trace = trace[::-1].astype(np.int16)
        self.move = move[::-1].astype(np.int16)

        self.signal = signal[len(signal) - 10 * tracelen:].astype(np.float64)
        self.fastq = fastq

        fastq_list = self.fastq.split('\n')
        self.sequence = fastq_list[1].replace('U', 'T')
        self.qscore = np.array([ascii_score_dict[symbol] for symbol in fastq_list[3]], dtype=np.int16)
        self.mean_qscore = sum(self.qscore) / len(self.qscore)

        self.normalized_signal = None
        # used for trim and normalize
        self.trimIdxbyHMM = 0
        self.trimIdxbyMapping = 0
        self.normalizeDelta = 0
        self.trimmedSignal = []
        self.normalizemed = 0
        #
        self.formatSignal = []

        #trim
        self.trimSuccess = False
        self.filterFlg = 0

    def setrefgenome(self,seq):
        self.refgenome = seq
