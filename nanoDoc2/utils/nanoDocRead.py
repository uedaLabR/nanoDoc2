import numpy as np

ascii_order = '!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
ascii_score_dict = {ascii_order[k]: k for k in range(len(ascii_order))}


class nanoDocRead():

    #read.read_id, chrom, strand, r_st, r_en, q_st, q_en, cigar_str, fastq, trace, move,row_data
    def __init__(self,read_id,chrom, strand, r_st, r_en, q_st, q_en, cigar_str,fastq,trace,move,signal):

        self.read_id = read_id
        self.chrom = chrom
        self.strand = strand
        self.r_st = r_st
        self.r_en = r_en
        self.q_st = q_st
        self.q_en = q_en
        self.cigar_str = cigar_str

        tracelen = len(trace)
        adoptorlen = (len(signal) - (10 * tracelen))
        self.signal = signal[adoptorlen:].astype(np.float64)
        self.normSignal = None

        if strand == 1:

            self.trace = trace[::-1].astype(np.int16)
            self.move = move[::-1].astype(np.int16)
            self.signal = self.signal[::-1]

        else:
            self.trace = trace.astype(np.int16)
            self.move = move.astype(np.int16)


        self.fastq = fastq

        fastq_list = self.fastq.split('\n')
        self.sequence = fastq_list[1].replace('U', 'T')
        self.qscore = np.array([ascii_score_dict[symbol] for symbol in fastq_list[3]], dtype=np.int16)
        self.mean_qscore = sum(self.qscore) / len(self.qscore)

        self.refgenome = ""
        self.traceboundary = None


    def setrefgenome(self,seq):
        self.refgenome = seq

    def settraceboundary(self,traceboundary):
        self.traceboundary = traceboundary
