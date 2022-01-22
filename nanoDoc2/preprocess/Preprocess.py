import nanoDoc2.preprocess.ViterbiSegmentation  as vs
import nanoDoc2.preprocess.SignalNormalization as ss


def preprocess(read,fmercurrent):

    #Viterbi segmentation
    seq, cigar, left, traceboundary, frombasecaller_idx,possiblemove_idx = vs.flipplopViterbiEach(read)
    read.cigar_str = cigar
    read.settraceboundary(traceboundary)
    #Signal Normalize
    read.signal = ss.normalizeSignal(read,traceboundary,fmercurrent)

    return read

