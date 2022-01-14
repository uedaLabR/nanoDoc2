import mappy as mp
import tRex.genomemap.AbstractGenomeMap as AM
from tRex.TRexMappingResult import MappingResult
import time
import multiprocessing
from multiprocessing import Pool
import time
import random
THRES_EDIT_DIST=50


class MappyGenomeMap(AM.AbstractGenomeMap):

    def resolveCons(self,temp,nameset):
        temp2 = []
        conspost = "_CONS"
        for hd in temp:
            tRNAName = hd.tRNA_species
            if conspost in tRNAName:
                tRNAorg = tRNAName.replace(conspost, "")
                if tRNAorg not in nameset:
                    hd.tRNA_species = tRNAorg
                    temp2.append(hd)
            else:
                temp2.append(hd)

        return temp2

    def genomeMap(self,reference_holder,reads,mapOption,MAXCORE):
        print('Mapping by minimap2')
        start_time = time.time()
        self.match = mapOption['match']
        self.mismatch = mapOption['mismatch']
        self.gapopen = mapOption['gapopen']
        self.gapex = mapOption['gapex']
        self.maxstartpos = mapOption['maxstartpos']
        self.maxrefstart = mapOption['maxrefstart']
        self.minrefmatch = mapOption['minrefmatch']
        self.minreadlen = mapOption['minreadlen']
        self.maxreadlen = mapOption['maxreadlen']
        self.k = mapOption['k']
        self.w = mapOption['w']
        self.bw = mapOption['bw']
        self.best_n = mapOption['best_n']
        self.min_cnt = mapOption['min_cnt']
        self.min_chain_score = mapOption['min_chain_score']
        self.min_dp_score = mapOption['min_dp_score']

        scoretp = (self.match,self.mismatch,self.gapopen,self.gapex)
        unmod_ref_path = reference_holder.get_full_fasta_path()
        a = mp.Aligner(unmod_ref_path,min_dp_score=self.min_dp_score,n_threads=MAXCORE,w=self.w,bw=self.bw,k=self.k,best_n=self.best_n,min_cnt=self.min_cnt,min_chain_score=self.min_chain_score, scoring = scoretp)  # load or build index
        if not a: raise Exception("ERROR: failed to load/build index")

        rcnt = 0
        mappedcnt = 0
        for read in reads:
            if len(read.sequence) < self.minreadlen or len(read.sequence) > self.maxreadlen:
                continue
            temp = []
            temps = set()
            rcnt += 1
            if rcnt%1000==0:print("total",rcnt,"mapped",mappedcnt)
            for hit in a.map(read.sequence,MD=True):  # traverse alignments
                hitcound1 = hit.q_st < self.maxstartpos and hit.r_st < self.maxrefstart
                hitcound2 = (hit.r_en - hit.r_st) > self.minrefmatch

                score = countMatch(hit.cigar_str,hit.NM)
                if hitcound1 and hitcound2: # else something wrong fot trna may be tandem read
                    result = MappingResult(tRNA_species=hit.ctg, r_st=hit.r_st, r_en=hit.r_en,
                                            q_st=hit.q_st, q_en=hit.q_en, cigar_str=hit.cigar_str, MD=hit.MD,score = score)

                    if score > 0:
                        temps.add(hit.ctg)
                        temp.append(result)


            temp = self.resolveCons(temp,temps)
            if len(temp) > 1:
                temp = filterCand(temp,THRES_EDIT_DIST)

            if len(temp) > 0:
                mappedcnt += 1

            for result in temp:
                read.add_mapping_result(result)

        mapped_reads = [read for read in reads if len(read.mapping_results) >0]
        print('Finish. {}reads {}s\n'.format(len(mapped_reads),time.time()-start_time))
        return mapped_reads

def filterCand(templist,thres_edit_distancefromtophit):

    maxscore = max(list(map(lambda x:x.score,templist)))
    # print(maxscore)
    temp = list(filter(lambda x: (maxscore - x.score) <= thres_edit_distancefromtophit ,templist))
    temp = sorted(temp, key=lambda x: x.score, reverse=True)
    temp2 =[]
    hitref = set()
    for mapresult in temp:
        if mapresult.tRNA_species not in hitref:
            temp2.append(mapresult)
        hitref.add(mapresult.tRNA_species)
    # print(list(map(lambda x:x.score,temp)))
    return temp2

import pysam
def countMatch(cigar,nm):

    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    match = 0
    for cigaroprator, cigarlen in a.cigar:
        if cigaroprator == 0:  # match
            match += cigarlen
    return match - nm

def random_assign(dlist):

    holderlen = len(dlist)
    if holderlen ==0:
        return None
    idx = random.randrange(holderlen)
    return dlist[idx]

class MappyGenomeMapMT(AM.AbstractGenomeMap):

    #Experimantal
    # If multi thread version following change is required
    # for now we skip implementation
    # # load/build the index before all threads:
    # a = mp.Aligner("ref.fa")
    #
    # # then in each thread
    # thrbuf = mp.ThreadBuffer()
    # for seq in seq_array:
    #     for hit in a.map(seq, buf=thrbuf):
    #
    # # blabla

    def genomeMap(self,reference_holder,reads,mapOption,MAXCORE):
        print('Mapping by minimap2')
        start_time = time.time()
        self.match = mapOption['match']
        self.mismatch = mapOption['mismatch']
        self.gapopen = mapOption['gapopen']
        self.gapex = mapOption['gapex']
        self.maxstartpos = mapOption['maxstartpos']
        self.maxrefstart = mapOption['maxrefstart']
        self.minrefmatch = mapOption['minrefmatch']
        self.minreadlen = mapOption['minreadlen']
        self.maxreadlen = mapOption['maxreadlen']
        self.k = mapOption['k']
        self.best_n = mapOption['best_n']
        self.min_cnt = mapOption['min_cnt']
        self.min_chain_score = mapOption['min_chain_score']

        scoretp = (self.match,self.mismatch,self.gapopen,self.gapex)
        unmod_ref_path = reference_holder.get_full_fasta_path()

        ncore = multiprocessing.cpu_count()
        if ncore > MAXCORE:
            ncore = MAXCORE

        self.thrbuf = mp.ThreadBuffer()

        self.alignerholder = []
        for n in range(ncore):
            a = mp.Aligner(unmod_ref_path,k=self.k,best_n=self.best_n,min_cnt=self.min_cnt,min_chain_score=self.min_chain_score, scoring = scoretp)  # load or build index
            if not a: raise Exception("ERROR: failed to load/build index")
            self.alignerholder.append(a)

        with Pool(ncore) as p:
            mapped_reads = p.map(self.apply_mappy_mapping, reads)

        mapped_reads = [read for read in reads if len(read.mapping_results) >0]
        print('Finish. {}reads {}s\n'.format(len(mapped_reads),time.time()-start_time))
        return mapped_reads



    def apply_mappy_mapping(self, read):

        aligner = random_assign(self.alignerholder)
        if not aligner:
            return read #something wrong
        temp = []
        for hit in aligner.map(read.sequence,buf=self.thrbuf,MD=True):  # traverse alignments
            hitcound1 = hit.q_st < self.maxstartpos and hit.r_st < self.maxrefstart
            hitcound2 = (hit.r_en - hit.r_st) > self.minrefmatch
            score = countMatch(hit.cigar_str, hit.NM)

            if hitcound1 and hitcound2: # else something wrong fot trna may be tandem read
                result = MappingResult(tRNA_species=hit.ctg, r_st=hit.r_st, r_en=hit.r_en,
                                        q_st=hit.q_st, q_en=hit.q_en, cigar_str=hit.cigar_str, MD=hit.MD,score = score)
                temp.append(result)

        if len(temp) > 1:
            temp = filterCand(temp,THRES_EDIT_DIST)
        for result in temp:
            read.add_mapping_result(result)

        return read


