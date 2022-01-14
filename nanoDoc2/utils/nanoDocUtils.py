

def get_fast5_reads_dirs(directories:list,MAX_CORE:int,readmax = -1):
    """
    return the list of reads from fast5 files in the directory

    Args:
        directory (str): path to the directory containing fast5 files
        MAX_CORE (int): maximum number of cores

    Returns:
        list: list of Read intances

    """
    f5list = []
    for directory in directories:
        print('load fast5 reads from {}'.format(directory))
        f5list.extend(get_fast5_files_in_dir(directory))

    #print('fast5 list {}'.format(f5list))
    if len(f5list) == 1:
       reads = get_fast5_reads_from_file(f5list[0])
       print('Finish. 1 reads are loaded\n'.format(len(reads)))
       return reads
    if readmax > 0:
        upto = min(readmax,len(f5list)-1)
        f5list = f5list[0:upto]
    ncore = get_number_of_core(MAX_CORE=MAX_CORE)
    with Pool(ncore) as p:
        reads = p.map(get_fast5_reads_from_file, f5list)
        reads = list(itertools.chain.from_iterable(reads)) # flatteining
    print('Finish. {}reads are loaded\n'.format(len(reads)))
    return reads

def get_fast5_reads_from_file(fast5_filepath:str):
    reads = []
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            readid = read.read_id

            row = read.handle["Raw"]
            signal = row["Signal"][()]
            channel_info = read.get_channel_info()
            digitisation = channel_info['digitisation']
            offset = channel_info['offset']
            range_value = channel_info['range']
            pA_signal = (signal + offset) * range_value / digitisation
            read_info = read.handle[read.raw_dataset_group_name].attrs
            duration = read_info['duration']
            basecall_run = read.get_latest_analysis("Basecall_1D")
            fastq = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Fastq")
            trace = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Trace")
            move = read.get_analysis_dataset(basecall_run, "BaseCalled_template/Move")

            #read_id, signal, tracelen, fastq
            if len(trace) >0:
                #read_id, signal, trace, move, fastq, duration):
                read = Read(read_id=readid,signal=pA_signal,trace=trace,move=move,fastq=fastq,duration=duration)
                reads.append(read)

        return reads
