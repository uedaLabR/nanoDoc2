import ont_fast5_api
from ont_fast5_api.conversion_tools.multi_to_single_fast5 import convert_multi_to_single, try_multi_to_single_conversion
from ont_fast5_api.conversion_tools.single_to_multi_fast5 import batch_convert_single_to_multi, get_fast5_file_list, \
    create_multi_read_file


input_path="/data/nanopore/IVT/koreaIVT/singleFast5"
save_path="/data/nanopore/IVT/koreaIVT/multifast5"
filename_base="batch"
batch_size=4000
threads=24
recursive=True


batch_convert_single_to_multi(input_path,
                              save_path,
                              filename_base,
                              batch_size,
                              threads,
                              recursive,
                              follow_symlinks=None,
                              target_compression=None)