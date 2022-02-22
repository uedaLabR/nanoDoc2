import pandas as pd
import shutil
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq

def writeToPqWithIdx(df,output_file):

    pschema = schema(
        [
            ('read_no', uint32()),
            ('read_id', string()),
            ('chr', string()),
            ('strand', bool_()),
            ('start', uint32()),
            ('end', uint32()),
            ('cigar', string()),
            ('genome', string()),
            ('fastq', string()),
            ('offset', uint16()),
            ('traceintervals', list_(uint16())),
            ('trace', list_(uint8()))
        ]
    )

    pyarrow_table = Table.from_pandas(df, pschema)
    pq.write_table(
        pyarrow_table,
        output_file,
        row_group_size = 4000,
        compression='snappy',
        flavor=['spark'],
    )

def writeToPq(df,output_file):

    pschema = schema(
        [
            ('read_id', string()),
            ('chr', string()),
            ('strand', bool_()),
            ('start', uint32()),
            ('end', uint32()),
            ('cigar', string()),
            ('genome', string()),
            ('fastq', string()),
            ('offset', uint16()),
            ('traceintervals', list_(uint16())),
            ('trace', list_(uint8()))
        ]
    )

    pyarrow_table = Table.from_pandas(df, pschema)
    pq.write_table(
        pyarrow_table,
        output_file,
        row_group_size = 4000,
        compression='snappy',
        flavor=['spark'],
    )