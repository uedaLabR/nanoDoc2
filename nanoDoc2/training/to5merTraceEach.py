import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import pandas as pd
from tensorflow.keras.callbacks import ModelCheckpoint
import pandas as pd
from Bio import SeqIO
from operator import itemgetter, attrgetter
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import numpy as np

def loadpq(path, samplesize):

    df = pq.read_table(path).to_pandas()
    return df

def toEach(path,pathout):

    samplesize = 5000


    df = loadpq(path, samplesize)
    labelidx = df["fmer"].unique().tolist()
    for label in labelidx:

        output_file = pathout+"/"+label+".pq"
        subdf = df.query('fmer == "'+label+'"')
        subdf = subdf.sample(n=100)

        pschema = schema(
            [
                ('flg', uint32()),
                ('fmer', string()),
                ('traces', list_(uint8()))
            ]
        )
        subdf = pd.DataFrame(subdf,
                          columns=['flg', 'fmer', 'traces'])
        pd.set_option('display.max_rows', None)

        pyarrow_table = Table.from_pandas(subdf, pschema)
        pq.write_table(
            pyarrow_table,
            output_file,
            row_group_size=5000,
            compression='snappy',
            flavor=['spark'],
        )


if __name__ == "__main__":

    path = "/data/nanopore/nanoDoc2/5000traceeachFix.pq"
    pathout = "/data/nanopore/nanoDoc2/100traceeach"
    toEach(path,pathout)
