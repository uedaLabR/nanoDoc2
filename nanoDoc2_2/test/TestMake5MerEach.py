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
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int32,int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import numpy as np

def loadpq(path, samplesize):

    df = pq.read_table(path).to_pandas()
    return df

def toEach(path,pathout):

    samplesize = 5000


    df = loadpq(path, samplesize)
    print(df.columns)
    print(df)
    print(df.head())
    labelidx = df["fmer"].unique().tolist()
    print(labelidx)
    print(len(labelidx))

    for label in labelidx:

        output_file = pathout+"/"+label+".pq"
        subdf = df.query('fmer == "'+label+'"')
        subdf = subdf.sample(n=5000)
        # print(subdf["flg"].head(10))
        # print(subdf["fmer"].head(10))
        # print(subdf["trace"].head(10))
        # print(subdf["signal"].head(10))

        pschema = schema(
            [
                ('flg', uint32()),
                ('fmer', string()),
                ('trace', list_(int32())),
                ('signal', list_(float32()))
            ]
        )
        subdf = pd.DataFrame(subdf,
                          columns=['flg', 'fmer', 'trace', 'signal'])
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

    path = "/data/nanoDoc2_2/varidate/5000signal.pq"
    pathout = "/data/nanoDoc2_2/5000traceeach"
    toEach(path,pathout)
