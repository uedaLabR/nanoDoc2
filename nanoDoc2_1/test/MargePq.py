import os
import glob
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import pandas as pd
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq

if __name__ == "__main__":

    path = "/data/nanopore/nanoDoc2_1/10000signal/"
    dirs = os.listdir(path=path)

    for dir in dirs:

        fout = path +"/" + dir+".pq"
        f = path +"/" + dir
        files = glob.glob(f + "/*.pq")
        dftotal = None
        for file in files:


            df = pq.read_table(file).to_pandas()
            if dftotal is None:
                dftotal = df
            else:
                dftotal = pd.concat([dftotal, df])

        pschema = schema(
            [
                ('flg', uint32()),
                ('fmer', string()),
                ('trace', list_(uint16())),
                ('signal', list_(float32()))
            ]
        )


        pyarrow_table = Table.from_pandas(dftotal, pschema)
        pq.write_table(
            pyarrow_table,
            fout,
            row_group_size=4000,
            compression='snappy',
            flavor=['spark'],
        )

