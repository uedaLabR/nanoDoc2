from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2.utils.PqFileReader import PqReader
from matplotlib import pyplot as plt
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import pandas as pd


def getPlot(data,pos):


    plt.suptitle(str(pos))
    fig = plt.plot(data, linewidth=1)
    return fig

if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    fr = PqReader(path,100)

    pathout = '/data/nanopore/nanoDoc2/testfig/test.pdf'
    output_file = '/data/nanopore/nanoDoc2/nanodoctest.pq'

    gm = GraphManager(pathout)
    datab =[]
    for n in range(2):

        pos = 1500+n
        data,cnt = fr.getRowData("cc6m_2459_t7_ecorv",True,pos)
        datab.extend(data)

    print(datab)
    pschema = schema(
            [
                ('key', string()),
                ('signal', list_(float32()))
            ]
    )
    r = []
    tp = ("ATGCA", datab)
    r.append(tp)
    print(r)
    df = pd.DataFrame(r,
                      columns=['key', 'signal'])
    pyarrow_table = Table.from_pandas(df, pschema)
    pq.write_table(
        pyarrow_table,
        output_file,
        row_group_size=4000,
        compression='snappy',
        flavor=['spark'],
    )



    gm.save()