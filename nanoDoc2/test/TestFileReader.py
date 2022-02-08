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

    #path = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    path = '/data/nanopore/nanoDoc2/testSARSCOV2'
    fr = PqReader(path,4000)

    # pathout = '/data/nanopore/nanoDoc2/testfig/test.pdf'
    # output_file = '/data/nanopore/nanoDoc2/nanodoctest.pq'

    #gm = GraphManager(pathout)
    datab =[]
    depth = fr.getRowData("hCoV-19_South_Korea", True, 10001)
    for n in range(2):

        pos = 10713+n
        depth = fr.getRowData("hCoV-19_South_Korea",True,pos)
        print(depth)
        data,cnt = fr.getRowData("hCoV-19_South_Korea",True,pos,takecnt=100)
        print(cnt)

    # print(datab)
    # pschema = schema(
    #         [
    #             ('key', string()),
    #             ('signal', list_(float32()))
    #         ]
    # )
    # r = []
    # tp = ("ATGCA", datab)
    # r.append(tp)
    # print(r)
    # df = pd.DataFrame(r,
    #                   columns=['key', 'signal'])
    # pyarrow_table = Table.from_pandas(df, pschema)
    # pq.write_table(
    #     pyarrow_table,
    #     output_file,
    #     row_group_size=4000,
    #     compression='snappy',
    #     flavor=['spark'],
    # )



    #gm.save()