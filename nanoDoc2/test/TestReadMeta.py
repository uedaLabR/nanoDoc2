
import pyarrow.parquet as pq


def readmeta(pqf):

    parquet_file = pq.ParquetFile(pqf)
    for n in range(11):
        print(parquet_file.metadata.row_group(0).column(n).statistics)

if __name__ == "__main__":


    pqf = "/data/nanopore/nanoDoc2/testCurlcakeIVT/cc6m_2244_t7_ecorv_1_0/batch_0_pre.pq"
    readmeta(pqf)