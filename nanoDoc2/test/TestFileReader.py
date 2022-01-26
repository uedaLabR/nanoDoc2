from nanoDoc2.utils.PqFileReader import PqReader

if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2/testCurlcakeIVT_old'
    fr = PqReader(path,4000)


    for n in range(1000):

        print(1100+n)
        fr.getRowData("cc6m_2459_t7_ecorv",True,1100+n)