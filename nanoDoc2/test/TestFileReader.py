from nanoDoc2.utils.PqFileReader import PqReader

if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    fr = PqReader(path,4000)

    fr.loadpq("cc6m_2459_t7_ecorv",1100,True)