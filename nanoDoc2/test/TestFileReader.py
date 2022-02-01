from nanoDoc2.graph.GraphManager import GraphManager
from nanoDoc2.utils.PqFileReader import PqReader
from matplotlib import pyplot as plt

def getPlot(data,pos):


    plt.suptitle(str(pos))
    fig = plt.plot(data, linewidth=1)
    return fig

if __name__ == "__main__":

    path = '/data/nanopore/nanoDoc2/testCurlcakeIVT'
    fr = PqReader(path,100)

    pathout = '/data/nanopore/nanoDoc2/testfig/test.pdf'
    gm = GraphManager(pathout)
    for n in range(100):

        pos = 1500+n
        data,cnt = fr.getRowData("cc6m_2459_t7_ecorv",True,pos)
        print(pos,data)
        # for n in range(10):
        #     d = data[n]
        #     print(d[1])
        #     plt.plot(d[1], linewidth=1)
        # plt.savefig('/data/nanopore/nanoDoc2/testfig/'+str(pos)+"_"+str(n)+".png")
            #gm.add_figure(fig)


    gm.save()