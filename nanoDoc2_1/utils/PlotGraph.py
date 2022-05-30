from Bio import SeqIO


def getData(p):

    f = open(p)
    lines = f.readlines()
    data = []
    data.append(0)
    data.append(0)
    data.append(0)
    data.append(0)
    for line in lines:
       if line.startswith("#"):
           continue
       line = line.split("\t")
       score = float(line[4])
       data.append(score)

    return data

def getAns(p):

    f = open(p)
    lines = f.readlines()
    data = []
    for line in lines:
       line = line.split("\t")
       pos = float(line[0])
       data.append(pos)

    return data

import matplotlib.pyplot as plt
import numpy as np
def plotGraph(fin,out,ans,figsize=(23,4)):

    answer = getAns(ans)
    scores = getData(fin)
    plt.figure(figsize=figsize)
    plt.plot(scores,lw=2)
    for xc in answer:
        plt.axvline(xc, color='red', alpha=.25,lw=2)

    plt.savefig(out)
    plt.ylim(-0.1, 1)


    for n in range(len(scores)):

        if n in answer:
            #print(n,max(scores[n-2],scores[n-1],scores[n],scores[n+1],scores[+2]))
            sudscore = scores[n-4:n+3]
            max_index = np.argmax(sudscore) -4
            print(n, max(sudscore),max_index)


if __name__ == "__main__":


    # fin = '/data/nanopore/nanoDoc2_1/varidate/16S_score.txt'
    # ans = '/data/nanopore/nanoDoc2_1/testrun/16sans.txt'
    # out = '/data/nanopore/nanoDoc2_1/varidate/16sout2.png'
    # figsize = (16, 4)
    fin = '/data/nanopore/nanoDoc2_1/varidate/23S_score.txt'
    ans = '/data/nanopore/nanoDoc2_1/testrun/23sans.txt'
    out = '/data/nanopore/nanoDoc2_1/varidate/23sout2.png'
    figsize = (23, 4)
    plotGraph(fin,out,ans,figsize)
