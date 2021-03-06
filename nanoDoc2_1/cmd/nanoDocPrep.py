import click

import os
import nanoDoc2_1.training.initialTraning as initialTraning
import nanoDoc2_1.training.initialTrainingDec as initialTrainingDec
import nanoDoc2_1.training.docTraining6mer as docTraining6mer
import itertools
import os


@click.group()
def cmd():
    pass


import nanoDoc2_1.trainprep.make6merParquet as make6merParquet

@cmd.command()
@click.option('-r', '--refs',multiple=True)
@click.option('-p', '--pqs',multiple=True)
@click.option('-o', '--out')
@click.option('-j', '--join',default=False)
@click.option('-takeCnt', '--samplesize',default=10000)
def make6mer(refs,pqs,out,takeCnt,join):


    make6merParquet.makeSamplePlan(refs, pqs, out, takeCnt,join)


@cmd.command()
@click.option('-in', '--in6mmer')
@click.option('-o', '--outdir')
@click.option('-epochs', '--epochs',default=200)
@click.option('-device', '--device')
def traincnn(in6mmer,outdir,samplesize,epochs,device='/GPU:0'):

    click.echo('trainCNN')
    initialTraning.main(in6mmer,outdir,epochs,device)

@cmd.command()
@click.option('-in', '--in6mmer')
@click.option('-o', '--outdir')
@click.option('-epochs', '--epochs',default=50)
@click.option('-device', '--device')
def traincnnAdd(in5mmer,outdir,epochs,device='/GPU:0'):

    click.echo('trainCNN Dec')
    initialTrainingDec.main(in5mmer,outdir,epochs,device)


@cmd.command()
@click.option('-d1', '--data1')
@click.option('-d2', '--data2')
@click.option('-o', '--outdir')
@click.option('-in', '--weightdir')
@click.option('-ssize', '--samplesize',default=12750)
@click.option('-epochs', '--epochs',default=3)
@click.option('-device', '--device')
def traindoc(data1, data2,outdir,weightdir,samplesize,epochs,device='/GPU:0'):

    click.echo('trainDoc')
    nucs = ('A','T','C','G')
    for n1,n2,n3,n4,n5,n6 in itertools.product(nucs, nucs, nucs,nucs, nucs, nucs):

        nuc = n1+n2+n3+n4+n5+n6
        print('training doc',nuc)
        docTraining6mer.train(data1, data2, outdir, nuc, weightdir, samplesize, epochs,device)



def main():
    cmd()


if __name__ == '__main__':
    main()