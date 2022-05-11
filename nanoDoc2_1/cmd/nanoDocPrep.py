import click
import SamplingPlan
import to5merpq
import os
import initialtrainingCNN
import make5merwight
import itertools
import os
import makeIndex

@click.group()
def cmd():
    pass


import nanoDoc2_1.trainprep.make6merParquet as make6merParquet

@cmd.command()
@click.option('-r', '--refs')
@click.option('-p', '--pqs')
@click.option('-o', '--out')
@click.option('-j', '--join',default=False)
@click.option('-takeCnt', '--samplesize',default=10000)
def make6mer(refs,pqs,out,takeCnt,join):

    make6merParquet.makeSamplePlan(refs, pqs, out, takeCnt,join)


@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-o', '--outwight')
@click.option('-ssize', '--samplesize',default=1200)
@click.option('-epochs', '--epochs',default=500)
def traincnn(in5mmer,outwight,samplesize,epochs):

    click.echo('trainCNN')
    initialtrainingCNN.main(in5mmer,outwight,samplesize,epochs)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-o', '--outwight')
@click.option('-ssize', '--samplesize',default=1200)
@click.option('-epochs', '--epochs',default=500)
def traincnnAdd(in5mmer,outwight,samplesize,epochs):

    click.echo('trainCNN')
    initialtrainingCNN.main(in5mmer,outwight,samplesize,epochs)

@cmd.command()
@click.option('-in', '--in5mmer')
@click.option('-o', '--outwight')
@click.option('-wight', '--bestwight')
@click.option('-ssize', '--samplesize',default=12000)
@click.option('-epochs', '--epochs',default=3)
def traindoc(in5mmer,outwight,bestwight,samplesize,epochs):

    click.echo('trainDoc')
    nucs = ('A','T','C','G')
    for n1,n2,n3,n4,n5,n6 in itertools.product(nucs, nucs, nucs,nucs, nucs, nucs):

        nuc = n1+n2+n3+n4+n5+n6
        print('training doc',nuc)
        make5merwight.train(in5mmer,outwight,nuc,bestwight,samplesize,epochs)

def main():
    cmd()


if __name__ == '__main__':
    main()