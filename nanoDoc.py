import click
import os
import nanoDoc2_1.analysis.comparisonAnalysisKmean as nanoDocAnalysis
import pathlib
import nanoDoc2_1.preprocess.fast5ToProcessedPq as fast5ToProcessedPq
import nanoDoc2_1.utils.PlotGraph as PlotGraph

@click.group()
def cmd():
    pass


@cmd.command()
@click.option('-i', '--input',required='True')
@click.option('-o', '--output',required='True')
@click.option('-r', '--ref',required='True')
@click.option('-fm', '--fmercurrent',required='True')
@click.option('-t', '--thread', default=12)
@click.option('-qv', '--qvaluethres', default=5)
@click.option('-mo', '--mappyoption', nargs=4, type=click.Tuple([int,int,int,int]), default=(12, 10, 30, 20))
def fast5ToReSegmentedPq(input,output,ref,fmercurrent,thread,qvaluethres,mappyoption):

    print(input,output)
    click.echo('make resegmanted pq file')
    MAX_CORE = thread
    fast5ToProcessedPq.h5tosegmantedPq(input,output,ref,MAX_CORE,qvaluethres,fmercurrent,mappyoption)


@cmd.command()
@click.option('-w', '--weight',required='True')
@click.option('-r', '--ref',required='True')
@click.option('-rpq', '--rpq',required='True')
@click.option('-tgpq', '--tgpq',required='True')
@click.option('-o', '--output',required='True')
@click.option('-tsid', '--transcriptid',default="")
@click.option('-s', '--start',default=1)
@click.option('-e', '--end',default=-1)
@click.option('-minreadlen', '--minreadlen',default=500)
def analysis(weight,ref,rpq,tgpq,output,transcriptid,start,end,minreadlen):

    click.echo('modification call')
    p_sub = pathlib.Path(output)
    if not os.path.exists(p_sub.parent):
        os.mkdir(p_sub.parent)
    print(rpq)
    print(tgpq)
    nanoDocAnalysis.modCall(weight, ref, rpq,tgpq, output, transcriptid, transcriptid, start, end,minreadlen)

@cmd.command()
@click.option('-f', '--fin',required='True')
@click.option('-o', '--fout')
@click.option('-a', '--ans')
@click.option('-s', '--figsize', nargs=2, type=click.Tuple([int,int]), default=(16, 4))
def plotGraph(fin, fout, ans, figsize):

    if fout is None:
        fout = os.path.dirname(fin) +"/result.png"

    PlotGraph.plotGraph(fin, fout, ans, figsize)



def main():
    cmd()


if __name__ == '__main__':
    main()