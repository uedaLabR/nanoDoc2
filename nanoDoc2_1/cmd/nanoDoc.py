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
@click.option('-i', '--input')
@click.option('-o', '--output')
@click.option('-r', '--ref')
@click.option('-fm', '--fmercurrent')
@click.option('-t', '--thread', default=12)
@click.option('-qv', '--qvalueThres', default=5)
def fast5ToReSegmentedPq(input,output,ref,fmercurrent,thread,qvalueThres):

    print(input,output)
    click.echo('make resegmanted pq file')
    MAX_CORE = thread
    fast5ToProcessedPq.h5tosegmantedPq(input,output,ref,MAX_CORE,qvalueThres,fmercurrent)


@cmd.command()
@click.option('-w', '--wight')
@click.option('-p', '--param')
@click.option('-r', '--ref')
@click.option('-rraw', '--refraw')
@click.option('-traw', '--tgraw')
@click.option('-o', '--output')
@click.option('-chrom', '--chrom',default="")
@click.option('-s', '--start',default=1)
@click.option('-e', '--end',default=-1)
@click.option('-st', '--strand',default="+")
@click.option('-minreadlen', '--minreadlen',default=200)
def analysis(wight,param,ref,refraw,tgraw,output,chrom,start,end,strand,minreadlen):

    click.echo('modification call')
    p_sub = pathlib.Path(output)
    if not os.path.exists(p_sub.parent):
        os.mkdir(p_sub.parent)
    print(refraw)
    print(tgraw)
    nanoDocAnalysis.modCall(wight,param, ref, refraw,tgraw, output, chrom, chrom, start, end, strand, minreadlen)

@cmd.command()
@click.option('-f', '-fin')
@click.option('-o', '-fout')
@click.option('-a', '--ans')
@click.option('-s', '--figsize', nargs=2, type=int, default=(16, 4))
def plotGraph(fin, fout, ans, figsize):

    PlotGraph.plotGraph(fin, fout, ans, figsize)



def main():
    cmd()


if __name__ == '__main__':
    main()