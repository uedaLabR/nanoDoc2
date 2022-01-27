from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt


class GraphManager():
    def __init__(self, output_file_name):
        name = output_file_name + '.pdf'
        print("output at %s" % name)
        self.pdf = PdfPages(name)

    def add_figure(self, figure):
        self.pdf.savefig(figure)
        figure.clf()
        plt.close(figure)

    def save(self):
        self.pdf.close()
