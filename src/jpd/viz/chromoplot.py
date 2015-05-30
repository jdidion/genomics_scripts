#!/usr/bin/env python
# encoding: utf-8

"""
Plot one or more data series on chromosome axes. Input can be in a single file
or one file per chromosome. The single input data file must be of the format: 
chr, x, y1, y2..., where y1, y2 are the data series. Per-chromosome files are
assumed not to have the chr field unless the --with_chr option is set. There is
one horizontal axis for each chromosome, but the x-axis tics, labels, etc only
appear at the bottom of the image. The x and y scales are kept the same for all
chromosomes. 

Series can be placed above/below the chromosome axis by making values positive
or negative. Options can also be used to specify the figure size, DPI, and 
title, and data series colors and lablels.
"""

from colors import rainbow
import csv
import fileinput
from math import ceil
import os
import os.path
import sys

import matplotlib
import matplotlib.colors
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, subplots_adjust
from matplotlib.font_manager import FontProperties
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Wedge
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator, MultipleLocator, FixedLocator

sys.path.append("../../lib/python")
import util.cl
from util.misc import compare_mixed

AUTOSOME = [str(i) for i in range(1, 20)]
SEX      = ['X', 'Y']
MT       = 'M'

class Chromoplot:
    def __init__(self, options=None, load=True):
        self.options = options
        if load: self.load_data()
        
    def load_data(self, src=self.options.data_file, opt=self.options):
        self.data = {}
        self.stats = { 'max_x': 0, 'max_pos_y': 0, 'max_neg_y': 0 }
        self.scale = opt.output_scale / opt.input_scale if opt.output_scale else opt.input_scale
        
        if opt.series_indices:
            self.nseries = len(opt.series_indices)
            self.pos_indices = []
            self.neg_indices = []
            for i in opt.series_indices:
                (self.pos_indices if i > 0 else self.neg_indices).append(abs(i)-1)
            self._init_series()
        else:
            self.pos_indices = self.neg_indices = None
            
        if os.path.isfile(src):
            self._load_data_file(src, pos_indices, neg_indices)
        else:
            for chrm in opt.chromosomes:
                self._load_data_file(os.path.join(src, 
                    self.options.file_name_format % chrm), chrm)
            
    def _init_series(self, opt=self.options):
        colors = opt.series_colors or rainbow(self.nseries)
        self.legend_colors = list(colors)
        # recycle colors if necessary
        while len(colors) < self.nseries:
            colors += colors[0:min(self.nseries-len(colors), len(colors))]
        self.pos_colors = [colors[i] for i in self.pos_indices]
        self.neg_colors = [colors[i] for i in self.neg_indices]
        
        labels = opt.series_labels or ["Series %s" % n for n in range(1, self.nseries)]
        self.legend_labels = list(labels)
        # recycle labels if necessary
        while len(labels) < self.nseries:
            labels += labels[0:min(self.nseries-len(labels), len(labels))]
        self.pos_labels = [labels[i] for i in self.pos_indices]
        self.neg_labels = [labels[i] for i in self.neg_indices]
                
    def _load_data_file(self, infile, chrm=None, opt=self.options):
        for row in csv.reader(open(infile)):
            if not chrm or opt.with_chr:
                chrm = row.pop(0)
            if chrm not in self.data:
                self.data[chrm] = { 'x': [] }    

            x = float(row[0])
            self.data[chrm]['x'].append(x)
            self.stats['max_x'] = max(self.stats['max_x'], x)
            
            y = [0 if n == '' or n == '\N' else float(n) for n in row[1:]]
            if self.pos_indices is None and self.neg_indices is None:
                self.nseries = len(y)
                self.pos_indices = range(0, self.nseries)
                self._init_series()
            if self.pos_indices:
                self._load_y(chrm, y, self.pos_indices, 'pos_y')
            if self.neg_indices:
                self._load_y(chrm, y, self.neg_indices, 'neg_y', True)
    
    def _load_y(self, chrm, y, indices, key, negate=False):
        elts = elements(y, indices)
        if negate:
            elts = [e * -1 for e in elts]

        if key not in self.data[chrm]:
            self.data[chrm][key] = [[] for i in range(0, len(elts))]
        matrix_append(self.data[chrm][key], elts)
        
    def plot(self, figure_file=self.options.figure_file, opt=self.options):
        def stacked_bar_vlines(fig, ax, x_range, data, opt):
            # TODO: deal with self.scale
            # TODO: return max_y, min_y
            x = data['x']
            
            if 'pos_y' in data:
                ymin = ymax = [0] * len(x)
                for y, c, l in zip(data['pos_y'], self.pos_colors, self.pos_labels):
                    ymin = ymax
                    ymax = y if ymin is None else sum_arrays(ymin, y)
                    ax.vlines(x, ymin, ymax, color=c, label=l)
                    
            if 'neg_y' in data:
                ymax = ymin = [0] * len(x)
                for y, c, l in zip(data['neg_y'], self.neg_colors, self.neg_labels):
                    ymax = ymin
                    ymin = y if ymax is None else sum_arrays(ymax, y)
                    ax.vlines(x, ymin, ymax, color=c, label=l)
        
        def stacked_bar_hist(fig, ax, x_range, data, opt):
            max_y = min_y = 0
            
            if 'pos_y' in data:
                x = zip(*[data['x']] * len(data['pos_y']))
                weights = zip(*data['pos_y'])
                n, bins, patches = ax.hist(x, bins=x_range[1] / self.scale, 
                    histtype='barstacked', range=x_range, weights=weights,
                    label=self.pos_labels)
                color_patches(self.pos_colors, patches)
                max_y = max([sum(l) for l in zip(*n)])
                
            if 'neg_y' in data:
                x = zip(*[data['x']] * len(data['neg_y']))
                weights = zip(*data['neg_y'])
                n, bins, patches = ax.hist(x, bins=x_range[1] / self.scale, 
                    histtype='barstacked', range=x_range, weights=weights,
                    label=self.neg_labels, bottom=-0.01)
                color_patches(self.neg_colors, patches)
                min_y = min([sum(l) for l in zip(*n)])

            return max_y, min_y
            
        def color_patches(colors, patches):
            if len(colors) == 1:
                for p in patches:
                    if p.get_height() != 0:
                        p.set_color(colors[0])
            else:
                for c, p in zip(colors, patches):
                    for r in p:
                        if r.get_height() == 0:
                            r.set_visible(False)
                        else:
                            r.set_color(c)
            
        def make_legend_patches(colors):
            patches = []
            for c in colors:
                patches.append(Rectangle((0,0), 0.01, 0.01, color=c, fill=True))
            return patches
                
        data = self.data
        stats = self.stats
        graph_fn = stacked_bar_hist if opt.bar_graph_mode == 'hist' else stacked_bar_vlines
        
        # intersection
        chromosomes = sorted(filter(lambda x:x in self.chromosomes, self.data.keys()), compare_mixed)
        nchr = len(chromosomes)
        chrm_sizes = parse_properties(opt.chromosome_sizes_file, fn=float)
        max_x=self.stats['max_x'] * 1.01
        
        fig = figure(figsize=(opt.image_width, opt.image_height), dpi=opt.dpi, facecolor='w')
        fig.suptitle(opt.title, fontsize='x-large')
        #TODO: adjust margins, label and legend positions based on figsize
        subplots_adjust(left=.06, bottom=.10, right=.95, top=.95, hspace=.25, wspace=0)
        axes = []
        
        for c, chrm in enumerate(self.chromosomes, start=1):
            print "Plotting chromosome %s" % chrm

            first = c == 1
            last = c == nchr
            ax = fig.add_subplot(nchr, 1, c, axis_bgcolor='0.9')
            axes.append(ax)
            x_range = (0, chrm_sizes[chrm] / opt.input_scale)

            # plot data
            max_y, min_y = graph_fn(fig, ax, x_range, data[chrm], opt)
            stats['max_pos_y'] = max(max_y, stats['max_pos_y'])
            stats['max_neg_y'] = min(min_y, stats['max_neg_y'])
            
            # adjust spines
            for loc, spine in ax.spines.iteritems():
                if first and loc == 'right':
                    spine.set_position(('outward', 10.5))
                elif last and loc == 'bottom':
                    spine.set_position(('outward', 10))
                else:
                    spine.set_color('none') # don't draw spine

            # format y axis
            ax.yaxis.set_label_text(chrm)
            ax.yaxis.set_label_coords(0.01, 0.5)
            ax.yaxis.get_label().set_rotation(0)
            if first:
                ax.yaxis.set_ticks_position('right')
                ax.yaxis.set_minor_locator(FixedLocator([0]))
            else:
                ax.yaxis.set_ticks([])
            
            # format x axis
            ax.set_xlim(xmin=0, xmax=max_x)
            if last:
                ax.xaxis.set_ticks_position('bottom')
                ax.xaxis.set_major_locator(MultipleLocator(50000 / opt.input_scale))
                ax.xaxis.set_minor_locator(MultipleLocator(10000 / opt.input_scale))
                ax.xaxis.set_label_text(opt.x_axis_label % opt.output_scale,
                    fontproperties=FontProperties(size='large'))
            else:
                ax.xaxis.set_ticks([])
                    
            # draw chromosome axis
            ax.add_line(Line2D(x_range, (0, 0), color='k'))
            
        max_y = stats['max_pos_y']
        min_y = stats['max_neg_y']
        grid_height = max_y + abs(min_y) + 10
        yax_radius = (float(max_y + abs(min_y)) / 2)
        yax_center = max_y - yax_radius
        grid_loc = axes[-1].xaxis.get_ticklocs(True)[1:]
        for i, ax in enumerate(axes):
            # set y scale
            ax.set_ylim(min_y, max_y)
            # add embelishments
            #ax.add_line(Line2D((-1 * yax_radius, 0), (0, 0), color='k',
            #    clip_on=False))
            ax.add_patch(Wedge((0, yax_center), yax_radius-1, 90, 270,
                color='0.9', fill=True, clip_on=False))
            ax.add_patch(Wedge((max_x, yax_center), yax_radius-1, 270, 90, 
                color='0.9', fill=True, clip_on=False))
            # draw x grid
            if opt.grid:
                height = max_y if i == 0 else grid_height
                for loc in grid_loc:
                    ax.add_line(Line2D((loc, loc), (0, height), color='0.5', 
                        linestyle=':', clip_on=False))
        min_tick = min(min_y, 0)
        max_tick = max(max_y, 0)
        axes[0].yaxis.set_ticks([min_tick, max_tick])
        axes[0].yaxis.set_ticklabels([str(abs(int(min_tick))), str(int(max_tick))])
        
        # add legend
        cols = int(ceil(float(self.nseries)/2))
        fig.legend(make_legend_patches(self.legend_colors), 
            self.legend_labels, loc='lower center', ncol=cols, 
            bbox_to_anchor=(0, 0, 1, 1), bbox_transform=fig.transFigure,
            prop=FontProperties(size='medium'))
        
        # add some figure labels
        fig.text(.01, .5, 'Chromosome', rotation=90, verticalalignment='center',
            fontproperties=FontProperties(size='large'))
        if opt.scale_title: 
            fig.text(.97, .5, opt.scale_title, rotation=270, 
                verticalalignment='center', 
                fontproperties=FontProperties(size='large'))
        
        fig.savefig(figure_file)

def main(argv=None):
    def add_opts(parser):
        parser.add_argument("-c", "--series_colors", type="str_list", metavar="LIST", default=None,
            help="Colors for data series, in the same order as series indices.")
        parser.add_argument("-g", "--grid", action="store_true", default=False, 
            help="Draw x-axis grid.")
        parser.add_argument("-i", "--series_indices", type="int_list", metavar="LIST", default=None, 
            help="Indices of data series, where index 1 is the first y value. Negative indices are "\
                "plotted below the chromosome axis.")
        parser.add_argument("-l", "--series_labels", metavar="LIST", type="str_list", default=None, 
            help="Labels for data series, in the same order as series-indices.")
        parser.add_argument("-p", "--file_name_format", metavar="FORMAT", default="chr%s.csv", 
            help="Format for per-chromosome data files if the supplied data source is a directory.")
        parser.add_argument("-r", "--chromosomes", metavar="LIST", action="overwrite", default=["N"], 
            help="Chromosomes to plot; * = all, A = autosome, S = sex, AX = autosome + X, "\
                "N = nuclear (all except M), or a comma-delimited list")
        parser.add_argument("-R", "--chromosome_sizes_file", metavar="FILE", type="readable_file",
            help="File with chromosome sizes (chr=size)")
        parser.add_argument("-s", "--input_scale", type=int, metavar="SCALE", default=100, 
            help="Scale of input data, in Kb.")
        parser.add_argument("-S", "--output_scale", type=int, metavar="SCALE",
            help="Scale of output data, in Kb.")
        parser.add_argument("-t", "--title", metavar="TITLE",
            help="Figure title.")
        parser.add_argument("-T", "--scale_title", metavar="TITLE",
            help="Label on right side to describe the scale.")
        parser.add_argument("-x", "--x_axis_label", metavar="LABEL", 
            default="Genomic Position (%iKb Windows)", 
            help="X-axis label.")
        parser.add_argument("--image_height", type=float, metavar="INCHES", default=11, 
            help="Image height in inches.")
        parser.add_argument("--image_width", type=float, metavar="INCHES", default=8.5, 
            help="Image width in inches.")
        parser.add_argument("--dpi", type=int, metavar="DPI", default=100, 
            help="Figure resolution.")
        parser.add_argument("--bar_graph_mode", metavar="MODE", default="hist", 
            help="Graph drawing mode; can be 'hist' or 'vlines'.")
        parser.add_argument("--with_chr", action="store_true", default=False,
            help="The input file has a chromosome field.")
        parser.add_argument("data_source", type="readable_path", metavar="FILE",
            help="Input data file.")
        parser.add_argument("figure_file", type="writeable_file", metavar="FILE",
            help="Output file.")

    args = parseopt.parse_options(add_opts, desc="Plot data on chromosomes.", version=1, args=argv)
    Chromoplot(args).plot()

if __name__ == '__main__':
    main()    