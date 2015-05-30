import Image, ImageDraw, ImageFont, ImageOps
import math
import numpy as np
from util.collections import which
from util.io import abspath
from util.misc import cumsum

DEFAULT_CMAP = [
    (  0,   0,    0),    # 0.0=black
    (  0,   0,   255),   # 0.25 = blue
    (  0, 255, 255),     # 0.5 = cyan
    (255, 255,   0),     # 0.75 = yellow
    (255,   0,   0)      # 1.0 = red
]
    
def heatmap(val, cmap):
    last = len(cmap)-1
    val = last*val
    index = int(val)
    if (index >= last):
        return cmap[last]
    elif (index < 0):
        return cmap[0]
    else:
        t = val - index
        s = 1.0 - t
        r = int(cmap[index][0]*s + cmap[index+1][0]*t + 0.5)
        g = int(cmap[index][1]*s + cmap[index+1][1]*t + 0.5)
        b = int(cmap[index][2]*s + cmap[index+1][2]*t + 0.5)
        return (r, g, b)

# diagonal space reserved at the bottom for labels
def rotate(x): return math.ceil(math.sqrt(math.pow(x, 2) / 2))
OFFSETS = (0, 190, 200, 400, 410)
X_OFFSETS = map(lambda x: int(rotate(x)), OFFSETS)
Y_OFFSETS = map(lambda x: int(rotate(OFFSETS[-1] - x)), OFFSETS)

def plot_ld_heatmap(mat_file, bin_file, outfile, window_size=500000, summary_stat="unlinked_pct", 
        cmap=DEFAULT_CMAP, missing_color=(255, 255, 255), bg_color=(255, 255, 255), chrm_color=(0, 0, 0),
        tics=False, font_file="/Library/Fonts/Microsoft/Arial.ttf", font_size=80):
    """
    Create a triangle heatmap plot from a numpy matrix of pairwise-bin values. Each cel of the
    matrix is a summary statistic (e.g. mean or 95th percentile) of the LD between the two bins.
    The matrix can either be square or upper-right-hand. bins is a three-column matrix:
    chrmosome, start position, summary value.
    
    TODO: support creating diamond heatmaps from a square matrix.
    """
    
    import csv
    def load_matrix(x, first):
        with open(abspath(x), "rU") as i:
            r = csv.reader(i, delimiter=",")
            h = r.next()[first:]
            m = np.array(list(map(lambda x: float("nan" if x == "NA" else x), row[first:]) for row in r), dtype=np.float64)
            return (h, m)

    mat = load_matrix(mat_file, 1)[1]
    #mat = mat[0:5063,0:5063]
    bin_cols, bins = load_matrix(bin_file, 0)
    #bins = bins[0:5063,]
    chrm_col = which("chrm", bin_cols)
    stat_col = which(summary_stat, bin_cols)

    chrms = bins[:,chrm_col]
    chrm_idxs = np.unique(chrms, return_index=True)[1]
    chrms = [int(chrms[index]) for index in sorted(chrm_idxs)]

    # TODO: this is a hack. Fix the process_ld_file function to write blank lines for missing windows
    row = 0
    nan = float("nan")
    blank = [nan] * (bins.shape[1]-3)
    for chrm in chrms:
        pos = 0
        while row < bins.shape[0] and int(bins[row, chrm_col]) == chrm:
            if int(bins[row, 1]) != pos:
                bins = np.insert(bins, row, [chrm, pos, 0] + blank, axis=0)
                mat = np.insert(mat, row, nan, axis=0)
                mat = np.insert(mat, row, nan, axis=1)
            row += 1
            pos += window_size

    mat_size = mat.shape[0]
    assert mat_size == mat.shape[1], "Matrix is not square"
    assert mat_size == len(bins), "Matrix and bins are of different size"
            
    def convert_chrm(c):
        c = int(c)
        return "X" if c == 23 else str(c)
    chrms = map(convert_chrm, bins[...,chrm_col])

    im = Image.new("RGB", (mat_size + X_OFFSETS[4], mat_size + Y_OFFSETS[0]), bg_color)
    pixels = im.load()
    draw = ImageDraw.Draw(im)
    font = ImageFont.truetype(font_file, font_size)

    chrm_id = None
    chrm_end = 0
    tic_interval = 10000000 / window_size
    
    for row, summary in enumerate(bins[...,stat_col]):
        chrm = chrms[row]

        def get_line_coords(start, end):
            x1 = row + X_OFFSETS[start]
            y1 = row + Y_OFFSETS[start]
            x2 = row + X_OFFSETS[end]
            y2 = row + Y_OFFSETS[end]
            return ((x1, y1), (x2, y2))
            
        # draw summary line
        color = heatmap(summary, cmap) if not np.isnan(summary) and summary > 0 else 0
        draw.line(get_line_coords(2, 3), fill=color)
        
        for col in xrange(row, mat_size):
            x = col + X_OFFSETS[4]
            val = mat[row,col]
            if np.isnan(val) or val < 0:
                pixels[x, row] = missing_color
            else:
                pixels[x, row] = heatmap(val, cmap)
        
        # Label chromosome boundaries
        if chrm != chrm_id or row == (mat_size - 1):
            draw.line(get_line_coords(0, 2), fill=chrm_color, width=5)
            
            if chrm_id is not None:
                # first draw text in a blank image, then rotate and paste into the main image
                txt = str(chrm_id)
                txt_size = draw.textsize(txt, font=font)
                txt_img = Image.new("L", txt_size)
                txt_draw = ImageDraw.Draw(txt_img)
                txt_draw.text((0, 0), txt, font=font, fill=255)
                txt_img = txt_img.rotate(-45, expand=1)
                offset = max(txt_size[0], txt_size[1])
                lft_indent = int(math.floor((txt_img.size[0] - offset) / 2))
                top_indent = int(math.floor((txt_img.size[1] - offset) / 2))
                txt_img = txt_img.crop((lft_indent, top_indent, 
                    txt_img.size[0] - lft_indent, txt_img.size[1] - top_indent))
                
                text_row = int(chrm_end + ((row - chrm_end) / 2))
                text_x = int(text_row + X_OFFSETS[0] + ((X_OFFSETS[2] - X_OFFSETS[0]) / 2) - (txt_img.size[0] / 2))
                text_y = int(text_row + Y_OFFSETS[0] + ((Y_OFFSETS[2] - Y_OFFSETS[0]) / 2) - (txt_img.size[1] / 2))
                im.paste(ImageOps.colorize(txt_img, bg_color, chrm_color), (text_x, text_y))
            
            chrm_id = chrm
            chrm_end = row
        
        # add tics
        if tics and ((row - chrm_end) % tic_interval) == 0:
            draw.line(get_line_coords(1, 2), fill=chrm_color, width=3)
    
    # convert to RGBA (required for background recoloring)
    im2 = im.convert("RGBA")
    
    # rotate
    im2 = im2.rotate(45, expand=True)
    
    # recolor background
    bg = Image.new("RGBA", im2.size, (255,)*4)
    im2 = Image.composite(im2, bg, im2)
    
    # convert back to RGB
    im = im2.convert(im.mode)
    
    # crop
    final_height = int(math.ceil(math.sqrt(2 * math.pow(mat_size, 2)) / 2)) + OFFSETS[-1]
    im = im.crop((0, 0, im.size[0], final_height))
    
    # save
    im.save(outfile)
