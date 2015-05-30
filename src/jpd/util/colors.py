"""
Colors
Several color-related functions borrowed from Sage.
TODO: rewrite this to use webcolors
"""

#*****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from math import floor, modf
import collections
from colorsys import hsv_to_rgb

def mod_one(x):
    """
    Reduce a number modulo 1.

    INPUT:

    - ``x`` - an instance of Integer, int, RealNumber, etc.; the
      number to reduce

    OUTPUT:

    - a float

    EXAMPLES::

        sage: from sage.plot.colors import mod_one
        sage: mod_one(1)
        1.0
        sage: mod_one(7.0)
        0.0
        sage: mod_one(-11/7)
        0.4285714285714286
        sage: mod_one(pi) + mod_one(-pi)
        1.0
    """
    x = float(x)
    if x != 1:
        x = modf(x)[0]
        if x < 0:
            x += 1
    return x

def float_to_html(r, g, b):
    """
    Converts a Red-Green-Blue (RGB) color tuple to a HTML hex color.
    Each input value should be in the interval [0.0, 1.0]; otherwise,
    the values are first reduced modulo one (see :func:`mod_one`).

    INPUT:

    - ``r`` - a number; the RGB color's "red" intensity

    - ``g`` - a number; the RGB color's "green" intensity

    - ``b`` - a number; the RGB color's "blue" intensity

    OUTPUT:

    - a string of length 7, starting with '#'

    EXAMPLES::

        sage: from sage.plot.colors import float_to_html
        sage: float_to_html(1.,1.,0.)
        '#ffff00'
        sage: float_to_html(.03,.06,.02)
        '#070f05'
        sage: float_to_html(*Color('brown').rgb())
        '#a52a2a'
        sage: float_to_html((0.2, 0.6, 0.8))
        Traceback (most recent call last):
        ...
        TypeError: float_to_html() takes exactly 3 arguments (1 given)
    """
    r, g, b = map(mod_one, (r, g, b))
    rr = "%x" % int(floor(r * 255))
    gg = "%x" % int(floor(g * 255))
    bb = "%x" % int(floor(b * 255))

    rr = '0' * (2 - len(rr)) + rr
    gg = '0' * (2 - len(gg)) + gg
    bb = '0' * (2 - len(bb)) + bb

    return '#' + rr + gg + bb


def rainbow(n, format='hex'):
    """
    Returns a list of colors sampled at equal intervals over the
    spectrum, from Hue-Saturation-Value (HSV) coordinates (0, 1, 1) to
    (1, 1, 1).  This range is red at the extremes, but it covers
    orange, yellow, green, cyan, blue, violet, and many other hues in
    between.  This function is particularly useful for representing
    vertex partitions on graphs.

    INPUT:

    - ``n`` - a number; the length of the list

    - ``format`` - a string (default: 'hex'); the output format for
      each color in the list; the other choice is 'rgbtuple'

    OUTPUT:

    - a list of strings or RGB 3-tuples of floats in the interval
      [0.0, 1.0]

    EXAMPLES::

        sage: from sage.plot.colors import rainbow
        sage: rainbow(7)
        ['#ff0000', '#ffda00', '#48ff00', '#00ff91', '#0091ff', '#4800ff', '#ff00da']
        sage: rainbow(7, 'rgbtuple')
        [(1.0, 0.0, 0.0), (1.0, 0.8571428571428571, 0.0), (0.28571428571428581, 1.0, 0.0), (0.0, 1.0, 0.57142857142857117), (0.0, 0.57142857142857162, 1.0), (0.2857142857142847, 0.0, 1.0), (1.0, 0.0, 0.85714285714285765)]

    AUTHORS:

    - Robert L. Miller

    - Karl-Dieter Crisman (directly use :func:`hsv_to_rgb` for hues)
    """
    R = []

    for i in range(n):
        R.append(tuple(map(float, hsv_to_rgb(float(i) / n, 1, 1))))

    if format == 'rgbtuple':
        return R
    elif format == 'hex':
        for j in range(len(R)):
            R[j] = float_to_html(*R[j])
        return R
