#!/usr/bin/env python

# Split a file into multiple files on a line-by-line basis. The simplest usage is to specify
# one or more pattern=output (-p) option, where the pattern is a regular expression and output
# is the output file. The output file can have placeholders ({0}, {1}, etc) that are filled in
# by values of corresponding capture groups from the regular expression. If more sophisticated
# processing is required (e.g. transforming each line), a python source file can be provided
# that contains a module variable called handlers, which is a list of objects that implement
# the LineHandler interface (see util.splitlib).

import os
import sys

sys.path.append(os.path.join(os.environ['LAB_HOME'], 'lib/python'))
from util.cl import parse
from util.io import DefaultFileHandler, split
from util.misc import load_module_from_file

def main(argv=None):
    def add_args(parser):
        parser.add_argument('-a', '--append', action='store_true', default=False)
        parser.add_argument('-h', '--header', default=0,
            help="Either an integer that is the number of header lines in the input file "\
                 "(use negative numbers to suppress copying of headers to output files) or "\
                 "a string that should be used as the header.")
        parser.add_argument('-P', '--pattern_file', type='readable_file', metavar="FILE")
        parser.add_argument('-p', '--pattern', type='mapping', action='append', 
            metavar="PATTERN=OUTPUT")
        parser.add_argument('-u', '--unmatched_file', type='writeable_file', metavar="FILE",
            help="File in which to write lines that do not match any pattern.")
        parser.add_argument('file', type='readable_file', metavar="FILE")
    
    ns = parse(add_args, args=argv)
    
    header = None
    if isinstance(ns.header, str):
        if ns.header.isdigit():
            header = int(ns.header)
        else:
            header = [ns.header]
    
    if ns.pattern_file:
        mod = load_module_from_file(ns.pattern_file)
        if not hasattr(mod, 'handlers'):
            raise Exception("Invalid pattern file: %s" % ns.pattern_file)
        else:
            handlers = mod.handlers
    else:
        handlers = (DefaultFileHandler(m[0], m[1], ns.append) for m in ns.pattern)

    if handlers:
        unmatched = None
        if ns.unmatched_file:
            unmatched = open(ns.unmatched_file, 'w')
        try:
            split(ns.file, handlers, header, unmatched)
        finally:
            if unmatched: unmatched.close()
    else:
        print "No handlers specified; doing nothing."
        
if __name__ == '__main__':
    main()