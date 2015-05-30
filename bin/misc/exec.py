#!/usr/bin/env python
"""
Executes a command-line program with arguments that are interpolated from a set of variables.
Variables can have multiple values, in which case the program will be executed once with each
possible combination of values.

Example::

    exec.py -v foo=1,2 -v bar=A,B myprog.sh -f '%(foo)s' -b 'bar%(bar)s'

    Calls the program myprog.sh four times:

    myprog.sh -f 1 -b barA
    myprog.sh -f 1 -b barB
    myprog.sh -f 2 -b barA
    myprog.sh -f 2 -b barB
"""

from ConfigParser import SafeConfigParser
from glob import glob
import os
import sys

sys.path.insert(0, os.path.join(os.environ['SCI_HOME'], "common", "lib"))
from util.cl import parse
from util.collections import VarArgGenerator
from util.fork import *
from util.io import safe_read_file_array

# handlers may be called after each execution or a single time with all results, depending on
# whether map_async or apply_async is called

class FileHandler(object):
    def __init__(self, filename=None):
        self.filename = filename

    def __enter__(self):
        self.results = []

    def __call__(self, result):
        if isinstance(result, list):
            self.results = result
        else:
            self.results.append(result)

    def __exit__(self, exc_type, exc_value, traceback):
        if self.results:
            if self.filename:
                out = open(self.filename, 'w')
            else:
                out = sys.stdout
            with out:
                write_table(self.results, out)

class FilePatternHandler(object):
    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, result):
        def _append(r):
            filename = self.pattern.format(**r[0])
            with open(filename, 'w') as outfile:
                outfile.write(r[1])

        if isinstance(result, list):
            for r in result:
                _append(r)
        elif isinstance(result, 'tuple'):
            _append(r)

def write_table(results, out):
    def esc(x): return '"%s"' % str(x).replace('"', '\\"')
    out.write(",".join(map(esc, results[0][0].keys()) + ['Result']))
    out.write("\n")
    for r in results:
        out.write(",".join(map(esc, r[0].values() + [r[1]])))
        out.write("\n")

if __name__ == '__main__':
    # Find the first non-option, which will be the program to execute
    for i in xrange(1, len(sys.argv), 2):
        if not sys.argv[i].startswith('-'):
            break

    pipe_help = """
    Some process don't write to stdout, so if the output of a process needs to be piped to another
    command, it is necessary to use a FIFO. The --pipe_file_pattern option specifies the file that
    should be created as a FIFO, so that when the process opens and writes to it, it can be read
    and handled by the command specified in --pipe_command. The --pipe_command argument must be
    a command that can read from the FIFO as it is written to and process the entire output,
    otherwise the process will block and/or die. Note that each instance of pipe_command needs its
    own thread, and so 1) -t must be an even number >= 2, and 2) only -t/2 threads will actually be
    executing the process.
    """

    def add_arguments(parser):
        parser.add_argument("-e", "--error_mode", choices=("exit", "log", "ignore"), default="ignore",
            help="What to do with errors.")
        parser.add_argument("-f", "--fork_mode", choices=("test","serial","thread","lsf"), default="thread",
            help="Execution mode. Set to 'test' for command logging only.")
        parser.add_argument("-g", "--glob", metavar="NAME=GLOB", action="dict", type="mapping",
            help="Define a file glob that will be expanded.")
        parser.add_argument("-k", "--fork_opts", action="extend_dict", type="mapping_list",
            help="Options specific to the fork mode.")
        parser.add_argument("-o", "--outfile_pattern", metavar="PATTERN", default=None,
            help="Pattern from which log files (containing stdout and stderr from each process) "\
                "are created by interpolation with variable values.")
        parser.add_argument("-s", "--sequence", metavar="NAME=LOW,HIGH,STEP", action="dict",
            type="delimited_mapping", help="Define a sequence variable.")
        parser.add_argument("-v", "--var", action="dict", type="delimited_macro_mapping",
            metavar="NAME=VALUE",
            help="Define a variable. VALUE can be a single value, comma-delimited list or range.")
        parser.add_argument("-V", "--var_file", type="readable_file", metavar="FILE",
            help="File from which to read variables. Each line is NAME=VALUE.")
        parser.add_argument("-F", "--file_var", action="dict", type="mapping", metavar="NAME=FILE",
            help="Define a variable whose values come frome FILE")
        parser.add_argument("-w", "--sliding_window", metavar="NAME=LOW,HIGH,STEP", action="dict",
            type="delimited_mapping", help="Define a sliding window variable.")
        pipe_group = parser.add_argument_group("pipes", pipe_help)
        pipe_group.add_argument("--pipe_file_pattern", metavar="PATTERN", default=None,
            help="Process output file to send to pipe_command.")
        pipe_group.add_argument("--pipe_command_pattern", metavar="COMMAND", default=None,
            help="Command to run on piped process output.")
        output_group = parser.add_mutually_exclusive_group()
        output_group.add_argument("-r", "--result_file", type="writeable_file", metavar="FILE",
            help="File where all results are written in CSV format (one column for each variable "\
                "value followed by a column with the result). All fields are quoted.")
        output_group.add_argument("-R", "--result_file_pattern", metavar="PATTERN",
            help="Pattern from which output file is created by interpolation with variable values.")

    ns = parse(add_arguments, args=sys.argv[1:i])
    prog = sys.argv[i]
    args = sys.argv[i+1:]

    argvars = VarArgGenerator()

    if ns.var_file:
        config = SafeConfigParser()
        config.read(ns.var_file)

        if config.has_section('constants'):
            argvars.update(config.items('constants'))

        if config.has_section['variables']:
            argvars.update(dict(parse_vars(m) for m in config.items('variables')))

    if ns.var:
        argvars.update(ns.var)

    if ns.glob:
        for name, globstr in ns.glob.iteritems():
            argvars.variables[name] = glob(globstr)
            argvars.add_oper(name, "{0}_base".format(name), lambda x: os.path.basename(x))
            argvars.add_oper(name, "{0}_noext".format(name), lambda x: os.path.splitext(x)[0])
            argvars.add_oper(name, "{0}_base_noext".format(name), lambda x: os.path.basename(os.path.splitext(x)[0]))

    if ns.file_var:
        for name, path in ns.file_var.iteritems():
            argvars.variables[name] = safe_read_file_array(path)

    if ns.sequence:
        for name, seq in ns.sequence.iteritems():
            s = map(int, seq)
            r = range(*s)
            if r[-1] < s[1]:
                r += [s[1]]
            argvars.variables[name] = r

    if ns.sliding_window:
        for name, seq in ns.sliding_window.iteritems():
            s = map(int, seq)
            r = range(*s)
            if r[-1] < s[1]:
                r += [s[1]]
            argvars.variables[name] = zip(r[0:(len(r)-1)], r[1:len(r)])

    fork_mode = ns.fork_mode
    fork_opts = ns.fork_opts or {}
    threads = int(fork_opts["threads"]) if "threads" in fork_opts else None
    threaded = threads is not None and threads > 1
    
    if threaded:
        if ns.pipe_file_pattern:
            threads = min(threads, len(argvars) * 2)
            assert threads >= 2
        else:
            threads = min(threads, len(argvars))
        fork_opts["threads"] = threads
    elif fork_mode == "thread":
        fork_mode = "serial"

    executor = get_executor(fork_mode, fork_opts)

    result_handler = None
    if ns.result_file:
        result_handler = FileHandler(ns.result_file)
    elif ns.result_file_pattern:
        result_handler = FilePatternHandler(ns.result_file_pattern)

    error_handler = None
    if ns.error_mode == 'exit':
        error_handler = reraise_error
    elif ns.error_mode == 'log':
        error_handler = log_error

    if ns.pipe_file_pattern is not None:
        cmd_iter = command_iter_with_pipe(prog, args, argvars,
            ns.pipe_file_pattern, ns.pipe_command_pattern, safe=threaded)
    else:
        cmd_iter = simple_command_iter(prog, args, argvars, ns.outfile_pattern, safe=threaded)

    exec_shell(cmd_iter, executor, result_handler, error_handler)
