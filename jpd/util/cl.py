"""
Simplified and enhanced interface to optparse, as well as validation methods.
"""
from argparse import ArgumentParser, FileType, ArgumentError, _AppendAction
import operator
import os
import re
from ConfigParser import ConfigParser

from .collections import charrange, PropDict
from .io import abspath, check_path, resolve_path, safe_read_file_array, parse_properties, csv_to_table
from .log import *

MACROS = dict(
    chrm={ "A" : range(1,20), "S" : ("X","Y"), "F": ("A","X"), "N": ("A","S"), "*" : ("N","M") },
    chrm_plink={ "A" : range(1,20), "X" : (23,), "Y" : (24,), "M" : (26,), "S" : ("X", "Y"), "F": ("A", "X"), "N": ("A","S"), "*" : ("N","M") },
    strain={
        "cccommon": ("129S1", "A_J", "CAST", "NOD", "NZO", "PWK", "WSB"),
        "cc": ("cccommon", "B6"),
        "common": ("cccommon","AKR","BALB","C3H","C57BL","CBA","DBA","LP_J"),
        "sanger": ("common","129P2","129S5","SPRET")
    }
)

class MyArgumentParser(ArgumentParser):
    def _check_value(self, action, value):
        # converted value must be one of the choices (if specified)
        if action.choices is not None:
            if isinstance(value, list):
                for v in value:
                    self._check_value(action, v)
            elif value not in action.choices:
                tup = value, ", ".join(map(repr, action.choices))
                msg = ("invalid choice: %r (choose from %s)") % tup
                import traceback
                print "".join(traceback.format_stack()[:-2])
                raise ArgumentError(action, msg)

def parse(add_func=None, version=None, unprocessed=None, auto_logging=True, args=None, actions=None,
          types=None, namespace=PropDict(), **kwargs):
    """Parses command line arguments using :module:`argparse`.

    This method optionally adds several arguments: for logging, version information, help,
    and for handing unprocessed arguments. Arguments can take advantage of the validation
    functions provided in this module.

    Returns a :class:`argparse.Namespace` object from which argument values can be retrieved

    Args:
        add_func (function): Called with the constructed ArgumentParser after default arguments
            are added.
        version: Program version used for printing help documentation. If specified, a --version
            argument is added.
        unprocessed (string): The name of the attribute in the returned Namespace to hold
            unprocessed arguments. If None, unprocessed arguments are not allowed and an
            exception is thrown if any are encountered.
        auto_logging (bool): If True, --log_level and --log_file options are added to the
            parser automatically and logging is configured (using util.log).
        args (iterable): Command line arguments to process. Defaults to sys.argv[1:].
        actions (dict): Additional actions to add to the parser.
        types (dict): Additional types to add to the parser.
        kwargs: any keyword arguments to pass to the ArgumentParser constructor.

    Example::
        def add_args(parser):
            parser.add_argument("-d", metavar="DIR", type="string")
        ns = parse(add_args, args=["-d", "/my/genotype/dir"])
        print ns.directory
    """

    parser = MyArgumentParser(**kwargs)
    _register_actions(parser, actions)
    _register_types(parser, types)

    if auto_logging:
        extend_logging()
        parser.add_argument("--log_level", metavar="LEVEL", default="WARNING",
            choices=get_level_names(), help="Logger level (default=WARNING)")
        parser.add_argument("--log_file", type="writeable_file", metavar="FILE",
            help="File where log messages should be written (default=stdout)")

    if version is not None:
        parser.add_argument("--version", action="version", version=version)

    # Let the caller add additional options. We do this before adding the unprocessed option
    # in case the caller wants to consume any positional args.
    if callable(add_func): add_func(parser)

    if unprocessed is not None:
        parser.add_argument(unprocessed, nargs="*")

    ns = parser.parse_args(args, namespace)

    if auto_logging:
        conf(ns.log_file, ns.log_level)

    return ns

def parse_with_config(*args, **kwargs):
    """
    Parse arguments that have an option named 'config' for a config file. Returns the merge of
    the config dict and the remainder of the command line args.
    """
    args = parse(*args, **kwargs)
    if args.config:
        conf = args.config
        del args.config
        conf.update(args)
        return conf
    else:
        return args

# Actions

class OverwriteAction(_AppendAction):
    """Like append, but overwrites (rather than appends to) default values."""
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest) == self.default:
            setattr(namespace, self.dest, None)
        _AppendAction.__call__(self, parser, namespace, values, option_string)

class ExtendAction(_AppendAction):
    """Like append, but uses extend rather than append when the value is a list."""
    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest, None) is None:
            items = []
        else:
            import copy
            items = copy.copy(getattr(namespace, self.dest))
        if isinstance(values, list):
            items.extend(values)
        else:
            items.append(values)
        setattr(namespace, self.dest, items)

class ExtendOverwriteAction(ExtendAction):
    """Like extend, but overwrites (rather than appends to) default values."""
    def __call__(self, parser, namespace, values, option_string=None):
        default = self.default
        if (callable(self.type)):
            default = self.type(default)
        if getattr(namespace, self.dest) == default:
            setattr(namespace, self.dest, None)
        ExtendAction.__call__(self, parser, namespace, values, option_string)

class DictAction(_AppendAction):
    """
    Like append, but assumes each value will be a tuple. Creates a dict at dest in which each
    key is value[0] and each value is value[1].
    """
    def __call__(self, parser, namespace, values, option_string=None):
        dest = getattr(namespace, self.dest)
        if dest == self.default:
            dest = {}
            setattr(namespace, self.dest, dest)
        assert isinstance(values, tuple), "DictAction only works with type=mapping"
        dest[values[0]] = values[1]

class ExtendDictAction(_AppendAction):
    """Like DictAction, but expects an interable of tuples rather than a single tuple."""
    def __call__(self, parser, namespace, values, option_string=None):
        dest = getattr(namespace, self.dest)
        if dest == self.default:
            dest = {}
            setattr(namespace, self.dest, dest)
        for k,v in values:
            dest[k] = v

_actions = dict(
    overwrite=OverwriteAction,
    extend=ExtendAction,
    extend_overwrite=ExtendOverwriteAction,
    dict=DictAction,
    extend_dict=ExtendDictAction
)
"""A dict of name-to-class mappings for all custom actions defined in this module."""

def _register_actions(parser, user_defined=None):
    """Add all actions defined in ``_actions`` to the parser."""
    actions = _actions
    if user_defined:
        actions = dict(actions)
        actions.update(user_defined)
    for k,v in actions.iteritems(): parser.register("action", k, v)

# Types

class CompositeType(object):
    def __init__(self, *types):
        self.types = types

    def __call__(self, string):
        result = string
        for t in self.types:
            result = t(result)
        return result

class TypeWithArgs(object):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs

    def __call__(self, string):
        return self._do_call(string, *self.args, **self.kwargs) or string

class ListWrapper(object):
    def __init__(self, data_type, collapse=True):
        self.data_type = data_type
        self.collapse = collapse

    def __call__(self, arg):
        if isinstance(arg, list):
            result = []
            for i in arg:
                v = self.data_type(i)
                if isinstance(v, list) and self.collapse:
                    result.extend(v)
                else:
                    result.append(v)
            return result
        else:
            return self.data_type(arg)

class delimited(TypeWithArgs):
    """Splits a string argument using a delimiter."""
    def _do_call(self, string, delim=",", range_delim="-", data_type=None, choices=None):
        if isinstance(string, str):
            vals = string.split(delim) if delim else (string,)
        else:
            vals = string
            
        if data_type is not None and not callable(data_type):
            data_type = _types[data_type]
        
        if vals[0] == "*" and choices is not None:
            vals = choices
        
        val_list = []
        for val in vals:
            if range_delim is None:
                val = (val,)
            else:
                val = val.split(range_delim)
            
            for v in _convert_range(val, data_type):
                if choices is not None and v not in choices:
                    raise ArgumentError(self, "invalid choice: {0} (choose from {1})".format(
                        v, ", ".join(map(repr, choices))))
                val_list.append(v)
        
        return val_list

class mapping(TypeWithArgs):
    def _do_call(self, string, delim="="):
        return tuple(string.split(delim))

delimited_mapping = CompositeType(mapping(), lambda x: (x[0], delimited()(x[1])))

mapping_list = CompositeType(delimited(range_delim=None), lambda l: tuple(mapping()(x) for x in l))

# Macro types

def expand_macro(name, value, macros):
    if name in macros:
        value = expand(value, macros[name])
    return value

def expand(value, value_dict):
    i = 0
    while i < len(value):
        if value[i] in value_dict:
            newv = value_dict[value[i]]
            value[i:i+1] = newv
        else:
            i += 1
    return value

delimited_macro_mapping = CompositeType(mapping(),
    lambda m: (m[0], expand_macro(m[0], delimited()(m[1]), MACROS)))

class delimited_macro(TypeWithArgs):
    def _do_call(self, value, macro, macros=MACROS):
        return expand_macro(macro, delimited()(value), macros)

# Validation types

class ComparisonValidator(TypeWithArgs):
    def _do_call(self, lhs, rhs, oper, expected=True):
        assert oper(lhs, rhs) == expected

def opened_file(mode="r", bufsize=None):
    return CompositeType(readable_file, FileType(mode, bufsize))

def existing_path(path):
    """Test that a path exists."""
    if path == "-":
        return path
    return resolve_path(path)

ACCESS = dict(r=os.R_OK, rU=os.R_OK, rb=os.R_OK, w=os.W_OK, wb=os.W_OK, x=os.X_OK)
class accessible_path(TypeWithArgs):
    """Test that a path is accessible"""
    def _do_call(self, path, type_, mode):
        if type_ == "f" and path == "-":
            return path
        return check_path(path, type_, ACCESS[mode])

class readable_file_group(TypeWithArgs):
    """Appends one or more extensions to a path prefix and checks that all files are accessible."""
    def _do_call(self, path, exts, mode="rU"):
        acc = accessible_path("f", mode)
        ret = []
        for e in exts:
            f = "{0}.{1}".format(path, e)
            f = existing_path(f)
            f = acc(f)
            ret.append(f)
        return ret

class writeable_path(object):
    def __init__(self, kind=None):
        self.kind = kind

    def __call__(self, p):
        if p == "-":
            return "-"
        p = abspath(p)
        try:
            path = resolve_path(p)
            check_path(path, self.kind, os.W_OK)
        except IOError:
            dirpath = os.path.dirname(p)
            if os.path.exists(dirpath):
                check_path(dirpath, "d", os.W_OK)
            else:
                os.makedirs(dirpath)
            path = os.path.join(dirpath, os.path.basename(p))
        return path

readable_dir = CompositeType(existing_path, accessible_path("d", "r"))
"""Test that a directory exists and is readable."""

readable_file = CompositeType(existing_path, accessible_path("f", "r"))
"""Test that a file exists and is readable."""

readable_path = CompositeType(existing_path, accessible_path(None, "r"))
"""Test that a file OR directory exists and is readable."""

writeable_file = writeable_path("f")
"""Test that a file exists and is writable, or doesn"t exist and its parent is writable."""

writeable_dir = writeable_path("d")
"""Test that a dir exists and is writable, or doesn"t exist and its parent is writable."""

readwriteable_dir = CompositeType(existing_path, accessible_path("d", "r"), writeable_path("d"))
"""Test that a dir exists and is readable and writeable."""

readwriteable_file = CompositeType(existing_path, accessible_path("f", "r"), writeable_path("f"))
"""Test that a file exists and is readable and writeable."""

executable_file=CompositeType(existing_path, accessible_path("f", "x"))
"""Test that a file exists and is executable."""

def file_glob(g):
    import glob
    return glob.glob(g)

def read_file(f):
    """Read contents of a file and return as an array of lines."""
    return safe_read_file_array(readable_file(f))

def read_csv(f):
    return csv_to_table(readable_file(f))

#def read_property_file(f, ordered=False):
#    """Read contents of a property file and return as a dict."""
#    config = parse_properties(readable_file(f), ordered=ordered)
#    config["__file__"] = f
#    return config

def read_property_file(f, defaults=None):
    config = ConfigParser(defaults, PropDict)
    config.read(readable_file(f))
    return config

int_fmt_re = re.compile("([\d\.]+)([KkMmGg]?)")
def int_fmt(s):
    """Similar to int(), but accepts K, M, and G abbreviations"""
    m = int_fmt_re.match(s)
    num, mult = m.groups()
    num = float(num)
    if mult is not None:
        mult = mult.lower()
        if mult == "k":
            num *= 1000
        elif mult == "m":
            num *= 1000000
        elif mult == "g":
            num *= 1000000000
        else:
            raise Exception("Unsupported multiplier {0}".format(mult))
    return int(num)

def positive(inclusive=False):
    """Test that a number is greater than (or equal to, if ``inclusive=True``) zero."""
    return ge(0) if inclusive else gt(0)

def gt(x):
    """Test that a number is greater than another number."""
    return CompositeType(int, ComparisonValidator(x, operator.gt))

def ge(x):
    """Test that a number is greater than another number."""
    return CompositeType(int, ComparisonValidator(x, operator.ge))

class list_apply(TypeWithArgs):
    def _do_call(self, l, func):
        return [func(e) for e in l]

_types = dict(
    infile=opened_file("r"),
    outfile=opened_file("w"),
    readable_dir=readable_dir,
    readable_path=readable_path,
    writeable_dir=writeable_dir,
    readable_file=readable_file,
    writeable_file=writeable_file,
    writeable_path=writeable_path(),
    readwriteable_dir=readwriteable_dir,
    readwriteable_file=readwriteable_file,
    file_glob=file_glob,
    read_file=read_file,
    read_csv=read_csv,
    property_file=read_property_file,
    pos_int=positive(),
    non_neg_int=positive(True),
    int_fmt=int_fmt,
    # these require the use of the dict action
    mapping=mapping(),
    delimited_mapping=delimited_mapping,
    delimited_macro_mapping=delimited_macro_mapping,
    # these require the use of the extend_dict action
    mapping_list=mapping_list,
    # these require the use of extend action (otherwise you get nested lists)
    str_list=delimited(data_type=str),
    int_list=delimited(data_type=int),
    float_list=delimited(data_type=float),
    readable_file_list=delimited(
        data_type=CompositeType(file_glob, ListWrapper(readable_file)), range_delim=None)
)
"""A dict of name-to-type mappings for custom types defined in this module."""

def _convert_range(rng, data_type=None):
    assert 1 <= len(rng) <= 2
    rng = map(str, rng)
    chars = all(len(s) == 1 for s in rng)

    if data_type is None:
        try:
            rng = map(int, rng)
            data_type = int
        except:
            if chars:
                data_type = str
            else:
                return rng
    elif data_type != str:
        rng = map(data_type, rng)

    if len(rng) == 1:
        return rng

    if data_type == int:
        return xrange(rng[0], rng[1] + 1)
    else:
        assert chars
        return charrange(rng[0], rng[1])

def _register_types(parser, user_defined):
    """Add all types defined in ``_types`` to the parser."""
    types = _types
    if user_defined:
        types = dict(types)
        types.update(user_defined)
    for k,v in types.iteritems(): parser.register("type", k, v)
