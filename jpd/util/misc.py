"""
Miscelaneous util methods.
"""

from contextlib import contextmanager
from functools import wraps
import logging
import math
import os
import sys
import types

# bitwise functions

def shift_sequence(n):
    """
    A generator that produces a sequence of integers where each is bit shifted one more position
    to the left.

    EXAMPLE::

        list(shift_sequence(5)) # => [1,2,4,8,16]
    """
    for i in range(0, n):
        yield 1 << i

def shift_to_int(s):
    """Returns the number of places a bit has been shifted left."""
    return int(math.log(s, 2))

def int_to_bitset(i):
    """Convert an integer to an array of bits."""
    bs = []
    while i >= 1:
        bs.append(i % 2)
        i = i >> 1
    bs.reverse()
    return bs

def bitset_to_int(b):
    """Convert a bitset to an integer."""
    b = b[b.index(1):len(b)]
    i = 0
    for flag in b:
        i = (i << 1) + flag
    return i

def or_list(l):
    """``or`` together all elements of a list."""
    if l is None or len(l) == 0: return 0
    return reduce(lambda x,y: x|y, l)


# math

def cumsum(seq):
    cum = 0
    for val in seq:
        cum += val
        yield cum

# object testing

def ismodule(obj):
    return type(obj) is types.ModuleType

def isfunction(obj):
    return type(obj) is types.FunctionType

def ismethod(obj):
    return type(obj) is types.MethodType

def isbuiltin(obj):
    return type(obj) is types.BuiltinFunctionType

def iscode(obj):
    return type(obj) is types.CodeType

def iterable(obj):
    return hasattr(obj, '__iter__')

def measurable(obj):
    return hasattr(obj, '__len__')

def ifnone(obj, default):
    return obj if obj is not None else default

# decorators

def accessor(func):
    """Decorates a property factory function.

    A property factory function returns a dict with args to property(). For example::

        class Foo:
            @prop
            def bar():
                def getf(self):
                    print "getf"
                    return self.__bar
                def setf(self, value):
                    print "setf"
                    self.__bar = value
                return dict(get(fget=getf, fset=setf))
        f=Foo()
        f.bar = 5
        print f.bar

    prints::

        getf
        setf
        5
    """
    return property(**func())

def attrs(func, **kwargs):
    """Decorator that adds attributes to a function."""
    for k,v in kwargs.iteritems():
        setattr(func, k, v)
    return func

def overrides(super_method):
    """
    Decorator for a method in a sub-class that overrides a method in its super-class.
    Automatically calls the super method before the overriding method.
    """
    def decorator(f):
        @wraps(f)
        def wrapper(self, *args, **kwargs):
            super_method(self, *args, **kwargs)
            f(*args, **kwargs)
        return wrapper
    return decorator

def compare_mixed(a, b, numbers_first=True):
    """Compare mixed alphanumeric values. Numbers come before characters by default."""
    if (type(a) == type(b)):
        return cmp(a,b)
    elif isinstance(a, int):
        return -1 if numbers_first else 1
    else:
        return 1 if numbers_first else -1

# other

def load_module_from_file(path):
    import imp
    mod_dir = os.path.dirname(path)
    mod_name = os.path.splitext(os.path.basename(path))[0]
    fp, mod_path, desc = imp.find_module(mod_name, [mod_dir])
    try:
        return imp.load_module(mod_name, fp, mod_path, desc)
    finally:
        if fp: fp.close()

@contextmanager
def buf(value=None):
    """Make creation and usage of a string buffer compatible with ``wait``."""
    from cStringIO import StringIO
    b = StringIO(value)
    try:
        yield b
    finally:
        b.close()

@contextmanager
def timeit(logfn=lambda x: sys.stdout.write(x), msg="Execution took {0:.2} seconds"):
    import time
    start = time.time()
    try:
        yield
    finally:
        finish = time.time()
        logfn(msg.format(finish - start))

class Label(object):
    """ContextManager that can be used to break out of nested loops.

    with Label() as search:
        for ...
            for ...
                search.break_()
    """
    class __break(Exception):
        def __init__(self, ctx):
            self.ctx = ctx

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return isinstance(exc_val, self.__break) and exc_val.ctx is self

    def break_(self):
        raise self.__break(self)

def bash(command, catch=True, safe=False, **popen_kwargs):
    """Execute a command using bash. Warning: this is vulnerable to shell command injection."""

    if not isinstance(command, str):
        command = ' '.join(command)
    logging.debug("Executing <{0}>".format(command))
    popen_kwargs.setdefault("shell", True)
    popen_kwargs.setdefault("executable", "/bin/bash")

    from subprocess import check_output
    try:
        return check_output(command, **popen_kwargs).strip()
    except:
        if catch:
            return ExcInfo(safe=safe)
        else:
            raise
    finally:
        if "stdout" in popen_kwargs and hasattr(popen_kwargs["stdout"], "close"):
            popen_kwargs["stdout"].close()
        if "stderr" in popen_kwargs and hasattr(popen_kwargs["stderr"], "close"):
            popen_kwargs["stderr"].close()

class Struct:
    """Dynamically creates a simple class that has only attributes and an __init__."""
    __slots__ = [ 'class_' ]

    def __init__(self, name, slots=[], kwslots={}, types={}, doc=None):
        self.class_ = type(name, (object,), {})
        self.class_.__slots__ = slots + kwslots.keys()
        def _init(self, *args, **kwargs):
            def convert(name, value):
                if name in types:
                    t = types[name]
                    if callable(t):
                        value = t(value)
                    elif hasattr(t, "__contains__"):
                        assert value in t, "Invalid value for %s" % name
                    else:
                        raise Exception("Invalid type: %s" % t)
                return value

            for i in range(0, len(args)):
                assert i < len(slots)
                setattr(self, slots[i], convert(slots[i], args[i]))
            for k in kwargs:
                assert k in kwslots
                setattr(self, k, convert(k, kwargs[k]))
            for k in set(kwargs.keys()) ^ set(kwslots.keys()):
                default = kwslots[k]
                if callable(default):
                    default = default()
                setattr(self, k, default)
        self.class_.__init__ = _init
        self.class_.__doc__ = doc
        self.__repr = "Struct('%s', %s, %s)" % (name, slots, kwslots)

    def __call__(self, *args, **kwargs):
        return self.class_(*args, **kwargs)

    def __repr__(self):
        return self.__repr

class ExcInfo(object):
    def __init__(self, exc_info=None, safe=False):
        if exc_info is None:
            exc_info = sys.exc_info()
        self.exc_info = exc_info
        # exception tracebacks are not pickle-able, so if we want this ExcInfo to
        # be safe (i.e. pickleable), we need to delete the traceback
        if safe:
            self.exc_info.trace = None

    def __repr__(self):
        return str(self.exc_info)

    @property
    def exc_class(self):
        return self.exc_info[0]

    @property
    def exc_obj(self):
        return self.exc_info[1]

    @property
    def exc_trace(self):
        return self.exc_info[2]

    def reraise(self):
        raise self.exc_info[0], self.exc_info[1], self.exc_info[2]
