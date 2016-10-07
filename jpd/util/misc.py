"""
Miscelaneous util methods.
"""
from contextlib import contextmanager
from functools import wraps
import logging
import math
import os
import sys

def cumsum(seq):
    cum = 0
    for val in seq:
        cum += val
        yield cum

@contextmanager
def timeit(logfn=lambda x: sys.stdout.write(x), msg="Execution took {0:.2} seconds"):
    import time
    start = time.time()
    try:
        yield
    finally:
        finish = time.time()
        logfn(msg.format(finish - start))

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
