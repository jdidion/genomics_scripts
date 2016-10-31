"""High-level framework for distributed computing."""

import logging
import sys
from .misc import bash

def make_pool(threads=None):
    from multiprocessing import Pool, cpu_count
    from signal import signal, SIGINT, SIG_IGN
    if threads: threads = min(threads, cpu_count())
    # make worker threads ignore keyboard interrupts so they get propagated to the main thread
    return Pool(threads, lambda: signal(SIGINT, SIG_IGN))

## Executing functions in parallel ##

def distribute(func, iterable, mode=None, result_handler=None, **kwargs):
    """Dispatches to a specific distributed processing function based on the ``mode`` parameter.

    Args:
        func (callable): The function to distribute
        iterable: The arguments to ``func``. One element is passed to each execution of ``func``.
        mode (str or callable): typically the name of one of the implementation functions in this
            module (e.g. 'thread'). Can also be a callable.
        result_handler (callable): Does something with the result of executing func.
        kwargs: Any specific parameters to pass along to the implementation method.

    Returns:
        A list of results.

    Raises:
        AttributeError
    """
    if mode == "test":
        func_name = func.func_name
        for i in iterable:
            logging.info("Calling {0} with argument {1}".format(func_name, str(i)))
    elif mode is None or mode == "serial":
        results = [func(i) for i in iterable]
        return [result_handler(r) for r in results]
    elif callable(mode):
        return mode(func, iterable, result_handler, **kwargs)
    else:
        self = sys.modules['util.fork']
        if mode in dir(self):
            dist_handler = getattr(self, mode)
            if callable(dist_handler):
                return dist_handler(func, iterable, result_handler, **kwargs)

    raise AttributeError("Invalid mode %s" % mode)

def thread(func, iterable, result_handler=None, threads=None):
    """Perform a single task over an iterable using a thread pool.

    NOTE: Many exceptions are not pickle-able, which results in uncaught exceptions that cause
    worker threads to hang. Therefore, func should always catch and log all exceptions and only
    raise exceptions that are pickle-able.
    """
    pool = make_pool(threads)
    try:
        results = pool.map(func, iterable)
        if results and result_handler is not None:
            results = [result_handler(r) for r in results]
        return results
    except BaseException as e:
        pool.terminate()
        if not isinstance(e, KeyboardInterrupt):
            raise
    finally:
        pool.close()
        pool.join()

def thread_async(func, iterable, result_handler, threads=None, timeout=sys.maxint, batch=True):
    """Perform a single task over an iterable using a thread pool. The timing of the call to
    `result_handler` depends on `batch`: when `batch` is True, all results will be passed in a
    single call to `result_handler` when all tasks are complete, otherwise `result_handler` is
    called with each result immediately after each task completes.

    NOTE: Setting `batch`=False will likely result in decreased performance.

    NOTE: Many exceptions are not pickle-able, which results in uncaught exceptions that cause
    worker threads to hang. Therefore, func should always catch and log all exceptions and only
    raise exceptions that are pickle-able.
    """
    pool = make_pool(threads)
    try:
        if batch:
            pool.map_async(func, iterable, callback=result_handler).get(timeout)
        else:
            results = tuple(pool.apply_async(func, args, callback=result_handler) for args in iterable)
            for r in results:
                r.get(timeout)
    except BaseException as e:
        pool.terminate()
        if not isinstance(e, KeyboardInterrupt):
            raise
    finally:
        pool.close()
        pool.join()

## Executing shell commands in parallel ##

class Dispatcher(object):
    """ContextManager that dispatches results of parallel operations."""

    def __init__(self, result_handler=None, error_handler=None):
        self.result_handler = result_handler
        self.error_handler = error_handler

    def __call__(self, result):
        if isinstance(result[1], ExcInfo):
            err = result[1]
            if err.exc_class ==  SystemExit or err.exc_class == KeyboardInterrupt:
                err.reraise()
            elif self.error_handler:
                self.error_handler(result)
        elif self.result_handler:
            self.result_handler(result)

    def __enter__(self):
        if hasattr(self.result_handler, '__enter__'):
            self.result_handler.__enter__()
        if hasattr(self.error_handler, '__enter__'):
            self.error_handler.__enter__()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if hasattr(self.result_handler, '__exit__'):
            self.result_handler.__exit__(exc_type, exc_value, traceback)
        if hasattr(self.error_handler, '__exit__'):
            self.error_handler.__exit__(exc_type, exc_value, traceback)

# error handlers

def reraise_error(result):
    """Error handler that throws exceptions."""
    result[1].reraise()

def log_error(result):
    """Error handler that logs errors."""
    logging.error('Error for args %s' % str(result[0]), exc_info=result[1].exc_info)

class ExcInfo(object):
    """Wrapper class for exception info tuple."""
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

# executors for command-line processes

def log_commands(cmd_iter, output_handler):
    for args in cmd_iter:
        command = args[1]
        if not isinstance(command, str):
            command = ' '.join(command)
        logging.info(command)

def exec_serially(cmd_iter, dispatcher):
    for args in cmd_iter:
        dispatcher(_exec(args))

class exec_threaded(object):
    def __init__(self, kwargs):
        self.threads = int(kwargs.get("threads", 1))
        self.timeout = kwargs.get("timeout", None)
        if isinstance(self.timeout, str):
            try:
                self.timeout = int(self.timeout)
            except:
                self.timeout = int(self.timeout, 16)
        if self.timeout and self.timeout < 0:
            self.timeout = None

    def __call__(self, cmd_iter, dispatcher):
        thread_async(_exec, cmd_iter, dispatcher, self.threads, self.timeout)

class exec_lsf(object):
    def __init__(self, kwargs):
        self.queue = kwargs.get("queue", "week")
        self.output_file = kwargs.get("output_file", None)
        self.send_mail = kwargs.get("output_file", self.output_file is None)
        self.resources = kwargs.get("resources", None)

    def __call__(self, cmd_iter, dispatcher):
        for args in cmd_iter:
            bsub = "bsub -q {0} ".format(self.queue)
            if self.send_mail:
                bsub += "-N "
            if self.output_file:
                out_file = self.output_file.format(args)
                bsub += "-o {0} ".format(out_file)
            if self.resources:
                for r in self.resources:
                    bsub += "-R {0} ".format(r)

            argvars, command, popen_kwargs = args
            if not isinstance(command, str):
                command = ' '.join(command)

            bash(bsub + command, **popen_kwargs)

def _exec(args):
    argvars, command, popen_kwargs = args
    try:
        result = bash(command, **popen_kwargs)
        # If this method is called within a worker thread, then both argvars and result
        # must be pickle-able. If not, an uncaught exception will result in the worker
        # thread (and thus the main process) hanging. There are two work-arounds:
        # 1. popen_kwargs["catch"] = False: causes bash() to raise exceptions rather than return
        # a (potentially unpickle-able) ExcInfo
        # 2. popen_kwargs["safe"] = True: causes bash to null out ExcInfo traceback to make it
        # pickleable.
        return (argvars, result)
    except BaseException as e:
        raise Exception(e.message)

def get_executor(name, kwargs):
    if name == "test":
        return log_commands
    elif name == "serial" or (name == "thread" and kwargs.get("threads") == 1):
        return exec_serially
    elif name == "thread":
        return exec_threaded(kwargs)
    elif name == "lsf":
        return exec_lsf(kwargs)

# command iterators

def string_formatter(command):
    return lambda av: command.format(**av)

def command_iter(formatter, argvars, outfile_pattern=None, **popen_kwargs):
    for av in argvars:
        command = formatter(av)
        if outfile_pattern:
            from subprocess import STDOUT
            popen_kwargs['stdout'] = open(outfile_pattern.format(**av), 'w')
            popen_kwargs['stderr'] = STDOUT
        yield (av, command, popen_kwargs)

def simple_command_iter(prog, args, argvars, outfile_pattern=None, **popen_kwargs):
    def formatter(av):
        return [prog] + [arg.format(**av) for arg in args]
    return command_iter(formatter, argvars, outfile_pattern, **popen_kwargs)

def command_iter_with_pipe(prog, args, argvars, pipe_file_pattern, pipe_command_pattern, **popen_kwargs):
    """ Special command iterator that injects a pipe redirection command before every process."""
    for av in argvars:
        # first yield the pipe commands
        pipe_file = pipe_file_pattern.format(**av)
        pipe_command = pipe_command_pattern.format(**av)
        pipe_script = "; ".join(("mkfifo -m 777 {0}".format(pipe_file), pipe_command, "rm -f {0}".format(pipe_file)))
        yield (av, pipe_script, popen_kwargs)

        # now yield the main command
        wait_command = "while [ ! -p {0} ] ; do sleep 1; done;".format(pipe_file)
        command = " ".join([wait_command, prog] + [arg.format(**av) for arg in args])
        yield (av, command, popen_kwargs)

# simple interface

def exec_shell(cmd_iter, executor=log_commands, result_handler=None, error_handler=None):
    with Dispatcher(result_handler, error_handler) as d:
        try:
            executor(cmd_iter, d)
        except (SystemExit, KeyboardInterrupt):
            logging.warn("Exiting program before completion")
            sys.exit()
