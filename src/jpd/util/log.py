from functools import partial
import logging

class MyLogger(logging.Logger):
    def __init__(self, name, level=logging.NOTSET):
        logging.Logger.__init__(self, name, level)
        def method_factory(v):
            return lambda: v > 0 and self.isEnabledFor(v)
        for k,v in logging._levelNames.items():
            if isinstance(k, str): 
                setattr(self, "is_%s" % k.lower(), method_factory(v))

def extend_logging():
    """
    Add a TRACE level (before DEBUG) and methods to the logging class that check for whether a 
    certain logging state is enabled.
    """
    logging.setLoggerClass(MyLogger)
    logging.addLevelName(5, 'TRACE')
    setattr(logging, 'trace', lambda msg, *args, **kwargs: logging.log(5, msg, *args, **kwargs))
    
    def method_factory(v):
        return lambda: v > 0 and logging.root.isEnabledFor(v)
    for k,v in logging._levelNames.items():
        if isinstance(k, str): 
            setattr(logging, "is_%s" % k.lower(), method_factory(v))

def is_extended():
    return 'TRACE' in get_level_names()

def get_level_names():
    return logging._levelNames

def conf(filename=None, level=None, extended=True):
    if extended and not is_extended():
        extend_logging()
    if isinstance(level, str):
        level = get_level_names()[level]
    logging.basicConfig(filename=filename, level=level)
    
def handle_error(msg, format_args=None):
    """
    An error-handler decorator that logs any exception thrown by a function call. The msg parameter 
    specifies the error message to be displayed. The optional format_args parameter allows you to 
    specify which args (as a tuple of integers) passed to the decorated command are to be used in 
    formatting the message. 
    
    Example:
    @handle_error("Your %s is dumb", (0,))
    def foo(relation):
       raise Exception("So's your mom")
    foo('cousin')
    """
    return partial(ErrorHandler, msg=msg, format_args=format_args)

class ErrorHandler:
    def __init__(self, command, msg, format_args=None):
        self.command = command
        self.msg = msg
        self.format_args = format_args
    
    def __call__(self, *args, **kwargs):
        try:
            return self.command(*args, **kwargs)
        except:
            format_args = []
            if self.format_args is not None:
                for a in self.format_args:
                    if instanceof(a, int):
                        format_args.append(args[a])
                    else:
                        format_args.append(kwargs[a])
            # Don't return the top frame of the traceback
            e = sys.exc_info()
            logging.error(self.msg, *format_args, exc_info=(e[0], e[1], e[2].tb_next))