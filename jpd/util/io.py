"""
Utility i/o methods, mostly for reading and writing files.
"""
from __future__ import absolute_import
from collections import OrderedDict
from contextlib import contextmanager
from csv import reader, writer
import errno
import fileinput
from functools import partial
import gzip
from itertools import cycle, izip
import logging
import os
import os.path
import shutil
import sys
from .collections import PropDict
from .misc import iterable

# Paths

def abspath(path):
    return os.path.abspath(os.path.expanduser(path))

def filename(path):
    return os.path.splitext(os.path.basename(path))[0]
    
def resolve_path(path, parent=None):
    """Resolves the absolute path of the specified file.
    
    Args:
        path (str): Path to resolve.
        parent (str): The directory containing ``path`` if ``path`` is relative.
    
    Returns:
        The absolute path.
    
    Raises:
        IOError: if the path does not exist.
    """
    apath = abspath(path)
    if not os.path.exists(apath) and parent is not None:
        apath = abspath(os.path.join(parent, path))
    if not os.path.exists(apath):
        raise IOError(errno.ENOENT, "%s does not exist" % apath, apath)
    return apath

def check_path(path, ptype=None, access=None):
    """Checks that a path exists, is of the specified type, and allows the specified access.
    
    Args:
        ptype: 'f' for file or 'd' for directory.
        access (int): One of the access values from :module:`os`
    
    Raises:
        IOError if the path does not exist, is not of the specified type, or doesn't allow the
        specified access.
    """
    if ptype == 'f' and not os.path.isfile(path):
        raise IOError(errno.EISDIR, "%s is not a file" % path, path)
    elif ptype == 'd' and not os.path.isdir(path): 
        raise IOError(errno.ENOTDIR, "%s is not a directory" % path, path)
    elif not os.path.exists(path):
        raise IOError(errno.ENOENT, "%s does not exist" % path, path)
    if access is not None and not os.access(path, access): 
        raise IOError(errno.EACCES, "%s is not accessable" % path, path)
    return path

def mkdir(path, subdir=None, create=True, overwrite=False, assert_exists=False):
    if subdir:
        path = os.path.join(path, subdir)
    
    if assert_exists:
        assert os.path.isdir(path), "Directory does not exist: {0}".format(path)
    
    if create:
        if os.path.isdir(path) and overwrite:
            shutil.rmtree(path)
        if not os.path.isdir(path):
            os.makedirs(path)
    
    return path

def find(root, pattern, types="f", recursive=True):
    """
    Find all files under 'root' that match 'pattern' (a ``RegexObject``). Depending on the
    value of 'types', return files ("f"), directories ("d") or both ("fd").
    """
    found = []
    for root, dirs, files in os.walk(root):
        if types != "f":
            for d in filter(lambda x: pattern.match(x), dirs):
                found.append(os.path.join(root, d))
        if types != "d":
            for f in filter(lambda x: pattern.match(x), files):
                found.append(os.path.join(root, f))
    return found

# Reading/writing files

def open_stream(path, mode):
    """Opens a file or returns stdin/stdout if path == '-'"""
    if path == "-":
        return sys.stdin if mode[0] == "r" else sys.stdout
    elif os.path.splitext(path)[1] == ".gz":
        if "r" in mode:
            mode = "rb"
        else:
            mode = "wb"
        return gzip.open(path, mode)
    else:
        return open(path, mode)

def linecount(filename):
    """
    Fastest way to count the lines in a file. !Only works with files that use \\n as the
    line delimiter (adding univeral newline support slows this function by 3x)!
    """
    lines = 0
    buf_size = 1024 * 1024
    with open_stream(filename) as f:
        read_f = f.read # loop optimization
        buf = read_f(buf_size)
        while buf:
            lines += buf.count('\n')
            buf = read_f(buf_size)
    return lines

def safe_read_file(path):
    """Returns the contents of a file, or None if the file does not exist."""
    
    if path is None: return None
    if not os.path.isfile(path) or not os.access(path, os.R_OK): 
        logging.warn("%s is not a readable file", path)
        return None
    with open_stream(path, 'rU') as f: return f.read()
    
def safe_read_file_array(path, convert=lambda x:x):
    """Read a file into an array. Returns None if the file doesn't exist or is not readable."""
    
    if path is None: return None
    if not os.path.isfile(path) or not os.access(path, os.R_OK): 
        logging.debug("%s is not a readable file", path)
        return None
    return [ convert(s.rstrip()) for s in fileinput.input(path) ]

def read_chunked(path, chunksize=1024):
    """Generator that reads a binary file in chunks of ``chunksize``."""
    with open(path, 'rb') as infile:
        while True:
            data = infile.read(chunksize)
            if data:
                yield data
            else:
                break

def read_pickled(binfile, compressed=True):
    import cPickle
    if compressed:
        import zlib
    with open(binfile, 'rb') as f:
        data = f.read()
        if compressed:
            data = zlib.decompress(data)
        return cPickle.loads(data)
        
def write_file(path, string):
    """Write the contents of ``string`` to the specified file."""
    if os.path.isfile(path) and not os.access(path, os.W_OK):
        raise Exception("%s is not a writable file" % path)
    with open_stream(path, 'w') as f: f.write(string)

def write_array(path, a):
    write_file(path, "\n".join(a))
    
def write_dict(path, d):
    """Write a dict to a file as name=value lines."""
    write_file(path, "\n".join(["{0}={1}".format(k,v) for k,v in d.iteritems()]))

class FileWrapper(object): 
    __slots__ = ['_file'] 
    
    def __init__(self, fh): 
        object.__setattr__(self, '_file', fh) 
    
    def __getattr__(self, name): 
        return getattr(self._file, name) 
    
    def __setattr__(self, name, value): 
        setattr(self._file, name, value)

def move_on_close(fh, dest):
    class FileMover(FileWrapper):
        def close(self):
            import shutil
            self._file.close()
            shutil.move(self._file.name, dest)
    return FileMover(fh)

def del_on_close(fh):
    class FileDeleter(FileWrapper):
        def close(self):
            self._file.close()
            os.remove(self._file.name)
    return FileDeleter(fh)

def csv_to_table(csv_file, delim=",", convert=lambda x:x):
    """Parses a CSV file and returns a tuple of rows.
    
    Returns:
        A tuple of rows from the CSV file, i.e.::
        
            ((row_1_col_1, row_1_col_2, ...), (row_2_col_1, row_2_col_2, ...), ...)
    """
    if csv_file is None:
        return None
    with open_stream(csv_file, 'rU') as f:
        return tuple(convert(row) for row in reader(f, delimiter=delim))

def csv_to_dict(csv_file, delim=",", key=0, convert=lambda x:x, skip_blank=True):
    """Parses a CSV file and returns a dict of rows.
    """
    
    if csv_file is None:
        return None
    if isinstance(key, int):
        keyfn = lambda row: row[key]
    else:
        keyfn = key
    d = {}
    with open_stream(csv_file, 'rU') as f:
        for row in reader(f, delimiter=delim):
            if len(row) == 0 and skip_blank: continue
            v = convert(row)
            k = keyfn(v)
            d[k] = v
    return d

def table_to_csv(rows, csv_file, delim=","):
    with open_stream(csv_file, "w") as o:
        w = writer(o, delimiter=delim)
        w.writerows(rows)

def parse_properties(source=None, fn=None, delim="=", ordered=False):
    """
    Read properties from a file (one per line) and/or a list. If ``fn`` is specified, apply 
    the function to each value. Result is an ordered dict if 'ordered' is True, otherwise a 
    PropDict, which allows property-style indexing.
    """
    if isinstance(source, str):
        source = safe_read_file_array(source)
    props = OrderedDict() if ordered else PropDict()
    for line in source:
        line = line.strip()
        if len(line) == 0 or line[0] == "#":
            continue
        key, value = line.split(delim)
        if fn: value = fn(value)
        props[key] = value
    return props

@contextmanager
def open_mgr(f, mode='rU'):
    """Context manager that frees you from checking if an argument is a path or a file object."""
    if isinstance(f, str):
        with open_stream(f, mode) as fp:
            yield fp
    else:
        yield f

class FileCloser(object):
    def __init__(self):
        self.files = {}
    
    def __getitem__(self, key):
        return self.files[key]
    
    def __contains__(self, key):
        return key in self.files
        
    def add(self, fh, key=None, mode="rU"):
        if isinstance(fh, str):
            fh = open_stream(fh, mode)
        if key is None:
            key = fh.name
        self.files[key] = fh
        return fh

    def close(self):
        for fh in self.files.values():
            fh.close()        

def decoding_reader(f, converters=str, *args, **kwargs):
    if not iterable(converters):
        assert callable(converters)
        converters = cycle(converters)
    with open_mgr(f) as fp:
        for row in reader(fp, *args, **kwargs):
            yield [fn(x) for fn,x in izip(converters, row)]

# Compression

def compress_file(path, compressed_path=None, ctype="gz"):
    """Compress a file.
    
    Args:
        path (str): The file to compress.
        compressed_path (str): The compressed file. If None, the file is compressed in place.
        type: Either 'gz' or 'bz2'.
    """
    inplace = compressed_path is None
    if inplace:
        import tempfile
        tmp, compressed_path = tempfile.mkstemp()
        tmp.close()
    
    if ctype == 'gz':
        import gzip
        cfile = gzip.GzipFile(compressed_path, 'wb')
    elif ctype == 'bz2':
        import bz2
        cfile = bz2.BZ2File(compressed_path, 'w')
    else:
        raise Exception("Invalid compression type %s" % ctype)
        
    # Perform sequential compression as the source file might be quite large
    for bytes in read_chunked(path): cfile.write(bytes)
    
    cfile.close()
    
    if inplace: 
        import shutil
        shutil.move(compressed_path, path)
        
class ArchiveWriter(object):
    def __init__(self, init):
        self.archive = init
    def open(self):
        assert isinstance(self.archive, partial)
        self.archive = self.archive()
    def close(self):
        self.archive.close()
    def __enter__(self):
        self.open()
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        
class ZipWriter(ArchiveWriter):
    """Simple, write-only interface to a zip file."""
    def __init__(self, path, compression=True, **kwargs):
        import zipfile
        compression = zipfile.ZIP_DEFLATED if compression else zipfile.ZIP_STORED
        ArchiveWriter.__init__(self, partial(zipfile.Zipfile, path, 'w', compression, **kwargs))
    def writefile(self, path, arcname):
        self.archive.write(path, arcname)
    def writestr(self, string, arcname):
        self.archive.writestr(arcname, string)

class TarWriter(ArchiveWriter):
    """Simple, write-only interface to a tar file."""
    def __init__(self, path, compression='gz', **kwargs):
        import tarfile
        ArchiveWriter.__init__(self, partial(tarfile.TarFile, path, 'w', **kwargs))
        self.compression = 'gz' if compression is True else compression
    def writefile(self, path, arcname):
        self.archive.add(path, arcname)
    def writestr(self, string, arcname):
        ti = tarfile.TarInfo(arcname)
        ti.frombuf(string)
        self.archive.addfile(ti)
    def close(self):
        ArchiveWriter.close(self)
        if self.compression:
            compress_file(self.archive.name, ctype=self.compression)
        
def write_compressed(path, contents, container_type="zip", compression=True, **kwargs):
    """Write entries to a compressed container file.
    
    Args:
        path (str): Path of the zip file to write.
        contents (iterable): An iterable of (name,content) tuples. A content item can be a path
            to a system file that should be added to the container.
        container_type: Either 'zip' or 'tgz'.
        compression: Either a bool or, if ``container_type`` is tar, 'gz' (the default) or 'bz2'.
        kwargs: Additional args to pass to the container constructor.
    """
    if container_type == 'zip':
        c = ZipWriter(path, compression, **kwargs)
    elif container_type == 'tar':
        c = TarWriter(path, compression, **kwargs)
    else:
        raise Exception("Invalid container type %s" % container_type)
    
    with c.open():
        for name,content in contents.items():
            if os.path.isfile(content):
                c.writefile(content, name)
            else:
                c.writestr(content, name)

# File splitting

class LineHandler(object):
    """Interface that must be implemented by objects supplied as handlers to the split() function."""
    
    def match(self, line):
        """
        If the given line should be handled by this LineHandler, returns an array of values
        that can be substituted into an output file template (or the empty array). Otherwise
        returns None.
        """
        return []
    
    def transform(self, line):
        """
        Applies any necessary transformation to line and returns the result. Any non-string
        return value will be interpreted as the empty string.
        """
        return line
    
    def output(self, headers, line, match_values):
        """Send line to a desired output. Default action is to do nothing."""
        pass

class REMatchMixin(object):
    """Mixin that implements the match function using regular expressions."""
    def __init__(self, pattern):
        import re
        self.pattern = re.compile(pattern)
        
    def match(self, line):
        m = self.pattern(line)
        return m.groups if m else None
        
class FileOutputMixin(object):
    """Mixin that implements the output function to write output to a file."""
    def __init__(self, path, append=False):
        self.template = path
        self.append = append
        self.fileobjs = {}
        
        def close_on_exit(self):
            for f in self.fileobjs.values():
                try:
                    f.close()
                except:
                    pass
        
        import atexit
        atexit.register(close_on_exit, self)
        
    def output(self, headers, line, match_values):
        """Write line to a file."""
        path = self.template
        if match_values:
            path = path.format(*match_values)
        if path in self.fileobjs:
            fileobj = self.fileobjs[path]
        else:
            mode = 'a' if self.append else 'w'
            fileobj = open_stream(path, mode)
            self.fileobjs[path] = fileobj
            if headers:
                for h in headers:
                    fileobj.write(h)
                    fileobj.write("\n")
        fileobj.write(line)
        fileobj.write("\n")

class DefaultFileHandler(LineHandler, REMatchMixin, FileOutputMixin):
    def __init__(self, pattern, path, append=False, autoclose=True):
        REMatchMixin.__init__(self, pattern)
        FileOutputMixin.__init__(self, path, append, autoclose)

def split(infile, handlers, header_lines=0, unmatched_file=None):
    import fileinput
    
    headers = None
    if isinstance(header_lines, int):
        if header_lines > 0:
            headers = []
        else:
            header_lines = abs(header_lines)
    else:
        if iterable(header_lines):
            headers = header_lines
        header_lines = 0
        
    for line in fileinput.input(infile):
        line = line.rstrip()
        
        if header_lines > 0:
            headers.append(line)
            header_lines -= 1
        else:
            for h in handlers:
                match_values = h.match(line)
                if match_values:
                    line = h.transform(line)
                    if line:
                        h.output(headers, line, match_values)
                    break
            else:
                if unmatched_file: unmatched_file.write(line)
