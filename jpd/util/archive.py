# -*- coding: utf-8 -*-
"""Interface to archive formats.
"""
from contextlib import contextmanager
import os
import re
from xphyle.compression import *
from xphyle.paths import splitext

archive_formats = {}
"""Dict of registered archive formats"""

def register_archive_format(format_class : 'class'):
    """Register a new compression format.
    
    Args:
        ``format_class`` -- a subclass of ArchiveFormat
    """
    aliases = set(format_class.exts)
    aliases.add(format_class.lib_name)
    for alias in aliases:
        # TODO: warn about overriding existing format?
        archive_formats[alias] = format_class

def get_archive_format(name : 'str') -> 'ArchiveWriter':
    """Get the ArchiveWriter class for the given name.
    """
    if name in archive_formats:
        return archive_formats[name]
    raise ValueError("Unsupported archive format: {}".format(name))

def guess_archive_format(name : 'str') -> '(str,str)':
    """Guess the archive format from the file extension.
    
    Returns:
        (archive_format, compression_format), or ``None`` if the format
        can't be determined.
    """
    if name in archive_formats:
        return (name, None)
    file_parts = splitext(name, False)
    if len(file_parts) < 2:
        return None
    if len(file_parts) == 2 and file_parts[1] in archive_formats:
        return (file_parts[1], None)
    # it might be a compressed archive
    if (guess_compression_format(file_parts[-1]) and
            file_parts[-2] in archive_formats):
        return (file_parts[-2], file_parts[-1])
    return None

def _parse_pattern(pattern):
    if callable(pattern):
        return pattern
    if any(isinstance(pattern, type_) for type_ in (list, tuple, set)):
        return _SetMatcher(set(pattern))
    if isinstance(pattern, str):
        pattern = re.compile(pattern)
    return _PatternMatcher(pattern)

class _PatternMatcher(object):
    def __init__(self, p):
        self.p = p
    def __call__(self, entry_name):
        return self.p.match(entry_name)

class _SetMatcher(object):
    def __init__(self, s):
        self.s = s
    def __call__(self, entry_name):
        return entry_name in self.s

def open_archive(path : 'str', atype : 'str' = None, ctype : 'str' = None,
                 **kwargs) -> 'ArchiveWriter':
    """Open a writer for an archive file.
    
    Args:
        path: The path to the archive file
        atype: The archive type, or None if should be guessed from the file name
        ctype: The compression type, or None if it should be guessed from the
            file name
        kwargs: Addtional keyword arguments to pass to the ArchiveWriter
            constructor
    
    Returns:
        An open ArchiveWriter
    """
    if atype is None:
        guess = guess_archive_format(path)
        if guess:
            atype = guess[0]
            if ctype is None:
                ctype = guess[1]
    
    if atype is None:
        raise Exception("Invalid archive type {}".format(atype))
    
    archive_format = get_archive_format(atype)
    return archive_format(path, ctype, **kwargs)

class ArchiveReader(object):
    """Mixin for ``ArchiveFormat``s that are readable.
    """
    def get_entry(self, archive_name : 'str') -> 'bytes':
        """Reads an archive entry:
        
        Args:
            archive_name: The name of the entry
        
        Returns:
            The entry contents as bytes
        """
        raise NotImplemented()
    
    def get_string_entry(self, archive_name : 'str',
                         encoding : 'str' = 'utf-8') -> 'str':
        """Reads an archive entry as a string:
        
        Args:
            archive_name: The name of the entry
            encoding: The byte encoding to use
        
        Returns:
            The entry contents as a string
        """
        return self.get_entry(archive_name).decode(encoding)
    
    def extract_file(self, archive_name : 'str', dest_dir : 'str' = ".",
                     dest_file : 'str' = None) -> 'string':
        """Extract the archive entry to a file.
        
        Args:
            archive_name: The name of the entry
            dest_dir: The root directory where the file will be extracted
            dest_file: The relative or absolute path to the destination file. If
                ``None``, ``archive_name`` is used. If an absolute path,
                ``dest_dir`` is ignored.
        
        Returns:
            The fully resolved path to the file
        
        Raises:
            IOError if the destination file is not writeable
        """
        if dest_file is None:
            dest_file = archive_name
        if not os.path.isabs(dest_file):
            dest_file = os.path.join(dest_dir, dest_file)
        dest_file = check_writeable_file(dest_file)
        return self._extract(archive_name, dest_file)
    
    def extract_all_files(self, dest_dir : 'str', pattern = None) -> 'int':
        """Extract all files in the archive (optionally matching a specific
        pattern).
        
        Args:
            dest_dir: Root directory where files will be extracted
            pattern: Regular expression (as a string or compiled regular
                expression), tuple/list of entries to extract, or a callable
                that returns a True value for entries that should be extracted.
        
        Returns:
            The number of files that were extracted
        """
        if pattern:
            pattern = _parse_pattern(pattern)
        return self._extract_all(dest_dir, pattern)

class ArchiveWriter(object):
    """Mixin for ``ArchiveFormat``s that are writeable.
    """
    def add_entry(self, bytes : 'bytes', archive_name : 'str'):
        """Write the bytes into the archive.
        
        Args:
            bytes: The bytes to write
            archive_name: The name to assign the string in the archive
        """
        raise NotImplemented()
    
    def add_string_entry(self, string : 'str', archive_name : 'str',
                         encoding : 'str' = 'utf-8'):
        """Write the contents of a string into the archive.
        
        Args:
            string: The string to write
            archive_name: The name to assign the string in the archive
        """
        self.write(string.encode(encoding), archive_name)
    
    def add_file(self, path : 'str', archive_name : 'str' = None):
        """Copy a file into the archive.
        
        Args:
            path: Path to the file
            archive_name: Name to give the file within the archive. Defaults
                to ``path`` if None.
        """
        if archive_name is None:
            archive_name = path
        path = check_readable_file(path)
        self._write_file(path, archive_name)
    
    def add_all_files(self, src_dir : 'str', depth=None, pattern=None) -> 'int':
        """Copy all files in src_dir (optionally matching a pattern) into
        the archive.
        
        Args:
            src_dir: The directory containing the files to add
            depth: The maximum number of levels to recurse; no limit if None
            pattern: Regular expression (as a string or compiled regular
                expression), tuple/list of entries to extract, or a callable
                that returns a True value for files that should be added.
        """
        if pattern:
            pattern = _parse_pattern(pattern)
        
    
class ArchiveFormat(FileFormat, ArchiveReader, ArchiveWriter):
    """Base class for writers that store arbitrary data and files within
    archives (e.g. tar, zip).
    
    Args:
        path: Path to the archive file
        mode: File open mode
        compression: If None or False, do not compress the archive. If True,
            compress the archive using the default format. Otherwise specifies
            the name of compression format to use.
        create_opened: Whether the archive should be opened immediately. If
            False, the caller must explicitly call ``open``.
        kwargs: Additional arguments to pass to the open method
    """
    def __init__(self, path : 'str', mode : 'str' = 'mode',
                 compression : 'bool' = False, create_opened : 'bool' = True,
                 **kwargs):
        self.path = path
        self.mode = mode
        self.compression = compression
        self.archive = None
        if create_opened:
            self.open(**kwargs)
    
    def close(self):
        self.archive.close()
    
    def __enter__(self):
        if self.archive is None:
            self.archive = self.open()
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

class CompressableArchiveFormat(ArchiveFormat):
    """Base class for an ArchiveWriter that is compressable by an external
    executable (e.g. tar.gz) rather than having its own internal compression
    scheme (e.g. zip).
    """
    def close(self):
        super(CompressedArchiveWriter, self).close()
        if self.compression:
            ctype = self.compression
            if ctype is True:
                ctype = self.default_compression
            compress_file(self.path, ctype=ctype)

class ZipWriter(ArchiveFormat):
    """Simple, write-only interface to a zip file.
    """
    exts = ('zip',)
    lib_name = 'zipfile'
    
    def open(self, **kwargs):
        self.archive = self.lib.Zipfile(
            self.path, 'w',
            self.lib.ZIP_DEFLATED if self.compression else self.lib.ZIP_STORED,
            **kwargs)
    
    def _write_file(self, path, archive_name):
        self.archive.write(path, archive_name)
    
    def write(self, string, archive_name):
        self.archive.writestr(archive_name, string)
register_archive_format(ZipWriter)

class TarWriter(CompressableArchiveFormat):
    """Simple, write-only interface to a tar file.
    """
    exts = ('tar',)
    lib_name = 'tarfile'
    default_compression = 'gzip'
    
    def open(self, **kwargs):
        self.archive = self.lib.TarFile(self.path, 'w', **kwargs)
    
    def _write_file(self, path, archive_name):
        self.archive.add(path, archive_name)
    
    def write(self, string, archive_name):
        ti = tarfile.TarInfo(arcname)
        ti.frombuf(string)
        self.archive.addfile(ti)
register_archive_format(TarWriter)
