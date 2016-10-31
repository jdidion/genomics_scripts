from unittest import TestCase, skipIf
import os
from xphyle.archive import *
from xphyle.compression import *
from xphyle.paths import get_executable_path
from . import *

def get_format(filename):
    fmt = get_archive_format(filename)
    if fmt:
        return get_archive_format(fmt[0])
    raise ValueError("Invalid archive format: {}".format(filename))

class ArchiveTests(TestCase):
    def test_guess_format(self):
        self.assertTupleEqual(('tar', None), guess_archive_format('tar'))
        self.assertTupleEqual(('tar', None), guess_archive_format('.tar'))
        self.assertTupleEqual(('tar', None), guess_archive_format('foo.tar'))
        self.assertTupleEqual(('tar', 'gz'), guess_archive_format('tar.gz'))
        self.assertTupleEqual(('tar', 'gz'), guess_archive_format('.tar.gz'))
        self.assertTupleEqual(('tar', 'gz'), guess_archive_format('foo.tar.gz'))
    
    def test_invalid_format(self):
        self.assertIsNone(guess_archive_format('foo'))
        self.assertIsNone(guess_archive_format('foo.bar.baz'))
        with self.assertRaises(ValueError):
            get_archive_format('foo')
    
    def write_read_strings(self, filename):
        strings = dict(('text{}'.format(i), random_text()) for i in range(3))
        with open_archive_writer(filename) as w:
            for name, text in strings.items():
                w.write(text, name)
    
    def write_read_files(self, filename):
        descriptors = [FileDescriptor(contents=random_text()) for i in range(3)]
        with make_files(*descriptors) as paths:
            with open_archive_writer(filename) as w:
                for path in paths:
                    w.write_file(path)
