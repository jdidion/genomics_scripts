"""
Module for working with data tables (square two-dimensional arrays).
"""

from __future__ import print_function
import pandas as p

# I/O

def from_text(path, rownames=True, header=True, delimiter="\t"):
    df = p.io.parsers.read_table(path, sep=delimiter, 
        header=0 if header is True else None, 
        index_col=0 if rownames is True else None)
    return Table(df)

def to_csv(table, csvfile, rownames=True, header=True, **kwargs):
    with open(csvfile, 'w') as f:
        w = csv.writer(f, **kwargs)
        for r in table.rows():
            w.writerow(r)

# Display

def head(table, rows=10, cols=10, dest=sys.stdout, sep="\t"):
    rowrng = xrange(0, rows)
    colrng = xrange(0, cols)
    
    if table.has_rownames():
        print(*table.rownames(rowrng), sep=sep, file=dest)
    
    cellitr = table.cells(rowrng, colrng)
    colnames = table.colnames(colrng) if table.has_colnames() else None
    for r in rowrng:
        if colnames is not None:
            print(colnames[r], sep, file=dest)
            
        for c in colrng:
            print(cellitr.next(), file=dest)
            if c == (cols - 1):
                print("\n", file=dest)
            else:
                print(sep, file=dest)
    
# Transformation

def transpose(table):
    """Returns a new table with rows and columns transposed."""
    
    if table.byrow():
        return ColTable(table.colnames(), table.rownames(), table.rows())
    else:
        return RowTable(table.colnames(), table.rownames(), table.cols())

# Model

class NameMap(object):
    def __init__(self, names=None):
        self._idx_map = bidict()
        self._str_map = bidict()
        if names is not None:
            self.add_all(names)
    
    def __len__(self):
        return len(self._str_map)
    
    def add(self, name):
        self._str_map[len(self)] = name
    
    def add_all(self, names):
        for n in names: self.add(n)
           
    def swap(self, a, b):
        """Swapping two rows/columns just adds/updates entry in _idx_map."""
        
        a = self.resolve(a)
        a_new = a
        if a in ~self._idx_map:
            a_new = (~self._idx_map)[a]
            
        b = self.resolve(b)
        b_new = b
        if b in ~self._idx_map:
            b_new = (~self._idx_map)[b]
            
        self._idx_map[a_new] = b
        self._idx_map[b_new] = a
        
    def names(self, idxs=None):
        if idxs is None:
            idxs = sorted(self._str_map.keys())
        if len(self._idx_map) > 0:
            idxs = [self._idx_map.get(i, i) for i in idxs]
        return [self._str_map[i] for i in idxs]
    
    def resolve(self, name):
        if name in ~self._str_map:
            return (~self._str_map)[name]
        elif name in self._idx_map:
            return self._idx_map[name]
        else:
            return name
    
    def resolve_all(self, names):
        if names is not None:
            if isinstance(names, list):
                names = [self.resolve(n) for n in names]
            else:
                names = [self.resolve(names)]
        return names

class BaseTable(object):
    """A data table (square two-dimensional array). Data can be stored by row or by column."""
    
    def __init__(self, rownames=True, colnames=True, data=None):
        self._data = []
        self._nrow = None
        self._ncol = None
        
        if util.misc.iterable(rownames):
            self._rows = NameMap(rownames)
        else:
            self._rows = NameMap()
        self._name_rows = (rownames is True)
        
        if util.misc.iterable(colnames):
            self._cols = NameMap(colnames)
        else:
            self._cols = NameMap()
        self._name_cols = (colnames is True)

        if data is not None:
            for d in data:
                self.append(d)
                
    @property
    def nrow(self):
        return self._nrow

    @property
    def ncol(self):
        return self._ncol

    def has_rownames(self):
        return len(self._rows) > 0

    def rownames(self, rows=None):
        return self._rows.names(rows)

    def rows(self, rows=None):
        if rows is None:
            rows = xrange(0, self.nrow())
        return self._iter_rows(self._rows.resolve_all(rows))

    def has_colnames(self):
        return len (self._cols) > 0

    def colnames(self, cols=None):
        return self._cols.names(cols)

    def cols(self, cols=None):
        if cols is None:
            cols = xrange(0, self.ncol())
        return self._iter_cols(self._cols.resolve_all(cols))
    
    def cells(self, rows=None, cols=None):
        return self._iter_cells(self._rows.resolve_all(rows), self._cols.resolve_all(cols))
    
class RowTable(BaseTable):
    @property
    def byrow(self):
        return True
        
    def append(self, row):
        if self._name_cols and len(self._cols) == 0
            self._cols.add_all(values)
            values = None
        elif self._name_rows:
            self._rows.add(values.pop(0))
            
        if self._ncol is None:
            self._nrow = 0
            self._ncol = len(row)
        else:
            assert len(row) == self._ncol, "Row must be of length %i" % self._ncol
        
        if values is not None:          
            self._data.append(row)
            self._nrow += 1
    
    def insert(self, row, at=0):
        self.append(row)
        self.swap(self._nrow - 1, at)
    
    def swap(self, a, b):
        self._rows.swap(a, b)
    
    def _iter_rows(self, rows):
        for r in rows:
            yield self._data[r]
    
    def _iter_cols(self, cols):
        for c in cols:
            yield [r[c] for r in self._data]
    
    def _iter_cells(self, rows, cols):
        for r in rows:
            for c in cols:
                yield self._data[r][c]

class ColTable(BaseTable):
    @property
    def byrow(self):
        return False
        
    def append(self, col):
        if self._name_rows and len(self._rows) == 0:
            self._rows.add_all(values)
            values = None
        elif self._name_cols:
            self._cols.add(values.pop(0))
        
        if self._nrow is None:
            self._nrow = len(col)
            self._ncol = 0
        else:
            assert len(col) == self._nrow, "Column must be of length %i" % self._nrow
        
        if values is not None:
            self._data.append(col)
            self._ncol += 1
    
    def insert(self, col, at=0):
        self.append(col)
        self.swap(self._ncol - 1, at)

    def swap(self, a, b):
        self._cols.swap(a, b)
                
    def _iter_rows(self, rows):
        for r in rows:
            yield [c[r] for c in self._data]
    
    def _iter_cols(self, cols):
        for c in cols:
            yield self._data[c]
    
    def _iter_cells(self, rows, cols):
        for r in rows:
            for c in cols:
                yield self._data[c][r]
