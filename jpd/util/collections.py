"""
Utility methods related to lists and dicts, and custom collection classes.
"""
from __future__ import absolute_import
from collections import OrderedDict, Sequence, defaultdict
import sys
from UserDict import UserDict, DictMixin
from .misc import accessor, buf, measurable, overrides

def all_equal(l):
    """Returns True if all items in a list are the same."""
    assert len(l) > 0
    return len(set(l)) == 1

def none_equal(l):
    """Returns True if all items in ``l`` are unique."""
    return len(l) == len(set(l))

def indices(l, value):
    """Returns a list of indicies of elements of ``l`` that equal ``value``."""
    return [i for i in range(0, len(l)) if l[i] == value]

def elements(l, indices):
    """Return a new list consisiting of elements of ``l`` at the specified ``indices``."""
    return [l[i] for i in indices]

def first(fn, seq):
    """Find the first element for which ``fn`` returns True."""
    for elt in seq:
        if fn(elt): return elt
            
def first_index(fn, seq):
    """Finds the index of the first element for which ``fn`` returns True (or -1 if none)."""
    for i, elt in enumerate(seq):
        if fn(elt): return i
    return -1

def which(val, seq):
    """Finds the index of the first element in seq that is equal to val."""
    return first_index(lambda x: x==val, seq)

def and_or(l1, l2):
    """Concatenate two lists, where either or both might be ``None``."""
    return (l1 or []) + (l2 or [])

def sum_arrays(a1, a2):
    """Sum the corresponding elements of two arrays."""
    return [sum(x) for x in zip(a1, a2)]

def add_array(a1, a2):
    """For each element in a1, add the corresponding element in a2."""
    for i, x in enumerate(a2): a1[i] += x
    return a1 

def matrix_append(mat, a):
    """Add a new column to a matrix."""
    for i, x in enumerate(a): mat[i].append(x)
    return mat

def flatten(seq, remove_none=True):
    """Flatten an iterable potentially composted of nested iterables.
    
    Returns:
        A list that contains all elements retrieved from the iterable and all recursively contained 
        sub-iterables. None values are removed if remove_none is True.
    
    Examples::
        >>> [1, 2, [3,4], (5,6)]
        [1, 2, [3, 4], (5, 6)]
        >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)], False)
        [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]
        >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)], True)
        [1, 2, 3, 42, 4, 5, 6, 7, 8, 9, 10]
    """
    result = []
    for el in seq:
        if isinstance(el, list) or isinstance(el, tuple):
            result.extend(flatten(el))
        elif not remove_none or el is not None:
            result.append(el)
    return result

def extremes(iterable, key=None):
    """Returns the min and max of iterable."""
    l = list(iterable)
    l.sort(key=key)
    return (l[0],l[-1])

def charrange(start='A', end='Z', step=1):
    """Generator function for a range of characters."""
    for c in xrange(ord(start), ord(end)+1, step):
        yield chr(c)

def zipiter(iterables):
    iterators = map(iter, iterables)
    while iterators:
        yield map(next, iterators)

def index_map(l):
    im = OrderedDict()
    for i,val in enumerate(l):
        im[val] = i
    return im

class wrapiter(object):
    """
    Wrap an interable (potentially None) so that next() always returns None 
    when the iterable is exhausted.
    """
    def __init__(self, iterable=None):
        if iterable:
            self.iterable = iter(iterable)
        else:
            self.curr = None

    def next(self):
        if not hasattr(self, 'curr') or self.curr:
            try:
                self.curr = self.iterable.next()
            except:
                self.curr = None
        return self.curr
            
# dict methods

def merge(d1, d2, fn):
    """Update dict ``d1`` with dict ``d2`` using ``fn`` to resolve conflicts.
    
    Like update, but instead of overwriting the value in d1 with d2, any time a key exists in both 
    hashes ``fn`` is called and the result is the new value of the key.
    
    Args:
        d1 (dict): destination
        d2 (dict): source
        fn (function): used to resolve conflicts when a key appears in both tables
    
    Returns:
        d1
    """
    for key, value in d2.iteritems():
        if key in d1:
            d1[key] = fn(d1[key], value)
        else:
            d1[key] = value
    return d1

# miscelaneous collection classes

class Cache(object):
    """Caches arbitrary values. Can be extended to define initializer functions.
    class MyCache(Cache):
        def get_foo(self):
            return "bar"
    
    for i in xrange(10):
        # get_foo is only called in the first iteration; after that
        # the cached value is returned
        print MyCache()["foo"]
    """
    def __init__(self):
        self.cache = {}
        
    def __getitem__(self, name):
        return self.get(name)
    
    def get(self, name, factory=None):
        if name not in self.cache:
            if factory is None:
                factory_name = "get_{0}".format(name)
                if hasattr(self, factory_name):
                    factory = getattr(self, factory_name)
            if factory is not None:
                self.cache[name] = factory()
            elif hasattr(self, "default"):
                factory = getattr(self, "default")
                self.cache[name] = factory(name)
            else:
                raise KeyError(name)
        return self.cache[name]
    
class Range(Sequence):
    def __init__(self, start, end, incl=(True,False)):
        self.start = start
        self.end = end
        self.incl = incl
        self._first = start + (not incl[0])
        self._last = end - (not incl[1])
        if self._first > self._last:
            raise ValueError("start must be <= end")
        
    def __len__(self):
        return self._last - self._first + 1

    def __getitem__(self, i):
        if abs(i) >= len(self):
            raise IndexError("Index out of bounds: %d" % i)
        
        if i < 0: 
            return i + self._last + 1
        else:
            return i + self._first

    def __iter__(self):
        return iter(xrange(self._first, self._last + 1))
    
    def __repr__(self):
        return "{0}{1}-{2}{3}".format('[' if self.incl[0] else '(', 
            self.start, self.end, ']' if self.incl[1] else ')')

    def __cmp__(self, other):
        """
        Two ranges are equal if they intersect, even if their start and end points are 
        not the same. A Range is equal to an integer if it contains that integer.
        """
        if self.intersects(other):
            return 0
        else:
            if isinstance(other, Range):
                other = other._first
            return -1 if self._first < other else 1
    
    def slice(self, l):
        return l[self.first:(self.last+1)]
    
    def identical(self, other):
        if not isinstance(other, Range):
            return False
        return (self._first == other._first and self._last == other.last)
        
    def contains(self, other):
        if isinstance(other, int):
            return other >= self._first and other <= self._last
        else:
            return self._first <= other._first and self._last >= other._last
        
    def intersects(self, other):
        if isinstance(other, int):
            return self.contains(other)
        else:
            return self._first <= other._last and self._last >= other._first

class View(Sequence):
    def __init__(self, alist, rng):
        self.rng = rng
        self.alist = alist

    def __len__(self):
        return len(self.rng)

    def __getitem__(self, i):
        return self.alist[self.rng[i]]
    
class SparseTable(object):
    """A two-dimensional table (rows and columns)."""
    
    def __init__(self, seq=None, row_sort=sorted, col_sort=sorted, index_name=u"", missing=u"", delim=","):
        """Create a SparseTable.
        
        Args:
            seq (iterable): Initial values to add to table.
            row_sort (callable or list): Either a function that sorts row labels into the order 
                they will be printed (top to bottom), or a list of row labels in sorted order.
            col_sort (callable or list): Either a function that sorts col labels into the order 
                they will be printed (left to right), or a list of column labels in sorted order.
            index_name: Name of the first (index) column.
            missing: Value used when a cell is empty.
            delim: Default delimiter to use when printing a table.
        """
        self.__d = {}
        self.__colnames = set()
        self.row_sort = row_sort
        self.col_sort = col_sort
        self.index_name = index_name
        self.missing = missing
        self.delim = delim
        if seq is not None: 
            self.update_cells(seq)
    
    def __repr__(self):
        return "SparseTable([%s])" % ",".join([repr(t) for t in self.itercells()])
    
    def __str__(self):
        with buf() as s:
            self.print_table(s)
            return s.getvalue()

    # container interface

    def __len__(self):
        """Returns the number of rows in this table."""
        return len(self.__d)
    
    def __contains__(self, key):
        """Tests if a row or cell exists in this table.
        
        If ``key`` is a (rowname,colname) tuple, test whether that cell exists (regardless of 
        whether it has a value). Otherwise ``key`` is treated as a rowname and this method tests
        for existence of that row.
        """
        if isinstance(key, tuple) and len(key) == 2:
            return key[0] in self.__d and key[1] in self.__colnames
        else:
            return key in self.__d
        
    def __getitem__(self, key):
        """Get a row or cell.
        
        If ``key`` is a (rowname,colname) tuple, return the value of that cell. Otherwise ``key`` 
        is treated as a rowname and the corresponding row tuple is returned.
        """
        if isinstance(key, tuple) and len(key) == 2:
            return self.get_cell(key[0], key[1])
        else:
            return self.get_row(key)
    
    def __setitem__(self, key, value):
        """Update a row or cell.
        
        If ``key`` is a (rowname,colname) tuple, update the value of that cell. Otherwise ``key`` 
        is treated as a rowname and the corresponding row is updated.
        """
        if isinstance(key, tuple) and len(key) == 2:
            return self.update_cell(key[0], key[1], value)
        else:
            return self.update_row(key, value)
    
    def __delitem__(self, key):
        """Delete a row or cell.
        
        If ``key`` is a (rowname,colname) tuple, that cell is deleted. Otherwise ``key`` 
        is treated as a rowname and the corresponding row is deleted.
        """
        if isinstance(key, tuple) and len(key) == 2:
            return self.delete_cell(key[0], key[1])
        else:
            return self.delete_row(key)
    
    def __iter__(self):
        """Returns a row iterator."""
        return self.iterrows()
        
    @property
    def rownames(self):
        return tuple(self._sort(self.__d.iterkeys(), self.row_sort)) # TODO: cache?
    
    @property
    def colnames(self):
        return tuple(self._sort(self.__colnames, self.col_sort)) # TODO: cache?
    
    @staticmethod
    def _sort(l, fn_or_list):
        if callable(fn_or_list):
            return fn_or_list(l)
        elif isinstance(fn_or_list, list):
            return fn_or_list
        else:
            return l
    
    def get_row(self, rowname, insert_missing=False):
        """Get the column values for a given row.
        
        Args:
            rowname: Row to return.
            insert_missing (bool): Whether missing columns should be added.
        
        Returns:
            A view (i.e. non-modifiable) of the specified row dict.
        """
        assert rowname in self.__d
        row = self.__d[rowname]
        if insert_missing and len(self.__d) != len(self.__colnames):
            row = dict(row)
            for cname in self.__colnames - set(row.keys()):
                row[cname] = self.missing
        return DictView(row)
    
    def get_cell(self, rowname, colname, **kwargs):
        """Returns the value of the cell at the given row and column."""
        if 'missing' in kwargs:
            missing = kwargs['missing']
        else:
            missing = self.missing
        return self.__d[rowname].get(colname, missing) if rowname in self.__d else missing
            
    def iterrows(self, header=True):
        """Iterate over the table one row at a time.
        
        Args:
            header (bool): Whether to emit a header row.
        
        Yields:
            A tuple of column values for each row.
        """
        if header:
            i = self.index_name
            if isinstance(i, str):
                i = (i,)
            yield i + tuple(self.colnames)
        for rname in self.rownames:
            yield (rname,) + tuple(self.get_cell(rname, cname) for cname in self.colnames)
        
    def itercells(self, emit_missing=False):
        """Iterate over the table one cell at a time.
        
        Args:
            emit_missing (bool): Whether cells with None value should be emitted.
            
        Yields:
            A tuple (rowname,colname,cellvalue) for each cell.
        """
        from itertools import product
        for rname, cname in product(self.rownames, self.colnames):
            value = self.get(rname, cname, None)
            if value is None:
                if emit_missing:
                    value = self.missing
                else:
                    continue
            yield (rname, cname, value)
        
    def update_cells(self, seq):
        """Update the table from an iterable of (row,col,value) tuples.
        
        Args:
            seq: An iterable of (rowname,colname,value) tuples.
        
        Returns:
            A list of previous values for each cell in the seq.
        """
        return [self.update_cell(rname, cname, value) for rname, cname, value in seq]

    def update_row(self, rowname, colvalues):
        """Update a row in the table from a dict of column values.
        
        Args:
            rowname (str): row to update
            colvalues (dict): key is column name
        
        Returns:
            A dict of the row's previous values, or None if the row was added.
        """
        old = None
        if rowname in self.__d:
            old = dict(self.__d[rowname])
        for colname, value in colvalues.iteritems():
            self.update_cell(rowname, colname, value)
        return old

    def update_cell(self, rowname, colname, value):
        """Update the value at the specified row and column.
        
        Returns:
            The previous value at the specified cell, if any.
        """
        old = None
        if rowname in self.__d:
            row = self.__d[rowname]
        else:
            row = {}
            self.__d[rowname] = row
        old = row.get(colname, None)
        row[colname] = value
        self.__colnames.add(colname)
        return old
    
    def delete_row(self, rowname):
        assert rowname in self.__d
        del self.__d[rowname]
        
    def delete_cell(self, rowname, colname):
        assert rowname in self.__d
        assert colname in self.__colnames
        row = self.__d[rowname]
        if colname in row:
            del row[colname]
        
    def print_table(self, to=sys.stdout, header=True, delim=None, newline="\n"):
        if delim is None:
            delim = self.delim
        for row in self.iterrows(header):
            to.write(delim.join(map(str, flatten(row))))
            to.write(newline)

class PropDict(DictMixin):
    """An dict-like objects whose elements can be accessed either as properties or indexes."""
    
    __hash__ = None
    
    def __init__(self, *args, **kwargs):
        if len(args) > 0:
            assert len(args) == 1 and hasattr(args, "__getitem__")
            self.__dict__.update(args[0])
        self.__dict__.update(kwargs)
    
    def __getitem__(self, key):
        if not hasattr(self, key):
            raise KeyError("Invalid key {0}".format(key))
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        if not hasattr(self, key):
            raise KeyError("Invalid key {0}".format(key))
        delattr(self, key)
    
    def __eq__(self, other):
        return vars(self) == vars(other)
    
    def __ne__(self, other):
        return not (self == other)
    
    def __repr__(self):
        return "PropDict({0})".format(
            ",".join("{0}={1}".format(name, value) for name, value in sorted(self.__dict__.items())))
    
    def __contains__(self, key):
        return key in self.__dict__
    
    def keys(self):
        return self.__dict__.keys()
    
    def copy(self):
        return PropDict(self)

class Lockable(object):
    """Mixin that manages a single variable: locked.
    
    ``locked`` is initially False and can only be set to True. When ``locked`` is set to True,
    the ``_do_lock`` method is called.
    """
    @accessor
    def locked():
        def fget(self):
            return getattr(self, '__locked', False)
        def fset(self, value):
            assert isinstance(value, bool)
            if getattr(self, '__locked', False):
                assert value, "Cannot unlock a locked %s"
            setattr(self, '__locked', value)
            self._do_lock()
        return dict(fget=fget, fset=fset)
    
    def _do_lock(self):
        """Override this method to perform custom locking behavior."""
        pass

class DictView(UserDict):
    """Wraps a dict and disallows modification."""
    def __init__(self, d):
        self.data = d

    def __setitem__(self, key, val):
        raise Exception("This dictionary view cannot be modified")

    def __delitem__(self, key):
        raise Exception("This dictionary view cannot be modified")
            
class LockableDict(dict, Lockable):
    """A dict that, when locked, prevents addition or deletion of items."""
    
    def __setitem__(self, key, val):
        assert not self.locked, "This dict has been locked"
        dict.__setitem__(self, key, val)
    
    def __delitem__(self, key):
        assert not self.locked, "This dict has been locked"
        dict.__delitem__(self, key)
        
class VarArgGenerator(Lockable):
    """Generates arg lists for every combination of a set of variables.
    
    Manages two dicts: one for constants and one for variables. Generates arg lists for every
    combination of variables, with constants added in to each list.
    """
    
    __slots__ = ['variables']
    
    def __init__(self, variables=None):
        self.variables = LockableDict()
        if variables:
            self.variables.update(variables)
        self.opers = defaultdict(lambda: {})
        self.__locked = False

    def __repr__(self):
        return "VarArgGenerator(%s)" % repr(self.variables)
    
    @overrides(Lockable._do_lock)
    def _do_lock(self):
        self.variables.lock()
    
    def add(self, key, value):
        self.variables[key] = value
    
    def add_oper(self, src_key, dest_key, fn):
        self.opers[src_key][dest_key] = fn
    
    def copy(self, locked=False):
        """Create a copy of this VarArgGenerator.
        
        Args:
            locked (bool): If None, the current locked status is maintained.
        """
        new = VarArgGenerator(self.variables)
        new.locked = locked if locked is not None else self.locked
        return new
    
    def update(self, other):
        assert not self.locked, "Cannot update a locked VarArgGenerator"
        if isinstance(other, VarArgGenerator):
            self.variables.update(other.variables)
        elif isinstance(other, dict):
            self.variables.update(other)
        else:
            raise TypeError("Cannot update from %s" % type(other))
    
    def __len__(self):
        """Returns the number of combinations that an iterator will produce."""
        size = 1
        for name, value in self.variables.iteritems():
            if measurable(value) and len(value) > 1:
                size *= len(value)
        return size
                
    def __iter__(self):
        args = {}
        variables = []
        divisors = []
        num_combos = 1
        
        for name, value in self.variables.iteritems():
            if not measurable(value):
                args[name] = value
            elif len(value) <= 1:
                args[name] = value[0]
            else:
                size = len(value)
                variables.append((name,value))
                divisors.append(size)
                num_combos *= size
            
        for i in range(0, num_combos):
            iargs = dict(args)
            n = i
            for v in range(0, len(variables)):
                (name,value) = variables[v]
                div = divisors[v]
                iargs[name] = value[n % div]
                n = n / div
            for src_key, val in self.opers.iteritems():
                for dest_key, fn in val.iteritems():
                    iargs[dest_key] = fn(iargs[src_key])
            yield iargs

class BitfieldEnum(object):
    """An enumeration where each value is associated with a power of two."""
    
    def __init__(self, *names):
        self.i_to_n = {}
        self.n_to_i = {}
        for x,n in enumerate(names):
            i = pow(2, x)
            self.i_to_n[i] = n
            self.n_to_i[n] = i
    
    def names(self):
        return self.n_to_i.keys()
        
    def add_or(self, name, *names):
        i = self.or_names(*names)
        self.i_to_n[i] = name
        self.n_to_i[name] = i
        return i
    
    def get_or(self, *names):
        return self.i_to_n[self.or_names(*names)]
    
    def or_names(self, *names):
        i = 0
        for n in names:
            i = i | self.n_to_i[n]
        return i

# Tree functions/classes

class TreeNode(OrderedDict):
    def __init__(self, name, value=None, key=None):
        OrderedDict.__init__(self)
        self.key = key
        self.name = name
        self.value = value
        
    @property
    def is_leaf(self):
        return len(self) == 0
    
    def add(self, child):
        if isinstance(child, str):
            child = TreeNode(child)
        self[child.name] = child
        return child
    
    def sorted_items(self):
        return sorted(self.iteritems())
    
    def __repr__(self):
        sup = dict.__repr__(self)
        if self.value and sup:
            return "<%s>%s" % (self.value, sup)
        else:
            return sup
            
def build_tree(rows, keys, xform=None, root_name='root', node_class=TreeNode):
    """Builds a tree from an iterable of dicts. 

    The depth of the tree is determined by the length of keys. A key can be a string or a tuple
    that defines a compound key. The value stored for a given row is by default the row itself,
    but this can be altered by passing xform, which can be either a function or a list of keys
    to keep from the row."""

    tree = node_class(root_name)
    for r in rows:
        node = tree
        for key in keys:
            if isinstance(key, tuple):
                val = tuple(r[k] for k in key)
            else:
                val = r[key]
            if not val in node:
                node.add(node_class(val, key=key))
            node = node[val]
        
        if xform is not None:
            if callable(xform):
                r = xform(r)
            else:
                r = dict((k,r[k]) for k in xform)
        
        node.value = r
    return tree

def rollup(tree, columns):
    """Rolls up column(s) at each level of a tree.

    Columns are specified as a dict with key being a dict key or attribute name for accessing 
    the column in leaf nodes, and value being a function that accepts an iterable and returns 
    the rollup value. The value of each non-leaf node is overwritten with a dict containing the
    rollup values.
    """
    assert isinstance(tree, TreeNode)
    values = dict((k,[]) for k in columns)
    for node in tree:
        child = tree[node]
        if not child.is_leaf:
            rollup(child, columns)
        for col,fn in columns.iteritems(): 
            values[col].append(
                child.value[col] if isinstance(child.value, dict) else getattr(child.value, col))
    tree.value = dict((k, columns[k](values[k])) for k in columns)

def walk_tree(tree, fn, leaves_only=False):
    """Visits each node of the tree and calls 'fn' with arguments ('parent', 'node'). Returns a
    tuple of all non-None results returned from 'fn'.
    """
    assert isinstance(tree, TreeNode)
    
    results = []
    
    def _walk(node, path, fn, leaves_only):
        if node.is_leaf or not leaves_only:
            result = fn(tree, path, node)
            if result:
                results.append(result)
        new_path = path + (node.name,)
        for child in node.values():
            _walk(child, new_path, fn, leaves_only)
    
    _walk(tree, (), fn, leaves_only)
    
    return results
