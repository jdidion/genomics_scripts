
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
