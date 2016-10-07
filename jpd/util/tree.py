
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
