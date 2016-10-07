
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
