#!/usr/bin/env python

#!/usr/bin/env python

from collections import defaultdict
from csv import reader, writer
import sys

def main():
    junctions = defaultdict(lambda: [])
    w = writer(sys.stdout, delimiter="\t")
    for row in reader(sys.stdin, delimiter="\t"):
        key = "{0}:{1}-{2}:{3}".format(row[0], row[1], row[2], row[5])
        junctions[key].append(row)
    for row_set in junctions.values():
        if len(row_set) == 1:
            w.writerow(row_set[0])
        else:
            row = row_set[0]
            row[4] = sum(map(lambda r: int(r[4]), row_set))
            w.writerow(row)

if __name__ == "__main__":
    main()