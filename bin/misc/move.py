#!/usr/bin/env python
import os
import sys
import glob
import shutil

for f in glob.glob(sys.argv[1]):
    o = os.path.join(sys.argv[2], os.path.basename(f))
    shutil.move(f, o)