#!/usr/bin/env python

from __future__ import print_function
import re, sys, subprocess

if len(sys.argv) < 2:
    print("Usage: {} <path to ldmd2 executable>".format(sys.argv[0]))
    sys.exit(1)

ldc = sys.argv[1].replace("ldmd2", "ldc2")
ldc_output = subprocess.check_output([ldc, "-version"])
version_re = r"""^.+\((?P<LDC>[^\)]+)\):\n\s*based on DMD (?P<DMD>\S+) and LLVM (?P<LLVM>\S+)\n\s*built with (?P<BOOTSTRAP>.*)\n"""
match = re.match(version_re, ldc_output, re.MULTILINE)

print("module utils.ldc_version_info_;")
for component, version in match.groupdict().items():
    print("immutable {}_VERSION_STRING = \"{}\";".format(component, version))
