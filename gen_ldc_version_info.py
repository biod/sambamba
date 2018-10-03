#!/usr/bin/env python

from __future__ import print_function
import re, sys, subprocess

if len(sys.argv) < 2:
    print("Usage: {0} <path to ldmd2 executable>".format(sys.argv[0]))
    sys.exit(1)

ldc = sys.argv[1].replace("ldmd2", "ldc2")
ldc_output = subprocess.Popen([ldc, '-version'], stdout=subprocess.PIPE).communicate()[0]
version_re = r"""^.+\((?P<LDC>[^\)]+)\):\n\s*based on DMD (?P<DMD>\S+) and LLVM (?P<LLVM>\S+)\n(?:\s*built with (?P<BOOTSTRAP>.*)\n)?"""
match = re.match(version_re, ldc_output.decode("utf-8") , re.MULTILINE)

if not match:
    sys.exit("ERROR: failed to generated LDC version information")

print("module utils.ldc_version_info_;")
for component, version in match.groupdict().items():
    if version is None:
        version = "version not available"
    print("immutable {0}_VERSION_STRING = \"{1}\";".format(component, version))
