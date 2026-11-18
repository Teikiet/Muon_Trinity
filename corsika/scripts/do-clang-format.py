#!/usr/bin/env python3

#
# (c) Copyright 2020 CORSIKA Project, corsika8@kit.edu
#
# This software is distributed under the terms of the 3-clause BSD license.
# See file LICENSE for a full version of the license.
#

"""
Run clang-format with the style file in the CORSIKA repository.

By default it finds new files and files with modifications with respect 
to the current main branch and prints the filenames which need clang-formatting. 
Returns 1 if there are files which need modifications and 0 otherwise, 
so it can be used as a test.
"""

import argparse
import subprocess as subp
import os
import sys
import re

debug = False
do_progress = False
try:    
    from progress.bar import ChargingBar
    do_progress = True
except ImportError:
    do_progress = False
    # no progress bar

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--apply', action="store_true",
    help="Apply clang-format to files which need changes.")
parser.add_argument("--all", action="store_true",
    help="Check all files below current path instead of new/modified.")
parser.add_argument("--docker", action="store_true",
    help="Use corsika/devel:clang-8 container to run clang-format. Make sure you are in group \"docker\". ")

args = parser.parse_args()

excludeDirs = [r"^(\./)?modules/", r"^(\./)?externals/", r"^(\./)?build", r"^(\./)?install",
               r"(\./)?\.git", r"^(\./)?corsika/framework/units",
               r"^(\./)?Random123/"]

filelist = []
if args.all:
    for dirpath, dirnames, filenames in os.walk("."):
        excl = False
        for excl_dir in excludeDirs:
            if re.findall(excl_dir, dirpath):
                excl = True
                break
        if excl:
            continue
        for f in filenames:
            if (f.endswith(".hpp") or f.endswith(".cpp") or f.endswith(".inl")):
                filename = os.path.join(dirpath, f)
                if not os.path.islink(filename):
                    filelist.append(filename)
    if not filelist:
        raise SystemExit("Error: You specified --all, but file list is empty. "
                         "Did you run from the build directory?")
else:
    cmd = "git diff main --name-status"
    for line in subp.check_output(cmd, shell=True).decode("utf8").strip().split("\n"):
        if line.startswith("D"): continue
        if line.startswith("R"):
            filelist.append(line.split()[-1])
        else:
            filelist.append(line[1:].lstrip())

    cmd = "git ls-files --exclude-standard --others"
    filelist2 = subp.check_output(cmd, shell=True).decode("utf8").strip().split("\n")
    filelist += filelist2
    # some cleanup
    filelist_clean = []
    for f in filelist:
        if not (f.endswith(".hpp") or f.endswith(".cpp") or f.endswith(".inl")):
            continue
        if os.path.islink(f):
            continue
        excl = False
        for excl_dir in excludeDirs:
            if re.findall(excl_dir, f):
                excl = True
                break
        if excl:
            continue
        filelist_clean.append(f)
    filelist = filelist_clean

if debug:
    print ("filelist: ", filelist)
    
cmd = os.getenv("CLANG_FORMAT", "clang-format")

if args.docker: 
  USER=os.environ["USER"]
  UID=os.getuid()
  GID=os.getgid()
  PWD=os.getcwd()
  # note, currently in container it is clang-8 
  cmd = "docker container run --rm -v {}:/corsika -w /corsika -u {}:{} corsika/devel:clang-8 clang-format-8".format(PWD,UID,GID)
cmd += " -style=file"

version = subp.check_output(cmd.split() + ["--version"]).decode("utf-8")
print (version)
print ("Note: the clang-format version has an impact on the result. Make sure you are consistent with current CI. Consider \'--docker\' option.")


bar = None
if do_progress:
    bar = ChargingBar('Processing', max=len(filelist))

if args.apply:
    changed = []
    for filename in filelist:        
        if bar: 
            bar.next()
        a = open(filename, "rb").read()
        subp.check_call(cmd.split() + ["-i", filename])
        b = open(filename, "rb").read()
        if a != b:
            changed.append(filename)
    if bar: 
        bar.finish()
    if debug:
        print ("changed: ", changed)
    
else:
    # only print files which need formatting
    files_need_formatting = 0
    for filename in filelist:
        if bar:
            bar.next()
        a = open(filename, "rb").read()
        b = subp.check_output(cmd.split() + [filename])
        if a != b:
            files_need_formatting += 1
            print(filename)
    if bar:
        bar.finish()            
    sys.exit(1 if files_need_formatting > 0 else 0)
