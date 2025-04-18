#!/usr/bin/env python3

import getopt
import glob
import os
import sys

import gi
gi.require_version('Gdk', '3.0')
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk as gtk
from gi.repository import Gdk as gdk

#from gi.repository import Gtk as gtk 
#import gtk.gdk

import diamond.debug as debug
import diamond.schema as schema

def Help():
  debug.dprint("Usage: update_options [OPTIONS] ... [FILES]\n" + \
               "\n" + \
               "Updates flml, bml, swml and adml files. If FILES is not specified, all .xml and mpml files in\n" + \
               "tests/*/., tests/*/*/. and examples/*/. will be updated. Options:\n" + \
               "\n" + \
               "-h  Display this help\n" + \
               "-v  Verbose mode", 0)

  return

try:
  opts, args = getopt.getopt(sys.argv[1:], "hv")
except getopt.getoptError:
  Help()
  sys.exit(-1)

if ("-h", "") in opts:
  Help()
  sys.exit(0)

if not ("-v", "") in opts:
  debug.SetDebugLevel(0)

rootDir = os.path.join(os.path.dirname(__file__), os.path.pardir)
testDir = os.path.join(rootDir, "tests")
longtestDir = os.path.join(rootDir, "longtests")
examplesDir = os.path.join(rootDir, "examples")

extdict = {"xml"  : "test_options.rng",
           "mpml" : "multiphase.rng"}

# cache parsed schema files
schemadict = {}
for k,v in list(extdict.items()):
  schemadict[k] = schema.Schema(os.path.join(rootDir, "schemas", v))

filenames = args
if len(filenames) == 0:
  filenames = []
  for k,v in list(extdict.items()):
    filenames += glob.glob(os.path.join(testDir, "*", "*."+k))
    filenames += glob.glob(os.path.join(testDir, "*", "*", "*."+k))
    filenames += glob.glob(os.path.join(examplesDir, "*", "*."+k))

invalidFiles = []
updated = 0
for filename in filenames:
  debug.dprint("Processing " + str(filename), 1)
  
  ext = filename.split(".")[-1]
  sch = schemadict[ext]
 
  # Read the file and check that either the file is valid, or diamond.schema
  # can make the file valid by adding in the missing elements
  optionsTree = sch.read(filename)
  lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
  if len(lost_eles) + len(lost_attrs) > 0 or not optionsTree.valid:
    debug.deprint(str(filename) + ": Invalid", 0)
    debug.deprint(str(filename) + " errors: " + str((lost_eles, added_eles, lost_attrs, added_attrs)), 1)
    invalidFiles.append(filename)
    continue
  
  # Write out the updated options file
  optionsTree.write(filename)
  debug.dprint(str(filename) + ": Updated", 0)
  updated += 1

debug.dprint("Summary:", 0)
debug.dprint("Invalid options files:", 0)
for filename in invalidFiles:
  debug.dprint(filename, 0)
debug.dprint("Invalid: " + str(len(invalidFiles)), 0)
debug.dprint("Updated: " + str(updated), 0)
