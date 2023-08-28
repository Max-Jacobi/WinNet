#!/usr/bin/env python

t.checklist = { \
   'flow/flow_0010.dat'    : { 'method':'flowcompare', 'tolerance':2.0e-2, 'lowerlimit':1.0e-6 }, \
   'flow/flow_0020.dat'    : { 'method':'flowcompare', 'tolerance':2.0e-2, 'lowerlimit':1.0e-6 }, \
   'flow/flow_0030.dat'    : { 'method':'flowcompare', 'tolerance':2.0e-2, 'lowerlimit':1.0e-6 }, \
   'flow/flow_0040.dat'    : { 'method':'flowcompare', 'tolerance':2.0e-2, 'lowerlimit':1.0e-6 }, \
   'flow/flow_0050.dat'    : { 'method':'flowcompare', 'tolerance':2.0e-2, 'lowerlimit':1.0e-6 }, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
