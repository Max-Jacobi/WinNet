#!/usr/bin/env python

# The reference was produced with the same trajectory, but with the neutrino
# source given as luminosity instead of number flux (see Notes in the archive).
t.checklist = { \
   'finab.dat'    : { 'method':'default', 'tolerance':1e-9 }, \
   'mainout.dat'  : { 'method':'listcompare', 'tolerance':[1e-8], 'x_column':1, 'y_column':[4]}, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
