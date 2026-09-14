#!/usr/bin/env python

# The reference mainout.dat is the one of the read_trajectory test, i.e., the
# same trajectory in cgs units. Comparison is a linear interpolation, hence the
# relatively large tolerance.
t.checklist = { \
   'mainout.dat'  : { 'method':'listcompare', 'tolerance':[1e-10,1e-10] ,'x_column':1, 'y_column':[2,3]}, \
}
t.program = t.basedir + "/bin/winnet"
t.testdir = t.basedir + "/test/" + t.testname
t.logfile = t.testdir + ".log"
t.arcfile = t.testdir + ".tar.gz"
