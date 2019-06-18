#!/usr/bin/python

#---------------------
# (c) Copyright, 2011 by the Regents of the University of California.
#Author: Ji Qiang.
#Description: Energy scan of driven phase (output file engout, col. 5 vs. 1).
#

import os, sys, string, pdb

#
# Note: python's index start from 0
# On Windows, this script has to be run under Cygwin
# The Cygwin can downloaded from: http://cygwin.com/install.html
#

#location of the variable in the "test.in" to be scanned.
#new_row is the line number -1 of that phase variable in the file.
#new_col is the location of phase variable in that line - 1. 
new_row = 40 #41 row
new_col = 6 #7th column

valueini = 208.0 #initial phase
delv = 0.1 #phase scan step size
valueend = 210.0 #end phase

newValueIndex = valueini
while newValueIndex <= valueend:
    test5file = open('test.in','r')
    lines = test5file.readlines()
    test5file.close()

    # modifiy array[2][4]
    new_line = string.split(lines[new_row])
    new_line[new_col] = str(newValueIndex)
    new_line = string.join(new_line) + '\n'

    lines[new_row] = new_line
    test5file2 = open('test.in', 'w')
    test5file2.writelines(lines)
    test5file2.close()
       
    os.system('./ImpactZ.exe > a')
    os.system('tail -1 fort.18 >> tmpeng')
    os.system('echo '+str(newValueIndex)+' >> tmpphase')
    print 'index: ',newValueIndex
    os.system('tail -1 fort.18')
    newValueIndex = newValueIndex+delv

os.system('paste tmpphase tmpeng > engout')
os.system('rm tmpeng tmpphase')
