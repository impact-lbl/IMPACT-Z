#!/usr/bin/python
import os, sys, string, pdb, math, fileinput

#
# Note: python's index start from 0
# On Windows, this script has to be run under Cygwin
# The Cygwin can downloaded from: http://cygwin.com/install.html
#
#---------------------
# (c) Copyright, 2012 by the Regents of the University of California.
#Authors: T. Wilson and J. Qiang.
#Description: PhaseOptz.py beta v. 1.0: 
#             Find the RF driven phase used in the test.in given the
#             design phase initially specified in the test.in.
#             The initial test.in is copied into file test_backup.in.
#             The found RF driven phases are put into the test.in file.
#             <NOTE>: Here, the executable file name of the IMPACTz code is IMPACTexe.
#------------------------------------------------------
# specify initial parameters for optimization
#------------------------------------------------------
# starting and ending row location of lattice in test.in
# start reading test.in at row
row_start = 16
row_end = 48
# restarting function switch row location
row_restart = 3
# reference frequency row location
row_rffreq = 15
#------------------------------------------------------

# phase angle scan specifications
angle_in = 0
angle_end = 360
delv = 10

# column indices to read input file test.in
# 3: ele. id; 5: freq; 6: ang
col_el = 3
col_freq = 5
col_ang = 6

# column indices to read output file
# 4: eng
col_eng = 4

# fractional precision of phase angle to be determined
tol = 0.0001

# max number of iterations in Brent's Method before quitting
it_max = 20

#-----------------------------------------------------

def main():

   # remove previous data files
   os.system('rm -f ang_* eng_* data_*')

   # create backup of test.in
   os.system('cp test.in test_backup.in')

   # from row_start to row_end, optimize beam elements one at a time
   row_end = fileLength('test.in')

   # default values for variables
   zedge = -1
   max_ang = -1
   final_ang = -1
   length = -1
   freq = -1
   prediction = -1

   # global values from initial parameter specification
   # Angle_in and angle_end will be global variables to store the start and end points of the initial scan.
   # However, for later scans, they are changed if we can make a good prediction of the angle
   # (if beta is close enough to 1). When this happens, to not lose the initial values angle_in and angle_end,
   # we store them in angle_inVal and angle_endVal.
   global angle_in
   global angle_end
   angle_inVal = angle_in
   angle_endVal = angle_end

   test5file = open('test.in','r')
   lines = test5file.readlines()
   test5file.close()

   reffreq = float(string.split(lines[row_rffreq])[4])
   omega = 2*3.1415926*reffreq

   #find the number of RF cavities
   nRF = 0
   for i in range(row_start, row_end):
      # read from test.in
      file = open('test.in', 'r')
      lines = file.readlines()
      file.close()

      if lines[i][0] != '!':
         if int(string.split(lines[i])[col_el]) > 100 and float((string.split(lines[i])[col_freq]).replace('d','E')) !=  0:
           nRF = nRF + 1
         else:
           pass

   print 'nRF:', nRF

   # scan through the rows of beam elements in test.in
   # row_start and row_end are acceptable here because the -7 and -99 lines are added and removed WITHIN each
   # iteration of the loop.
   zedge0 = 0.0
   length0 = 0.0
   for i in range(row_start, row_end+nRF):

      # set angle_in and angle_end back to default values stored in angle_inVal and angle_endVal
      # This is reset for every row, just in case our prediction is bad (cannot find a max within +/- 20 degrees
      # of the prediction), so we are forced to scan the entire default range to find a max. 
      angle_in = angle_inVal
      angle_end = angle_endVal

      # read from test.in
      file = open('test.in', 'r')
      lines = file.readlines()
      file.close()

      # if not a commented line
      if lines[i][0] != '!':
         
         length0 = float(string.split(lines[i])[0])

         # if bpm > 100 and rf field != 0
         if int(string.split(lines[i])[col_el]) > 100 and float((string.split(lines[i])[col_freq]).replace('d','E')) !=  0:

            # read in the user specified angle relative to zero phase and frequency
            spec_ang = string.split(lines[i])[col_ang]
            prev_freq = freq
            freq = float(string.split(lines[i])[col_freq].replace('d','E'))

            print 'finding max for beam element in line '+str(i)+'...'


            # insert -99 and -7 line after line i

            # 'x' is the empty string that will be the string version of test.in. We insert -7 and -99 lines
            # after line i.
            x = ''
            # makes sure that we only insert one set of -7 and -99 lines (this could have happened if we have
            # duplicate beam elements)
            inserted = False

            for line in fileinput.input('test.in'):
               if line == lines[i] and inserted == False:
                  # makes sure that we only insert one set of -7 and -99 lines
                  # more than one can be inserted if we have duplicate beam elements
                  inserted = True

                  # set position parameters for -7 and -99 line
                  length_str = string.split(lines[i])[0]
                  prev_length = length
                  # convert scientific notation from fortran 'd' to python 'E'
                  length = float(length_str.replace('d','E'))
                  #zedge_str = string.split(lines[i])[4]
                  prev_zedge = zedge
                  #prev_zedge = zedge - length
                  #zedge = float(zedge_str.replace('d','E'))
                  zedge = zedge0 
                  cav_end = zedge + length
                  # create -7 and -99 lines
                  #tmp3line = '0 0 51 -7 '+str(cav_end)+' '+str(cav_end)+' '+str(cav_end)+' / \n'
                  #tmp99line = '0 0 0 -99 '+str(cav_end)+' '+str(cav_end)+' '+str(cav_end)+' / \n'
                  tmp3line = '0 0 51 -7 '+' / \n'
                  tmp99line = '0 0 0 -99 '+' / \n'
                  x = x+line
                  x = x+tmp3line+tmp99line

               else:
                  x = x+line


            # overwrite test.in with -7 and -99 lines
            file = open('test.in', 'w')
            file.write(x)
            file.close()

            # predict next phase angle if RstartFlag = 1 (if this isn't the first beam element)
            if int(string.split(lines[row_restart])[1]) == 1:
               prediction = predict(zedge, prev_zedge, max_ang, prev_freq, prev_length,omega)
            
            # determine phase angle for which energy is a maximum
            max_ang = findMax(i, cav_end, prediction)
            
            # write user specified angle relative to zero phase (max angle) into test.in
            file = open('test.in','r')
            lines = file.readlines()
            file.close()
            # define user specified angle
            new_line = string.split(lines[i])
            final_ang = max_ang+float(spec_ang)
            # confine angle to 0-360 degrees
            final_ang = angMod(final_ang)
            # write user specified angle into test.in
            new_line[col_ang] = str(final_ang)
            new_line = string.join(new_line)+'\n'
            lines[i] = new_line
            file = open('test.in','w')
            file.writelines(lines)
            file.close()
            
            file = open('test.in','r')
            lines = file.readlines()
            file.close()
            # if RstartFlag = 1 (if this isn't the first beam element)
            if int(string.split(lines[row_restart])[1]) == 1:
               # copies previous fort.51 file to be used to restart (more information on this in the returnEng() function)
               os.system('cp temp_fort.51 fort.51')

            # return energy of specified angle (this makes the last line in 'engout_spec' have the data we want)
            # this keeps 'engout''s last line as still being for the max energy
            # run ImpactZ simulation and output eng data to a file
            os.system('./IMPACTexe > a')
            os.system('tail -1 fort.18 >> tmpeng')
            os.system('echo '+str(final_ang)+' >> tmpphase')
            os.system('paste tmpphase tmpeng > engout_spec')
            os.system('rm tmpeng tmpphase')
            file = open('engout_spec', 'r')
            line = file.readline()
            file.close
            eng = float(string.split(line)[col_eng].replace('d','E'))
            os.system(' echo ' + str(eng) + ' >> eng_data_' + str(i))
            end_pos = float(string.split(line)[1])
            print 'Optimized values'
            print 'angle: '+str(final_ang)
            print 'energy: '+str(eng)
            print '----------'

            # if we just optimized for the first element, set Rstartflag to 1
            # Rstartflag set to 1 makes the simulation restart from the last paused position; default set to 0.
            file = open('test.in', 'r')
            lines = file.readlines()
            file.close()
            if int(string.split(lines[row_restart])[1]) == 0:
               new_line = string.split(lines[row_restart])
               new_line[1] = str(1)
               new_line = string.join(new_line)+'\n'
               lines[row_restart] = new_line
               file = open('test.in', 'w')
               file.writelines(lines)
               file.close

            # remove -7 and -99 line from test.in
            file = open('test.in', 'r')
            lines = file.readlines()
            file.close()
            x = ''
            for line in fileinput.input('test.in'):
               #if line == tmp3line or line == tmp99line:
               if line == tmp99line:
                  pass
               else:
                  x = x+line
            file = open('test.in', 'w')
            file.write(x)
            file.close()

            # copy fort.51 into a temporary file temp_fort.51.
            # We will need this data for the next beam element because we restart simulation right after
            # current beam element.
            os.system('cp fort.51 temp_fort.51')

         zedge0 = zedge0 + length0
 

   # remove -7 and -99 line from test.in
   file = open('test.in', 'r')
   lines = file.readlines()
   file.close()
   x = ''
   for line in fileinput.input('test.in'):
      #if line == tmp3line or line == tmp99line:
      if line == tmp3line:
                  pass
      else:
                  x = x+line
   file = open('test.in', 'w')
   file.write(x)
   file.close()

   # After looping over all lines in test.in, we're done.
   # Set Rstartflag back to 0 for next run.
   file = open('test.in', 'r')
   lines = file.readlines()
   file.close()
   new_line = string.split(lines[row_restart])
   new_line[1] = str(0)
   new_line = string.join(new_line)+'\n'
   lines[row_restart] = new_line
   file = open('test.in', 'w')
   file.writelines(lines)
   file.close()

#------------------------------------------------------

# returns the number of lines in a file
def fileLength(filename):
   with open(filename) as f:
      i = 0
      for i, l in enumerate(f):
         pass
   return i + 1

#------------------------------------------------------

# confines angles to a range of 0-360 degrees; not necessary but makes things cleaner
def angMod(angle):
   while angle > 360:
      angle -= 360
   while angle < 0:
      angle += 360
   return angle

#-----------------------------------------------------

def predict(zedge, prev_zedge, max_ang, prev_freq, prev_length, omega):

   # Set default diff1 and diff2 to 1. diff1 and diff2 are difference values that tell us how close
   # we are to an edge of a certain cavity (see the 'for loops' later). We desire the smallest diff1
   # and diff2 we can get, so setting diff1 and diff2 to 1 by default is fine, as 1 is a suitably
   # 'large' difference.
   diff1 = 1
   diff2 = 1
   
   file = open('test.in', 'r')
   lines = file.readlines()
   file.close()

 
   # make prediction for phase angle to make angular scan range smaller, reduce running time

   # First, we need to find the elapsed time between the time when the beam entered and exited
   # the previous cavity. This will be found by looking at fort.18_max because it is not easily
   # calculable (the beam is being accelerated).
   
   # read in time from fort.18_max
   file = open('fort.18_max', 'r')
   lines18 = file.readlines()
   file.close

   #t1 = 0.0
   #t2 = 0.0
   for ind in range(len(lines18)):

      #print prev_zedge,prev_length,float(string.split(lines18[ind])[0]),ind

      # Check to see if the position of the beam in fort.18_max in a certain line is within 0.01
      # of the position of the leading edge of the previous cavity.
      if math.fabs(float(string.split(lines18[ind])[0]) - prev_zedge) < 0.01:
         # Save this difference; it will be useful for comparison, since we want the time value
         # for a position as close as we can get to the leading edge of the previous cavity
         diff_new1 = math.fabs(float(string.split(lines18[ind])[0]) - prev_zedge)
         # Save the time value
         t1 = float(string.split(lines18[ind])[1])/omega
         # If there is a line that provides a more accurate estimate, use it:
         if diff_new1 < diff1:
            diff1 = diff_new1
            t1 = float(string.split(lines18[ind])[1])/omega
      # This section of code is identical to the last, except we are findig the time at the back edge
      # instead of the leading edge of the previous cavity.
      if math.fabs(float(string.split(lines18[ind])[0]) - (prev_zedge + prev_length)) < 0.01:
         diff_new2 = math.fabs(float(string.split(lines18[ind])[0]) - (prev_zedge + prev_length))
         t2 = float(string.split(lines18[ind])[1])/omega
         if diff_new2 < diff2:
            diff2 = diff_new2
            t2 = float(string.split(lines18[ind])[1])/omega

   #print t1,t2,prev_zedge,prev_length
   # delt is the time elapsed between the beam entering and exiting the previous cavity
   delt = t2 - t1
   #delt = 0.0

   # Find out how fast the beam was traveling at the end of the previous cavity
   file = open('engout', 'r')
   line = file.readline()
   file.close()
   beta = float(string.split(line)[5])
   # 0.9975 is an empirical result based on accuracy of predictions for a range of beta values.
   # For beta values > 0.9975, the predictions were accurate to within about 20 degrees.
   if beta > 0.9975:
      print 'Estimating angle range:'
      # d is the distance between the end of the previous cavity and the beginning to the current cavity
      d = (zedge - prev_zedge)-prev_length
      # theta0 is the angle for which there was an energy maximum for the previous cavity
      theta0 = max_ang
      v = beta*3*10**8
      t = d/v
      # total_t is the predicted time from the zedge of previous cavity to zedge of current cavity
      total_t = t + delt
      # period is the number of periods elapsed during this time
      period = (total_t)*prev_freq
      # Calculate the fractional part of the period; tells us what part of the period the cavity is at
      # after the predicted elapsed time
      while period > 1:
         period -= 1

      # theta1 is the predicted angle in degrees
      theta1 = theta0 - 360*period
      theta1 = angMod(theta1)

      return theta1

   else:
      return -1

#------------------------------------------------------

def findMax(row, cav_end, prediction):
   # set default values to -1
   eng_max = -1
   ang_max = -1

   theta1 = prediction
   
   if theta1 != -1:
      # Set new angular scan range to prediction +/- 20 degrees
      global angle_in
      global angle_end
      angle_in = theta1 - 20
      angle_end = theta1 + 20

   # Sets angular scan range either to default values or predicted values.
   # range() arguments in Python must be integers.
   ang_in = int(angle_in)
   ang_end = int(angle_end)

   # Find phase angle for which energy is a maximum from angle_in to angle_end in increments of delv.
   print 'Initial scan from '+str(ang_in)+' to '+str(ang_end)+' by '+str(delv)+'...'
   for i in range(ang_in, ang_end+1, delv):
      # Find maximum energy and corresponding angle
      eng = returnEng(i, row, cav_end)
      if (eng > eng_max):
         eng_max = eng
         ang_max = i
      
   # If the maximum turns out to be either the max or min of the angular range, we have reason
   # to be suspicious that our initial prediction range doesn't contain the maximum.
   # This could be true either if:
   # 1) energies for angles less than ang_in are greater
   # 2) energies for angles greater than ang_end are greater
   if ((ang_max == ang_in and returnEng(ang_in-0.001, row, cav_end) > returnEng(ang_in, row, cav_end)) or  (ang_max == ang_end and returnEng(ang_end+0.001, row, cav_end) > returnEng(ang_in, row, cav_end))):
      print 'Try scan from 0 to 360 by 10...'
      for i in range(0, 360+1, 10):
         eng = returnEng(i, row, cav_end)
         if (eng > eng_max):
            eng_max = eng
            ang_max = i         

   print 'initial max angle: '+str(ang_max)
   print 'initial max energy: '+str(eng_max)
   
   # Set bracketed maximum triplet for Brent's method, along with corresponding energy values
   bracket_low = ang_max - 10
   bracket_high = ang_max + 10
   eng_max = returnEng(ang_max, row, cav_end)
   eng_low = returnEng(bracket_low, row, cav_end)
   eng_high = returnEng(bracket_high, row, cav_end)

   # Global max energy angle is determined by searching using Brent's method.
   global_max = brents(bracket_low, ang_max, bracket_high, eng_low, eng_max, eng_high, tol, it_max, row, cav_end)
   return global_max

#-----------------------------------------------------------

# Objective function; returns the energy for a specified phase angle
def returnEng(ang, row, cav_end):

   # confines angle to 0-360 degrees
   ang = angMod(ang)

   # We want fort.51 to contain the simulation information up to the previous beam element--this makes it so that 
   # we don't have to re-do that part of the simualtion. However, because of the way IMPACT-Z runs, every time it
   # encounters a -7 line, it outputs data to fort.51. This is a problem because when we are testing multiple phase
   # angles for a single beam element, the fort.51 file that we are using to restart will contain information from
   # the last simulation run--the previous angle for the same element--but we want the simulation data for the
   # optimized phase of the previous element in fort.51. So, once we get the optimized angle, we copy this data into
   # temp_fort.51, and rewrite this data into fort.51 before each new run.
   file = open('test.in','r')
   lines = file.readlines()
   file.close()
   # if RstartFlag = 1 (if this isn't the first beam element)
   if int(string.split(lines[row_restart])[1]) == 1:
      os.system('cp temp_fort.51 fort.51')

   # outputs ang data to a file
   os.system(' echo '+str(ang)+' >> ang_data_'+str(row))

   # open test.in and change phase angle
   file = open('test.in','r')
   lines = file.readlines()
   file.close()
   new_line = string.split(lines[row])
   new_line[col_ang] = str(ang)
   new_line = string.join(new_line)+'\n'
   lines[row] = new_line
   file = open('test.in', 'w')
   file.writelines(lines)
   file.close()

   # run ImpactT simulation and output eng data to a file
   os.system('./IMPACTexe > a')
   # The following line is included so that we can easily access a fort.18 file that contains the data for the
   # maximum energy. This will be useful when we make our predictions for the optimal angle of the next beam element.
   # We need to create fort.18_max because if the user specifies a non-zero angle from the zero phase, the latest fort.18
   # file will be for that angle, and not the max energy.
   os.system('cp fort.18 fort.18_max')
   os.system('tail -1 fort.18 >> tmpeng')
   os.system('echo '+str(ang)+' >> tmpphase')
   os.system('paste tmpphase tmpeng > engout')
   os.system('rm tmpeng tmpphase')
   file = open('engout', 'r')
   line = file.readline()
   file.close
   eng = string.split(line)[col_eng]
   os.system(' echo ' + str(eng) + ' >> eng_data_' + str(row))
   end_pos = float(string.split(line)[1])

   #print end_pos,cav_end
   # if eng is NaN or if final z position is not within 2*dz (dz = cdt; dt = 4*10**-12) of end of cavity, set eng to 0
   if math.isnan(float(eng)) or math.fabs(cav_end - end_pos) > 2*1.2*10**-3:
      eng = 0.0

   print ang
   print eng
   print '----------'

   return float(eng)
   
#------------------------------------------------------

# Brent's method of parabolic interpolation finds maximum given a bracketed triplet
def brents(ax, bx, cx, fa, fb, fc, tol, it_max, row, cav_end):
    
    print "Brent's method to find maximum with fractional tolerance " + str(tol) + "..."
    cgold = 0.3819660
    zeps = 0.0000000001

    # a and b must be in ascending order
    a = min(ax, cx)
    b = max(ax, cx)

    v = bx
    w = v
    x = v
    e = 0
    fx = returnEng(x, row, cav_end)
    fv = fx
    fw = fx

    # main function loop
    for i in range(it_max):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2. * tol1

        # check if done
        if (abs(x - xm) <= (tol2 - 0.5 * (b - a))):
           break # to line 309

        if (abs(e) > tol1):
           # construct a parabolic fit
           r = (x - w) * (fx - fv)
           q = (x - v) * (fx - fw)
           p = (x - v) * q - (x - w) * r
           q = 2. * (q - r)
           if (q > 0):
              p = -p
           q = abs(q)
           etemp = e
           
           # This line is in the example code for the Brent's method in 'Numerical Recipes in Fortran', but I'm not
           # sure what purpose it serves. It seems invalid since there is no previous declaration of 'd' in the first
           # iteration of this for loop.
           #e = d  

           # This is a clumsy workaround in trying to adapt code from 'Numerical Recipes in Fortran', in which
           # they use 'goto' statements, which don't exist in Python. I was able to recreate the basic structure
           # without too much repetitious code by having a single iteration 'for loop' that I could 'break' out of.
           for i in range(1):
              
              # condition to determine acceptability of parabolic fit
              if (not ((abs(p) >= abs(0.5 * q * etemp)) or (p <= q * (a - x)) or (p >= q * (b - x)))):
                 # parabolic fit works
                 d = p / q
                 u = x + d
                 
                 if (((u - a) < tol2) or ((b - u) < tol2)):
                    # skip over golden section step
                    break
                 
                 else:
                    # do golden section step
                    if (x >= xm):
                       e = a - x      
                    else:
                       e = b - x
                       d = cgold * e
                       break

              else:
                 # parabolic fit does not work
                 # do golden section step
                 if (x >= xm): 
                    e = a - x
                 else:
                    e = b - x
                 d = cgold * e

                 break

           # arrive here with d computed from parabolic fit OR golden section step
           if (abs(d) >= tol1):
              u = x + d
           else:
              u = x + math.copysign(tol1,d)

           # objective function evaluation
           fu = returnEng(u, row, cav_end)
           if (fu >= fx):
              if(u >= x):
                 a = x
              else:
                 b = x
              v = w
              fv = fw
              w = x
              fw = fx
              x = u
              fx = fu
           else:
              if (u < x):
                 a = u
              else:
                 b = u
              if ((fu >= fw) or (w == x)):
                 v = w
                 fv = fw
                 w = u
                 fw = fu
              elif ((fu >= fv) or (v == x) or (v == w)):
                 v = u
                 fv = fu

        # if the condition in line 187 is false, we still need to do the golden section step and everything
        # else that follows
        else:
           if (x >= xm):  
              e = a - x
           else:
              e = b - x
              d = cgold * e
              if (abs(d) >= tol1):
                 u = x + d
              else:
                 u = x + math.copysign(tol1, d)
                 fu = returnEng(u, row, cav_end)
                 if (fu >= fx):
                    if (u >= e):
                       a = x
                    else:
                       b = x
                       v = w
                       fv = fw
                       w = x
                       fw = fx
                       x = u
                       fx = fu
                 else:
                    if (u < x):  
                       a = u
                    else:    
                       b = u
                       if ((fu >= fw) or (w == x)):
                          v = w
                          fv = fw
                          w = u
                          fw = fu
                       elif ((fu >= fv) or (v == x) or (v == w)):
                          v = u
                          fv = fu

    # break in line 185 also brings us here
    xmax = x
    brent = fx

    print 'Zero phase (maximum)'
    print 'angle:  '+str(xmax)
    print 'energy: '+str(brent)

    return xmax

#-----------------------------------------------------------

main()
