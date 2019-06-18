#!/usr/bin/python3
import os,  re, sys,  pdb, math, fileinput
from glob import glob
from shutil import copyfile
#
#--- (c) Copyright, 2014 by the Regents of the University of California.
#This PreProcessing.py is based on the PhaseOpt.py beta v. 1.1: 
#Authors: Zhicong Liu, T. Wilson,  and J. Qiang.
#Description: 
#             Find the RF driven phase used in the ImpactT.in given the
#             design phase initially specified in the ImpactT.in.
#             The initial ImpactT.in is copied into file ImpactT_backup.in.
#             The found RF driven phases are put into the ImpactT.in file.
#             <NOTE>: Here, we assume that there is NO OVERLAP of RF cavities
#             and the executable file of the IMPACT-T code is ImpactT.
#------------------------------------------------------
# specify initial parameters for optimization
#------------------------------------------------------
# starting and ending row location of lattice in ImpactT.in


# row location of restart switch in the ImpactT.in
row_restart = 7
#------------------------------------------------------
_IMPACT_T_NAME ='ImpactTexe'
# phase angle scan specifications
angle_in = 0
angle_end = 360
delv = 20
#
# column indices to read input file ImpactT.in
# 3: bpm; 7: ang
col_el = 3
col_freq = 6
col_ang = 7

# column indices to read output file
# 4: eng
col_eng = 4

# fractional precision of phase angle to be determined
tol = 0.0001

# max number of iterations in Brent's Method before quitting
it_max = 20

    
def process(parent):
    # remove previous data files
    files = glob(os.path.join('.' + "/ang_*"))
    files += glob(os.path.join('.' + "/eng_*"))
    files += glob(os.path.join('.' + "/data_*"))
    for f in files:
        os.remove(f)
    
    # create backup of ImpactT.in
    copyfile('ImpactT.in','ImpactT_backup.in')
    
    # from row_start to row_end, optimize beam elements one at a time
    row_start = 14
    row_end = len(open('ImpactT.in').readlines())

    # default values for variables
    zedge = -1
    max_ang = -1
    final_ang = -1
    length = -1
    freq = -1
    prediction = -1
    
    global angle_in
    global angle_end
    global delv
    global col_el
    global col_freq
    global col_ang
    global tol
    global it_max
    angle_in = 0
    angle_end = 360
    delv = 20
    col_el = 3
    col_freq = 6
    col_ang = 7
    tol = 0.0001
    it_max = 20


    # global values from initial parameter specification
    # Angle_in and angle_end will be global variables to store the start and end points of the initial scan.
    # However, for later scans, they are changed if we can make a good prediction of the angle
    # (if beta is close enough to 1). When this happens, to not lose the initial values angle_in and angle_end,
    # we store them in angle_inVal and angle_endVal.
    angle_inVal     = angle_in
    angle_endVal    = angle_end

    # scan through the rows of beam elements in ImpactT.in
    # row_start and row_end are acceptable here because the -3 and -99 lines are added and removed WITHIN each
    # iteration of the loop.
    for i in range(row_start, row_end):
        # set angle_in and angle_end back to default values stored in angle_inVal and angle_endVal
        # This is reset for every row, just in case our prediction is bad (cannot find a max within +/- 20 degrees
        # of the prediction), so we are forced to scan the entire default range to find a max. 
        angle_in  = angle_inVal
        angle_end = angle_endVal

        # read from ImpactT.in
        file = open('ImpactT.in', 'r')
        lines = file.readlines()
        file.close()
        
        print(lines[i].split(),col_el,col_freq)


        # if bpm > 100 and rf field != 0
        if int(lines[i].split()[col_el]) > 100 and float((lines[i].split()[col_freq]).replace('d','E')) !=  0:

            # read in the user specified angle relative to zero phase and frequency
            spec_ang = (lines[i].split())[col_ang]
            prev_freq = freq
            freq = float(lines[i].split()[col_freq].replace('d','E'))

            print('finding max for beam element in line '+str(i)+'...')

            # insert -99 and -3 line after line i
            # 'x' is the empty  that will be the version of ImpactT.in. We insert -3 and -99 lines
            # after line i.
            x = ''
            # makes sure that we only insert one set of -3 and -99 lines (this could have happened if we have
            # duplicate beam elements)
            inserted = False

            for line in fileinput.input('ImpactT.in'):
                if line == lines[i] and inserted == False:
                    # makes sure that we only insert one set of -3 and -99 lines
                    # more than one can be inserted if we have duplicate beam elements
                    inserted = True
                    
                    # set position parameters for -3 and -99 line
                    length_str = (lines[i].split())[0]
                    prev_length = length
                    length = float(length_str.replace('d','E'))
                    
                    zedge_str = lines[i].split()[4]
                    prev_zedge = zedge
                    zedge = float(zedge_str.replace('d','E'))
                    
                    cav_end = zedge + length
                    # create -3 and -99 lines
                    tmp3line = '0 0 80 -3 '+str(cav_end)+' '+str(cav_end)+' '+str(cav_end)+' / \n'
                    tmp99line = '0 0 0 -99 '+str(cav_end)+' '+str(cav_end)+' '+str(cav_end)+' / \n'
                    x = x+line
                    x = x+tmp3line+tmp99line
                    
                    print(zedge,length,cav_end)

                else:
                    x = x+line


            # overwrite ImpactT.in with -3 and -99 lines
            file = open('ImpactT.in', 'w')
            file.write(x)
            file.close()

            # predict next phase angle if RstartFlag = 1 (if this isn't the first beam element)
            if int((lines[row_restart].split())[1]) == 1:
               prediction = predict(zedge, prev_zedge, max_ang, prev_freq, prev_length)
            
            # determine phase angle for which energy is a maximum
            max_ang = findMax(i, cav_end, prediction)
            
            # write user specified angle relative to zero phase (max angle) into ImpactT.in
            file = open('ImpactT.in','r')
            lines = file.readlines()
            file.close()
            # define user specified angle
            new_line = lines[i].split()
            final_ang = max_ang+float(spec_ang)
            # confine angle to 0-360 degrees
            final_ang = angMod(final_ang)
            # write user specified angle into ImpactT.in
            new_line[col_ang] = str(final_ang)
            new_line = ' '.join(new_line)+'\n'
            lines[i] = new_line
            file = open('ImpactT.in','w')
            file.writelines(lines)
            file.close()
            
            file = open('ImpactT.in','r')
            lines = file.readlines()
            file.close()
            # if RstartFlag = 1 (if this isn't the first beam element)
            if int(lines[row_restart].split()[1]) == 1:
               # copies previous fort.80 file to be used to restart (more information on this in the returnEng() function)
                copyfile('temp_fort.80','fort.80')

            # return energy of specified angle (this makes the last line in 'engout_spec' have the data we want)
            # this keeps 'engout''s last line as still being for the max energy
            # run ImpactT simulation and output eng data to a file
            ImpactExe = os.path.join(sys.path[0],'src',_IMPACT_T_NAME)
            os.system(ImpactExe+' > a')
            #os.system('tail -1 fort.18 >> tmpeng')
            tailAppend('fort.18','tmpeng')

            os.system('echo '+str(final_ang)+' >> tmpphase')

            #os.system('paste tmpphase tmpeng > engout_spec')
            pasteL('tmpphase','tmpeng','engout_spec')

            os.remove('tmpeng')
            os.remove('tmpphase')
            file = open('engout_spec', 'r')
            line = file.readline()
            file.close
            eng = float(line.split()[col_eng].replace('d','E'))
            os.system(' echo ' + str(eng) + ' >> eng_data_' + str(i))
            end_pos = float(line.split()[2])
            print('Optimized values')
            print('angle: '+str(final_ang))
            print('energy: '+str(eng))
            print('----------')

            # if we just optimized for the first element, set Rstartflag to 1
            # Rstartflag set to 1 makes the simulation restart from the last paused position; default set to 0.
            file = open('ImpactT.in', 'r')
            lines = file.readlines()
            file.close()
            if int(lines[row_restart].split()[1]) == 0:
               new_line = lines[row_restart].split()
               new_line[1] = str(1)
               new_line = ' '.join(new_line)+'\n'
               lines[row_restart] = new_line
               file = open('ImpactT.in', 'w')
               file.writelines(lines)
               file.close

            # remove -3 and -99 line from ImpactT.in
            file = open('ImpactT.in', 'r')
            lines = file.readlines()
            file.close()
            x = ''
            for line in fileinput.input('ImpactT.in'):
               if line == tmp3line or line == tmp99line:
                  pass
               else:
                  x = x+line
            file = open('ImpactT.in', 'w')
            file.write(x)
            file.close()

            # copy fort.80 into a temporary file temp_fort.80.
            # We will need this data for the next beam element because we restart simulation right after
            # current beam element.
            copyfile('fort.80','temp_fort.80')
 

   # After looping over all lines in ImpactT.in, we're done.
   # Set Rstartflag back to 0 for next run.
    file = open('ImpactT.in', 'r')
    lines = file.readlines()
    file.close()
    new_line = lines[row_restart].split()
    new_line[1] = str(0)
    new_line = ' '.join(new_line)+'\n'
    lines[row_restart] = new_line
    file = open('ImpactT.in', 'w')
    file.writelines(lines)
    file.close()

# returns the number of lines in a file
def fileLength(filename):
   with open(filename) as f:
      i = 0
      for i, l in enumerate(f):
         pass
   return i + 1

# confines angles to a range of 0-360 degrees; not necessary but makes things cleaner
def angMod(angle):
   while angle > 360:
      angle -= 360
   while angle < 0:
      angle += 360
   return angle

def predict(zedge, prev_zedge, max_ang, prev_freq, prev_length):

   # Set default diff1 and diff2 to 1. diff1 and diff2 are difference values that tell us how close
   # we are to an edge of a certain cavity (see the 'for loops' later). We desire the smallest diff1
   # and diff2 we can get, so setting diff1 and diff2 to 1 by default is fine, as 1 is a suitably
   # 'large' difference.
   diff1 = 1
   diff2 = 1
   
   file = open('ImpactT.in', 'r')
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
   #t1=0.0
   #t2=0.0
   for ind in range(len(lines18)):
      # Check to see if the position of the beam in fort.18_max in a certain line is within 0.01
      # of the position of the leading edge of the previous cavity.
      print(float((lines18[ind].split())[1]),prev_zedge)
      if math.fabs(float((lines18[ind].split())[1]) - prev_zedge) < 0.01:
         # Save this difference; it will be useful for comparison, since we want the time value
         # for a position as close as we can get to the leading edge of the previous cavity
         diff_new1 = math.fabs(float(lines18[ind].split()[1]) - prev_zedge)
         # Save the time value
         t1 = float(lines18[ind].split()[0])
         # If there is a line that provides a more accurate estimate, use it:
         if diff_new1 < diff1:
            diff1 = diff_new1
            t1 = float(lines18[ind].split()[0])
      # This section of code is identical to the last, except we are findig the time at the back edge
      # instead of the leading edge of the previous cavity.
      if math.fabs(float(lines18[ind].split()[1]) - (prev_zedge + prev_length)) < 0.01:
         diff_new2 = math.fabs(float(lines18[ind].split()[1]) - (prev_zedge + prev_length))
         t2 = float(lines18[ind].split()[0])
         if diff_new2 < diff2:
            diff2 = diff_new2
            t2 = float(lines18[ind].split()[0])
   
   # delt is the time elapsed between the beam entering and exiting the previous cavity
   print(float((lines18[ind].split())[1]),prev_zedge,prev_length)
   try:
       delt = t2 - t1
   except:
       print("check code!!Especially the space!!!")
       sys.exit()

   # Find out how fast the beam was traveling at the end of the previous cavity
   file = open('engout', 'r')
   line = file.readline()
   file.close()
   beta = float(line.split()[5])
   # 0.9975 is an empirical result based on accuracy of predictions for a range of beta values.
   # For beta values > 0.9975, the predictions were accurate to within about 20 degrees.
   if beta > 0.9975:
      print('Estimating angle range:')
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
   print('Initial scan from '+str(ang_in)+' to '+str(ang_end)+' by '+str(delv)+'...')
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
      print('Try scan from 0 to 360 by 10...')
      for i in range(0, 360+1, 10):
         eng = returnEng(i, row, cav_end)
         print('i: ',i,eng)
         if (eng > eng_max):
            eng_max = eng
            ang_max = i         

   print('initial max angle: '+str(ang_max))
   print('initial max energy: '+str(eng_max))
   
   # Set bracketed maximum triplet for Brent's method, along with corresponding energy values
   bracket_low = ang_max - 10
   bracket_high = ang_max + 10
   eng_max = returnEng(ang_max, row, cav_end)
   eng_low = returnEng(bracket_low, row, cav_end)
   eng_high = returnEng(bracket_high, row, cav_end)

   # Global max energy angle is determined by searching using Brent's method.
   global_max = brents(bracket_low, ang_max, bracket_high, eng_low, eng_max, eng_high, tol, it_max, row, cav_end)
   return global_max

# Objective function; returns the energy for a specified phase angle
def returnEng(ang, row, cav_end):

   # confines angle to 0-360 degrees
   ang = angMod(ang)

   # We want fort.80 to contain the simulation information up to the previous beam element--this makes it so that 
   # we don't have to re-do that part of the simualtion. However, because of the way IMPACT-T runs, every time it
   # encounters a -3 line, it outputs data to fort.80. This is a problem because when we are testing multiple phase
   # angles for a single beam element, the fort.80 file that we are using to restart will contain information from
   # the last simulation run--the previous angle for the same element--but we want the simulation data for the
   # optimized phase of the previous element in fort.80. So, once we get the optimized angle, we copy this data into
   # temp_fort.80, and rewrite this data into fort.80 before each new run.
   file = open('ImpactT.in','r')
   lines = file.readlines()
   file.close()
   # if RstartFlag = 1 (if this isn't the first beam element)
   if int(lines[row_restart].split()[1]) == 1:
      copyfile('temp_fort.80','fort.80')

   # outputs ang data to a file
   os.system('echo '+str(ang)+' >> ang_data_'+str(row))

   # open ImpactT.in and change phase angle
   file = open('ImpactT.in','r')
   lines = file.readlines()
   file.close()
   new_line = lines[row].split()
   new_line[col_ang] = str(ang)
   new_line = ' '.join(new_line)+'\n'
   lines[row] = new_line
   file = open('ImpactT.in', 'w')
   file.writelines(lines)
   file.close()

   # run ImpactT simulation and output eng data to a file
   ImpactExe = os.path.join(sys.path[0],'src',_IMPACT_T_NAME)
   os.system(ImpactExe+' > a')

   # The following line is included so that we can easily access a fort.18 file that contains the data for the
   # maximum energy. This will be useful when we make our predictions for the optimal angle of the next beam element.
   # We need to create fort.18_max because if the user specifies a non-zero angle from the zero phase, the latest fort.18
   # file will be for that angle, and not the max energy.

   copyfile('fort.18','fort.18_max')
   
   #os.system('tail -1 fort.18 >> tmpeng')
   tailAppend('fort.18','tmpeng')
   
   os.system('echo '+str(ang)+' >> tmpphase')
   
   
   #os.system('paste tmpphase tmpeng > engout')
   pasteL('tmpphase','tmpeng','engout')
   os.remove('tmpeng')
   os.remove('tmpphase')
   file = open('engout', 'r')
   line = file.readline()
   file.close
   eng = line.split()[col_eng]
   os.system(' echo ' + str(eng) + ' >> eng_data_' + str(row))
   end_pos = float(line.split()[2])

   # if eng is NaN or if final z position is not within 2*dz (dz = cdt; dt = 4*10**-12) of end of cavity, set eng to 0
   if math.isnan(float(eng)) or math.fabs(cav_end - end_pos) > 4*1.2*10**-3:
      eng = 0.0

   print(ang,cav_end,end_pos)
   print(eng)
   print('----------')

   return float(eng)
   
# Brent's method of parabolic interpolation finds maximum given a bracketed triplet
def brents(ax, bx, cx, fa, fb, fc, tol, it_max, row, cav_end):
    
    print("Brent's method to find maximum with fractional tolerance " + str(tol) + "...")
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

    print('Zero phase (maximum)')
    print('angle:  '+str(xmax))
    print('energy: '+str(brent))

    return xmax

def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))

def tailAppend(src,dest):
    f1 = open(src,'r')
    li = f1.readlines()
    f2 = open(dest,'a')
    f2.write(li[-1])
    f1.close()
    f2.close()
    
def pasteL(src1,src2,dest):
    f1 = open(src1,'r')
    f2 = open(src2,'r')
    l1 = f1.readlines()
    l2 = f2.readlines()

    fdest = open(dest,'w')
    minline = min(len(l1),len(l2))
    maxline = max(len(l1),len(l2))
    for i in range(0,minline):
        fdest.write(l1[i].rstrip()+l2[i])
    for i in range(minline,maxline):
        try:
            fdest.write(l1[i])
        except:
            fdest.write(l2[i])
    f1.close()
    f2.close()
    fdest.close()
