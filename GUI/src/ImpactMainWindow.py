#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#GUI for Impact code V 1.0
#by Zhicong Liu and Ji Qiang, 2018

import os, sys, threading, subprocess
import math
from shutil import copyfile

if sys.version_info[0] < 3:
    print("Error: need python version 3.x!")
    exit(0)

import tkinter as tk
from tkinter import messagebox,filedialog,ttk

import numpy as np

import ImpactFile
import PreProcessing
import ConvertFunc
import LatticeFrame
import ImpactTSet
import ImpactTPlot
import ImpactZSet
import ImpactZPlot

_MPINAME   ='mpirun'
_IMPACT_T_NAME = 'ImpactT executable not defined'
_IMPACT_Z_NAME = 'ImpactZ executable not defined'
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    _IMPACT_Z_NAME = 'ImpactZ.serial.linux'
# elif platform == "darwin":
    # OS X
elif platform == "win32":
    # Windows...
    _IMPACT_T_NAME = 'ImpactTexe'
    _IMPACT_Z_NAME = 'ImpactZexe'

_width = 560
_height= 633


PARTICLE_TYPE       = {'Electron'   :'510998.9461 -1.0',
                       'Proton'     :'938272081.3 1.0',
                       'Positron'   :'510998.9461 1.0',
                       'Antiproton' :'938272081.3 -1.0',
                       'Other'      :'Other_NONE'}

DISTRIBUTION_T_TYPE = {'Uniform'   :'1',
                     'Gauss'     :'2',
                     'WaterBag'  :'3',
                     'SemiGauss' :'4',
                     'KV'        :'5',
                     'Read'      :'16',
##                     'Read Parmela'   :'24',
##                     'Read Elegant'   :'25',
                     'CylcoldZSob'    :'27'}

DISTRIBUTION_Z_TYPE = {'Uniform'   :'1',
                     'Gauss'     :'2',
                     'WaterBag'  :'3',
                     'SemiGauss' :'4',
                     'KV'        :'5',
                     'Read'      :'19',
                     'Multi-charge-state WaterBag'   :'16',
                     'Multi-charge-state Gaussian'   :'17'}

DIAGNOSTIC_TYPE   = {'At given time'    :'1',
                     'At bunch centroid':'2',
                     'No output'        :'3'}

OUTPUT_Z_TYPE        = {'Standard'          :'1',
                        '95% Emittance'     :'2'}

BOUNDARY_TYPE        = {'Trans:open,  Longi:open'     :'1',
                        'Trans:open,  Longi:period'   :'2',
                        'Trans-Round, Longi-open'     :'3',
                        'Trans-Round, Longi-perod'    :'4',
                        'Trans:Rect,  Longi-open'     :'5',
                        'Trans-Rect,  Longi-perod'    :'6'}
INTEGRATOR_TYPE      = {'Linear'    :'1',
                        'Non Linear':'2'}

class ImpactMainWindow(tk.Tk):
    LABEL_TEXT =[
        "This is a beta version of the IMPACT user interface..."
        ]

    PLOTTYPE = {'Centriod location' :2,
                'Rms size'          :3,
                'Centriod momentum' :4,
                'Rms momentum'      :5,
                'Twiss'             :6,
                'Emittance'         :7}

    count = 0
    AccKernel = 'ImpactT'
    ImpactThread=threading.Thread()
    ImpactThread.setDaemon(True)
    def __init__(self, master=None):  
        tk.Tk.__init__(self, master)
        self.title("Impact")
        self.createWidgets(master)
        #self.master.iconbitmap()
          
        for item in ImpactMainWindow.LABEL_TEXT:
            print(item)
            
    def createWidgets(self, master):
        self.t= startWindow(self)
        w1  = 400
        h1  = 300
        ws1 = self.t.winfo_screenwidth()  # width of the screen
        hs1 = self.t.winfo_screenheight() # height of the screen
        x1 = (ws1/2) - (w1/2)
        y1 = (hs1/2) - (h1/2)
        self.t.overrideredirect(0)
        self.t.geometry('%dx%d+%d+%d' % (w1, h1, x1, y1))
        self.frame_left = tk.Frame(self, height =_height, width = _width)
        self.frame_left.grid(row=0,column=0)
        
        self.frame_logo = tk.Frame(self.frame_left, height =_height/10, width = _width)
        self.frame_logo.pack(side = 'top')
        #print(resource_path(os.path.join('icon','ImpactT.gif')))
        try:
            self.logo_ImpactT = tk.PhotoImage(file = resource_path(os.path.join('icon','ImpactT.gif')), format="gif")
            self.logo_ImpactZ = tk.PhotoImage(file = resource_path(os.path.join('icon','ImpactZ.gif')), format="gif")
            self.label_logo  = tk.Label(self.frame_logo,image = self.logo_ImpactT)
            self.label_logo.pack(fill = 'y',side = 'left')
        except:
            self.label_logo  = tk.Label(self.frame_logo,bitmap='error')
            self.label_logo.pack(fill = 'y',side = 'left')
        
        '''
        LARGE_FONT= ("Helvetica", 20,'italic')
        self.button_switch = tk.Button(self.frame_logo)
        self.button_switch["text"]        = "Switch Kernel"
        #self.run["font"]        = LARGE_FONT
        self.button_switch["command"]     = lambda: self.switch()
        self.button_switch.pack(fill = 'both',expand =1,side = 'top',padx = 150)
        '''
        """Frame1: CPU GPU information"""
        self.frame_input1 = tk.LabelFrame(self.frame_logo, height =_height/10, width = _width, 
                                                text="Configuration")
        #self.frame_input1.grid(row=1,column=0,columnspan=2,sticky='W')
        self.frame_input1.pack(side = 'right')
        
        self.frame_CPU = tk.Frame(self.frame_input1, height =_height/10, width = _width)
        self.frame_CPU.pack(side = 'top')
        
        vcmd = (self.register(self.validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.label_noc1 = tk.Label(self.frame_CPU, text="# of cores at Y",pady=1)
        self.entry_noc1 = tk.Entry(self.frame_CPU, 
                                   validate = 'key',
                                   validatecommand = vcmd,
                                   width=4)
        self.entry_noc1.insert(0, '1')
        self.label_noc1.pack(side='left')
        self.entry_noc1.pack(side='left')

        self.label_noc2 = tk.Label(self.frame_CPU, text="# of cores at Z",pady=1)
        self.entry_noc2 = tk.Entry(self.frame_CPU, 
                                   validate = 'key',
                                   validatecommand = vcmd,
                                   width=4)
        self.entry_noc2.insert(0, '1')
        self.label_noc2.pack(side='left')
        self.entry_noc2.pack(side='left')
        
        self.GPUflag   = tk.IntVar()
        self.check_GPU  = tk.Checkbutton(self.frame_CPU, text="GPU", variable=self.GPUflag)
        self.check_GPU.pack(side='left')
        '''
        self.label_dic = tk.Label(self.frame_input1, text="Work Dictionary",pady=1)
        self.label_dic.pack(side='left')
        '''
        self.entry_dic = tk.Entry(self.frame_input1, width=45)
        #self.entry_dic.insert(0, os.path.dirname(os.path.abspath(__file__)))
        self.entry_dic.insert(0, sys.path[0])
        self.entry_dic.pack(side='left')
        
        self.button_dic = tk.Button(self.frame_input1)
        self.button_dic["text"]        = "..."
        self.button_dic["command"]     = lambda: self.changeDic()
        self.button_dic.pack(side = 'left')
        
        self.button_dic.bind("<Enter>", lambda event, h=self.button_dic: h.configure(bg="#00CD00"))
        self.button_dic.bind("<Leave>", lambda event, h=self.button_dic: h.configure(bg="#FFFFFF"))
        
        """Frame2: Time step"""
        self.frame1 = tk.Frame(self.frame_left, 
                               height =_height/2, width = _width)
        #self.frame_input2.grid(row=2,column=0,sticky='W')
        self.frame1.pack(fill = 'x',side = 'top')
        
        self.frame_input2 = tk.LabelFrame(self.frame1, 
                                          height =_height/2, width = _width, 
                                          text="Numerical parameters")
        #self.frame_input2.grid(row=2,column=0,sticky='W')
        self.frame_input2.pack(fill = 'x',side = 'left')
        rowtemp=1
        frame2width = 7
        
        """dt"""
        self.label_dt    = tk.Label(self.frame_input2, text="Timestep")
        self.entry_dt    = tk.Entry(self.frame_input2, 
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_dt.insert(0, '1.0e-9')
        self.label_dt.grid(row=rowtemp,sticky='W')
        self.entry_dt.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1

        """Nstep"""
        self.label_Nstep = tk.Label(self.frame_input2, text="# of step")
        self.entry_Nstep = tk.Entry(self.frame_input2,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Nstep.insert(0, '100')
        self.label_Nstep.grid(row=rowtemp,sticky='W')
        self.entry_Nstep.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1

        """Np"""
        self.label_Np    = tk.Label(self.frame_input2, text="Particle number")
        self.entry_Np    = tk.Entry(self.frame_input2,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Np.insert(0, '1.6e5')
        self.label_Np.grid(row=rowtemp,sticky='W')
        self.entry_Np.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Grid X"""
        self.label_Ngx   = tk.Label(self.frame_input2, text="Grid number X")
        self.entry_Ngx   = tk.Entry(self.frame_input2,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Ngx.insert(0, '64')
        self.label_Ngx.grid(row=rowtemp,sticky='W')
        self.entry_Ngx.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Grid Y"""
        self.label_Ngy   = tk.Label(self.frame_input2, text="Grid number Y")
        self.entry_Ngy   = tk.Entry(self.frame_input2,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Ngy.insert(0, '64')
        self.label_Ngy.grid(row=rowtemp,sticky='W')
        self.entry_Ngy.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Grid Z"""
        self.label_Ngz   = tk.Label(self.frame_input2, text="Grid number Z")
        self.entry_Ngz   = tk.Entry(self.frame_input2,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Ngz.insert(0, '64')
        self.label_Ngz.grid(row=rowtemp,sticky='W')
        self.entry_Ngz.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1

        
        """Twiss"""
        twisswidth = [9,11,11]
        self.frame_Twiss = tk.LabelFrame(self.frame1, 
                                         height =_height/2, width = _width, text="Beam Twiss")
        self.frame_Twiss.pack(fill = 'y',side = 'left')
            
        self.twiss_t = []
        self.twiss_Tchara = ["alpha","beta","emittance"]
        for column in range(3):
            label = tk.Label(self.frame_Twiss, text=self.twiss_Tchara[column], 
                            borderwidth=0,
                            width=twisswidth[column])
            label.grid(row=0, column=column+1, sticky="ns", padx=1, pady=1)
            self.twiss_t.append(label)
        
        self.twiss_x = []
        self.twiss_xchara = ["X","Y","Z"]
        for row in range(3):
            label = tk.Label(self.frame_Twiss, text=self.twiss_xchara[row], 
                                 borderwidth=0, width=2)
            label.grid(row=row+1, column=0, sticky="ns", padx=1, pady=1)
            self.twiss_x.append(label)
        
        '''Twiss and Sigma string'''
        self.twiss_param = [["3.84562e-4","0.001" ,"0.0", "1.0","1.0","0.0","0.0"],
                            ["3.84562e-4","0.001" ,"0.0", "1.0","1.0","0.0","0.0"],
                            ["1.18618e-6","0.0014","0.0", "1.0","1.0","0.5","0.002"]]
        self.string_sigma = []
        for row in range(3):
            current_row = []
            for column in range(7):
                t = tk.StringVar(value=self.twiss_param[row][column])
                current_row.append(t)
            self.string_sigma.append(current_row)
        
        self.string_twiss = []
        for row in range(3):
            current_row = []
            for column in range(3):
                t = tk.StringVar(value='0')
                current_row.append(t)
            self.string_twiss.append(current_row)

        '''Twiss Entry'''   
        self.entry_twiss = []
        for row in range(3):
            current_row = []
            for column in range(3):
                label = tk.Entry(self.frame_Twiss, 
                                 textvariable=self.string_twiss[row][column], 
                                 borderwidth=0,
                                 width=twisswidth[column])
                label.grid(row=row+1, column=column+1, sticky="ns", padx=2, pady=1)
                current_row.append(label)
            self.entry_twiss.append(current_row)
        
        """Distribution"""
        self.frame_Dist = tk.Frame(self.frame_Twiss,
                                   height =_height/2, width = _width)
        self.frame_Dist.grid(row=4, column=0, rowspan=1, columnspan=4,sticky="nsew", padx=1, pady=1)
        
        self.label_dist    = tk.Label(self.frame_Dist, text="Distribution")
        self.label_dist.grid(row=0, column=0,
                             sticky="nw", padx=1, pady=1)
        self.distTypeComx = tk.StringVar(self.frame_Dist,'WaterBag')
        self.distTypeNumb = tk.StringVar(self.frame_Dist,'3')
        self.distType = ttk.Combobox(self.frame_Dist,text=self.distTypeComx,
                                     values=['Uniform', 'Gauss', 'SemiGauss',
                                             'WaterBag','KV', 'Read',
                                             'CylcoldZSob','Other'])
        self.distType.grid(row=0, column=1,
                             sticky="nsew", padx=1, pady=1)
        
        """Particle Type"""
        self.label_ptcType    = tk.Label(self.frame_Dist, text="Particle")
        self.label_ptcType.grid(row=1, column=0,
                                sticky="nw", padx=1, pady=1)
        self.ptcTypeComx = tk.StringVar(self.frame_Dist,'Electron')
        self.ptcMass     = tk.StringVar(self.frame_Dist,'510998.9461')
        self.ptcCharge   = tk.StringVar(self.frame_Dist,'-1.0')
        self.ptcType = ttk.Combobox(self.frame_Dist,text=self.ptcTypeComx,
                                     values=['Electron', 'Proton', 'Positron', 'Antiproton','Other'])
        self.ptcType.grid(row=1, column=1,
                          sticky="nsew", padx=1, pady=1)
        
        """Particle"""
        self.frame_Particle = tk.LabelFrame(self.frame1, 
                                            height =_height/2, width = _width, 
                                            text="Beam")
        self.frame_Particle.pack(fill = 'both',side = 'top')
        rowtemp=0
                
        """Current"""
        self.label_cur   = tk.Label(self.frame_Particle, text="Current(A)")
        self.entry_cur   = tk.Entry(self.frame_Particle,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_cur.insert(0, '0.13')
        self.label_cur.grid(row=rowtemp,sticky='W')
        self.entry_cur.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Energy"""
        self.label_Ek    = tk.Label(self.frame_Particle, text="Energy(eV)")
        self.entry_Ek    = tk.Entry(self.frame_Particle,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_Ek.insert(0, '1.0')
        self.label_Ek.grid(row=rowtemp,sticky='W')
        self.entry_Ek.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Frequency"""
        self.label_frq   = tk.Label(self.frame_Particle, text="Frequency(Hz)")
        self.entry_frq   = tk.Entry(self.frame_Particle,
                                    validate = 'key',
                                    validatecommand = vcmd,
                                    width=frame2width)
        self.entry_frq.insert(0, '1.3e9')
        self.label_frq.grid(row=rowtemp,sticky='W')
        self.entry_frq.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        '''Advanced Setting'''
        self.MPI_EXE      =tk.StringVar(value=_MPINAME)
        self.IMPACT_T_EXE =tk.StringVar(value=os.path.join(sys.path[0],'src',_IMPACT_T_NAME))
        self.IMPACT_Z_EXE =tk.StringVar(value=os.path.join(sys.path[0],'src',_IMPACT_Z_NAME))
        
        self.Nbunch      = tk.StringVar(value='1')
        self.Dim         = tk.StringVar(value='6')
        self.Flagmap     = tk.StringVar(value='Linear')
        self.Flagerr     = tk.StringVar(value='0')
        self.Flagdiag    = tk.StringVar(value='At bunch centroid')
        self.Flagimag    = tk.StringVar(value='0')
        self.Zimage      = tk.StringVar(value='0.02')
        self.Flagbc      = tk.StringVar(value='Trans:open,  Longi:open')
        
        self.Xrad        = tk.StringVar(value='0.15')
        self.Yrad        = tk.StringVar(value='0.15')
        self.Zrad        = tk.StringVar(value='1.0e5')
        
        self.FlagRestart  = tk.StringVar(value='0')
        
        self.Nemission   = tk.StringVar(value='1')
        self.Temission   = tk.StringVar(value='8.0e-11')
        self.Tinitial    = tk.StringVar(value='0.0')
        
        '''Advanced Set for ImpactZ'''
        self.FlagOutput_Z    = tk.StringVar(value='Standard')
        self.ZperiodSize     = tk.StringVar(value='0.1')
        self.FlagSubcycle    = tk.StringVar(value='0')
        self.NpList          = tk.StringVar(value='16000')
        self.CurrentList     = tk.StringVar(value='0.0')
        self.SchargeList     = tk.StringVar(value='1.06577993775845e-09')
           
        self.button_AdvancedControl = tk.Button(self.frame1,text='Advanced Setting',
                                                command = self.makeAdvancedSet)
        self.button_AdvancedControl.pack(fill='both',expand=1,padx = 15, pady=20)
        LARGE_FONT= ("Helvetica", 10)
        self.button_AdvancedControl["font"]        = LARGE_FONT
        
        self.button_AdvancedControl.bind("<Enter>", lambda event, h=self.button_AdvancedControl: h.configure(bg="#00CD00"))
        self.button_AdvancedControl.bind("<Leave>", lambda event, h=self.button_AdvancedControl: h.configure(bg="#FFFFFF"))
        
        """Lattice"""
        self.frame_input3 = tk.LabelFrame(self.frame_left, 
                                          height =_height/6, width = _width,
                                          text="Lattice")
        self.frame_input3.pack(fill="both", expand=1, side=tk.TOP)
        
        self.lattice = LatticeFrame.LatticeFrameC(self.frame_input3, height = _height/6,width = _width)
        self.lattice.pack(fill="both", expand=1, side=tk.TOP)
              
        """Console"""
        self.frame_console = tk.LabelFrame(self.frame_left, 
                                           height =_height/5, width = _width,
                                           text="Console")
        self.frame_console.pack(fill="both", expand=1, side=tk.TOP)
        
        self.console_sv = tk.Scrollbar(self.frame_console, orient=tk.VERTICAL)
        self.console = LatticeFrame.ConsoleText(self.frame_console,
                                   width = LatticeFrame._TextWidth, height=LatticeFrame._TextHeight,
                                   bg='black',fg='green', yscrollcommand=self.console_sv.set)
        self.console_sv.config(command=self.console.yview)  
        self.console_sv.pack(fill="y", expand=0, side=tk.RIGHT, anchor=tk.N)  
        self.console.pack(fill="both", expand=1, side=tk.TOP)
        self.console.start()
        
        """Control"""
        
        self.frame_control = tk.Frame(self.frame_left, height =_height/5, width = _width)
        self.frame_control.pack(fill="both", expand=1, side=tk.TOP)
        '''
        self.frame_output = tk.LabelFrame(self.frame_control,
                                          height =_height/5, width = _width,
                                          text="Plot")
        self.frame_output.pack(side='left')
        
        self.frame_plotControl = PlotControlFrame(self.frame_output,
                                             height =_height/5, width = _width)
        self.frame_plotControl.pack(side='left')
        '''
        
        '''Final button'''
        LARGE_FONT= ("Helvetica", 20,'italic')                
        self.button_initial = tk.Button(self.frame_control)
        self.button_initial["text"]        = "Pre-Process"
        self.button_initial["font"]        = LARGE_FONT
        self.button_initial["command"]     = lambda: self.thread_it(self.preprocessing)
        self.button_initial.pack(fill = 'both',expand =1,side = 'left')
        
        self.button_initial.bind("<Enter>", lambda event, h=self.button_initial: h.configure(bg="#00CD00"))#green
        self.button_initial.bind("<Leave>", lambda event, h=self.button_initial: h.configure(bg="#FFFFFF"))
        
        self.button_run = tk.Button(self.frame_control)
        self.button_run["text"]         = "Run"
        self.button_run["font"]         = LARGE_FONT
        self.button_run["command"]      = lambda: self.thread_it(self.startImpact)
        #self.button_run["command"]      = lambda: self.startImpact()
        self.button_run.pack(fill = 'both',expand =1,side = 'left')
        self.run_lock = threading.RLock()
        
        self.button_run.bind("<Enter>", lambda event, h=self.button_run: h.configure(bg="#00CD00"))#green
        self.button_run.bind("<Leave>", lambda event, h=self.button_run: h.configure(bg="#FFFFFF"))
        
        self.button_plot = tk.Button(self.frame_control)
        self.button_plot["text"]        = "Post-Process"
        self.button_plot["font"]        = LARGE_FONT
        self.button_plot["command"]     = lambda: self.makeAdvancedPlot()
        self.button_plot.pack(fill = 'both',expand =1,side = 'left')
        
        self.button_plot.bind("<Enter>", lambda event, h=self.button_plot: h.configure(bg="#00CD00"))#green
        self.button_plot.bind("<Leave>", lambda event, h=self.button_plot: h.configure(bg="#FFFFFF"))
        
        
        '''
        self.test = tk.Button(self.frame_left)
        self.test["text"] = "debug"
        self.test["command"] = lambda: self.debug()
        self.test.pack(fill = 'both',expand =1,side = 'top')
        
        self.test2 = tk.Button(self.frame_left)
        self.test2["text"] = "debug2"
        self.test2["command"] = lambda: self.debug2()
        self.test2.pack(fill = 'both',expand =1,side = 'top')
        '''
        self.updateTwissLock    =0
        
        for i in range(3):    
            for j in range(3):
                self.string_sigma[i][j].trace('w',lambda a,b,c,h=i: self.updateTwiss(h))
                self.string_twiss[i][j].trace('w',lambda a,b,c,h=i: self.updateSigma(h))
        for row in range(3):
            for column in range(3):
                self.string_sigma[row][column].set(self.twiss_param[row][column])
                pass
        self.updatePtcTypeLock  =0
        self.ptcTypeComx.trace('w', self.updatePtc)
        self.ptcMass.trace('w', self.updatePtcType)
        self.ptcCharge.trace('w', self.updatePtcType)
        
        self.updateDistTypeLock  =0
        self.distTypeComx.trace('w', self.updateDist)
        self.distTypeNumb.trace('w', self.updateDistType)

        self.update()
        self.winwidth     = self.winfo_width()
        self.winheight    = self.winfo_height()
        self.screenwidth  = self.winfo_screenwidth()  # width of the screen
        self.screenheight = self.winfo_screenheight() # height of the screen
        self.winPosX = (self.screenwidth/2) - (self.winwidth/2)
        self.winPosY = (self.screenheight/2) - (self.winheight/2)
        self.withdraw()

        self.t.lift()
        '''degue'''
        #self.t.startImpactT(self)
        #self.makeAdvancedPlot()
    
    #(Kilean)#######################################################
    def getBeam4pImpact(self):
      self.beam = {
        'nCore_y':int(self.entry_noc1.get()),
        'nCore_z':int(self.entry_noc2.get()),
        'mass': float(self.ptcMass.get()), # eV/c^2
        'energy': float(self.entry_Ek.get()), # eV
        'n_particles': int(float(self.entry_Np.get())),
        'error study flag':int(float(self.Flagerr.get())),
        'restart flag':int(float(self.FlagRestart.get())),
        'standard output':int(OUTPUT_Z_TYPE[self.FlagOutput_Z.get()]),
        'current' : float(self.entry_cur.get()), # ampere
        'x00': float(self.string_sigma[0][0].get()),
        'x11': float(self.string_sigma[0][1].get()),
        'x01': float(self.string_sigma[0][2].get()),
        'y00': float(self.string_sigma[1][0].get()),
        'y11': float(self.string_sigma[1][1].get()),
        'y01': float(self.string_sigma[1][2].get()),
        'z00': float(self.string_sigma[2][0].get()),
        'z11': float(self.string_sigma[2][1].get()),
        'z01': float(self.string_sigma[2][2].get()),
        'frequency': float(self.entry_frq.get()), # Hz
        'phase': float(self.Tinitial.get()),   #radian
        'mesh_x' : int(self.entry_Ngx.get()),
        'mesh_y' : int(self.entry_Ngy.get()), 
        'mesh_z' : int(self.entry_Ngz.get())}
      if self.distType.get() == 'Other':
        self.beam['distribution id'] = int(self.distTypeNumb.get())
      else:
        self.beam['distribution id'] = int(DISTRIBUTION_Z_TYPE[self.distType.get()])
      self.beam['charge per mass'] = float(self.ptcCharge.get())/self.beam['mass']

    def getLattice4pImpact(self):
      '''
      get lattice data for pImpact
      '''
      def str2lattice(latticeStr):
        lattice = []
        for i in range(len(latticeStr)):
            if latticeStr[i] in ['\n','\r\n',''] or latticeStr[i][0]=='!':
                continue
            elem = str2elem(latticeStr[i])
            if elem : #check if elem is not empty
                lattice.append(elem)
        return lattice
      def str2elem(elemStr): 
        elemStr = elemStr.split()
        elemID=int(float(elemStr[3]))
        if elemID == 0:
            elemtDict = {'type':'drift',
                         'length': float(elemStr[0]),
                         'n_sckick': int(elemStr[1]), 
                         'n_map': int(elemStr[2]), 
                         'radius': float(elemStr[4])
                        }   
        elif elemID == 1:
            elemtDict = {'type':'quad',
                         'length': float(elemStr[0]),
                         'n_sckick': int(elemStr[1]), 
                         'n_map': int(elemStr[2]), 
                         'B1': float(elemStr[4]), 
                         'input file id': int(elemStr[5]), 
                         'radius': float(elemStr[6]) }
            if len(elemStr)>8:
                         elemtDict['dx']=float(elemStr[7])
            if len(elemStr)>9:
                         elemtDict['dy']=float(elemStr[8]) 
        elif elemID == 4:
            elemtDict = {'type':'bend',
                         'length': float(elemStr[0]),
                         'n_sckick': int(elemStr[1]), 
                         'n_map': int(elemStr[2]), 
                         'angle': float(elemStr[4]), 
                         'k1': float(elemStr[5]), 
                         'input switch': int(float(elemStr[6])),
                         'radius': float(elemStr[7]),
                         'entrance edge': float(elemStr[8]),
                         'exit edge': float(elemStr[9]),
                         'entrance curvature': float(elemStr[10]),
                         'exit curvature': float(elemStr[11]),
                         'FINT': float(elemStr[12]),
                        }                     
        elif elemID == 104:
            elemtDict= {'type':'scrf',
                        'length': float(elemStr[0]),
                        'n_sckick': int(elemStr[1]), 
                        'n_map': int(elemStr[2]), 
                        'scale': float(elemStr[4]), 
                        'frequency': float(elemStr[5]), 
                        'phase': float(elemStr[6]), 
                        'input file id': int(elemStr[7]), 
                        'radius': float(elemStr[8])
                       }
        elif elemID == -2:
            elemtDict= {'type':'write full',
                        'file id': int(elemStr[2])}                   
        elif elemID == -7:
            elemtDict= {'type':'restart'}
        elif elemID == -99:
            elemtDict= {'type':'halt'} 
        else :
            elemtDict= {}
        return elemtDict
      latticeStr = self.lattice.getHide().splitlines()
      self.lattice = str2lattice(latticeStr)
      
    def exit(self):
      self.getBeam4pImpact()
      self.getLattice4pImpact()
      self.destroy()
      
    #######################################################(Kilean)#
    def debug(self):
        self.lattice.convertNtoW(self.lattice.get('0.0', tk.END))
    
    def debug2(self):
        self.lattice.updateHide()
        
    def changeDic(self):
        filename = filedialog.askdirectory(parent=self)
        try:
            if filename=='':
                return
            os.chdir(filename)
        except:
            print("Cannot change the folder. Check the path Please.")
        
        self.entry_dic.delete(0, 'end')
        self.entry_dic.insert(0, filename)
        print(filename)
    
    def updatePtc(self,*args):
        if self.updatePtcTypeLock ==1:
            return
        self.updatePtcTypeLock = 1
        
        if self.ptcTypeComx.get() !='Other':
            try:
                ptc = PARTICLE_TYPE[self.ptcTypeComx.get()].split()
                self.ptcMass.set(ptc[0])
                self.ptcCharge.set(ptc[1])
            except:
                pass
        self.updatePtcTypeLock = 0
    
    def updatePtcType(self,*args):
        if self.updatePtcTypeLock ==1:
            return
        self.updatePtcTypeLock = 1
        
        invermap = dict(map(lambda t:(t[1],t[0]), PARTICLE_TYPE.items()))
        ptcFound = 0
        for key in invermap.keys():
            ptc = key.split()
            try:
                if math.isclose(float(ptc[0]), float(self.ptcMass.get()),rel_tol=1e-3):
                    if math.isclose(float(ptc[1]), float(self.ptcCharge.get()),rel_tol=1e-3):
                        self.ptcTypeComx.set(invermap[key])
                        ptcFound = 1
                        break
            except:
                pass
        if ptcFound==0:
            self.ptcTypeComx.set('Other')
        self.updatePtcTypeLock = 0

    def updateDist(self,*args):
        if self.updateDistTypeLock ==1:
            return
        self.updateDistTypeLock = 1
        
        if self.distTypeComx.get() != 'Other':
            try:
                if self.AccKernel=='ImpactT':
                    dist = DISTRIBUTION_T_TYPE[self.distTypeComx.get()]
                    self.distTypeNumb.set(dist)
                elif self.AccKernel=='ImpactZ':
                    dist = DISTRIBUTION_Z_TYPE[self.distTypeComx.get()]
                    self.distTypeNumb.set(dist)
                else:
                    print("Kernel type ERROR!")
                    sys.exit()
            except:
                pass
        self.updateDistTypeLock = 0
        
    def updateDistType(self,*args):
        if len(self.distTypeNumb.get())==0:
            return
        if self.updateDistTypeLock ==1:
            return
        self.updateDistTypeLock = 1
        
        invermap = {'':''}
        if self.AccKernel=='ImpactT':
            invermap = dict(map(lambda t:(t[1],t[0]), DISTRIBUTION_T_TYPE.items()))
        elif self.AccKernel=='ImpactZ':
            invermap = dict(map(lambda t:(t[1],t[0]), DISTRIBUTION_Z_TYPE.items()))
        else:
            print("Kernel type ERROR!")
            sys.exit()

        distFound = 0
        for key in invermap.keys():
            try:
                if key == self.distTypeNumb.get():
                    self.distTypeComx.set(invermap[key])
                    distFound = 1
                    break
            except:
                pass
        if distFound==0:
            self.distTypeComx.set('Other')
        self.updateDistTypeLock = 0
        
    def updateTwiss(self,i):
        if self.updateTwissLock ==1:
            return
        self.updateTwissLock = 1
        
        try:
            s1,s2,s3 = float(self.string_sigma[i][0].get()), float(self.string_sigma[i][1].get()), float(self.string_sigma[i][2].get())
            f, m, k  = float(self.entry_frq.get()),    float(self.ptcMass.get()),      float(self.entry_Ek.get())
            re=[0]*3
            if i==1 or i==0:
                re[0],re[1],re[2] = ConvertFunc.Sigma2Twiss(s1,s2,s3, f,m,k)
            if i==2:
                re[0],re[1],re[2] = ConvertFunc.Sigma2TwissZ(s1,s2,s3, f,m,k)    
            digits = [2,4,4]
            for j in range(3):
                #if not math.isclose(re[j], float(self.string_twiss[i][j].get()),rel_tol=1e-06):
                    self.string_twiss[i][j].set('{:.{}e}'.format(re[j], digits[j]))
        except:
            #print('Error')
            pass
        self.updateTwissLock = 0
        
    def updateSigma(self,i):
        if self.updateTwissLock == 1:
            return
        self.updateTwissLock = 1

        try:
            s1,s2,s3 = float(self.string_twiss[i][0].get()), float(self.string_twiss[i][1].get()), float(self.string_twiss[i][2].get())
            f, m, k  = float(self.entry_frq.get()),    float(self.ptcMass.get()),      float(self.entry_Ek.get())
            re=[0]*3
            if i==1 or i==0:
                re[0],re[1],re[2] = ConvertFunc.Twiss2Sigma(s1,s2,s3, f,m,k)
            if i==2:
                re[0],re[1],re[2] = ConvertFunc.Twiss2SigmaZ(s1,s2,s3, f,m,k)
            digits = [4,4,4]
            for j in range(3):
                #if not math.isclose(re[j], float(self.string_sigma[i][j].get()),rel_tol=1e-06):
                self.string_sigma[i][j].set('{:.{}e}'.format(re[j], digits[j]))
                #else:
                #    pass
        except:
            #print('Error')
            pass
        self.updateTwissLock = 0
        
    def preprocessing(self):
        try:
            os.chdir(self.entry_dic.get())
        except:
            print("Error: cannot get to the dictionary" + self.entry_dic.get())
            return
        
        if self.AccKernel=='ImpactT':
            np = self.save('ImpactT.in')
            self.button_run     .config(state='disable')
            self.button_plot    .config(state='disable')
            self.button_initial .config(state='disable')
            PreProcessing.process(self)
            self.load('ImpactT.in')
            self.button_run     .config(state='normal')
            self.button_plot    .config(state='normal')
            self.button_initial .config(state='normal')
        elif self.AccKernel=='ImpactZ':
            print('PreProcessing cannot use at Z code')
            sys.exit()
            #np = self.save('ImpactZ.in')
            #PreProcessing.main('Z')
        else:
            print('Cannot find kernel: '+self.AccKernel)
    
    def switch(self):
        if self.AccKernel=='ImpactT':
            self.switchToImpactZ()
        elif self.AccKernel == 'ImpactZ':
            self.switchToImpactT()
        print("Switch kernel to "+self.AccKernel)
        
    def switchToImpactZ(self):
        self.AccKernel='ImpactZ'
        try:
            self.AdvancedSet.destroy()
            self.AdvancedPlot.destroy()
        except:
            pass
        
        try:
            self.label_logo.config(image = self.logo_ImpactZ)
        except:
            self.label_logo.config(bitmap='error')
        self.entry_Nstep.config(state='disabled')
        self.entry_dt.config(state='disabled')
        
        distZ=['Uniform', 'Gauss', 'SemiGauss',
               'WaterBag','KV', 'Read',
               'Multi-charge-state WaterBag',
               'Multi-charge-state Gaussian',
               'Other']
        if self.distTypeComx.get() not in distZ:
            self.distTypeComx.set(distZ[0])
        self.distType['values'] =distZ
        self.lattice.titleZ()
        
    def switchToImpactT(self):
        self.AccKernel='ImpactT'
        try:
            self.AdvancedSet.destroy()
            self.AdvancedPlot.destroy()
        except:
            pass
        try:
            self.label_logo.config(image = self.logo_ImpactT)
        except:
            self.label_logo.config(bitmap='error')
        self.entry_Nstep.config(state='normal')
        self.entry_dt.config(state='normal')
        self.Flagmap.set('Linear')
        distT =['Uniform', 'Gauss', 'SemiGauss',
                'WaterBag','KV', 'Read',
                'CylcoldZSob','Other']
        self.distType['values'] =distT
        if self.distTypeComx.get() not in distT:
            self.distTypeComx.set('Other')
        self.lattice.titleT()
        
    def makeAdvancedSet(self):
        try:
            self.AdvancedSet.destroy()
        except:
            pass
        if self.AccKernel=='ImpactT':
            self.AdvancedSet = ImpactTSet.AdvancedSetFrame(self,self.AccKernel)
        elif self.AccKernel=='ImpactZ':
            self.AdvancedSet = ImpactZSet.AdvancedSetFrame(self,self.AccKernel)
        else:
            print('Cannot find kernel: '+self.AccKernel)
                    
    def makeAdvancedPlot(self):
        try:
            self.AdvancedPlot.destroy()
        except:
            pass
        if self.AccKernel=='ImpactT':
            self.AdvancedPlot = ImpactTPlot.AdvancedPlotControlFrame(self)  
        elif self.AccKernel=='ImpactZ':
            self.AdvancedPlot = ImpactZPlot.AdvancedPlotControlFrame(self)
        else:
            print('Cannot find kernel: '+self.AccKernel)
            
    def thread_it(self,func):
        self.button_run['state']='disabled'
        self.run_lock.acquire()
        ImpactThread=threading.Thread(target=func)
        ImpactThread.setDaemon(True)
        ImpactThread.start()
        self.run_lock.release()
        self.button_run['state']='normal'
        
    def startImpact(self):
        print("Start " + self.AccKernel)
        if self.AccKernel=='ImpactT':
            try:
                os.chdir(self.entry_dic.get())
            except:
                print("Error: cannot get to the dictionary" + self.entry_dic.get())
                return
            np = self.save('ImpactT.in')
            
            #ImpactExe = os.path.join(sys.path[0],'src',self.IMPACT_T_EXE)
            ImpactExe = self.IMPACT_T_EXE.get()
            
            if np==1:
                cmd = ImpactExe
            elif np>1:
                cmd = self.MPI_EXE.get()+' -n '+str(np)+' '+ImpactExe
            print(cmd)
            p=subprocess.Popen(cmd,stdout=subprocess.PIPE,bufsize=1)
            for line in iter(p.stdout.readline,b''):
                print(('>>{}'.format(line.rstrip())))
            p.stdout.close()
        elif self.AccKernel=='ImpactZ':
            try:
                os.chdir(self.entry_dic.get())
            except:
                print("Error: cannot get to the dictionary" + self.entry_dic.get())
                return
            np = self.save('ImpactZ.in')
            
            #ImpactExe = os.path.join(sys.path[0],'src',_IMPACT_Z_NAME)
            ImpactExe = self.IMPACT_Z_EXE.get()
            if np==1:
                cmd = ImpactExe
            elif np>1:
                cmd = self.MPI_EXE.get()+' -n '+str(np)+' '+ImpactExe
            print(cmd)
            p=subprocess.Popen(cmd,stdout=subprocess.PIPE,bufsize=1)
            for line in iter(p.stdout.readline,b''):
                print(('>>{}'.format(line.rstrip())))
            p.stdout.close()
        else:
            print('Cannot find kernel: '+self.AccKernel)
  
    def validate(self, action, index, value_if_allowed,
                       prior_value, text, validation_type, trigger_type, widget_name):
        for word in text:
            if word in 'eE0123456789.-+*/':
                try:
                    #float(value_if_allowed)
                    return True
                except ValueError:
                    self.bell()
                    return False
            else:
                self.bell()
                return False
            
    def save(self,fileName):
        if os.path.isfile(fileName) == True and os.path.isfile(fileName+'.old') == False:
            copyfile(fileName,fileName+'.old')
        if self.AccKernel=='ImpactT':
            return self.saveImpactT(fileName)
        elif self.AccKernel=='ImpactZ':
            return self.saveImpactZ(fileName)
        else:
            print('Cannot find kernel: '+self.AccKernel)

    def saveImpactT(self,fileName):
        try:
            ImpactInput = open(fileName,'w')
        except:
            print(( "  ERROR! Can't save file as '" + fileName + "'"))
            return False

        ImpactInput.writelines(["!This is the input file for ImpactT. It is generated by GUI.\n", 
                                 "!If you meet any bug, please contact jqiang@lbl.gov\n\n"])
        
        ImpactInput.write( str(int(self.entry_noc1.get()))+' '
                          +str(int(self.entry_noc2.get()))+' '
                          +('5\n' if self.GPUflag.get() else '0\n'))
        np = int(self.entry_noc1.get())*int(self.entry_noc2.get())
        
        ImpactInput.write( str(float(   self.entry_dt.get()))+' '
                          +str(int(     self.entry_Nstep.get())) + ' '
                          +str(int(     self.Nbunch.get()))+'\n')
        ImpactInput.write(str(int(      self.Dim.get()))+' '+
                          str(int(float(self.entry_Np.get())))+' '+
                          INTEGRATOR_TYPE[self.Flagmap.get()] +' '+
                          str(int(float(self.Flagerr.get()))) +' '+
                          DIAGNOSTIC_TYPE[  self.Flagdiag.get()]+ ' '+
                          str(int(      self.Flagimag.get()))+' '+
                          str(float(    self.Zimage.get())) +'\n')
        ImpactInput.write( str(int(     self.entry_Ngx.get()))+' '
                          +str(int(     self.entry_Ngy.get()))+' '
                          +str(int(     self.entry_Ngz.get()))+' '
                          #+BOUNDARY_TYPE[   self.Flagbc.get()]+' '
                          +'1 '
                          +str(float(   self.Xrad.get()))+' '
                          +str(float(   self.Yrad.get()))+' '
                          +str(float(   self.Zrad.get()))+'\n')
        if self.distTypeComx.get() == 'Other':
            ImpactInput.write(self.distTypeNumb.get()    +' '+
                  str(int(float(self.FlagRestart.get())))   +' '+
                  '0 '+    #Flagsbstp
                  str(int(float(self.Nemission.get())))     +' '+
                  str(float(self.Temission.get()))          +'\n')
        else:
            ImpactInput.write(DISTRIBUTION_T_TYPE[self.distTypeComx.get()]    +' '+
                              str(int(float(self.FlagRestart.get())))   +' '+
                              '0 '+    #Flagsbstp
                              str(int(float(self.Nemission.get())))     +' '+
                              str(float(self.Temission.get()))          +'\n')

        for row in range(3):
            for column in range(7):
                ImpactInput.write(str(float(self.string_sigma[row][column].get()))+' ')
            ImpactInput.write('\n')
        '''
        if self.ptcType.get()=='Other':
            pass
        else:
            self.ptcMass.set(PARTICLE_TYPE[self.ptcType.get()].split()[0])
            self.ptcCharge.set(PARTICLE_TYPE[self.ptcType.get()].split()[1])
        '''
        ImpactInput.write(str(float(self.entry_cur.get()))+' '
                          +str(float(self.entry_Ek.get()))+' '
                          +self.ptcMass.get()     + ' '
                          +self.ptcCharge.get()   + ' '
                          +str(float(self.entry_frq.get()))+' '
                          +str(float(self.Tinitial.get()))
                          +'\n\n')
        
        ImpactInput.write('!===========LATTICE===========\n')
        lattice = self.lattice.getHide().splitlines()
        print(lattice)
        for line in lattice:
            if line !='':
                ImpactInput.write(line+'\n')
        ImpactInput.close()
        return np
    
    def saveImpactZ(self,fileName):
        try:
            ImpactInput = open(fileName,'w')
        except:
            print(( "  ERROR! Can't save file as '" + fileName + "'"))
            return False    

        ImpactInput.writelines(["!This is the input file for ImpactZ. It is generated by GUI.\n", 
                                "!If you meet any bug about the GUI, please contact jqiang@lbl.gov\n\n"])
        
        ImpactInput.write( str(int(self.entry_noc1.get()))+' '
                          +str(int(self.entry_noc2.get()))+' '
                          +('5\n' if self.GPUflag.get() else '0\n'))
        np = int(self.entry_noc1.get())*int(self.entry_noc2.get())
        
        ImpactInput.write(str(int(              self.Dim.get()))+' '+
                          str(int(float(        self.entry_Np.get())))+' '+
                          INTEGRATOR_TYPE[      self.Flagmap.get()] +' '+
                          str(int(float(        self.Flagerr.get()))) +' '+
                          OUTPUT_Z_TYPE[        self.FlagOutput_Z.get()]+ '\n')
        
        ImpactInput.write( str(int(     self.entry_Ngx.get()))+' '
                          +str(int(     self.entry_Ngy.get()))+' '
                          +str(int(     self.entry_Ngz.get()))+' '
                          +BOUNDARY_TYPE[   self.Flagbc.get()]+' '
                          +str(float(   self.Xrad.get()))+' '
                          +str(float(   self.Yrad.get()))+' '
                          +str(float(   self.ZperiodSize.get()))+'\n')
        
        if self.distType.get() == 'Other':
            ImpactInput.write(self.distTypeNumb.get()    +' '+
                          str(int(float(self.FlagRestart.get())))   +' '+
                          str(int(float(self.FlagSubcycle.get())))  +' '+
                          str(int(float(self.Nbunch.get())))  +'\n')
        else:
            ImpactInput.write(DISTRIBUTION_Z_TYPE[self.distType.get()]    +' '+
                          str(int(float(self.FlagRestart.get())))   +' '+
                          str(int(float(self.FlagSubcycle.get())))  +' '+
                          str(int(float(self.Nbunch.get())))  +'\n')

        
        ImpactInput.write(self.NpList.get()+ '\n')
        ImpactInput.write(self.CurrentList.get()+ '\n')
        ImpactInput.write(self.SchargeList.get()+ '\n')

        for row in range(3):
            for column in range(7):
                ImpactInput.write(str(float(self.string_sigma[row][column].get()))+' ')
            ImpactInput.write('\n')
        
        ImpactInput.write(str(float(self.entry_cur.get()))+' '
                          +str(float(self.entry_Ek.get()))+' '
                          +self.ptcMass.get()     + ' '
                          +self.ptcCharge.get()   + ' '
                          +str(float(self.entry_frq.get()))+' '
                          +str(float(self.Tinitial.get()))
                          +'\n\n')
        
        ImpactInput.write('!===========LATTICE===========\n')
        lattice = self.lattice.getHide().splitlines()
        print(lattice)
        for line in lattice:
            if line !='':
                ImpactInput.write(line+'\n')
        ImpactInput.close()
        #ImpactInput.write('0.1  5 20 2 7.5 7.5 0.0 6.0d-3 /\n\n\n')
        ImpactInput.close()
        return np
    
    def load(self,inputFileName):
        
        if self.AccKernel=='ImpactT':
            self.loadImpactT(inputFileName)
        elif self.AccKernel=='ImpactZ':
            self.loadImpactZ(inputFileName)
        else:
            print('Cannot find kernel: '+self.AccKernel)
        path = os.path.dirname(os.path.abspath(inputFileName))
        os.chdir(path)
        print(path)
        self.entry_dic.delete(0, 'end')
        self.entry_dic.insert(0,path)

    def loadImpactT(self,inputFileName):
        dataList = ImpactFile.conciseReadInput(inputFileName)

        if dataList ==False:
            return
        '''CPU'''
        self.entry_noc1.delete(0,tk.END)
        self.entry_noc1.insert(0, dataList[0][0])
        self.entry_noc2.delete(0,tk.END)
        self.entry_noc2.insert(0, dataList[0][1])
        try:
            self.GPUflag.set(0 if int(dataList[0][2])==0 else 1)
        except:
            self.GPUflag.set(0)
        '''TimeStep'''
        self.entry_dt.delete(0,tk.END)
        self.entry_dt.insert(0, dataList[1][0])
        self.entry_Nstep.delete(0,tk.END)
        self.entry_Nstep.insert(0, dataList[1][1])
        self.Nbunch.set(dataList[1][2])
              
        '''Particle'''

        self.Dim.set(dataList[2][0])
        self.entry_Np.delete(0,tk.END)
        self.entry_Np.insert(0, dataList[2][1])

        #self.Flagmap.set(1)
        self.Flagerr.set(0 if int(dataList[2][3])==0 else 1)
        
        invermap = dict(map(lambda t:(t[1],t[0]), DIAGNOSTIC_TYPE.items()))
        if dataList[2][4] not in ['1','2','3']:
            dataList[2][4] = '3'
        self.Flagdiag.set(invermap[dataList[2][4]])
        
        self.Flagimag.set(0 if int(dataList[2][5])==0 else 1)
        self.Zimage.set(dataList[2][6])
        
        '''Grid'''
        self.entry_Ngx.delete(0,tk.END)
        self.entry_Ngx.insert(0, dataList[3][0])
        self.entry_Ngy.delete(0,tk.END)
        self.entry_Ngy.insert(0, dataList[3][1])
        self.entry_Ngz.delete(0,tk.END)
        self.entry_Ngz.insert(0, dataList[3][2])
        
        invermap = dict(map(lambda t:(t[1],t[0]), BOUNDARY_TYPE.items()))
        self.Flagbc.set(invermap[dataList[3][3]])
        
        self.Xrad.set(dataList[3][4])
        self.Yrad.set(dataList[3][5])
        self.Zrad.set(dataList[3][6])
        
        '''Distribution'''
        invermap = dict(map(lambda t:(t[1],t[0]), DISTRIBUTION_T_TYPE.items()))
        try:
            if dataList[4][0] in invermap.keys():
                self.distTypeComx.set(invermap[dataList[4][0]])
            else:
                self.distTypeComx.set('Other')
                self.distTypeNumb.set(dataList[4][0])
        except:
            self.distTypeComx.set(invermap['3'])
            print('Cannot recognize distribution type, Set as waterbag')

        self.FlagRestart.set(0 if int(dataList[4][1])==0 else 1)
        self.Nemission.set(dataList[4][3])
        self.Temission.set(dataList[4][4])

        '''Twiss'''
        for row in range(3):
            for column in range(7):
                self.string_sigma[row][column].set(dataList[5+row][column])
        
        '''Particle Type'''
        invermap = dict(map(lambda t:(t[1],t[0]), PARTICLE_TYPE.items()))
        ptcFound = 0
        for key in invermap.keys():
            ptc = key.split()
            try:
                if math.isclose(float(ptc[0]), float(dataList[8][2]),rel_tol=1e-3):
                    if math.isclose(float(ptc[1]), float(dataList[8][3]),rel_tol=1e-3):
                        self.ptcTypeComx.set(invermap[key])
                        ptcFound = 1
                        break
            except:
                pass
        if ptcFound==0:
            self.ptcTypeComx.set('Other')
        self.ptcMass    .set(dataList[8][2])
        self.ptcCharge  .set(dataList[8][3])
         
        """Current"""
        self.entry_cur.delete(0, tk.END)
        self.entry_cur.insert(0, dataList[8][0])
        
        """Energy"""
        self.entry_Ek.delete(0, tk.END)
        self.entry_Ek.insert(0, dataList[8][1])
        
        """Frequency"""
        self.entry_frq.delete(0, tk.END)
        self.entry_frq.insert(0, dataList[8][4])
        
        self.Tinitial.set(dataList[8][5]) 
        
        """Lattice"""
        self.lattice.latticeTextHide.delete(0.0, tk.END)
        for lineSet in dataList[9:]:
            for word in lineSet:
                self.lattice.latticeTextHide.insert(tk.END, word + ' ')
            self.lattice.latticeTextHide.insert(tk.END,'\n')
        self.lattice.update()
        
    def loadImpactZ(self,inputFileName):
        dataList = ImpactFile.conciseReadInput(inputFileName)
        
        if dataList ==False:
            return
        rowtemp=0
        '''CPU'''
        self.entry_noc1.delete(0,tk.END)
        self.entry_noc1.insert(0, dataList[rowtemp][0])
        self.entry_noc2.delete(0,tk.END)
        self.entry_noc2.insert(0, dataList[rowtemp][1])
        try:
            self.GPUflag.set(0 if int(dataList[rowtemp][2])==0 else 1)
        except:
            self.GPUflag.set(0)
        rowtemp+=1
              
        '''Particle'''
        self.Dim.set(dataList[rowtemp][0])
        self.entry_Np.delete(0,tk.END)
        self.entry_Np.insert(0, dataList[rowtemp][1])

        invermap = dict(map(lambda t:(t[1],t[0]), INTEGRATOR_TYPE.items()))
        self.Flagmap.set(invermap[dataList[rowtemp][2]])
        
        self.Flagerr.set(0 if int(dataList[rowtemp][3])==0 else 1)
        
        invermap = dict(map(lambda t:(t[1],t[0]), OUTPUT_Z_TYPE.items()))
        self.Flagdiag.set(invermap[dataList[rowtemp][4]])
        rowtemp+=1
        
        '''Grid'''
        self.entry_Ngx.delete(0,tk.END)
        self.entry_Ngx.insert(0, dataList[rowtemp][0])
        self.entry_Ngy.delete(0,tk.END)
        self.entry_Ngy.insert(0, dataList[rowtemp][1])
        self.entry_Ngz.delete(0,tk.END)
        self.entry_Ngz.insert(0, dataList[rowtemp][2])
        
        invermap = dict(map(lambda t:(t[1],t[0]), BOUNDARY_TYPE.items()))
        self.Flagbc.set(invermap[dataList[rowtemp][3]])
        
        self.Xrad.set(str(float(dataList[rowtemp][4])))
        self.Yrad.set(str(float(dataList[rowtemp][5])))
        self.ZperiodSize.set(str(float(dataList[rowtemp][6])))
        rowtemp+=1
        
        '''Distribution'''
        invermap = dict(map(lambda t:(t[1],t[0]), DISTRIBUTION_Z_TYPE.items()))
        
        try:
            if dataList[rowtemp][0] in invermap.keys():
                self.distTypeComx.set(invermap[dataList[rowtemp][0]])
            else:
                self.distTypeComx.set('Other')
                self.distTypeNumb.set(dataList[rowtemp][0])
        except:
            self.distTypeComx.set(invermap['3'])
            print('Cannot recognize distribution type, Set as waterbag')

        self.FlagRestart.set( 0 if int(dataList[rowtemp][1])==0 else 1)
        self.FlagSubcycle.set(0 if int(dataList[rowtemp][2])==0 else 1)
        self.Nbunch.set(dataList[rowtemp][3])
        rowtemp+=1
        
        '''Multiple Charge State'''
        self.NpList.set(' '.join(dataList[rowtemp])) 
        rowtemp+=1
        self.CurrentList.set(' '.join(dataList[rowtemp])) 
        rowtemp+=1
        self.SchargeList.set(' '.join(dataList[rowtemp])) 
        rowtemp+=1
        
        '''Twiss'''
        for row in range(3):
            for column in range(7):
                self.string_sigma[row][column].set(dataList[rowtemp][column])
            rowtemp+=1

        '''Particle Type'''
        #print(PARTICLE_TYPE)
        invermap = dict(map(lambda t:(t[1],t[0]), PARTICLE_TYPE.items()))
        #print(invermap)
        ptcFound = 0
        for key in invermap.keys():
            ptc = key.split()
            try:
                if math.isclose(float(ptc[0]), float(dataList[rowtemp][2]),rel_tol=1e-3):
                    if math.isclose(float(ptc[1]), float(dataList[rowtemp][3]),rel_tol=1e-3):
                        self.ptcTypeComx.set(invermap[key])
                        ptcFound = 1
                        break
            except:
                pass
        if ptcFound==0:
            self.ptcTypeComx.set('Other')
        
        self.ptcMass    .set(dataList[rowtemp][2])
        self.ptcCharge  .set(dataList[rowtemp][3])
        #print(rowtemp,dataList[rowtemp])
            
        """Current"""
        self.entry_cur.delete(0, tk.END)
        self.entry_cur.insert(0, dataList[rowtemp][0])
        
        """Energy"""
        self.entry_Ek.delete(0, tk.END)
        self.entry_Ek.insert(0, dataList[rowtemp][1])
        
        """Frequency"""
        self.entry_frq.delete(0, tk.END)
        self.entry_frq.insert(0, dataList[rowtemp][4])
        
        self.Tinitial.set(dataList[rowtemp][5]) 
        
        rowtemp+=1
        """Lattice"""
        self.lattice.latticeTextHide.delete(0.0, tk.END)
        for lineSet in dataList[rowtemp:]:
            for word in lineSet:
                self.lattice.latticeTextHide.insert(tk.END, word + ' ')
            self.lattice.latticeTextHide.insert(tk.END,'\n')
        self.lattice.update()

        
        
class PlotControlFrame(tk.Frame):
    """Output"""
    def __init__(self, master=None, cnf={}, **kw):
        tk.Frame.__init__(self, master, cnf, **kw)
        """Plot Control"""
        self.frame_plotButton = tk.Frame(self, height =_height/5, width = _width)
        self.frame_plotButton.grid(column=0, row = 0, pady=5 ,padx=10, sticky="e")
        
        self.frame_se = tk.Frame(self.frame_plotButton)
        self.frame_se.pack(side='left')
        self.frame_radio = tk.Frame(self.frame_se)
        self.frame_radio.pack(side='left')
        
        self.plotDirct = tk.IntVar()
        self.plotDirct.set(0)
        self.frame_radio.x = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="X", value=0)
        self.frame_radio.x.pack(side='left')
        self.frame_radio.y = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="Y", value=1)
        self.frame_radio.y.pack(side='left')
        self.frame_radio.z = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="Z", value=2)
        self.frame_radio.z.pack(side='left')
        
        self.plotTypeComx = tk.StringVar(self.frame_se,'Rms size')
        self.plotType = ttk.Combobox(self.frame_se,text=self.plotTypeComx,
                                     values=list(ImpactMainWindow.PLOTTYPE.keys()))
        self.plotType.pack(side = 'left')
        
        self.plot = tk.Button(self.frame_plotButton,text='plot',command=self.makePlot)
        self.plot.pack(fill = 'both',expand =1,side = 'left',padx=10)
        
        self.t = ttk.Separator(self, orient=tk.HORIZONTAL).grid(column=0, row = 1, sticky="we")

           
        self.frame2 = tk.Frame(self, height =_height/5, width = _width)
        self.frame2.grid(column=0, row = 2, pady=5 ,padx=10, sticky="we")
        
        self.ParticlePlot = tk.Button(self.frame2,text='Phase Space Plot',
                               command = self.makeParticlePlot)
        self.ParticlePlot.pack(fill = 'both',expand =1,side = 'left')
        
        self.button_AdvancedPlot = tk.Button(self.frame2,text='Advanced Plot',
                               command = self.makeAdvancedPlot)
        self.button_AdvancedPlot.pack(fill = 'both',expand =1,side = 'left')
    
    def makePlot(self):
        print((self.plotType))
        
        PlotFileName='fort.'+str(self.plotDirct.get()+24)
        yx=ImpactMainWindow.PLOTTYPE[self.plotType.get()]
        yl=yx if self.plotDirct.get()!=2 else yx-1

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=ImpactTPlot.PlotFrame(plotWindow,PlotFileName,0,yl)
        l.pack()
        
    def makeParticlePlot(self):
        print((self.plotType))
        filename = filedialog.askopenfilename(parent=self)
        try:
            t=open(filename)
            t.close()
        except:
            return
        
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Phase Space Plot')
        
        l=ImpactTPlot.PhaseSpaceFrame(plotWindow,filename)
        l.pack()            
    def makeAdvancedPlot(self):
        try:
            self.AdvancedPlot.destroy()
        except:
            pass
        if self.AccKernel=='ImpactT':
            self.AdvancedPlot = ImpactTPlot.AdvancedPlotControlFrame(self)  

class startWindow(tk.Toplevel):
    def __init__(self, master, cnf={}, **kw):
        tk.Toplevel.__init__(self, master, cnf, **kw)
        self.title('Start')
        self.focus_set()
        
        self.frame = tk.Frame(self,width = 300,height = 200)
        self.frame.pack(fill = 'both',expand =1,side = 'top')
        LARGE_FONT= ("Helvetica", 20,'italic')     
        
        
        self.button_ImpactT = tk.Button(self.frame,text='ImpactT',font = LARGE_FONT,
                                        command = lambda: self.startImpactT(master))
        self.button_ImpactT.pack(fill = 'both',expand =1,side = 'left')
        
        self.button_ImpactT.bind("<Enter>", lambda event, h=self.button_ImpactT: h.configure(bg="#00CD00"))
        self.button_ImpactT.bind("<Leave>", lambda event, h=self.button_ImpactT: h.configure(bg="#FFFFFF"))

        
        self.button_ImpactZ = tk.Button(self.frame,text='ImpactZ',font = LARGE_FONT,
                                        command = lambda: self.startImpactZ(master))
        self.button_ImpactZ.pack(fill = 'both',expand =1,side = 'left')
        
        self.button_ImpactZ.bind("<Enter>", lambda event, h=self.button_ImpactZ: h.configure(bg="#00CD00"))
        self.button_ImpactZ.bind("<Leave>", lambda event, h=self.button_ImpactZ: h.configure(bg="#FFFFFF"))
        
        self.button_close = tk.Button(self,text='QUIT',font = LARGE_FONT,
                                        command = lambda: master.quit())
        self.button_close.pack(fill = 'both',expand =0,side = 'bottom')
        
        self.button_close.bind("<Enter>", lambda event, h=self.button_close: h.configure(bg="#00CD00"))
        self.button_close.bind("<Leave>", lambda event, h=self.button_close: h.configure(bg="#FFFFFF"))
        
        self.protocol('WM_DELETE_WINDOW', lambda: master.quit())
        
        
    def startImpactT(self,master):
        master.switchToImpactT()
        master.deiconify()
        self.resize(master)
        self.destroy()
        
    def startImpactZ(self,master):
        master.switchToImpactZ()
        master.deiconify()
        self.resize(master)
        self.destroy()
        
    def closeWindow(self,master):
        master.quit()
        
    def resize(self,master):
        w  = master.winwidth
        h  = master.winheight

        
        x = master.winPosX
        y = master.winPosY
        print(w,h,x,y)
        master.geometry('%dx%d+%d+%d' % (w, h, x, y))
        
        
class MyMenu():
    def __init__(self, root):

        self.menubar = tk.Menu(root)
        
        filemenu = tk.Menu(self.menubar, tearoff=0)
        filemenu.add_command(label="Load", command=lambda: self.file_open(root))
        filemenu.add_command(label="Save", command=self.file_save)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=root.exit)
        
        controlMenu = tk.Menu(self.menubar, tearoff=0)
        controlMenu.add_command(label="Switch to ImpactT", command=lambda: self.control_switchToImpactT(root))
        controlMenu.add_command(label="Switch to ImpactZ", command=lambda: self.control_switchToImpactZ(root))

        
        helpmenu = tk.Menu(self.menubar, tearoff=0)
        helpmenu.add_command(label="Help", command=self.help_help)
        helpmenu.add_command(label="About", command=self.help_about)
        self.menubar.add_cascade(label="File", menu=filemenu)
        self.menubar.add_cascade(label="Switch", menu=controlMenu)
        self.menubar.add_cascade(label="Help", menu=helpmenu)
        
        root.config(menu=self.menubar)
      
    def file_new(self):
        messagebox.showinfo('about', 'GUI authors: Zhicong Liu and Ji Qiang \n verion 1.0 \n jqiang@lbl.gov ')
        pass
          
    def file_open(self,root):
        filename = filedialog.askopenfilename(parent=root)
        if filename=='':
            return
        root.load(filename)
        pass

    def file_save(self):
        filename = filedialog.asksaveasfilename(parent=root)
        if filename=='':
            return
        root.save(filename)
        pass
        
    def control_switchToImpactT(self,root):
        root.switchToImpactT()
        root.bell()
        pass
    
    def control_switchToImpactZ(self,root):
        root.switchToImpactZ()
        root.bell()
        pass
        
    def edit_copy(self):
        messagebox.showinfo('about', 'GUI authors: Zhicong Liu, Ji Qiang \n verion 1.0 \n jqiang@lbl.gov ')
        pass
        
    def edit_paste(self):
        messagebox.showinfo('about', 'GUI authors: Zhicong Liu, Ji Qiang \n verion 1.0 \n jqiang@lbl.gov ')
        pass
    
    def help_help(self):
        messagebox.showinfo('Help', 'Please refer "ImpactTv1.8.pdf" at this folder ')
        pass
    
    def help_about(self):
        messagebox.showinfo('About', ' Version 1.0 \n\n Kernel author: Ji Qiang \n jqiang@lbl.gov \n\n GUI authors: Zhicong Liu, Ji Qiang \n jqiang@lbl.gov ')
        pass
        
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
        print(1,base_path)
    except Exception:
        base_path = os.path.abspath(".")
    #print(base_path,relative_path,os.path.join(base_path, relative_path))
    return os.path.join(base_path, relative_path)


def quitConfirm():
    if messagebox.askokcancel("Quit", "Do you really wish to quit?"):
        root.destroy()


'''
            pyinsteller -F -add-data "icon\ImpactT.gif;icon" ImpactGUI.py
'''

'''
1, label at plot
2, unit
3, plot size/location
4, emit growth divided by zero
ok:
1,3,4
'''

'''
1, density plot
'''
