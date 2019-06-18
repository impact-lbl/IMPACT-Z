'''
'''

import tkinter as tk
from tkinter import filedialog,ttk
import os,sys

_height=200
_width =200
class AdvancedSetFrame(tk.Toplevel):            
    def __init__(self, master, kernel, cnf={}, **kw):
        tk.Toplevel.__init__(self, master, cnf, **kw)
        self.title('Set')
        self.focus_set()  
        #self.grab_set()
        rowtemp=0
        self.frame1 = tk.LabelFrame(self, height =_height/10, width = _width,
                                          text="Configuration")
        self.frame1.pack(side = 'top')
        
        """Exe path"""
        self.frame_exe = tk.LabelFrame(self.frame1)
        self.frame_exe.grid(row=rowtemp,column=0,rowspan=2,columnspan=2,sticky='w')
        rowtemp+=2
        
        self.frame_input1 = tk.Frame(self.frame_exe)
        self.frame_input1.grid(row=0,column=0,columnspan=2,sticky='w')
        
        self.label_exePat    = tk.Label(self.frame_input1, text="mpirun:")
        self.label_exePat.pack(side='left')
        self.entry_exePath = tk.Entry(self.frame_input1, width=45,textvariable=master.MPI_EXE)
        self.entry_exePath.pack(side='left')
        
        self.frame_input = tk.Frame(self.frame_exe,borderwidth=0, highlightthickness=0)
        self.frame_input.grid(row=1,column=0,columnspan=2,sticky='w')
        
        self.label_exePat    = tk.Label(self.frame_input, text="Exe:")
        self.label_exePat.pack(side='left')
        self.entry_exePath = tk.Entry(self.frame_input, width=45,textvariable=master.IMPACT_T_EXE)
        self.entry_exePath.pack(side='left')
        
        self.button_exePath = tk.Button(self.frame_input)
        self.button_exePath["text"]        = "..."
        self.button_exePath["command"]     = lambda: self.changeExe(master.IMPACT_T_EXE)
        self.button_exePath.pack(side = 'left')
        
        self.button_exePath.bind("<Enter>", lambda event, h=self.button_exePath: h.configure(bg="#00CD00"))
        self.button_exePath.bind("<Leave>", lambda event, h=self.button_exePath: h.configure(bg="#FFFFFF"))
        
        
        """Flag restart"""
        self.check_restart  = tk.Checkbutton(self.frame1, 
                                             text="Restart at previous check point", 
                                             variable=master.FlagRestart)
        self.check_restart  .grid(row=rowtemp,column=0,columnspan=2,sticky='w')
        rowtemp+=1
        
        """Flag error"""
        self.check_error    = tk.Checkbutton(self.frame1, text="Misalignmnet and rotation error", variable=master.Flagerr)
        self.check_error    .grid(row=rowtemp,column=0,columnspan=2,sticky='w')
        rowtemp+=1
        
        """Flag image"""
        self.check_image  = tk.Checkbutton(self.frame1, text="Image current", variable=master.Flagimag)
        self.check_image.grid(row=rowtemp,column=0,columnspan=2,sticky='w')
        if kernel!='ImpactT':
            self.check_image     .config(state='disabled')
        rowtemp+=1
        
        """Zimage"""
        self.label_dt5    = tk.Label(self.frame1, text="Z image")
        self.entry_dt5    = tk.Entry(self.frame1,textvariable=master.Zimage)
        self.label_dt5.grid(row=rowtemp,sticky='W')
        self.entry_dt5.grid(row=rowtemp,column=1,sticky='E')
        if kernel!='ImpactT':
            self.entry_dt5     .config(state='disabled')
        rowtemp+=1
        
        self.t = ttk.Separator(self.frame1, orient=tk.HORIZONTAL).grid(row = rowtemp,column=0,columnspan=2, sticky="we")
        rowtemp+=1
         
        self.label_mass    = tk.Label(self.frame1, text="Mass (eV)")
        self.entry_mass    = tk.Entry(self.frame1,textvariable=master.ptcMass)
        self.label_mass.grid(row=rowtemp,sticky='W')
        self.entry_mass.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        self.label_charge    = tk.Label(self.frame1, text="Charge (e)")
        self.entry_charge    = tk.Entry(self.frame1,textvariable=master.ptcCharge)
        self.label_charge.grid(row=rowtemp,sticky='W')
        self.entry_charge.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        """Nbunch"""
        self.label_dt0    = tk.Label(self.frame1, text="Nbunch")
        self.entry_dt0    = tk.Entry(self.frame1,textvariable=master.Nbunch)
        self.label_dt0.grid(row=rowtemp,sticky='W')
        self.entry_dt0.grid(row=rowtemp,column=1,sticky='E')
        self.entry_dt0.config(state='disabled')
        rowtemp+=1
        
        """Dimension"""
        self.label_dt1    = tk.Label(self.frame1, text="Dimension")
        self.entry_dt1    = tk.Entry(self.frame1,textvariable=master.Dim)
        self.label_dt1.grid(row=rowtemp,sticky='W')
        self.entry_dt1.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Dimension"""
        self.label_dist    = tk.Label(self.frame1, text="Distribution")
        self.entry_dist    = tk.Entry(self.frame1,textvariable=master.distTypeNumb)
        self.label_dist.grid(row=rowtemp,sticky='W')
        self.entry_dist.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        '''Type of integrator
        self.label_integrator       = tk.Label(self.frame1, text="Integrator")
        self.box_integrator         = ttk.Combobox(self.frame1,textvariable=master.Flagmap,
                                           values=['Linear', 'Non Linear'])
        self.label_integrator   .grid(row=rowtemp,sticky='W')
        self.box_integrator     .grid(row=rowtemp,column=1,sticky='E')
        if kernel=='ImpactT':
            self.box_integrator     .config(state='disabled')
        rowtemp+=1
        
        '''
        """Flag diagnostics"""
        self.label_diag     = tk.Label(self.frame1, text="Output")
        self.label_diag     .grid(row=rowtemp,sticky='W')
        self.box_diag       = ttk.Combobox(self.frame1,textvariable=master.Flagdiag,
                                           values=['At given time', 'At bunch centroid', 'No output'])
        self.box_diag       .grid(row=rowtemp, column=1, sticky='E')
        rowtemp+=1
        

        
        self.t = ttk.Separator(self.frame1, orient=tk.HORIZONTAL).grid(row = rowtemp,column=0,columnspan=2, sticky="we")
        rowtemp+=1
        
        """Open boundary condition
        self.label_bc       = tk.Label(self.frame1, text='Boundary')
        self.label_bc       .grid(row=rowtemp,sticky='W')
        self.box_bc         = ttk.Combobox(self.frame1,textvariable=master.Flagbc,
                                           values=['Open boundary condition'])
        self.box_bc         .grid(row=rowtemp, column=1, sticky='E')
        if kernel=='ImpactT':
            self.box_bc.config(state='disabled')
        rowtemp+=1
        """

        
        """Computational domain"""
        self.label_domain   = tk.Label(self.frame1, text='Computational Domain:')
        self.label_domain   .grid(row=rowtemp,column=0,columnspan=2,sticky='w')
        rowtemp+=1
        self.label_Xrad     = tk.Label(self.frame1, text='X:')
        self.entry_Xrad     = tk.Entry(self.frame1,textvariable=master.Xrad)
        self.label_Xrad .grid(row=rowtemp,sticky='E')
        self.entry_Xrad .grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        self.label_Yrad    = tk.Label(self.frame1, text='Y:')
        self.entry_Yrad    = tk.Entry(self.frame1,textvariable=master.Yrad)
        self.label_Yrad.grid(row=rowtemp,sticky='E')
        self.entry_Yrad.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        self.label_Zrad    = tk.Label(self.frame1, text='Z:')
        self.entry_Zrad    = tk.Entry(self.frame1,textvariable=master.Zrad)
        self.label_Zrad.grid(row=rowtemp,sticky='E')
        self.entry_Zrad.grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        self.t = ttk.Separator(self.frame1, orient=tk.HORIZONTAL).grid(row = rowtemp,column=0,columnspan=2, sticky="we")
        rowtemp+=1
        

        
        self.label_Nemission = tk.Label(self.frame1, text='Emission Step')
        self.entry_Nemission = tk.Entry(self.frame1,textvariable=master.Nemission)
        self.label_Nemission .grid(row=rowtemp,sticky='E')
        self.entry_Nemission .grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        self.label_Temission = tk.Label(self.frame1, text='Emission Time')
        self.entry_Temission = tk.Entry(self.frame1,textvariable=master.Temission)
        self.label_Temission .grid(row=rowtemp,sticky='E')
        self.entry_Temission .grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        self.label_Tinitial = tk.Label(self.frame1, text='Initial reference time')
        self.entry_Tinitial = tk.Entry(self.frame1,textvariable=master.Tinitial)
        self.label_Tinitial .grid(row=rowtemp,sticky='E')
        self.entry_Tinitial .grid(row=rowtemp,column=1,sticky='E')
        rowtemp+=1
        
        """Advanced distribution"""
        self.frame_Twiss = tk.LabelFrame(self, height =_height/10, width = _width,
                                          text="Initial Distribution")
        self.frame_Twiss.pack(side = 'top')
        twisswidth = 10

        self.twiss_s = []
        self.twiss_chara = ["sigma(m)","sigmaP","muxpx",'xScale','PxSacle','x shift','px shift']
        for column in range(7):
            label = tk.Label(self.frame_Twiss, text=self.twiss_chara[column], 
                            borderwidth=0,
                            width=twisswidth)
            label.grid(row=0, column=column+1, sticky="ns", padx=1, pady=1)
            self.twiss_s.append(label)
        
        self.twiss_x = []
        self.twiss_xchara = ["X","Y","Z"]
        for row in range(3):
            label = tk.Label(self.frame_Twiss, text=self.twiss_xchara[row], 
                                 borderwidth=0, width=2)
            label.grid(row=row+1, column=0, sticky="ns", padx=1, pady=1)
            self.twiss_x.append(label)
        
        self._twiss = []
        for row in range(3):
            current_row = []
            for column in range(7):
                label = tk.Entry(self.frame_Twiss, 
                                 textvariable=master.string_sigma[row][column], 
                                 borderwidth=0,
                                 width=twisswidth)
                label.grid(row=row+1, column=column+1, sticky="ns", padx=2, pady=1)
                current_row.append(label)
            self._twiss.append(current_row)
            
    def changeExe(self,exePath):
        filename = filedialog.askopenfilename(parent=self)
        if filename=='':
            return
        
        exePath.set(filename)
        
        print(exePath.get())
        
