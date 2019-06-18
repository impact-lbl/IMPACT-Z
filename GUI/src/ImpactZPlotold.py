#!/usr/bin/env python
#This code is to plot the result from ImpactZ
#Zhicong@21/10/2016
#Input : fort.xx
#Output: figures about beam size and emittance
# plots are saved at '/post'

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import ttk,filedialog
import time,os,sys
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

#import scipy.ndimage
from scipy.stats import gaussian_kde
import numpy as np
import ParticlePlot

_height=300
_width =200

IMPACT_Z_ADVANCED_PLOT_TYPE= {'Centriod location (mm)'    :1,
                     'Rms size (mm)'             :2,
                     'Centriod momentum (MC)'    :3,
                     'Rms momentum (MC)'         :4,
                     'Twiss'                     :5,
                     'Emittance (mm-mrad)'       :6}

IMPACT_Z_SciFormatter = FormatStrFormatter('%2.1E')
IMPACT_Z_sciMaxLimit  = 99999 *2
IMPACT_Z_sciMinLimit  = 0.0001*2

class AdvancedPlotControlFrame(tk.Toplevel):
    """Output"""
            
    def __init__(self, master=None, cnf={}, **kw):
        tk.Toplevel.__init__(self, master, cnf, **kw)
        self.title('ImpactZ Plot')
        self.focus_set()  
        """Plot Control"""
        self.frame_plotButton = tk.Frame(self)
        self.frame_plotButton.grid(column=0, row = 0, pady=5 ,padx=10, sticky="we")
        
        self.frame_radio = tk.Frame(self.frame_plotButton)
        self.frame_radio.pack(side='top')
        
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
        
        self.plotTypeComx = tk.StringVar(self.frame_plotButton,'Rms size (mm)')
        self.plotType = ttk.Combobox(self.frame_plotButton,text=self.plotTypeComx,
                                     width = 20,
                                     values=list(IMPACT_Z_ADVANCED_PLOT_TYPE.keys()))
        self.plotType.pack(side = 'top')
        self.plot = tk.Button(self.frame_plotButton,text='plot',command=self.makePlot)
        self.plot.pack(fill = 'both',expand =1,side = 'top',padx=10)
        
        self.t = ttk.Separator(self, orient=tk.HORIZONTAL).grid(column=0, row = 1, sticky="we")

           
        self.frame2 = tk.Frame(self, height =_height/5, width = _width)
        self.frame2.grid(column=0, row = 2, pady=5 ,padx=10, sticky="nswe")
        
        rowN=0
        
        self.button_overall = tk.Button(self.frame2,text='Overall',
                               command = self.overallPlot)
        self.button_overall.grid(row = rowN, column=0,  pady=5 ,padx=5, columnspan = 2, sticky="nswe")
        rowN+=1
        
        self.button_emitGrowth      = tk.Button(self.frame2,text='EmitGrowth',
                                                command = self.emitGrowthPlot)
        self.button_emitGrowth      .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_Ek              = tk.Button(self.frame2,text='Kinetic Energy',
                                                command = lambda: self.energyPlot(3,'Kinetic Energy (MeV)'))
        self.button_Ek              .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        '''
        self.button_beta            = tk.Button(self.frame2,text='Beta',
                                                command = lambda: self.energyPlot(4,'Beta'))
        self.button_beta            .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_gamma           = tk.Button(self.frame2,text='Gamma',
                                                command = lambda: self.energyPlot(2,'Gamma'))
        self.button_gamma           .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        '''
        self.button_rmax            = tk.Button(self.frame2,text='Rmax',
                                                command = lambda: self.energyPlot(5,'Rmax (mm)'))
        self.button_rmax            .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_dw              = tk.Button(self.frame2,text='Absolute phase',
                                                command = lambda: self.energyPlot(1,'Absolute phase (rad)'))
        self.button_dw              .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        
        self.button_Temperature         = tk.Button(self.frame2,text='Temperature Plot',
                                                    command = self.makeTemperaturePlot)
        self.button_Temperature         .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_Loss                = tk.Button(self.frame2,text='live Particle #',
                                                    command = self.liveParticlePlot)
        self.button_Loss                .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        
        self.t = ttk.Separator(self.frame2, orient=tk.HORIZONTAL).grid(column=0, row = rowN, columnspan=2,sticky="we")        
        rowN+=1
        
        self.max                        = tk.Button(self.frame2,text='Max amplitude',
                                                    command = self.maxPlot)
        self.max                        .grid(row = rowN, column=0,  pady=5 ,padx=5, columnspan=2,sticky="nswe")
        rowN+=1
        
        self.button_3order              = tk.Button(self.frame2,text='3 order parameter',
                                                    command = self.make3orderPlot)
        self.button_3order              .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_4order              = tk.Button(self.frame2,text='4 order parameter',
                                                    command = self.make4orderPlot)
        self.button_4order              .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        
        self.t = ttk.Separator(self.frame2, orient=tk.HORIZONTAL).grid(column=0, row = rowN, columnspan=2,sticky="we")        
        rowN+=1

        scaling = float(master.entry_frq.get())*2*3.1415926/299792458
        self.button_Particle            = tk.Button(self.frame2,text='Phase Space Plot',
                                                    command = lambda:self.ParticlePlot(scaling))
        self.button_Particle            .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_ParticleDesity1D    = tk.Button(self.frame2,text='Density1D',
                                                    command = lambda:self.ParticleDensityPlot1D(scaling))
        self.button_ParticleDesity1D    .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        
        self.button_ParticleDensity     = tk.Button(self.frame2,text='Density2D (by Grid)',
                                                    command = lambda:self.ParticleDensityPlot(scaling))
        self.button_ParticleDensity     .grid( row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_ParticleDensity2    = tk.Button(self.frame2,text='Density2D (by Ptc)',
                                                    command = lambda:self.ParticleDensityPlot2(scaling))
        self.button_ParticleDensity2    .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1

    def overallPlot(self):
        print(self.__class__.__name__)

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=OverallFrame(plotWindow)
        l.pack()    
        
    def energyPlot(self,y,yLabel):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title(sys._getframe().f_back.f_code.co_name)
        
        l=PlotFrame(plotWindow,'fort.18',0,y,yLabel)
        l.pack()
    
    def emitGrowthPlot(self):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=EmitGrowthFrame(plotWindow)
        l.pack()   
        
    def makeTemperaturePlot(self):
        print((self.plotType))

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=TemperatureFrame(plotWindow)
        l.pack()
        
    def liveParticlePlot(self):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title(sys._getframe().f_back.f_code.co_name)
        
        l=PlotFrame(plotWindow,'fort.28',0,3,'Live particle number')
        l.pack()
        
    def ParticlePlot(self,scaling):
        print(self.__class__.__name__)
        filename = filedialog.askopenfilename(parent=self)
        try:
            t=open(filename)
            t.close()
        except:
            return
        
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Phase Space Plot')
        
        l=ParticlePlot.ParticleFrame(plotWindow,filename,scaling,'ImpactZ')
        l.pack() 
                
    def ParticleDensityPlot(self,scaling):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=ParticlePlot.ParticleDensityFrame_weight2D(plotWindow,fileName,scaling,'ImpactZ')
        l.pack()
        
    def ParticleDensityPlot1D(self,scaling):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=ParticlePlot.ParticleDensityFrame_weight1D(plotWindow,fileName,scaling,'ImpactZ')
        l.pack()
                
    def ParticleDensityPlot2(self,scaling):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=ParticlePlot.ParticleDensityFrame2D_slow(plotWindow,fileName,scaling,'ImpactZ')
        l.pack()
        
    def makePlot(self):
        print(self.__class__.__name__)
        
        PlotFileName='fort.'+str(self.plotDirct.get()+24)        
        yx=IMPACT_Z_ADVANCED_PLOT_TYPE[self.plotType.get()]
        yl=yx if self.plotDirct.get()!=2 else yx-1

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')
        
        l=PlotFrame(plotWindow,PlotFileName,0,yl,self.plotType.get())
        l.pack()
        

    def maxPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.27'
        try:
            t=open(filename)
            t.close()
        except:
            return
        
        plotWindow = tk.Toplevel(self)
        plotWindow.title('maxPlot')
        
        l=PlotMaxFrame(plotWindow,filename)
        l.pack() 
    def make3orderPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.29'
        try:
            t=open(filename)
            t.close()
        except:
            return
        
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Cubic root of 3rd moment')
        
        l=PlotHighorderFrame(plotWindow,filename)
        l.pack() 
        
    def make4orderPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.30'
        try:
            t=open(filename)
            t.close()
        except:
            return
        
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Square root, square root of 4th moment')
        
        l=PlotHighorderFrame(plotWindow,filename)
        l.pack() 

class PlotBaseFrame(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.fig = Figure(figsize=(6,5), dpi=100)
        self.subfig = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
    
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        

        
class PlotFrame(tk.Frame):
    def __init__(self, parent,PlotFileName,xl,yl,labelY):
        tk.Frame.__init__(self, parent)
        #LARGE_FONT= ("Verdana", 12)
        #label = tk.Label(self, font=LARGE_FONT,
        #                 text='plot '+PlotFileName+
        #                 ' use '+str(xl)+':'+str(yl))
        #label.pack(pady=10,padx=10)

        try:
            fin = open(PlotFileName,'r')
        except:
            print(( "  ERRPR! Can't open file '" + PlotFileName + "'"))
        
        linesList  = fin.readlines()
        fin .close()
        linesList  = [line.split() for line in linesList ]
        x   = np.array([float(xrt[xl]) for xrt in linesList])
        y   = np.array([float(xrt[yl]) for xrt in linesList])
        
        if labelY in ['Centriod location (mm)','Rms size (mm)','Rmax (mm)']:
            y = y*1.0e3       # unit convert from m to mm
        elif labelY in ['Emittance (mm-mrad)']:
            y = y*1.0e6       # unit convert from (m-rad) to (mm-mrad)
        
        fig = Figure(figsize=(7,5), dpi=100)
        subfig = fig.add_subplot(111)
        subfig.plot(x,y)
        subfig.set_xlabel('Z (m)')
        subfig.set_ylabel(labelY)
        
        
        xMax = np.max(x)
        xMin = np.min(x)
        yMax = np.max(y)
        yMin = np.min(y)
        if (xMax-xMin)>IMPACT_Z_sciMaxLimit or (xMax-xMin)<IMPACT_Z_sciMinLimit:
            self.subfig.xaxis.set_major_formatter(IMPACT_Z_SciFormatter)
        if (yMax-yMin)>IMPACT_Z_sciMaxLimit or (yMax-yMin)<IMPACT_Z_sciMinLimit:
            self.subfig.yaxis.set_major_formatter(IMPACT_Z_SciFormatter)
        
        box = subfig.get_position()
        subfig.set_position([box.x0*1.3, box.y0*1.1, box.width, box.height])
        
        canvas = FigureCanvasTkAgg(fig, self) 
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def quit(self):
        self.destroy()
                
class OverallFrame(tk.Frame):
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.fig = Figure(figsize=(12,5), dpi=100)
        self.subfig = []
        self.subfig.append(self.fig.add_subplot(221))
        self.subfig.append(self.fig.add_subplot(222))
        self.subfig.append(self.fig.add_subplot(223))
        self.subfig.append(self.fig.add_subplot(224))

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
    
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.plot()
    def plot(self):
        picNum = 4
        fileList    = [[]*2]*picNum
        saveName    = []
        labelList   = [[]*2]*picNum
        xdataList   = [[]*2]*picNum
        ydataList   = [[]*2]*picNum
        xyLabelList = [[]*2]*picNum
        
        saveName.append('sizeX')
        fileList[0]     = ['fort.24','fort.27']
        labelList[0]    = ['rms.X','max.X']
        xdataList[0]    = [0,0]
        ydataList[0]    = [2,1]
        xyLabelList[0]  = ['z drection (m)','beam size in X (mm)']
        
        saveName.append('sizeY')
        fileList[1]     = ['fort.25','fort.27']
        labelList[1]    = ['rms.Y','max.Y']
        xdataList[1]    = [0,0]
        ydataList[1]    = [2,3]
        xyLabelList[1]  = ['z drection (m)','beam size in Y (mm)']
        
        saveName.append('sizeZ')
        fileList[2]     = ['fort.26','fort.27']
        labelList[2]    = ['rms.Z','max.Z']
        xdataList[2]    = [0,0]
        ydataList[2]    = [2,5]
        xyLabelList[2]  = ['z drection (m)','beam size in Z (mm)']
        
        saveName.append('emitXY')
        fileList[3]     = ['fort.24','fort.25']
        labelList[3]    = ['emit.nor.X','emit.nor.Y']
        xdataList[3]    = [0,0]
        ydataList[3]    = [6,6]
        xyLabelList[3]  = ['z drection (m)','emittance at X and Y (mm*mrad)']
        
        lineType = ['r-','b--']

        for i in range(0,picNum):
            for j in range(0,2):
                try:
                    fin = open(fileList[i][j],'r')
                except:
                    print("ERRPR Can't open file ' " + fileList[i][j] + "'")
                    return
                linesList  = fin.readlines()
                fin .close()
                linesList  = [line.split() for line in linesList ]
                xId = xdataList[i][j]
                yId = ydataList[i][j]
                x   = np.array([float(xrt[xId]) for xrt in linesList])
                y   = np.array([float(xrt[yId]) for xrt in linesList])
                if i in range(0,picNum-1):
                    y=y*1.0e3
                elif i == picNum-1:
                    y=y*1.0e6
                self.subfig[i].plot(x, y, lineType[j], linewidth=2, label=labelList[i][j])

            self.subfig[i].set_xlabel(xyLabelList[i][0])
            self.subfig[i].set_ylabel(xyLabelList[i][1])
            box = self.subfig[i].get_position()
            self.subfig[i].set_position([box.x0*1.1, box.y0*1.1, box.width, box.height *0.88])
            
            xMax = np.max(x)
            xMin = np.min(x)
            yMax = np.max(y)
            yMin = np.min(y)
            if (xMax-xMin)>IMPACT_Z_sciMaxLimit or (xMax-xMin)<IMPACT_Z_sciMinLimit:
                self.subfig[i].xaxis.set_major_formatter(IMPACT_Z_SciFormatter)
            if (yMax-yMin)>IMPACT_Z_sciMaxLimit or (yMax-yMin)<IMPACT_Z_sciMinLimit:
                self.subfig[i].yaxis.set_major_formatter(IMPACT_Z_SciFormatter)
            #xmajorFormatter = FormatStrFormatter('%2.2E')
            #self.subfig[i].yaxis.set_major_formatter(xmajorFormatter)  
            
            self.subfig[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.21),fancybox=True, shadow=True, ncol=5)

        self.canvas.draw()
        
class EmitGrowthFrame(PlotBaseFrame):
    def __init__(self, parent):
        PlotBaseFrame.__init__(self, parent)
        self.plot()
    def plot(self):        
        fileList        = ['fort.24','fort.25']
        xdataList       = [1,1]
        ydataList       = [7,7]
        xyLabelList     = ['Z (m)','Avg emit growth in X and Y']
        
        lineType = ['r-','b--']
        
        try:
            fin1 = open(fileList[0],'r')
        except:
            print("  ERRPR! Can't open file '" + fileList[0] + "'")
            return
        try:
            fin2 = open(fileList[1],'r')
        except:
            print("  ERRPR! Can't open file '" + fileList[1] + "'")
            return
        linesList1  = fin1.readlines()
        linesList2  = fin2.readlines()
        fin1 .close()
        fin2 .close()
        linesList1  = [line.split() for line in linesList1 ]
        linesList2  = [line.split() for line in linesList2 ]
        xId = xdataList[0]-1
        yId = ydataList[0]-1
        try:
            x   = [float(xrt[xId]) for xrt in linesList1]
            start = (float(linesList1[0][yId]) + float(linesList2[0][yId]))/2
            if start < 1.0e-16:
                start=1.0e-16
            y   = [(float(linesList1[k][yId]) + float(linesList2[k][yId]))/2 / start -1 for k in range(len(linesList1))]
        except:
            print("  ERRPR! Can't read data '" + fileList[1] + "'")
            
        self.subfig.cla()
        self.subfig.plot(x, y, lineType[0], linewidth=2, label='emit.growth')
        self.subfig.set_xlabel(xyLabelList[0])
        self.subfig.set_ylabel(xyLabelList[1])
        self.subfig.legend()
        
        self.canvas.draw()
        
class TemperatureFrame(PlotBaseFrame):
    def __init__(self, parent):
        PlotBaseFrame.__init__(self, parent)
        self.plot()
    def plot(self):
        arg=['ct','fort.24','fort.25','fort.26']
        labelList= ['X','Y','Z']
        lineType = ['-','--',':']
        col      = ['b','g','r']
        linew    = [2,2,3]
        picNum = len(arg) - 1
        plotPath = './post'
        if os.path.exists(plotPath) == False:
            os.makedirs(plotPath)
            
        self.subfig.cla()
        for i in range(1,picNum+1):
            try:
                fin = open(arg[i],'r')
            except:
                print( "  ERRPR! Can't open file '" + arg[i] + "'")
                return
    
            linesList  = fin.readlines()
            fin .close()
            linesList  = [line.split() for line in linesList ]
            x   = [float(xrt[0]) for xrt in linesList]
            yl=4
            y   = [float(xrt[yl])*float(xrt[yl]) for xrt in linesList]
            self.subfig.plot(x, y, color = col[(i-1)],linestyle=lineType[i-1], linewidth=linew[i-1],label=labelList[i-1])
            
        self.subfig.set_xlabel('T (s)')
        self.subfig.set_ylabel('Temperature')
        self.subfig.legend()
        
        self.canvas.draw()

class PlotHighOrderBaseFrame(tk.Frame):
    ParticleDirec = {'X (mm)'    :1,
                     'Px (MC)'   :2,
                     'Y (mm)'    :3,
                     'Py (MC)'   :4,
                     'Z (deg)'    :5,
                     'Pz (MeV)'   :6}
    data = np.array([])
    def __init__(self, parent, PlotFileName):
        tk.Frame.__init__(self, parent)
        try:
            self.data = np.loadtxt(PlotFileName)
        except:
            print(( "  ERROR! Can't open file '" + PlotFileName + "'"))
            return
        
        self.data = np.transpose(self.data)
        for i in range(0,4,2):
            self.data[i] = self.data[i] * 1e3  # from m to mm
            
        self.frame_PlotParticleControl = tk.Frame(self)
        self.frame_PlotParticleControl.pack()
        
        self.label_x    = tk.Label(self.frame_PlotParticleControl, text="Direction:")
        self.label_x.pack(side='left')

        self.ppc1Value  = tk.StringVar(self.frame_PlotParticleControl,'X (mm)')
        self.ppc1       = ttk.Combobox(self.frame_PlotParticleControl,text=self.ppc1Value,
                                       width=6,
                                       values=['X (mm)', 'Px (MC)', 'Y (mm)', 'Py (MC)','Z (mm)','Pz (MC)'])
        self.ppc1.pack(fill = 'both',expand =1,side = 'left')
        
        LARGE_FONT= ("Verdana", 12)
        self.button_ppc=tk.Button(self.frame_PlotParticleControl)
        self.button_ppc["text"] = "Plot"
        self.button_ppc["foreground"] = "blue"
        self.button_ppc["bg"] = "red"
        self.button_ppc["font"] = LARGE_FONT
        self.button_ppc["command"] = self.plot
        self.button_ppc.pack(fill = 'both',expand =1,side = 'left')

        x   = 0
        y   = self.ParticleDirec[self.ppc1.get()]
        
        self.fig = Figure(figsize=(6,5), dpi=100)
        self.subfig = self.fig.add_subplot(111)
        self.subfig.scatter(self.data[x],self.data[y],s=1)
        
        xmajorFormatter = FormatStrFormatter('%2.2E')
        self.subfig.yaxis.set_major_formatter(xmajorFormatter)
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.4, box.y0, box.width, box.height])

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.plot()
        
class PlotMaxFrame(PlotHighOrderBaseFrame):
    def __init__(self, parent,ifile):    
        PlotHighOrderBaseFrame.__init__(self, parent, ifile)
        
    def plot(self):
        y   = self.ParticleDirec[self.ppc1.get()]
        
        self.subfig.cla()
        self.subfig.plot(self.data[0],self.data[y])
    
        axis_format_Z(self.data[0],self.data[y], self.subfig)
        
        self.subfig.set_xlabel('Z (m)')
        if y%2==1:
            self.subfig.set_ylabel('Max '+ self.ppc1.get())
        else:
            self.subfig.set_ylabel('Max '+ self.ppc1.get())
        self.canvas.draw()
        
        
class PlotHighorderFrame(PlotHighOrderBaseFrame):
    def __init__(self, parent,ifile):    
        PlotHighOrderBaseFrame.__init__(self, parent, ifile)
        
    def plot(self):
        y   = self.ParticleDirec[self.ppc1.get()]
        
        self.subfig.cla()
        self.subfig.plot(self.data[0],self.data[y])
        
        xmajorFormatter = FormatStrFormatter('%2.2E')
        self.subfig.yaxis.set_major_formatter(xmajorFormatter)

        self.subfig.set_xlabel('time (secs)')
        if y==1 or y ==3:
            self.subfig.set_ylabel(self.ppc1.get())
        elif y==2 or y ==4:
            self.subfig.set_ylabel(self.ppc1.get())
        elif y==5:
            self.subfig.set_ylabel('phase (degree)')
        elif y==6:
            self.subfig.set_ylabel('Energy deviation (MeV)')
        self.canvas.draw()
        
def axis_format_Z(xData,yData,subfig):
    xMax = np.max(xData)
    xMin = np.min(xData)
    yMax = np.max(yData)
    yMin = np.min(yData)
    if (xMax-xMin)>IMPACT_Z_sciMaxLimit or (xMax-xMin)<IMPACT_Z_sciMinLimit:
        subfig.xaxis.set_major_formatter(IMPACT_Z_SciFormatter)
    if (yMax-yMin)>IMPACT_Z_sciMaxLimit or (yMax-yMin)<IMPACT_Z_sciMinLimit:
        subfig.yaxis.set_major_formatter(IMPACT_Z_SciFormatter)