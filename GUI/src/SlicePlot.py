import tkinter as tk
from tkinter import ttk,filedialog

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import time,os,sys
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 

import numpy as np

class SliceBaseFrame(tk.Frame):
    SlicePlotDirec={'Number of macroparticle per cell':             1,
                    'Current':                                      2,
                    'X slice emmittance (m-rad)':                   3,
                    'Y slice emmittance (m-rad)':                   4,
                    'Energy spread per cell - Correlated (eV)':     5,
                    'Energy spread per cell - Uncorrelated (eV)':   6}
    
    sciFormatter = FormatStrFormatter('%2.1E')
    sciMaxLimit  = 99999 *2
    sciMinLimit  = 0.0001*2
    data = np.array([])
    
    def __init__(self, parent, PlotFileName):
        tk.Frame.__init__(self, parent)
        try:
            self.data = np.loadtxt(PlotFileName)
        except:
            print( "ERROR! Can't open file '" + PlotFileName + "'")
            return
        
        self.data = np.transpose(self.data)

        self.frame_PlotControl = tk.Frame(self)
        self.frame_PlotControl.pack()
        
        self.label_y    = tk.Label(self.frame_PlotControl, text='Axi:')
        self.label_y.pack(side='left')
        self.ppc2Value  = tk.StringVar(self.frame_PlotControl,'Current')
        self.ppc2       = ttk.Combobox(self.frame_PlotControl,text=self.ppc2Value,
                                       width=40,
                                       values=['Number of macroparticle per cell',
                                               'Current',
                                               'X slice emittance (m-rad)',
                                               'Y slice emittance (m-rad)',
                                               'Energy spread per cell - Correlated (eV)',
                                               'Energy spread per cell - Uncorrelated (eV)'])
        self.ppc2.pack(fill = 'both',expand =1,side = 'left')
        
        LARGE_FONT= ("Verdana", 12)
        self.button_ppc=tk.Button(self.frame_PlotControl)
        self.button_ppc["text"]         = "Plot"
        self.button_ppc["foreground"]   = "#FF0000"
        self.button_ppc["bg"]           = "#FFFF00"
        self.button_ppc["font"]         = LARGE_FONT
        self.button_ppc["command"]      = self.plot
        self.button_ppc.pack(fill = 'both',expand =1,side = 'right')
        
        x   = 0
        try:
            y   = self.SlicePlotDirec[self.ppc2.get()]
        except:
            print("SlicePlot Direction Error!")
            print("No "+self.ppc2.get() + "or " + y +" colume doesn't exist")
            self.quit()
            
        
        self.fig = Figure(figsize=(7,6), dpi=100)
        self.subfig = self.fig.add_subplot(111)
        self.subfig.plot(self.data[x],self.data[y])
        
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.4, box.y0, box.width, box.height])
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        self.plot()
        
    def plot(self):
        xData   = self.data[0]
        yData   = self.data[self.SlicePlotDirec[self.ppc2.get()]]

        self.subfig.cla()
        
        self.subfig.plot(xData,yData)
        #self.subfig.relim()
        self.subfig.autoscale()
        
        xMax = np.max(xData)
        xMin = np.min(xData)
        yMax = np.max(yData)
        yMin = np.min(yData)
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        if (yMax-yMin)>self.sciMaxLimit or (yMax-yMin)<self.sciMinLimit:
            self.subfig.yaxis.set_major_formatter(self.sciFormatter)
        
        self.subfig.set_xlabel('bunch length (m)')
        self.subfig.set_ylabel(self.ppc2.get())
        self.canvas.draw()
        
    def quit(self):
        self.destroy()
        