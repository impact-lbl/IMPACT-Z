import time,os,sys


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.stats import gaussian_kde


import tkinter as tk
from tkinter import ttk,filedialog

import numpy as np

class ParticleBaseFrame(tk.Frame):
    ParticleDirecWithUnit_T = {'X (mm)'       :0,
                               'Px (MC)'      :1,
                               'Y (mm)'       :2,
                               'Py (MC)'      :3,
                               'Z (mm)'       :4,
                               'Pz (MC)'      :5}

    ParticleDirecWithUnit_Z = {'X (mm)'       :0,
                               'Px (MC)'      :1,
                               'Y (mm)'       :2,
                               'Py (MC)'      :3,
                               'Z (deg)'      :4,
                               'Pz (MC)'      :5}
    
    ParticleDirec = {'X'       :0,
                     'Px'      :1,
                     'Y'       :2,
                     'Py'      :3,
                     'Z'       :4,
                     'Pz'      :5}
    
    sciFormatter = FormatStrFormatter('%2.1E')
    sciMaxLimit  = 99999 *2
    sciMinLimit  = 0.0001*2
    DefaultUnit_T = ['mm','MC','mm','MC','mm','MC']
    DefaultUnit_Z = ['mm','MC','mm','MC','deg','MC']
    
    data = np.array([])
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        tk.Frame.__init__(self, parent)
        try:
            self.data = np.loadtxt(PlotFileName)
        except:
            print( "ERROR! Can't open file '" + PlotFileName + "'")
            return
        
        self.data = np.transpose(self.data)
        
        if TorZ == 'ImpactZ':
            print("The position of X and Y had been multiplied by omege/c to meet the unit Conversion from ImpactZ")
            try:
                self.data[0] = self.data[0] * 1000 * scaling
            except:
                print( "Warning: Can't read the first column @ '" + PlotFileName + "'")
            try:
                self.data[2] = self.data[2] * 1000 * scaling
            except:
                print( "Warning: Can't read the third column @ '" + PlotFileName + "'")
            try:
                self.data[4] = self.data[4] /(3.1415926)*180
            except:
                print( "Warning: Can't read the fifth column @ '" + PlotFileName + "'")
        elif TorZ=='ImpactZ':
            for i in range(0,6,2):
                try:
                    self.data[i] = self.data[i] * 1000 * scaling
                except:
                    print( "Warning: Can't read the column " + str(i)+" @ '" + PlotFileName + "'")
        else:
            print("Warning: cannot recognize T or Z.")
            
        self.frame_PlotParticleControl = tk.Frame(self)
        self.frame_PlotParticleControl.pack()
        
        self.label_scalingX        = tk.Label(self.frame_PlotParticleControl, text="ScalingX:")
        self.label_scalingX.pack(side='left')
        self.scalingX       = tk.Entry(self.frame_PlotParticleControl,  width=7)
        self.scalingX.insert(0, '1.0')
        self.scalingX.pack(fill = 'both',expand =1,side = 'left')
        
        self.label_scalingY        = tk.Label(self.frame_PlotParticleControl, text="ScalingY:")
        self.label_scalingY.pack(side='left')
        self.scalingY       = tk.Entry(self.frame_PlotParticleControl,  width=7)
        self.scalingY.insert(0, '1.0')
        self.scalingY.pack(fill = 'both',expand =1,side = 'left')
        
        self.label_unitX        = tk.Label(self.frame_PlotParticleControl, text="UnitAxi1:")
        self.label_unitX.pack(side='left')
        self.unitX       = tk.Entry(self.frame_PlotParticleControl,  width=6)
        self.unitX.insert(0, 'mm')
        self.unitX.pack(fill = 'both',expand =1,side = 'left')
        
        self.label_unitY        = tk.Label(self.frame_PlotParticleControl, text="UnitAxi2:")
        self.label_unitY.pack(side='left')
        self.unitY       = tk.Entry(self.frame_PlotParticleControl,  width=6)
        self.unitY.insert(0, 'MC')
        self.unitY.pack(fill = 'both',expand =1,side = 'left')
        
        self.label_x        = tk.Label(self.frame_PlotParticleControl, text="Axi1:")
        self.label_x.pack(side='left')

        self.ppc1Value  = tk.StringVar(self.frame_PlotParticleControl,'X')
        self.ppc1       = ttk.Combobox(self.frame_PlotParticleControl,text=self.ppc1Value,
                                       width=5,
                                       values=['X', 'Px', 'Y', 'Py','Z','Pz'])
        #                               values=['X (mm)', 'Px (MC)', 'Y (mm)', 'Py (MC)','Z (deg)','Pz (MC)'])
        self.ppc1.pack(fill = 'both',expand =1,side = 'left')
        
        self.label_y        = tk.Label(self.frame_PlotParticleControl, text="Axi2:")
        self.label_y.pack(side='left')
        self.ppc2Value  = tk.StringVar(self.frame_PlotParticleControl,'Px')
        self.ppc2       = ttk.Combobox(self.frame_PlotParticleControl,text=self.ppc2Value,
                                       width=5,
                                       values=['X', 'Px', 'Y', 'Py','Z','Pz'])
        #                               values=['X (mm)', 'Px (MC)', 'Y (mm)', 'Py (MC)','Z (deg)','Pz (MC)'])
        self.ppc2.pack(fill = 'both',expand =1,side = 'left')
        
        LARGE_FONT= ("Verdana", 12)
        self.button_ppc=tk.Button(self.frame_PlotParticleControl)
        self.button_ppc["text"]         = "Plot"
        self.button_ppc["foreground"]   = "#FF0000"
        self.button_ppc["bg"]           = "#FFFF00"
        self.button_ppc["font"]         = LARGE_FONT
        self.button_ppc["command"]      = self.plot
        self.button_ppc.pack(fill = 'both',expand =1,side = 'right')
        
        self.ppc1Value.trace('w',lambda a,b,c,direc='X': self.update(direc))
        self.ppc2Value.trace('w',lambda a,b,c,direc='Y': self.update(direc))

        x   = self.ParticleDirec[self.ppc1.get()]
        y   = self.ParticleDirec[self.ppc2.get()]
        
        self.fig = Figure(figsize=(7,6), dpi=100)
        self.subfig = self.fig.add_subplot(111)
        self.subfig.scatter(self.data[x],self.data[y],s=1)
        
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.4, box.y0, box.width, box.height])
        
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
    def update(self,direction):
        if direction == 'X':
            self.scalingX.delete(0, 'end')
            self.scalingX.insert(0, '1.0')
            self.unitX.delete(0, 'end')
            try:
                ind = self.ParticleDirec[self.ppc1.get()]
                if   TorZ=='ImpactT':
                    self.unitX.insert(0, self.DefaultUnit_T[ind])
                elif TorZ=='ImpactZ':
                    self.unitX.insert(0, self.DefaultUnit_Z[ind])
                else:
                    print("Warning: cannot recognize T or Z.")
            except:
                pass
        elif direction == 'Y':
            self.scalingY.delete(0, 'end')
            self.scalingY.insert(0, '1.0')
            self.unitY.delete(0, 'end')
            try:
                ind = self.ParticleDirec[self.ppc2.get()]
                if   TorZ=='ImpactT':
                    self.unitY.insert(0, self.DefaultUnit_T[ind])
                elif TorZ=='ImpactZ':
                    self.unitY.insert(0, self.DefaultUnit_Z[ind])
                else:
                    print("Warning: cannot recognize T or Z.")
            except:
                pass
        else:
            print("Warning: no this direction")

class ParticleFrame(ParticleBaseFrame):
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        ParticleBaseFrame.__init__(self, parent,PlotFileName,scaling,TorZ)
        self.plot()
        
    def plot(self):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        yData   = self.data[self.ParticleDirec[self.ppc2.get()]] * float(self.scalingY.get())
        
        #self.fig.clf()
        #self.subfig = self.fig.add_subplot(111)
        self.subfig.cla()
        
        self.subfig.scatter(xData,yData,s=1)
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
        
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel(self.ppc2.get()+' ('+self.unitY.get()+')')
        self.canvas.draw()
        
    def quit(self):
        self.destroy()
        


class ParticleDensityFrame_weight1D(ParticleBaseFrame):
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        ParticleBaseFrame.__init__(self, parent,PlotFileName,scaling,TorZ)
        self.ppc2.pack_forget()
        self.label_y.pack_forget()
        
        self.unitY.pack_forget()
        self.label_unitY.pack_forget()
        
        self.label_scalingY.pack_forget()
        self.scalingY.pack_forget()
        
        self.plot()
        
    def plot(self):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        
        self.subfig.cla()
        
        nx = 200
        
        xMax = np.max(xData)
        xMin = np.min(xData)
        
        hx = (xMax-xMin)/(nx-1)
        
        count = np.zeros(nx)
        tickx  = [xMin + i * (xMax-xMin)/(nx-1) for i in range(nx)]
        
        for i in range(0,len(xData)):
            ix = int((xData[i] - xMin)/hx)
            if ix<0:
                ix=0
                print("Error at density plot weight 1D! ix<0")
            if ix>=nx-1:
                ix=nx-2
            ab = (xData[i] - (xMin+ix*hx))/hx
            
            count[ix  ] += 1.0-ab
            count[ix+1] += ab
            pass
        count = count/np.max(count)
        self.subfig.fill_between(tickx,0,count)#,extent=(xMin,xMax,yMin,yMax))#plt.cm.ocean)
        #plt.colorbar()
        
        xMax = np.max(xData)
        xMin = np.min(xData)
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel('Density')

        self.canvas.draw()        
class ParticleDensityFrame_weight2D(ParticleBaseFrame):
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        ParticleBaseFrame.__init__(self, parent,PlotFileName,scaling,TorZ)
        
        self.label_gridSizeX        = tk.Label(self.frame_PlotParticleControl, text="GridSize:")
        self.label_gridSizeX.pack(side='left')
        self.gridSizeX       = tk.Entry(self.frame_PlotParticleControl,  width=5)
        self.gridSizeX.insert(0, '200')
        self.gridSizeX.pack(fill = 'both',expand =1,side = 'left')
        
        '''
        self.label_gridSizeY        = tk.Label(self.frame_PlotParticleControl, text="GridSizeY:")
        self.label_gridSizeY.pack(side='left')
        self.gridSizeY       = tk.Entry(self.frame_PlotParticleControl,  width=5)
        self.gridSizeY.insert(0, '100')
        self.gridSizeY.pack(fill = 'both',expand =1,side = 'left')
        '''
        
        '''
        self.button_ppc["text"] = "ContourPlot"
        LARGE_FONT= ("Verdana", 12)
        self.button_ppc1=tk.Button(self.frame_PlotParticleControl)
        self.button_ppc1["text"] = "gridDensity"
        self.button_ppc1["foreground"] = "red"
        self.button_ppc1["bg"] = "yellow"
        self.button_ppc1["font"] = LARGE_FONT
        self.button_ppc1["command"] = lambda:self.plot(flag = 'gridDensity')
        self.button_ppc1.pack(fill = 'both',expand =1,side = 'right')
        '''
        
        self.plot()
        
    def plot(self,flag='ContourPlot'):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        yData   = self.data[self.ParticleDirec[self.ppc2.get()]] * float(self.scalingY.get())
        
        self.subfig.cla()
        
        try:
            nx=int(self.gridSizeX.get())
            ny=int(self.gridSizeX.get())
        except:
            nx=200
            ny=200
            print("Warning: cannot get gridSizeX or gridSizeY, set to 100")
        if nx<10:
            nx=10
        if ny<10:
            ny=10
        xMax = np.max(xData)
        yMax = np.max(yData)
        xMin = np.min(xData)
        yMin = np.min(yData)
        
        hx = (xMax-xMin)/(nx-1)
        hy = (yMax-yMin)/(ny-1)
        
        count = np.zeros([ny,nx])
        
        for i in range(0,len(xData)):
            if xData[i] < xMin or xData[i] > xMax:
                continue
            if yData[i] < yMin or yData[i] > yMax:
                continue
            ix = int((xData[i] - xMin)/hx)
            iy = int((yData[i] - yMin)/hy)
            if ix<0:
                ix=0
                print("Error at density plot weight 2D! ix<0")
            if iy<0:
                iy=0
                print("Error at density plot weight 2D! iy<0")
            if ix>=nx-1:
                ix=nx-2
            if iy>=ny-1:
                iy=ny-2
            ab = (xData[i] - (xMin+ix*hx))/hx
            cd = (yData[i] - (yMin+iy*hy))/hy
            
            #iy=ny-iy-2
            count[iy  ,ix  ] += (1.0-ab) * (1.0-cd)
            count[iy+1,ix  ] += (    ab) * (1.0-cd) 
            count[iy  ,ix+1] += (1.0-ab) * (    cd) 
            count[iy+1,ix+1] += (    ab) * (    cd) 
            pass
        
        count[count == 0.0] = -0.0000001
        tmap = plt.cm.jet
        tmap.set_under('white',0.)
        #tmap.set_bad('white',0.)
        if flag=='ContourPlot':
            x = np.linspace(xMin, xMax, nx)
            y = np.linspace(yMin, yMax, ny)
            #count = scipy.ndimage.zoom(count, 3)
            self.msh = self.subfig.contourf(x, y, count,level=12,interpolation='gaussian',cmap =tmap , vmin=0.0001)
        else:
            self.msh = self.subfig.imshow(count, origin = "lower", interpolation='bilinear', 
                                      cmap=tmap,vmin=0.0000001,
                                      extent=(xMin,xMax,yMin,yMax),aspect="auto")#plt.cm.ocean)
        
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        if (yMax-yMin)>self.sciMaxLimit or (yMax-yMin)<self.sciMinLimit:
            self.subfig.yaxis.set_major_formatter(self.sciFormatter)
        '''
        ntick = 7
        tickx  = [i*(nx-0)/(ntick-1) for i in range(ntick)]
        labelx = ['{:2.2e}'.format(xMin+i*(xMax-xMin)/(ntick-1)) for i in range(ntick)]
        self.subfig.set_xticks(tickx)
        self.subfig.set_xticklabels(labelx)
        
        ticky  = [i*(ny-0)/(ntick-1) for i in range(ntick)]
        labely = ['{:2.2e}'.format(yMin+i*(yMax-yMin)/(ntick-1)) for i in range(ntick)]
        self.subfig.set_yticks(ticky)
        self.subfig.set_yticklabels(labely)
'''
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel(self.ppc2.get()+' ('+self.unitY.get()+')')

        self.canvas.draw()
        
class ParticleDensityFrame1D(ParticleBaseFrame):
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        ParticleBaseFrame.__init__(self, parent,PlotFileName,scaling,TorZ)
        self.ppc2.pack_forget()
        self.label_y.pack_forget()
        
        self.label_scalingY.pack_forget()
        self.scalingY.pack_forget()
        self.plot()
        
    def plot(self):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        
        self.subfig.cla()
        self.subfig.hist(xData,bins=100)

        xMax = np.max(xData)
        xMin = np.min(xData)
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel('Density')

        self.canvas.draw()
        
class ParticleDensityFrame2D(ParticleBaseFrame):
    def __init__(self, parent,ifile,scaling,TorZ):    
        ParticleBaseFrame.__init__(self, parent, ifile, scaling,TorZ)
        self.plot()
        
    def plot(self):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        yData   = self.data[self.ParticleDirec[self.ppc2.get()]] * float(self.scalingY.get())
        
        self.subfig.cla()
        self.subfig.hist2d(xData,yData,(100, 100),cmap = 'jet')

        xMax = np.max(xData)
        xMin = np.min(xData)
        yMax = np.max(yData)
        yMin = np.min(yData)
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        if (yMax-yMin)>self.sciMaxLimit or (yMax-yMin)<self.sciMinLimit:
            self.subfig.yaxis.set_major_formatter(self.sciFormatter)
        
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel(self.ppc2.get()+' ('+self.unitY.get()+')')

        self.canvas.draw()
        
class ParticleDensityFrame2D_slow(ParticleBaseFrame):
    def __init__(self, parent, PlotFileName,scaling,TorZ):
        ParticleBaseFrame.__init__(self, parent,PlotFileName,scaling,TorZ)
        self.plot()
        
    def plot(self):
        xData   = self.data[self.ParticleDirec[self.ppc1.get()]] * float(self.scalingX.get())
        yData   = self.data[self.ParticleDirec[self.ppc2.get()]] * float(self.scalingY.get())
        
        self.subfig.cla()       
        '''Calculate the point density'''
        xy = np.vstack([xData,yData])
        z = gaussian_kde(xy)(xy)
    
        '''Sort the points by density, so that the densest points are plotted last'''
        idx = z.argsort()
        x, y, z = xData[idx], yData[idx], z[idx]        
        
        self.subfig.scatter(x, y, c=z, s=10, edgecolor='')

        xMax = np.max(xData)
        xMin = np.min(xData)
        yMax = np.max(yData)
        yMin = np.min(yData)
        if (xMax-xMin)>self.sciMaxLimit or (xMax-xMin)<self.sciMinLimit:
            self.subfig.xaxis.set_major_formatter(self.sciFormatter)
        if (yMax-yMin)>self.sciMaxLimit or (yMax-yMin)<self.sciMinLimit:
            self.subfig.yaxis.set_major_formatter(self.sciFormatter)
        
        self.subfig.set_xlabel(self.ppc1.get()+' ('+self.unitX.get()+')')
        self.subfig.set_ylabel(self.ppc2.get()+' ('+self.unitY.get()+')')

        self.canvas.draw()