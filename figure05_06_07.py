
# Split by 2 figures: Phy O2 
# and DON PON 

'''
Created on 28. jun. 2017

@author: ELP
'''


import os

from netCDF4 import Dataset,num2date,date2num,date2index
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.dates as mdates
import datetime 

import tkinter as tk # python3
from tkinter.filedialog import askopenfilename,askdirectory   # python 3
root = tk.Tk()
root.withdraw()
#plt.style.use('ggplot')
#plt.style.use('bmh')

from matplotlib import rc
font = {'size' : 12}
rc('font', **font)
textsize = 13
x_text = 0.05
y_text = 1.2
labels = ('A) ','B) ','C) ')
import sys
from matplotlib.backends.backend_qt5agg import (
FigureCanvasQTAgg as FigureCanvas)
    
directory =  askdirectory() 

ice_fname = os.path.abspath(os.path.join(directory,'ice.nc')) 
water_fname = os.path.join(directory,'water.nc')
sediments_fname = os.path.join(directory,'sediments.nc')
   
fh_ice =  Dataset(ice_fname)  
fh_water =  Dataset(water_fname)  
fh_sediments =  Dataset(sediments_fname) 

depth = fh_ice.variables['z'][:] 
depth_faces = fh_ice.variables['z_faces'][:] 
   
max_faces = np.max(depth_faces)

min_ice = np.min(depth)
max_ice = np.max(depth)

ice_res = 6
depth_ice = (depth - max_ice)*-ice_res
depth_ice_faces = np.array((depth_faces - max_faces)*-ice_res)
depth_ice = np.array(depth_ice)

min_ice = np.amin(depth_ice_faces)
max_ice =  250 #np.amax(depth_ice_faces)

depth_water = np.array(fh_water.variables['z_faces'][:])
depth_sed = fh_sediments.variables['z_faces'][:] 

min_water = np.amin(depth_water)
max_water = np.amax(depth_water[:])

min_sed = np.amin(depth_sed)
max_sed = np.amax(depth_sed)
         
#########################
# Values for time axis  #
#########################

time = fh_ice.variables['time']      
time2 = fh_ice.variables['time'][:]
time_units = fh_ice.variables['time'].units
           
to_start = datetime.datetime(1982,1,1,12,0)
to_stop= datetime.datetime(1984,1,1,12,0)

start = date2index(to_start, time,#units = time_units,
                    calendar=None, select='nearest')
stop = date2index(to_stop, time,#units = time_units,
                    calendar=None, select='nearest')

#  

time = fh_ice.variables['time']      
time2 = fh_ice.variables['time'][start:stop]



def read_var(name): 
    var_ice =  ma.masked_invalid (
        np.array(fh_ice.variables[name][start:stop]).T )
    var_water = np.array(fh_water.variables[name][start:stop]).T 
    var_sediments = np.array(fh_sediments.variables[name][start:stop]).T                     
    data_units = fh_ice.variables[name].units                           
    return var_ice,var_water,var_sediments,data_units


def plot(var,axis0,axis1,axis2,tit): #ice,water,sed
    ticks_dict = {'B_BIO_DON':([0,1,2,3,4],[10,15,20,25,30],[0,10,20,30,40,50]),
                  'B_BIO_PON':([0,1,2,3,4],[0,1,2,3,4,5,6,10,15,20,25,50],[0,50,100,150,200]),
                  'B_BIO_O2':([0,10,20,30,40,50],[0,50,100,150,200,250,300,350],[0,25,50,75,100,125,150]),
                  'B_S_H2S':([0,0.1,0.2,0.3],[0,0.01,0.02,0.03,0.04],[0,0.01,0.02,0.03]),
                  'B_Fe_Fe2':([0,0.005,0.01,0.015,0.02],[0,0.01,0.02,0.03],[0,0.25,0.5,0.75,1]),
                  'B_Mn_Mn2':([0,0.01,0.02,0.03,0.04,0.05],[0,0.2,0.4,0.6,0.8,1],[0,2.5,5,7.5,10,12.5,15]),
                  'B_NUT_Si':([0,0.5,1,1.5,2],[5,10,15],[5,10,15]),
                  'B_NUT_NO3':([0,1,2,3,4,5],[0,5,10,15,20,25],[0,5,10,15,20,25]),
                  'B_NUT_NH4':([0,0.25,0.5,0.75,1,1.25,1.5],[0,1,2,3,4,5],[0,10,20,30,40,50]),
                  }

    ### Time ticks ### 
    from dateutil.relativedelta import relativedelta
    
    def add_colorbar(CS,axis,ma1,t):
        if ma1 > 10000 or ma1 < 0.001:
            cb = plt.colorbar(CS,ax = axis,pad=0.02,
                 aspect = 7,shrink = 0.9,format=ticker.FuncFormatter(fmt),ticks = t) 
        else: 
            cb = plt.colorbar(CS,ax = axis,pad=0.01,
                 aspect = 7,shrink = 0.9,ticks = t) 
        return cb    
    
    
    cmap = plt.get_cmap('viridis') 
    cmap_water = plt.get_cmap('CMRmap') 
    data = read_var(var)
    
    var_ice = data[0]
    var_water = data[1]

    var_sed = data[2]
    data_units = data[3]
                             
    X,Y = np.meshgrid(time2[:],depth_ice_faces)
    X  = num2date(X,units = time_units)  
    #start_f = num2date(time2[start],units = time_units) 
    #stop_f = num2date(time2[stop],units = time_units) 
    
    X_water,Y_water = np.meshgrid(time2[:],depth_water[5:])
    X_water = num2date(X_water,units = time_units)   
    
    X_sed, Y_sed = np.meshgrid(time2[:],depth_sed)
    X_sed = num2date(X_sed,units = time_units)   
    #X_sed = format_time
    #fh_ice.close()
    #fh_water.close()
    #fh_sediments.close()          
    CS1 = axis0.pcolor(X,Y,var_ice[:,:],cmap = cmap )       
    CS4 = axis1.pcolor(X_water,Y_water,var_water[5:,:],
                      cmap = cmap_water)
    CS7 = axis2.pcolor(X_sed,Y_sed,var_sed[:,:], cmap = cmap)  
    
    ma1 = ma.max(var_ice[:,:])  
     
    cb0 = add_colorbar(CS1,axis0,ma1,ticks_dict[var][0])
    cb1 = add_colorbar(CS4,axis1,ma1,ticks_dict[var][1])
    cb2 = add_colorbar(CS7,axis2,ma1,ticks_dict[var][2])
    
    axis2.axhline(max_water, color='w', linestyle = '--',linewidth = 1 ) 
    axis2.annotate('  Sediment Water Interface',
            xy =(start,max_water),
            xytext=(start,max_water-0.01),color = 'w')
    axis0.set_title(tit +', $\mu M$')
    axis0.set_ylim(min_ice,max_ice)
    
    axis1.set_ylim(max_water,min_water)
    axis2.set_ylim(max_sed,min_sed)  #

    # hide horizontal axis labels 
    axis0.set_xticklabels([])    
    axis1.set_xticklabels([])
    axis2.set_xticklabels([])    
    #letters = ['(a)','(b)','(c)']
    labels = ["Ice thickness, cm", "Depth, m","Depth, m" ]
    n = 1
    ice_ticks = np.arange(50,max_ice,50)
    
    axis0.set_yticks(ice_ticks)
    
    for axis in (axis1,axis2): 
        axis.yaxis.set_label_coords(-0.08, 0.45)
        axis.set_ylabel(labels[n])  
        n=n+1
    axis0.yaxis.set_label_coords(-0.08, 0.6)   
    axis0.set_ylabel(labels[0]) 
     
    axis2.annotate('  Sediment Water Interface',
            xy =(to_start,max_water),
            xytext=(to_start,max_water-0.01),color = 'w')   
    if var in ('B_BIO_O2','B_Mn_Mn2','B_NUT_NH4'):
        if (stop-start)>= 365*6:            
            axis2.xaxis.set_major_formatter(
                mdates.DateFormatter('%Y')) 
                #ticks = np.arange(time[start:stop],time[start:stop],50)
        elif (stop-start) > 367 and (stop-start) < 365*6:
            axis2.xaxis.set_major_formatter(
                mdates.DateFormatter("%b '%y"))   
        else : 
            axis2.xaxis.set_major_formatter(
                mdates.DateFormatter('%b'))    
        import matplotlib.ticker as ticker
        
        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)


def figure(numfig):
    if numfig == 5:
        figure = plt.figure(figsize=(8.27, 11.69), dpi=100,
                        facecolor='None',edgecolor='None')   
        gs0 = gridspec.GridSpec(3, 1)
        gs0.update(left=0.09, right= 1,top = 0.97,bottom = 0.05,
                           wspace=0.1,hspace=0.1)
        dy = 0.04
        
        gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0],hspace=dy)
        gs_1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],hspace=dy)
        gs_2 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2],hspace=dy)
        
        #add subplots

        ax0 = figure.add_subplot(gs[0]) # ice 
        ax1 = figure.add_subplot(gs[1]) # water
        ax2 = figure.add_subplot(gs[2]) # sed
        
        ax0_1 = figure.add_subplot(gs_1[0]) # ice       
        ax1_1 = figure.add_subplot(gs_1[1]) # water
        ax2_1 = figure.add_subplot(gs_1[2]) # sed
        
        ax0_2 = figure.add_subplot(gs_2[0]) # ice        
        ax1_2 = figure.add_subplot(gs_2[1]) # water
        ax2_2 = figure.add_subplot(gs_2[2]) # sed

        for i,axis in enumerate((ax0,ax0_1,ax0_2)):
            axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
                
        plot('B_BIO_DON',ax0,ax1,ax2,'DON')              
        plot('B_BIO_PON',ax0_1,ax1_1,ax2_1,'PON')  
        plot('B_BIO_O2',ax0_2,ax1_2,ax2_2,'O$_2$')
        #plt.show()
        plt.savefig('Figure05.pdf')
        plt.clf()
        
        
    elif numfig == 7:
        
        figure = plt.figure(figsize=(8.27, 11.69), dpi=100,
                        facecolor='None',edgecolor='None')   
        gs0 = gridspec.GridSpec(3, 1)
        gs0.update(left=0.09, right= 1,top = 0.97,bottom = 0.05,
                           wspace=0.1,hspace=0.1)
        dy = 0.04
        
        gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0],hspace=dy)
        gs_1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],hspace=dy)
        gs_2 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2],hspace=dy)
        
        #add subplots
        ax0 = figure.add_subplot(gs[0]) # ice 
        ax1 = figure.add_subplot(gs[1]) # water
        ax2 = figure.add_subplot(gs[2]) # sed
        
        ax0_1 = figure.add_subplot(gs_1[0]) # ice 
        ax1_1 = figure.add_subplot(gs_1[1]) # water
        ax2_1 = figure.add_subplot(gs_1[2]) # sed
        
        ax0_2 = figure.add_subplot(gs_2[0]) # ice 
        ax1_2 = figure.add_subplot(gs_2[1]) # water
        ax2_2 = figure.add_subplot(gs_2[2]) # sed
        
        for i,axis in enumerate((ax0,ax0_1,ax0_2)):
            axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')
                            
        plot('B_S_H2S',ax0,ax1,ax2,'H$_2$S')              
        plot('B_Fe_Fe2',ax0_1,ax1_1,ax2_1,'Fe(II)')  
        plot('B_Mn_Mn2',ax0_2,ax1_2,ax2_2,'Mn(II)')
        plt.savefig('Figure07.pdf')
        #plt.show()
        plt.clf()
        
    elif numfig == 6:
        
        figure = plt.figure(figsize=(8.27, 11.69), dpi=100,
                        facecolor='None',edgecolor='None')   
        gs0 = gridspec.GridSpec(3, 1)
        gs0.update(left=0.09, right= 1,top = 0.97,bottom = 0.05,
                           wspace=0.1,hspace=0.1)
        dy = 0.04
        
        gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[0],hspace=dy)
        gs_1 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[1],hspace=dy)
        gs_2 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs0[2],hspace=dy)
        
        #add subplots
        ax0 = figure.add_subplot(gs[0]) # ice 
        ax1 = figure.add_subplot(gs[1]) # water
        ax2 = figure.add_subplot(gs[2]) # sed
        
        ax0_1 = figure.add_subplot(gs_1[0]) # ice 
        ax1_1 = figure.add_subplot(gs_1[1]) # water
        ax2_1 = figure.add_subplot(gs_1[2]) # sed
        
        ax0_2 = figure.add_subplot(gs_2[0]) # ice 
        ax1_2 = figure.add_subplot(gs_2[1]) # water
        ax2_2 = figure.add_subplot(gs_2[2]) # sed
        for i,axis in enumerate((ax0,ax0_1,ax0_2)):
            axis.text(x_text, y_text, labels[i], transform=axis.transAxes,
                 fontsize=14, fontweight='bold', va='top', ha='right')                
        plot('B_NUT_Si',ax0,ax1,ax2,'Si')              
        plot('B_NUT_NO3',ax0_1,ax1_1,ax2_1,'NO$_3^-$')  
        plot('B_NUT_NH4',ax0_2,ax1_2,ax2_2,'NH$_4^+$')
        plt.savefig('Figure06.pdf')
        #plt.show()
        plt.clf()
        
        
        
        
#### specify figure here
figure(7)
figure(6)
figure(5)
#plt.show()
#
   

    