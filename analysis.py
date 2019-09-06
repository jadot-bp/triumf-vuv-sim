#   VUV RAYTRACING SCRIPT FOR TRIUMF
#   Contact: Ben Page (benjaminpage.acer@gmail.com)
#   Last Updated: Sep 4 2019

import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import pandas as p
import cPickle
import sys
import os
import re

try:
    from chroma.geometry import Solid
    from chroma.demo.optics import vacuum
    from chroma import make,view
    import builder
    import mfunctions
except:
    if os.path.exists('/singularity'):
        raise Exception('Chroma container is loaded, but it seems that builder.py and simulator.py cannot be loaded! Are they present?\nPlease consult the documentation.')
        exit()
    else:
        print 'Chroma container not loaded! To avoid errors, please enter the container before running the script!\n'
        expression = 'singularity shell --nv chroma.img'
        print 'REMINDER: The command is: ',expression
        exit()


def main(suppress=False):
    
    ##### Loading in Results #####
    print '-'*90
    print 'ANALYSIS SCRIPT FOR VUV SIMULATION'.center(90,'=')
    print '-'*90
    print 'Available files:\n'
    filenames = []
    display = ['#','Date:','Time:','n:',u'\u03bb (nm)'.encode('utf-8'),u'd (\u03bcm)'.encode('utf-8'),u's (\u03bcm)'.encode('utf-8')]
    print '{:<3}{:<12}{:<9}{:<9}{:<9}{:<9}{:<9}'.format(*display)
 
    for pos,filename in enumerate(os.listdir('results/')):
        filenames.append(filename)
        date = re.findall('\d+\-\d+\-\d+',filename)[0]
        time = re.findall('\d{2}\:\d{2}',filename)[0]
        args = re.findall('\d+\:\d+:\d+\.\d+\:\d+\.\d+',filename)[0]
        display = [str(pos),date,time] + args.split(':')
        print '{:<3}{:<12}{:<9}{:<9}{:<8}{:<8}{:<9}'.format(*display)
    print '(d - Slit Width, s - Beam Width)'
    print '='*90

    proceed = False
    while proceed == False:
        usrin = raw_input('Please select a file: ').strip()
        if usrin not in [str(i) for i in list(range(len(filenames)))]:
            print 'Bad input!'
        else:
            proceed = True

    selection = filenames[int(usrin)] 
    n,wavelength,width,beamsize = [float(i) for i in re.findall('\d+\:\d+:\d+\.\d+\:\d+\.\d+',selection)[0].split(':')]

    #print '\nSelected: '+ selection + '\n'

    ##### Data Loading and Results #####

    itm = 25.4 #Inches to millimetre conversion
    tol = 0.05  #Tolerance for selecting photons 
    res = 10 #Seed resolution
    seed = np.genfromtxt('config/seed.dat',delimiter=',') #display seed

    x = 0.01    #marker params
    h = x*np.sqrt(3)/2
    pyramid = make.linear_extrude([-x/2,0,x/2],[-h/2,h/2,-h/2],h,[0]*3,[0]*3)

    beg_marker = Solid(pyramid,vacuum,vacuum,color=0x32a8a2)    #photon start marker
    end_marker = Solid(pyramid,vacuum,vacuum,color=0xfc0303)    #photon end marker
    dir_marker = Solid(pyramid,vacuum,vacuum,color=0x00ff00) #Direction

    data = p.read_csv('results/'+selection,sep = '\t').values[:]
    
    photon_ids = data[:,0].astype(int)
    beg_pos = data[:,1:4]
    end_pos = data[:,4:7]
    wavelengths = data[:,7].astype(int)
    flags = data[:,8].astype(int)

    if suppress == False:
        if os.path.exists('config/geometry.pickle'):
            print 'Loading geometry from file...'
            geometry_path = open('config/geometry.pickle','rb')
            world = cPickle.load(geometry_path)
            world.add_solid(dir_marker,displacement = mfunctions.get_center('161')+np.array((0,1,0)))
            geometry_path.close()
        else:
            world = builder.build()   #Regenerate geometry
   
    sipm_ids = []
    pmt_ids = []
    for p_id in photon_ids:
        if (flags[p_id] & 0x1 << 2) == 4:
            if abs(end_pos[p_id,2]-9.08933)<tol:   #Photons detected at SiPM
                sipm_ids.append(p_id)
            elif abs(end_pos[p_id,0]-23.6555999)<tol:  #Photons detected at PMT
                pmt_ids.append(p_id)
        if suppress == False and seed[p_id] <= 1.0/res:
            world.add_solid(beg_marker,displacement=beg_pos[p_id])
            world.add_solid(end_marker,displacement=end_pos[p_id])
            
    pcount = len(pmt_ids)
    scount = len(sipm_ids)
        
    print 'Total photons:\t\t', n
    print 'Detections (PMT):\t', pcount
    print 'Detections (SiPM):\t',scount
    print 'Relative detection rate: {:.2f}%'.format(100*float(scount+pcount)/len(photon_ids))
    
    if suppress == False:
        view(world)
    elif suppress == True and scount == 0:
        print 'No photons detected!'
        exit()
    else:
        detectedx = np.asarray(list(end_pos[i,0]*itm for i in sipm_ids)) #Compiling points to plot and converting to mm
        detectedy = np.asarray(list(end_pos[i,1]*itm for i in sipm_ids))
        detectedwavs = np.array(list(wavelengths[i] for i in pmt_ids))
        
        meanx = np.mean(detectedx) #Mean detected positions
        meany = np.mean(detectedy)

        detectedx -= meanx #Centering plot about detected photons
        detectedy -= meany
      
        sampleset = np.sqrt(detectedx**2 + detectedy**2) #Heatmap info        
        binwidth = 2*ss.iqr(sampleset)/np.cbrt(len(sampleset)) #Freedman-Diaconis rule
        bins = (max(sampleset)-min(sampleset))/binwidth
        heatmap,xedges,yedges = np.histogram2d(detectedx,detectedy,bins=bins)
        extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
       
        fig = plt.figure()
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.axis('equal')
        ax2.axis('equal')
        fig.suptitle(u'n: {} | \u03bb: {}nm | d: {}\u03bcm | s: {}\u03bcm'.format(int(n),int(wavelength),width,beamsize))
        
        r = 0.059813*itm #Iris radius
        theta = np.linspace(0,2*np.pi,1000)
        iris_center = np.array((18.4429,-21.8872))*itm

        iris_x = r*np.cos(theta)+iris_center[0]-meanx
        iris_y = r*np.sin(theta)+iris_center[1]-meany

        ax1.scatter(detectedx,detectedy,s=0.3,c='k',label='Photon Hits')
        #ax1.plot(iris_x,iris_y,c='r',label='Iris Overlay')
        ax1.set_xlabel('x (mm)')
        ax1.set_ylabel('y (mm)')
        ax1.set_title('SiPM Hits')
        ax1.legend()

        im = ax2.imshow(heatmap.T,extent=extent,origin='lower',cmap='hot_r',interpolation='gaussian')
        ax2.set_xlabel('x (mm)')
        ax2.set_ylabel('y (mm)')
        ax2.set_title('SiPM Heat Map')
                
        fig.colorbar(im,ax=ax2) 
        plt.show()


if __name__ == '__main__': 
    try:
        flag = sys.argv[1]
    except:
        flag = None

    if flag == '--help':
        print "For scatter plot mode, the flag is '1'. To view a rendering of the results, the flag is '0'.\nFor more information, please consult the user manual"
        exit()
    elif flag in ['1','0']:
        main(bool(int(flag)))
    else:
        print 'Invalid Flag!'
        print "For scatter plot mode, the flag is '1'. To view a rendering of the results, the flag is '0'.\nFor more information, please consult the user manual."
        exit()
