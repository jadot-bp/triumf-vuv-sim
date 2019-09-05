#   VUV RAYTRACING SCRIPT FOR TRIUMF
#   Contact: Ben Page (benjaminpage.acer@gmail.com)
#   Last Updated: Aug 8 2019

"""
This script handles the creation and propagation of photons through the setup geometry defined by the builder script (builder.py). Please read the associated documentation for further info and/or contact the email address listed.
"""

from builder import build
from chroma.sim import Simulation
from chroma.event import Photons
from chroma.loader import load_bvh
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import mfunctions
import sys
import os
import cPickle


def main(n,wavelength,width,beamsize,regen=False):

    config_path = open('config/vars.txt','rb')
    test = config_path.readline()
    config_path.close()
    if not(abs(float(test)-width) < 0.1/25400):
        regen = True
    if os.path.exists('config/geometry.pickle') and regen == False:
        geometry_path = open('config/geometry.pickle','rb')
        world = cPickle.load(geometry_path)
        geometry_path.close()
    else:
         world = build(width)
    
    world.flatten()
    world.bvh = load_bvh(world) #Create bounding volume hierarchy
    sim = Simulation(world)

    def init_source(n,wavelength,width,beamsize):
            """Generates laser profile of n photons, with wavelength, emanating from position pos_offset with direction 'direction'."""
           
            pos, direction = mfunctions.get_source(n,width,beamsize,wavelength/25400000.0) 
            #Note: Chroma only natively handles integer wavelengths, so we convert to inches here instead of earlier in the code.
            source_center = mfunctions.get_center('161') + np.array((0,-0.1,0)) #Position source just in front of slits
            pos = pos + np.tile(source_center,(n,1))
            pol = np.cross(direction,(0,0,1))   #Polarisation
            wavelengths = np.repeat(wavelength,n)
            return Photons(pos,direction,pol,wavelengths)

    start = []
    end = []
    
    print 'Simulating photons...'

    for ev in sim.simulate([init_source(n,wavelength,width,beamsize)],keep_photons_beg=True,keep_photons_end=True,run_daq=False,max_steps=100):
        #print ev.photons_end.flags

        start.append(ev.photons_beg.pos)
        end.append(ev.photons_end.pos)
 
    print 'Saving data to file...'

    photon_id = list(range(n))
    flags = ev.photons_end.flags
    wavs = ev.photons_end.wavelengths    

    ##### Saving data to file #####
        
    current = dt.datetime.now() 
    filename = 'results/'+current.strftime('%Y-%m-%d_%H:%M')+'_%d:%d:%.2f:%.2f.dat' % (n,wavelength,width*25400,beamsize*25400)
    out_file = open(filename,'w')
    out_file.write('\t'.join(['ID','xi','yi','zi','xf','yf','zf','wavelength (nm)','flag\n']))
    
    for i in range(n):
        output = [photon_id[i]] + [item for item in start[0][i]] + [item for item in end[0][i]] + [int(wavs[i]),flags[i]]
        out_file.write('\t'.join([str(item) for item in output]) + '\n')
    
    out_file.close()
    
    print 'Generating seed file...'
    seed = np.random.uniform(size=n)    #This is the seed file for the analysis script
    np.savetxt('config/seed.dat',seed,delimiter=',')
    print 'Done.'

if __name__ == '__main__':
        
    try:
        n = int(sys.argv[1])
        wavelength = int(sys.argv[2])
        width = float(sys.argv[3])
        beamsize = float(sys.argv[4])
        regen = bool(sys.argv[5])
        main(n,wavelength,width,beamsize,regen)
    except IndexError:
        print 'Error! Arguments should be of the form:'
        print 'n - Number of photons (int)\nradius - Radius of beam profile (float)\nhwhm - HWHM of beam profile (float)\nregen - Regenerate geometry (bool)'

