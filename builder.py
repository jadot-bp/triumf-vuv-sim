#   VUV RAYTRACING SCRIPT FOR TRIUMF
#   Contact: Ben Page (benjaminpage.acer@gmail.com)
#   Last Updated: Aug 18 2019

"""
This script builds the setup geometry for use in the photon simulation script (simulator.py). Please read the associated documentation for further info and/or contact the email address listed.
"""

from chroma import make,view,mesh_from_stl
from chroma.geometry import Solid,Geometry
from chroma.demo.optics import water,vacuum,black_surface,r7081hqe_photocathode,shiny_surface
import numpy as np
import cPickle
import os
import re
import mfunctions
import sys

def build(slit_width):
    
    
    directory = 'stls/' #directory containing individual detector .stls
    
    
    ##### Generating World #####
    print 'Generating world...'
    size = (40,40,40) #Bounds for world box
    world = Geometry(vacuum)
    bounds = Solid(make.box(*size),vacuum,vacuum,color=0x33ffffff)
    world.add_solid(bounds,displacement=(18,-18,10))

    ##### Generating Non-Special Solids #####
    print 'Generating non-special solids...'
    blacklist = mfunctions.get_blacklist()  #Mesh read-in blacklist
    tol = 0.005 #Tolerance for painting meshes using conditionals

    setup_solid = Solid(make.cube(1),vacuum,vacuum,color=0xf0fc03) #origin

    for filename in os.listdir(directory):
        path = os.path.join(directory,filename)
        tmp_mesh = mesh_from_stl(path)
        mesh_no = re.findall('\d+',path)[0] #Select mesh number
                   
        if mesh_no in ['141','142']:    #Angled mirrors
            setup_surf, setup_col = mfunctions.mirror(mesh_no)
            """
        elif mesh_no in ['161','206']:  #Slit housing
            setup_surf = shiny_surface
            setup_col = 0x7d6b56 #grey
            """
        else:   #Non-optical components
            setup_surf = black_surface
            setup_col = 0x7d6b56 #grey
        #print 'Adding solid '+mesh_no+'/249' 
        if mesh_no not in blacklist:    #Remove blacklisted meshes from setup
            setup_solid += Solid(tmp_mesh,vacuum,vacuum,surface=setup_surf,color=setup_col)

    world.add_solid(setup_solid)

    ##### Generating Special Solids #####
    print 'Generating special solids...'
    """
    #Laser Absorb
    laser_absorb = Solid(make.segmented_cylinder(1,0.1),vacuum,vacuum,surface=black_surface,color=0xffff00) 
    cap_center = np.array((18.470600,-26.771999,19.469999))
    world.add_solid(laser_absorb,displacement=cap_center)
    """
    #SiPM Plate
    sipm_plate = Solid(mesh_from_stl(directory+'vuv_setup - Mesh 188.stl'),vacuum,vacuum,surface=black_surface,color=0xfc7f03)
    rotation_matrix = mfunctions.Ry(np.arctan(0.0755832))   #Correction for SiPM plate rotation in the STL
    sipm_center = np.array((18.48695,-22.89630,9.08933)) #Centre of SiPM plate
    correction = sipm_center-np.matmul(rotation_matrix,sipm_center) #Off-centre rotation induces displacement, this corrects for this
    world.add_solid(sipm_plate,rotation=rotation_matrix,displacement=correction+np.array((0,1,0)))

    #Detector Plate/SiPM
    detector = Solid(make.box(2,3,0.01),vacuum,vacuum,surface=mfunctions.get_sipm_surface(),color=0x0dff00) 
    world.add_solid(detector,displacement=sipm_center+np.array((0,0.8,0)))

    #Closed Mirror 
    #mirror_surf,mirror_col = mfunctions.mirror('250')
    mirror_surf = mfunctions.get_mirror_surface()
    mirror_col = 0xe6ad12 
    closed_mirror = Solid(mesh_from_stl(directory+'vuv_setup - Mesh 250.stl'),vacuum,vacuum,surface=mirror_surf,color=mirror_col)
    world.add_solid(closed_mirror,rotation=mfunctions.Ry(np.pi),displacement=np.array((18.470901,-21.916050,19.602501))+np.array((0,0,0.613293)))
    
    #Post-Monochromator Slit Doors
    basex = np.array((0.0,0.0,0.074744,0.6,0.6))#x-points for door base
    basey = np.array((0.0,0.00488564,0.0800017,0.0800017,0.0))#y-points
    height = 0.65
    #Max slit width: 1000micron
    door1 = make.linear_extrude(basex,basey,height)
    door2 = make.linear_extrude(-basex,basey,height)
    slit_center = mfunctions.get_center('161')
    door1_solid = Solid(door1,vacuum,vacuum,surface=shiny_surface,color=0xffff00)
    door2_solid = Solid(door2,vacuum,vacuum,surface=shiny_surface,color=0xffff00)
    world.add_solid(door1_solid,displacement=slit_center+np.array((0.5*slit_width,-0.144076,-0.05)))
    world.add_solid(door2_solid,displacement=slit_center+np.array((-0.5*slit_width,-0.144076,-0.05)))
    
    #PMT plate
    pmt_center = mfunctions.get_center('210')
    pmt_center += np.array((0.4,0,0))
    pmt_solid = Solid(make.segmented_cylinder(1.5,0.1),vacuum,vacuum,surface=mfunctions.get_pmt_surface(),color=0x0dff00)
    world.add_solid(pmt_solid,rotation=mfunctions.Rz(np.pi/2),displacement=pmt_center)

    ##### Saving Geometry to File #####

    config_path = open('config/geometry.pickle','wb')
    cPickle.dump(world,config_path)   
    config_path.close() 
    config_path = open('config/vars.txt','wb')
    config_path.write(str(slit_width)+'\n')
    config_path.write('This is the config file for geometry.pickle.\nAltering this file will cause geometry.pickle to regenerate!')
    config_path.close()
    print 'Done.'
    
    return world

if __name__ == '__main__':
    try: 
        slit_width = float(sys.argv[1])
    except:
        slit_width = raw_input('Please enter the slit width in inches: ')
    world = build(slit_width)
    view(world)
