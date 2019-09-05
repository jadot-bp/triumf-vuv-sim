#   Functions for simulation script
#   Contact: Ben Page (benjamnpage.acer@gmail.com)
#   Last Updated: Jul 31 2019

from chroma import mesh_from_stl
from chroma.transform import rotate
from chroma.geometry import Solid,Geometry,Surface
from chroma.demo.optics import water,vacuum,black_surface,shiny_surface
import numpy as np
import matplotlib.pyplot as plt
import cPickle
import sys
import os

def get_pmt_surface():
    data = np.genfromtxt('config/qe_curve.dat',delimiter=',')
    wavelength = data[:,0]
    efficiency = data[:,1]
    
    r8486_pmt = Surface('R8486 PMT')
    r8486_pmt.set('detect',efficiency,wavelength)
    return r8486_pmt

def get_sipm_surface():
    wavelength = np.arange(100,600,0.5)
    hamamatsu_sipm = Surface('Hamamatsu SiPM')
    hamamatsu_sipm.set('detect',np.ones(len(wavelength)),wavelength)
    return hamamatsu_sipm

def get_mirror_surface():
    data = np.genfromtxt('config/gold_reflectivity.dat',delimiter=',')
    wavelength = data[:,0]
    reflectivity = data[:,1]/100
    specular_percent = 0.995 #This is the % of reflectance that is specular
    diffuse_percent = 1-specular_percent
    gold_mirror = Surface('EO mirror')
    gold_mirror.set('absorb',1-reflectivity,wavelength)
    gold_mirror.set('reflect_specular',specular_percent*reflectivity,wavelength)
    gold_mirror.set('reflect_diffuse',diffuse_percent*reflectivity,wavelength)
    return gold_mirror

def get_blacklist():
    """Returns the blacklisted mesh numbers from blacklist.dat"""
    items = np.genfromtxt('config/blacklist.dat',dtype='string',delimiter=',',skip_header=True)
    blacklist = []
    for pos,item in enumerate(items):
        if item[3] == '1':
            blacklist.append(str(item[0]))
    return blacklist

def get_center(mesh_no):
    """Returns the cartesian centre of a mesh."""
    mesh = mesh_from_stl('stls/vuv_setup - Mesh %s.stl' % mesh_no)
    array = mesh.get_triangle_centers()
    x,y,z = np.hsplit(array,3)
    xc = (max(x)+min(x))/2
    yc = (max(y)+min(y))/2
    zc = (max(z)+min(z))/2
    return xc[0],yc[0],zc[0]

def get_source(n,width,beamsize,wavelength):
    """Creates the photon source array for the setup."""
    height = 0.5 #Height of slit (inches)
    
    x = np.random.uniform(-0.5,0.5,size=n)*width
    y = np.zeros(n)
    z = np.empty(n)

    i = 0
    while i < n:  #Clip photons outside of source region
        val = np.random.normal(0,scale=beamsize)
        if not(val > height/2 or val < -height/2):
            z[i] = val
            i += 1

    pos =  np.vstack((x,y,z)).transpose()
    
    theta = np.linspace(-np.pi/4,np.pi/4,100000) 
    prob = np.sinc(width*np.sin(theta)/wavelength)**2 #Diffraction pattern
    prob /= sum(prob) #Ensuring probabilities add to one

    pthetas = np.random.choice(theta,size=n,p=prob) #Thetas for each photon
    adj = 1.56407 #Mirror alignment correction
    tmp =  np.vstack((np.cos(pthetas),np.sin(pthetas),np.zeros(n))).transpose()
    
    dirn = np.empty((n,3))
    for i in range(n):
        tmpx = tmp[i,0]*np.cos(adj)-tmp[i,1]*np.sin(adj)
        tmpy = tmp[i,0]*np.sin(adj)+tmp[i,1]*np.cos(adj)
        tmpz = 0
        dirn[i] = tmpx,tmpy,tmpz
    
    return pos,dirn

def mirror(mesh_no,tol=0.005,regen=False):
    """Returns the surface and color arrays for the mesh of either the open (142) or closed (141) mirrors. 
    ====================
    Args:
    mesh_no - Mesh number for mirror (141 or 142)
    tol     - Tolerance for polygon-mirror approximation
    """
    gold_mirror = get_mirror_surface()
    mesh = mesh_from_stl('stls/vuv_setup - Mesh %s.stl' % mesh_no)
    mesh_resolution = len(mesh.get_triangle_centers()) #Get number of triangles in mesh
    config_path = 'config/mirror_properties_%s.pickle' % mesh_no
    remake_flag = False #Determines whether mirror properties file must be regenerated (takes a long time)
    if os.path.isfile(config_path):
        try:
            config_file = open(config_path,'rb')
            data = cPickle.load(config_file)
            config_file.close()
        except:
            data = []
        if len(data) == mesh_resolution: #Checks load discrepancy within tolerance
            mirror_surf = data[:,0]
            mirror_col = data[:,1]
            return mirror_surf,mirror_col
        else:
            remake_flag = True
    else:
        remake_flag = True
    
    if remake_flag == True:
        print 'Regenerating mirror data'
        filename = 'stls/vuv_setup - Mesh '+mesh_no+'.stl'
        mesh = mesh_from_stl(filename)

        centers = mesh.get_triangle_centers()

        x = centers[:,0]
        y = centers[:,1]
        z = centers[:,2]

        in_rad = 0.248  #Inner wall radius
        out_rad = 0.5   #Outer wall radius

        center_x = (max(x)+min(x))/2
        center_y = (max(y)+min(y))/2 

        mirror_surf = []    #Surface array
        mirror_col = []     #Color array

        #Constructing surface and color arrays for mesh
        for i in range(len(centers)):
            radius = ((x[i]-center_x)**2+(y[i]-center_y)**2)**0.5

            if int(mesh_no) == 141 or int(mesh_no) == 250:   #Closed Mirror
                condition = (abs(radius-out_rad) < tol) or z[i] == max(z)
            elif int(mesh_no) == 142: #Open Mirror
                condition = (abs(radius-in_rad) < tol) or (abs(radius-out_rad) < tol) or z[i] == min(z)
            
            elif int(mesh_no) == 221:   #Monochromator Mirror
                condition = False
            
            #Applying condition to mesh:
            if condition == True: 
                mirror_surf.append(black_surface)
                mirror_col.append(0x7d6b56) #grey
            elif condition == False:
                mirror_surf.append(gold_mirror)
                mirror_col.append(0xe6ad12) #gold
        
        config_file = open(config_path,'wb')
        cPickle.dump(np.vstack((mirror_surf,mirror_col)).transpose(),config_file)
        config_file.close()
        print 'New mirror properties config file created.'
        return mirror_surf, mirror_col

    raise RuntimeError, 'Mirror properties generation failed. mfunctions.mirror() is somehow out of loop! Have you broken something?'
    return None

def lorentz(x,g):
    """Return normalized Lorentz (Cauchy) profile with HWHM g
    CURRENTY DEPRECATED
    """
    f = (1/(np.pi))*(g/(x**2+g**2))
    return f/sum(f)

def lineshape(radius,hwhm,n=10000):
    """Returns Lorentzian lineshape position density in polar coordinates.
    CURRENTLY DEPRECATED
    """
    if hwhm > radius:
        raise Exception('HWHM cannot exceed radius of lineshape.')
    
    x = np.linspace(0,radius,n)

    theta = np.zeros(n)
    r = np.zeros(n)

    for i in range(len(r)):
        theta[i] = np.random.random()*np.pi*2 #Uniformly random polar angle
        r[i] = np.random.choice(x,p=lorentz(x,hwhm)) #Radius sampled from Lorentz distribution

    return r,theta 
   
def laser(radius,hwhm,direction=(0,0,1),size=None,dtype=np.double):
    """
    Creates a collimated laser source (positions) with Lorentzian (Cauchy) distribution.
    CURRENTLY DEPRECATED
    """
    r,theta = lineshape(radius,hwhm,size) #Theta aziumthal

    if np.equal(direction, (0,0,1)).all():
        rotation_axis=(0,0,1)
        rotation_angle = 0.0
    else:
        rotation_axis = np.cross((0,0,1),direction)
        rotation_angle = -np.arccos(np.dot(direction, (0,0,1))/np.linalg.norm(direction))


    if size == None:
        raise Exception("laser: size is None!")
    
    points = np.empty((size,3),dtype)

    points[:,0] = r*np.cos(theta)
    points[:,1] = r*np.sin(theta)
    points[:,2] = np.zeros(size)

    return rotate(points,rotation_angle,rotation_axis)

def Rx(t):
    """Rotation matrix for angle t about x-axis"""
    return np.array(((1,0,0),(0,np.cos(t),-np.sin(t)),(0,np.sin(t),np.cos(t)))) 

def Ry(t):
    """Rotation matrix for angle t about y-axis"""
    return np.array(((np.cos(t),0,np.sin(t)),(0,1,0),(-np.sin(t),0,np.cos(t)))) 

def Rz(t):
    """Rotation matrix for angle t about z-axis"""
    return np.array(((np.cos(t),-np.sin(t),0),(np.sin(t),np.cos(t),0),(0,0,1))) 


