#   VUV RAYTRACING SCRIPT FOR TRIUMF
#   Contact: Ben Page (benjaminpage.acer@gmail.com)
#   Last Updated: Sep 4 2019

"""
This script is the terminal frontend for the simulation code. It handles and filters user inputs. Please read the associated documentation for further info and/or contact the email address listed.
"""
import os
import subprocess
import sys
import numpy as np

try:
    import builder
    import simulator
except:
    if os.path.exists('/singularity'):
        raise Exception('Chroma container is loaded, but it seems that builder.py and simulator.py cannot be loaded! Are they present?\nPlease consult the documentation.')
        exit()
    else:
        print 'Chroma container not loaded! To avoid errors, please enter the container before running the script!\n'
        expression = 'singularity shell --nv chroma.img'
        print 'REMINDER: The command is: ',expression
        exit()

def main(path=None,regen=False):

    updated = 'Sep 4 2019'
    #ASCII graphics generation
    logo = open('config/logo.txt','r')
    print logo.read()
    logo.close()
    print "-"*90 
    print "CHROMA RAYTRACER FOR VUV SETUP".center(90,'=')
    print "-"*90+'\n'
    print "For help, please consult the user manual\nor contact Ben Page at benjaminpage.acer@gmail.com\n"
    print "Last Updated " + updated
    print "="*90

    if path!=None:
        print 'Loading in predefined parameters...'
        params = open(path,'r')
        expression = params.readline()
        print expression
        n,wavelength,width,beamsize = [float(i) for i in expression.split(':')] 
    else:
        print 'PLEASE SPECIFY THE SETUP PARAMETERS'
        print '[You can exit at any time by typing "exit"]\n'
        
        flag = False
        while flag == False:
            statements = [
                'Input number of photons: ',
                'Input wavelength (nm): ',
                u'Input slit width (\u03bcm): '.encode('utf-8'), 
                u'Input width of light source (\u03bcm): '.encode('utf-8')]
            values = []
            for pos, statement in enumerate(statements):
                while True:
                    usrin = raw_input(statement).strip()
                    if usrin in ['Exit','exit']:
                        exit()
                    else:
                        try:
                            values.append(float(usrin))
                            break
                        except:
                            print 'Bad input! Please enter a type int/float!'

            n,wavelength,width,beamsize = values
            if width > 1000:
                print u'Width cannot be greater than 1000\u03bcm'.encode('utf-8')
            else:
                print '\nSelected Inputs:\n'+'='*20
                print 'Number of photons = ',n
                print 'Wavelength = ',wavelength
                print 'Slit width = ',width
                print 'Beam Size = ',beamsize,'\n'
    
                while True:
                    usrin = raw_input('Proceed? [Y/N] : ').strip()
                    if usrin in ['yes','Yes','Y','y']:
                        flag = True
                        break
                    elif usrin in ['no','No','N','n']:
                        break
                    else: 
                        print 'Bad Input: Please confirm your selection.'
    
    #Cleaning inputs       
    n = int(n)
    wavelength = int(wavelength)
    width /= 25400.0 #Metric to Imperial conversion
    beamsize /= 25400.0
               
    simulator.main(n,wavelength,width,beamsize,regen)        
    
if __name__ == '__main__':

    if len(sys.argv) == 1:
        main(path=None)
    else:
        main(path=sys.argv[1])


