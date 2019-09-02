import numpy as np
import os
import string
import argparse

#densities = np.linspace(1e14, 1e15, 10)
#dens14 = np.linspace(1e14, 9e14, 9)
#dens15 = np.linspace(1e15, 9e15, 9)
#dens16 = np.linspace(1e16, 1e17, 10)
#densities = np.concatenate((dens14, dens15, dens16))

#print densities

# define constants
#m_e = 9.11e-31
#e = 1.6e-19
#e0 = 8.85e-12
#c = 3e8

# beam parameters
#sigr = 20e-6    # sigma r
#sigt = 200e-15  # bunch length (time)
#sigz = sigt * c
#Qb = 200e-12    # bunch charge
#N = Qb/e        # number of electrons

#for i in range(len(densities)):

#    ne = densities[i] * 1e6

#    omega_p = np.sqrt((ne * e**2)/(m_e * e0))
#    kp = omega_p/c
#    inv_kp = 1/kp

#    wbfield = m_e * c * omega_p / e

#    br = sigr * kp
#    bz = sigz * kp

#    beamdens = 1/(2*np.pi)**(1.5) * N/(sigr**2 * sigz)

#    I = 2 * np.pi * e * c * beamdens * sigr * sigz

#    normI = I/17e3

#    print '{:.2e} {:.2e} {:.2e} {:.2e} {:.2e} {:.2e}'.format(ne/1e6, inv_kp,
#  br,
 #                                                      bz, beamdens/1e6, normI)


def make_simulations(density, energy, energyspread, sigr, sigz,
                     charge, directory):

    # define constants
    m_e = 9.11e-31
    e = 1.6e-19
    e0 = 8.85e-12
    c = 3e8

    # calculate values needed in simulation
    ne = density * 1e6
    N = charge/e

    omega_p = np.sqrt((ne * e ** 2) / (m_e * e0))
    kp = omega_p / c
    inv_kp = 1 / kp

    wbfield = m_e * c * omega_p / e

    br = sigr * kp
    bz = sigz * kp

    beamdens = 1 / (2 * np.pi) ** (1.5) * N / (sigr ** 2 * sigz)

    I = 2 * np.pi * e * c * beamdens * sigr * sigz

    normI = -I / 17e3

    normE = 0.511e-3

    beamenergy = energy / normE

    beamenergyspread = beamenergy * energyspread

    # update the simulation configuration file

    print "Making simulation for density {:.2e}".format(density)
    sg = os.path.join(directory, 'lcode.cfg')
    fh = open(sg, "wb")

    beamlength_string = str(bz)
    lcodecfg1 = string.replace(lcodecfg, 'beamlengthval', beamlength_string)

    beamrad_string = str(br)
    lcodecfg2 = string.replace(lcodecfg1, 'beamradval', beamrad_string)

    beamenergy_string = str(beamenergy)
    lcodecfg3 = string.replace(lcodecfg2, 'beamenergyval', beamenergy_string)

    beamespread_string = str(beamenergyspread)
    lcodecfg4 = string.replace(lcodecfg3, 'beamespreadval', beamespread_string)

    beamcurrent_string = str(normI)
    lcodecfg5 = string.replace(lcodecfg4, 'beamcurrentval', beamcurrent_string)

    fh.write(lcodecfg5)
    fh.close()

    # update the submission file

    ss = os.path.join(directory, 'sub_script.bash')
    fs = open(ss, "wb")

    work_dir_string = directory
    subscript1 = string.replace(subscript, 'DEFworkdir', work_dir_string)
    subscript2 = string.replace(subscript1, 'DEFchangedir', work_dir_string)

    name_string = 'dens_' + format(density, '.2e')
    subscript3 = string.replace(subscript2, 'DEFname', name_string)

    fs.write(subscript3)
    fs.close()


lcodecfg = '''

# Simulation area:
geometry = c
window-width = 5;		r-step = 0.005
window-length = 10;		xi-step = 0.005
time-limit = 1.1;		time-step = 1

# Particle beams:
beam-current = beamcurrentval
rigid-beam = n
beam-substepping-energy = 810000
beam-particles-in-layer = 160
beam-profile = """
xishape=c, ampl=1.0, length=beamlengthval, rshape=g, radius=beamradval, 
angshape=l, angspread=4.5e-5, energy=beamenergyval, m/q=1, eshape=g, 
espread=beamespreadval
"""

# Plasma:
plasma-model = P # Plasma model (fluid/particles/newparticles, f/p/P)
plasma-particles-number = 40000
plasma-profile = 1 # Initial profile (uniform/stepwise/gaussian/arbitrary/channel)
plasma-width = 5
plasma-temperature = 0
ion-model = Y # Model of plasma ions (mobile/background/absent/equilibrium, Y/y/n/N)
ion-mass = 73382
substepping-depth = 3
substepping-sensivity = 0.2

# Every-time-step diagnostics:
indication-line-format = 1 # On-screen indication line format (eacht/eachxi)
output-Ez-minmax = y;		output-Ez-local = y
output-Phi-minmax = y;		output-Phi-local = y
output-lost-particles = y

# Periodical diagnostics:
output-time-period = 1

#  Colored maps: (Er,Ef,Ez,Phi,Bf,Bz,pr,pf,pz,pri,pfi,pzi
#                 nb,ne,ni,Wf,dW,SEB,Sf,Sf2,Sr,Sr2,dS,dS2):
colormaps-full = ""
colormaps-subwindow = "Er,Ez,Phi,Bf,ne,ni,nb,pz"
colormaps-type = n
drawn-portion = 1 # Drawn portion of the simulation window
subwindow-xi-from = 0;		subwindow-xi-to = -10
subwindow-r-from = 0;		subwindow-r-to = 5
output-reference-energy = 880000
output-merging-r = 1;		output-merging-z = 5
palette = d # Colormaps palette (default/greyscale/hue/bluewhitered, d/g/h/b)
                E-step = 0.059;	               nb-step = 0.00056
              Phi-step = 0.059;	               ne-step = 0.1
               Bf-step = 0.059;	               ni-step = 0.01
               Bz-step = 0.059;	             flux-step = 0.02
 electron-momenta-step = 0.1;	 r-corrected-flux-step = 0.02
      ion-momenta-step = 0.1;	           energy-step = 10

#  Output of various quantities as functions of xi:
#   (ne,nb,Ez,<Ez>,Bz,Phi,pz,emitt,dW,Wf,ni,pzi)
#   (nb2,Er,Ez2,Bf,Bz2,Fr,pr,pf,<rb>,dS,Sf,SEB,pri,pfi,Ef)
f(xi) = ne,nb,Ez,Eza,Phi
f(xi)-type = Y
axis-radius = 0;		auxillary-radius = 1
               E-scale = 0.59;	              nb-scale = 0.02
             Phi-scale = 0.59;	              ne-scale = 2
              Bz-scale = 0.59;	              ni-scale = 0.1
electron-momenta-scale = 0.5;	            flux-scale = 0.5
     ion-momenta-scale = 0.5;	          energy-scale = 1
     beam-radius-scale = 5;	       emittance-scale = 300

#  Beam particle information as pictures (r,pr,pz,M):
output-beam-particles = r,pr
draw-each = 1
beam-picture-height = 900
beam-pr-scale = 1000
beam-a-m-scale = 1000;		beam-pz-scale = 15000

# Output of beam characteristics in histogram form (r,z,M,a):
histogram-output = ""
histogram-output-accel = ""
histogram-type = y
histogram-bins = 300;		beam-angle-scale = 0.02

#  Trajectories of plasma particles:
trajectories-draw = n
trajectories-each = 1;		trajectories-spacing = 10
trajectories-min-energy = 1;	trajectories-energy-step = 0.5

# Saving run state periodically:
saving-period = 1000
save-beam = n
save-plasma = n
'''

subscript = '''

#!/bin/bash -l

#$ -S /bin/bash

module unload mpi
module load mpi/openmpi/3.0.0/intel-2017

#$ -l h_rt=48:00:00
#$ -l mem=0.5G
### -l tmpfs=15G
#$ -N DEFname
#$ -m be
#$ -M james.chappell.17@ucl.ac.uk
#$ -wd DEFworkdir
#$ -e ./logs/
#$ -o ./logs/
#$ -pe mpi 32

cd DEFchangedir

export PATH="/home/ucapjch/Binaries/:$PATH"

gerun /home/ucapjch/Binaries/lcode

'''



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    This script generates a number of LCODE simulations.""",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--energy', dest='energy', default=None,
                        help='''
    This is the value of the electron beam energy in GeV.
    
    e.g. --energy 1.0"''')

    parser.add_argument('--energyspread', dest='energyspread', default=None,
                        help='''
        This is the value of the relative electron beam energy spread.

        e.g. --energyspread 0.01 is a 1% energy spread"''')

    parser.add_argument('--length', dest='length', default=None,
                        help='''
            This is the bunch length, measured in seconds. 

            e.g. --length 200e-15 is 200fs bunch."''')

    parser.add_argument('--radius', dest='radius', default=None,
                        help='''
                This is the bunch radius, measured in metres. 

                e.g. --length 20e-6 is 20um radius bunch."''')

    parser.add_argument('--charge', dest='charge', default=None,
                        help='''
                    This is the bunch charge, measured in C. 

                    e.g. --charge 200e-12 is a 200pC bunch."''')

    arguments = parser.parse_args()

    energy = float(arguments.energy)
    energyspread = float(arguments.energyspread)
    sigz = float(arguments.length) * 3e8
    sigr = float(arguments.radius)
    charge = float(arguments.charge)

    dens14 = np.linspace(1e14, 9e14, 9)
    dens15 = np.linspace(1e15, 9e15, 9)
    dens16 = np.linspace(1e16, 1e17, 10)
    densities = np.concatenate((dens14, dens15, dens16))

    cwd = os.getcwd()

    for i in range(len(densities)):

        density = densities[i]
        os.chdir(cwd)
        print '{:.2e}'.format(density)
        directory = 'dens_' + format(density, '.2e')
        print "Making directory: ", directory
        if os.path.isdir(directory) is False:
            os.mkdir(directory)
        else:
            print "Directory %s already exists. Stopping." % (directory)
            break
        print "Copy template files: "

        make_simulations(density, energy, energyspread, sigr, sigz,
                         charge, cwd + "/" + directory)
        os.chdir(directory)
        os.mkdir('logs')
        run_command = "qsub sub_script.bash"
        print run_command
        os.system(run_command)
