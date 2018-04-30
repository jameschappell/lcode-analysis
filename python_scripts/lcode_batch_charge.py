import numpy as np
import os
import string
import argparse


lcode_cfg = '''
# Simulation area:
geometry = c
window-width = 6;       r-step = 0.01
window-length = 7;     xi-step = 0.01
time-limit = 1000000.1; time-step = 50

continuation = n # Plasma continuation (no/beam/longplasma, n/y/Y)
# Particle beams:
beam-current = -0.059
# 1 kA

rigid-beam = n
beam-substepping-energy = 2
focusing = n # External beam focusing (no/cosine/rectangular, n/c/r)
foc-period = 100
foc-strength = 0.01
beam-particles-in-layer = 500
beam-profile = """
xishape=l, ampl=8.71, length=0.60, rshape=g, radius=0.5, angshape=l, 
angspread=1e-16, energy=8.81e15, vshift=0, eshape=m, espread=0, m/q=1836e10
xishape=l, ampl=0.00, length=3.20, rshape=g, radius=1.000, angshape=l, angspread=1e-16, energy=1000, vshift=0, eshape=m, espread=0, m/q=1
xishape=c, electron_beam_charge, length=0.80, rshape=g, radius=0.034, angshape=l, angspread=3.2e-4, energy=1000, vshift=0, eshape=m, espread=0, m/q=1
"""
beam-tune-charge = n;       beam-regulate-layer = n
rng-seed = 1
# Laser beams:
laser = n
# Plasma:
plasma-model = P # Plasma model (fluid/particles/newparticles, f/p/P)
magnetic-field = 0
magnetic-field-type = c;    magnetic-field-period = 200
plasma-particles-number = 1000
plasma-profile = 1 # Initial profile (1-6, uniform/stepwise/gaussian/arbitrary/channel/subchannel)
plasma-zshape = ""
plasma-width = 5
plasma-width-2 = 1
plasma-density-2 = 0.5
plasma-temperature = 0
ion-model = y # Model of plasma ions (mobile/background/absent/equilibrium, Y/y/n/N)
ion-mass = 1836
substepping-depth = 3
substepping-sensivity = 0.2
trapped-path-limit = 0 # Path limit for trapped plasma electrons (0 to disable)
viscosity = 0 # Artificial viscosity
# Every-time-step diagnostics:
indication-line-format = 1 # On-screen indication line format (eacht/eachxi)
output-Ez-minmax = n;       output-Ez-local = n
output-Phi-minmax = n;      output-Phi-local = n
output-beam-survival = n;   output-survival-xi = "-2, -4, -6, -8, -10, -12"
output-beam-slices = n;     output-slices-xi = "-2, -4, -6, -8, -10, -12"
write-beam-particles = y;   write-beam-particles-each = 40
write-beam-particles-from = 0;  write-beam-particles-to = -5
output-lost-particles = n
write-beam-particles-q-m-from = 0.5;  write-beam-particles-q-m-to = 1.5
# Periodical diagnostics:
output-time-period = 20000
#  Colored maps: (Er,Ef,Ez,Phi,Bf,Bz,pr,pf,pz,pri,pfi,pzi
#                 nb,ne,ni,Wf,dW,SEB,Sf,Sf2,Sr,Sr2,dS,dS2):
colormaps-full = ""
colormaps-subwindow = "Ez,nb"
colormaps-type = n
drawn-portion = 1 # Drawn portion of the simulation window
subwindow-xi-from = -0;     subwindow-xi-to = -7
subwindow-r-from = 0;       subwindow-r-to = 5
output-reference-energy = 1000
output-merging-r = 1;       output-merging-z = 1
output-smoothing-r = 0;     output-smoothing-z = 0
palette = b # Colormaps palette (default/greyscale/hue/bluewhitered, d/g/h/b)
                E-step = 0.01;                  nb-step = 0.08
              Phi-step = 0.01;                  ne-step = 0.02
               Bf-step = 0.01;                  ni-step = 0.01
               Bz-step = 0.1;                flux-step = 0.01
 electron-momenta-step = 0.1;    r-corrected-flux-step = 0.01
      ion-momenta-step = 0.1;              energy-step = 0.01
#  Output of various quantities as functions of xi:
#   (ne,nb,Ez,<Ez>,Bz,Phi,pz,emitt,dW,Wf,ni,pzi)
#   (nb2,Er,Ez2,Bf,Bz2,Fr,pr,pf,<rb>,dS,Sf,SEB,pri,pfi,Ef)
f(xi) = 
f(xi)-type = 
axis-radius = 0;        auxillary-radius = 0.2
               E-scale = 0.4;                 nb-scale = 1
             Phi-scale = 0.4;               ne-scale = 2
              Bz-scale = 0.4;               ni-scale = 1
electron-momenta-scale = 1;             flux-scale = 1
     ion-momenta-scale = 1;           energy-scale = 1
     beam-radius-scale = 5;        emittance-scale = 5
#  Beam particle information as pictures (r,pr,pz,M):
output-beam-particles = ""
draw-each = 20
beam-picture-height = 300
beam-pr-scale = 100
beam-a-m-scale = 100;       beam-pz-scale = 20000
# Output of beam characteristics in histogram form (r,z,M,a):
histogram-output = ""
histogram-output-accel = ""
histogram-type = y
histogram-bins = 300;       beam-angle-scale = 0.1
#  Trajectories of plasma particles:
trajectories-draw = n
trajectories-each = 10;     trajectories-spacing = 10
trajectories-min-energy = 1;    trajectories-energy-step = 0.5
#  Detailed (substepped) plasma response:
substepping-output-depth = 4
substepping-output-map = n
substepping-output-f(xi) = n
substepping-output-particles = n
substepping-output-particles-area = f
# Saving run state periodically:
saving-period = 200000
save-beam = n
save-plasma = n
# Logging preferences (error/warning/info/debug, e/w/i/d):
log-stdout-level = d;       log-file-level = w
log-to-file = y;        log-filename = lcode.log
log-file-clean = y
save-config = y;        save-config-filename = lcode.runas.cfg

# Generated with LCODE version trunk/678

'''

subscript = '''
#!/bin/bash -l

#$ -S /bin/bash

module unload mpi
module load mpi/openmpi/3.0.0/intel-2017

#$ -l h_rt=02:30:00
#$ -l mem=4G
name_def
#$ -m be
#$ -M james.chappell.17@ucl.ac.uk
work_directory_def
#$ -e ./logs/
#$ -o ./logs/
#$ -pe mpi 250

change_dir_def

export PATH="/home/ucapjch/Binaries/:$PATH"

gerun /home/ucapjch/Binaries/lcode
'''


def current(n_0, sigma_x, sigma_y):

    c = 299792458   # m/s
    e = 1.6e-19

    current_val = 2 * np.pi * e * c * n_0 * sigma_x * sigma_y

    return current_val


def beam_dens(charge, sigma_x, sigma_y, sigma_z):

    # Calculate number of particles from beam charge [measured in pC]

    N = charge * 1e-12 / 1.6e-19

    beam_density = (1/(2*np.pi * np.sqrt(2*np.pi))) * (N / (sigma_x * sigma_y *
                                                         sigma_z))

    return beam_density


def make_environment(resultsdir, charge, sigma_x, sigma_y, sigma_z):

    density = beam_dens(charge, sigma_x, sigma_y, sigma_z)
    current_val = current(density, sigma_x, sigma_y)
    normalised_current = current_val / 17e3

    sg = os.path.join(resultsdir, 'lcode.cfg')
    fh = open(sg, "wb")

    beam_charge_string = 'ampl=' + format(normalised_current, '.3f')
    lcode_cfg1 = string.replace(lcode_cfg, 'electron_beam_charge',
                                beam_charge_string)

    fh.write(lcode_cfg1)
    fh.close()

    ss = os.path.join(resultsdir, 'sub_script.bash')
    fs = open(ss, "wb")

    name_string = '#$ -N charge_' + format(charge, 'g') + '_pC'
    subscript1 = string.replace(subscript, 'name_def', name_string)

    work_dir_string = '#$ -wd ' + os.getcwd() + '/' + res_dir
    subscript2 = string.replace(subscript1, 'work_directory_def',
                                work_dir_string)

    change_dir_string = 'cd ' + os.getcwd() + '/' + res_dir
    subscript3 = string.replace(subscript2, 'change_dir_def', change_dir_string)

    fs.write(subscript3)
    fs.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
        This script generates a number of input scripts for LCODE 
        simulations that investigate beam loading effects and dephasing 
        between a proton drive beam and electron witness beam. The proton 
        beam uses parameters from the LIU-SPS.""",
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)

    parser.add_argument('--min_charge', dest='min_charge', default=None,
                        help='''
        This is the minimum charge for the witness electron beam to be 
        scanned over. It is measured in pico-Coulombs.''')

    parser.add_argument('--max_charge', dest='max_charge', default=None,
                        help='''
            This is the maximum charge for the witness electron beam to be 
            scanned over. It is measured in pico-Coulombs.''')

    parser.add_argument('--charge_step', dest='charge_step', default=None,
                        help='''
            This is the size of the step in charge for the witness electron 
            beam to be scanned over. It is measured in pico-Coulombs.''')


    arguments = parser.parse_args()

    min_charge = float(arguments.min_charge)
    max_charge = float(arguments.max_charge)
    charge_step = float(arguments.charge_step)

    sigma_x = 6.74e-6
    sigma_y = 6.74e-6
    sigma_z = 60e-6

    charge_values = np.arange(min_charge, max_charge + 0.0001, charge_step)

    cwd = os.getcwd()

    for charge in charge_values:

        """Looping over the range of electron beam charges, creating 
        a new directory for each simulation."""

        os.chdir(cwd)
        print "Beam Charge: %s" % (charge)
        res_dir = 'charge_' + format(charge, 'g') + '_pC'
        print "Making directory: ", res_dir
        if os.path.isdir(res_dir) is False:
            os.mkdir(res_dir)
        else:
            print "Directory %s already exists. Stopping." % (res_dir)
            break

        make_environment(res_dir, charge, sigma_x, sigma_y, sigma_z)
        os.chdir(res_dir)
        os.mkdir('logs')
        run_command = "qsub sub_script.bash"
        print run_command
        os.system(run_command)
