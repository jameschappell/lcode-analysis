import numpy as np
import os

lcode_cfg = '''
# Simulation area:
geometry = c
window-width = 6;       r-step = 0.01
window-length = 7;     xi-step = 0.01
time-limit = 1000000.1; time-step = 25

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
xishape=c, ampl=0.50, length=0.80, rshape=g, radius=0.034, angshape=l, angspread=3.2e-4, energy=1000, vshift=0, eshape=m, espread=0, m/q=1
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
write-beam-particles = y;   write-beam-particles-each = 20
write-beam-particles-from = -2;  write-beam-particles-to = -5
output-lost-particles = n
write-beam-particles-q-m-from = 0;  write-beam-particles-q-m-to = 0
# Periodical diagnostics:
output-time-period = 10000
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
saving-period = 100000
save-beam = n
save-plasma = n
# Logging preferences (error/warning/info/debug, e/w/i/d):
log-stdout-level = d;       log-file-level = w
log-to-file = y;        log-filename = lcode.log
log-file-clean = y
save-config = y;        save-config-filename = lcode.runas.cfg

# Generated with LCODE version trunk/678

'''


def beam_dens_from_current(I_max, sigma):

    c = 299792458   # m/s
    e = 1.6e-19
    n = 2.8e21      # m^-3
    e_0 = 8.85e-12
    m_e = 9.11e-31

    w_p = np.sqrt((n * e**2)/(e_0 * m_e))
    sigma_x = sigma
    sigma_y = sigma

    peak_beam_density = I_max/(2 * np.pi * e * c * sigma_x * sigma_y)

    return peak_beam_density


def beam_dens(charge, sigma_x, sigma_y, sigma_z):

    # Calculate number of particles from beam charge [measured in pC]

    N = charge * 1e-12 / 1.6e-19

    beam_density = (1/(2*np.pi * np.sqrt(2*np.pi))) * (N / (sigma_x * sigma_y *
                                                         sigma_z))

    return beam_density

