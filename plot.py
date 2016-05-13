#!/usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as mp
import numpy as np
from astropy import constants as cst
from astropy import units as u
# from scipy import integrate.simps as simp

# wave bands are defined at www.astrophysicsspectator.com/topics/observation/MagnitudesAndColors.html

curve = lambda l,t: 2*cst.h*cst.c**2/l**5  *  1/(np.exp( cst.h*cst.c/(l*cst.k_B*t*u.K*u.m) ) - 1 )

def mags(wv1, wv2, temp):
    sum = curve(wv1, temp) + curve(wv2, temp)
    h = (wv2 - wv1) / 1000
    for i in range(1, 1000, 2):
        sum += 4 * curve(wv1 + i * h )
    for i in range(2, 999, 2):
        sum += 2 * curve(wv1 + i * h )
    sum = s * 1e-3 / 3

mu = np.linspace(1e-7, 5e-6, 1e3, endpoint=True)
temp = 5000
flux = curve(mu, temp)
print (flux)
xax = np.linspace(0, 5, 1e3, endpoint=True)
mp.plot(xax, flux)
mp.xlabel("Wavelength (um)")
mp.ylabel("Flux (J/s*m^3*sr)")
mp.show()

